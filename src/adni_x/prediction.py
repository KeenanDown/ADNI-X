from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, roc_auc_score
from sklearn.model_selection import train_test_split

from .analyses import context
from .numerics import (adjusted_features, apply_adjustment, normalise_tdmi, regional_sums,
                       hierarchy, correlation)
from .provenance import rng_for
from .results import InsufficientData, Result


def features(prepared, subjects):
    t = prepared.config.transform
    indices = subjects.tensor_row.to_numpy(int)
    if t.tdmi == "none":
        return np.asarray(prepared.array("raw_regional")[indices])
    raw = prepared.array("raw")
    return np.stack([regional_sums(normalise_tdmi(raw[i], t.epsilon, t.exclude_diagonal)[0],
                                  t.exclude_diagonal, t.regional_axis) for i in indices])


def split_indices(labels, settings, name):
    labels = np.asarray(labels)
    if len(set(labels)) < 2 or min(np.unique(labels, return_counts=True)[1]) < 2:
        raise InsufficientData("Stratified evaluation requires at least two participants per class")
    rng = rng_for(settings.seed, name)
    try:
        return train_test_split(np.arange(len(labels)), test_size=settings.test_fraction,
                                stratify=labels, random_state=int(rng.integers(2**31 - 1)))
    except ValueError as e:
        raise InsufficientData(str(e)) from e


def lda(prepared):
    df, atlas, settings = context(prepared, 4)
    if df.subject.duplicated().any():
        raise ValueError("Prediction input must have one scan per participant")
    raw = features(prepared, df).reshape(len(df), -1)
    binary = np.flatnonzero(df.analysis_group.isin([prepared.groups[0], prepared.groups[-1]]))
    train, test = split_indices(df.analysis_group.iloc[binary], settings, prepared.target("lda") + ".split")
    train, test = binary[train], binary[test]
    train_y, test_y, design, beta, nobs = adjusted_features(
        df.iloc[train], df.iloc[test], raw[train], raw[test], prepared.config.transform)
    keep = np.isfinite(train_y).all(axis=0) & (np.nanstd(train_y, axis=0) > 0)
    if not keep.any():
        raise InsufficientData("No finite, varying training features")
    model = LinearDiscriminantAnalysis(solver="svd")
    model.fit(train_y[:, keep], df.analysis_group.iloc[train])
    all_adjusted = apply_adjustment(raw, design.matrix(df), beta, design.reference, prepared.config.transform.adjustment)
    valid = np.isfinite(all_adjusted[:, keep]).all(axis=1)
    output = df[["scan_id", "subject", "analysis_group"]].rename(columns={"analysis_group": "group"})
    output["split"] = "projection_only"
    output.loc[train, "split"], output.loc[test, "split"] = "train", "test"
    output["prediction"] = pd.Series([None]*len(output), dtype="object")
    probability = f"probability_{prepared.groups[-1]}"
    auc = f"roc_auc_{prepared.groups[-1]}"
    output["decision"], output[probability] = np.nan, np.nan
    output["status"] = np.where(valid, "ok", "nonfinite_heldout_features")
    output.loc[valid, "prediction"] = model.predict(all_adjusted[valid][:, keep])
    output.loc[valid, "decision"] = model.decision_function(all_adjusted[valid][:, keep])
    probability_column = list(model.classes_).index(prepared.groups[-1])
    output.loc[valid, probability] = model.predict_proba(all_adjusted[valid][:, keep])[:, probability_column]
    held = output.loc[(output.split == "test") & valid]
    metrics = {"test_participants": len(test), "valid_test_participants": len(held), "accuracy": np.nan,
               "balanced_accuracy": np.nan, auc: np.nan,
               "status": "ok" if len(held) == len(test) else "incomplete_heldout_evaluation"}
    if len(held):
        metrics["accuracy"] = accuracy_score(held.group, held.prediction)
        metrics["balanced_accuracy"] = balanced_accuracy_score(held.group, held.prediction)
        if held.group.nunique() == 2:
            metrics[auc] = roc_auc_score(held.group == prepared.groups[-1], held[probability])
    matrix = confusion_matrix(held.group, held.prediction, labels=[prepared.groups[0], prepared.groups[-1]]) if len(held) else np.zeros((2, 2), int)
    result = Result(tables={"predictions": output, "metrics": pd.DataFrame([metrics]), "atlas": atlas},
                    arrays={"confusion": matrix, "feature_mask": keep, "confound_coefficients": beta,
                            "confound_nobs": nobs, "lda_coefficients": model.coef_, "lda_intercept": model.intercept_},
                    metadata={"cohort": "complete-case selected participants", "exploratory": True,
                              "confound_design": design.to_dict(), "classes": model.classes_.tolist(), "confusion_labels": [prepared.groups[0], prepared.groups[-1]],
                              "probability_class": prepared.groups[-1],
                              "fit_subjects": df.subject.iloc[train].tolist(), "test_subjects": df.subject.iloc[test].tolist(),
                              "feature_shape": [len(atlas), 16], "decision_positive_class": str(model.classes_[1])})
    result.figure("lda.confusion", "confusion")
    result.figure("lda.projections", "prediction_distribution", table="predictions", value="decision")
    return prepared.finish(result)


def heldout_hierarchy(prepared):
    df, atlas, settings = context(prepared, 6)
    raw = features(prepared, df).reshape(len(df), -1)
    train, test = split_indices(df.analysis_group, settings, prepared.target("heldout-hierarchy") + ".split")
    train_y, test_y, design, beta, nobs = adjusted_features(
        df.iloc[train], df.iloc[test], raw[train], raw[test], prepared.config.transform)
    train_y, test_y = train_y.reshape(-1, len(atlas), 16), test_y.reshape(-1, len(atlas), 16)
    cn = train_y[df.analysis_group.iloc[train].to_numpy() == prepared.groups[0]]
    if len(cn) < 2:
        raise InsufficientData(f"Held-out hierarchy requires at least two {prepared.groups[0]} training participants")
    reference = hierarchy(cn.mean(axis=0))
    if not np.isfinite(reference).all() or np.ptp(reference) == 0:
        raise InsufficientData("Training reference hierarchy is missing or constant")
    rows = []
    for indices, arr, split in ((train, train_y, "train"), (test, test_y, "test")):
        for index, values in zip(indices, arr):
            r = correlation(reference, hierarchy(values))
            rows.append({"scan_id": df.scan_id.iloc[index], "subject": df.subject.iloc[index],
                         "group": df.analysis_group.iloc[index], "split": split, "correlation": r,
                         "status": "ok" if np.isfinite(r) else "constant_or_missing_hierarchy"})
    result = Result(tables={"scores": pd.DataFrame(rows), "atlas": atlas},
                    arrays={"CN_training_hierarchy" if prepared.grouping == "diagnosis" else "reference_training_hierarchy": reference, "confound_coefficients": beta, "confound_nobs": nobs},
                    metadata={"cohort": "complete-case selected participants", "exploratory": True,
                              "confound_design": design.to_dict(), "fit_subjects": df.subject.iloc[train].tolist(),
                              "test_subjects": df.subject.iloc[test].tolist()})
    result.figure("heldout-hierarchy.correlations", "prediction_distribution", table="scores", value="correlation")
    return prepared.finish(result)

