from __future__ import annotations

from dataclasses import asdict, dataclass
import numpy as np
import pandas as pd
from scipy import stats


def connection_mask(parcels, exclude_diagonal=True, upper=False):
    if upper:
        return np.triu(np.ones((parcels, parcels), bool), k=1 if exclude_diagonal else 0)
    return ~np.eye(parcels, dtype=bool) if exclude_diagonal else np.ones((parcels, parcels), bool)


def normalise_tdmi(tensor, epsilon=1e-12, exclude_diagonal=True):
    a = np.asarray(tensor, float)
    tdmi = a.sum(axis=-1)
    valid = np.isfinite(tdmi) & (np.abs(tdmi) >= epsilon)
    out = np.full_like(a, np.nan)
    np.divide(a, tdmi[..., None], out=out, where=valid[..., None])
    mask = connection_mask(a.shape[0], exclude_diagonal)
    return out, int(np.sum(~valid & mask))


def regional_sums(tensor, exclude_diagonal=True, axis=1):
    a = np.array(tensor, dtype=float, copy=True)
    if exclude_diagonal:
        a[np.arange(len(a)), np.arange(len(a)), :] = 0
    return a.sum(axis=axis)


def network_means(tensor, network_ids, exclude_diagonal=True, edges="both"):
    network_ids = np.asarray(network_ids)
    values = []
    for net in sorted(set(network_ids)):
        idx = np.flatnonzero(network_ids == net)
        block = tensor[idx[:, None], idx, :]
        mask = connection_mask(len(idx), exclude_diagonal, upper=edges == "upper")
        values.append(block[mask].mean(axis=0) if mask.any() else np.full(tensor.shape[-1], np.nan))
    return np.asarray(values)


def global_tdmi(tensor, exclude_diagonal=True):
    mask = connection_mask(len(tensor), exclude_diagonal)
    return float(tensor[mask].sum(axis=-1).mean())


@dataclass
class Design:
    confounds: list[str]
    columns: list[str]
    medians: dict[str, float]
    reference: list[float]
    missing: str
    age: str

    @classmethod
    def fit(cls, df, transform):
        medians = {}
        columns = []
        for name in transform.confounds:
            if name == "sex":
                columns.append("sex_Male")
                if transform.missing == "median":
                    columns.append("sex_missing")
            else:
                source = "AGE" if name == "visit_age" and transform.age == "baseline" else name
                vals = pd.to_numeric(df[source], errors="coerce").replace([np.inf, -np.inf], np.nan)
                if not vals.notna().any():
                    raise ValueError(f"No observed values for confound {name}")
                medians[name] = float(vals.median())
                columns.append(name)
        model = cls(list(transform.confounds), columns, medians, [], transform.missing, transform.age)
        C = model.encode(df)
        model.reference = C.mean(axis=0).tolist()
        return model

    def encode(self, df):
        cols = []
        for name in self.confounds:
            if name == "sex":
                if self.missing == "complete" and df.sex.isna().any():
                    raise ValueError("Missing sex in complete-case design")
                cols.append(df.sex.eq("Male").fillna(False).to_numpy(dtype=float))
                if self.missing == "median":
                    cols.append(df.sex.isna().to_numpy(dtype=float))
            else:
                source = "AGE" if name == "visit_age" and self.age == "baseline" else name
                vals = pd.to_numeric(df[source], errors="coerce").replace([np.inf, -np.inf], np.nan)
                if self.missing == "median":
                    vals = vals.fillna(self.medians[name])
                if vals.isna().any():
                    raise ValueError(f"Missing {name} in complete-case design")
                cols.append(vals.to_numpy(dtype=float))
        return np.column_stack(cols) if cols else np.empty((len(df), 0))

    def matrix(self, df):
        # Fit with centred confounds for numerical stability; intercept is mean level.
        return np.column_stack([np.ones(len(df)), self.encode(df) - np.asarray(self.reference)])

    def to_dict(self):
        return asdict(self)


def complete_cases(df, transform):
    reasons = [[] for _ in range(len(df))]
    if transform.missing == "complete":
        for name in transform.confounds:
            source = "AGE" if name == "visit_age" and transform.age == "baseline" else name
            col = df[source] if source in df else pd.Series(np.nan, index=df.index)
            bad = col.isna().to_numpy() if name == "sex" else ~np.isfinite(pd.to_numeric(col, errors="coerce").to_numpy(dtype=float))
            for i in np.flatnonzero(bad):
                reasons[i].append(f"missing_{name}")
    return np.array([not r for r in reasons]), [";".join(r) for r in reasons]


def fit_coefficients(X, Y):
    """Multi-output regression grouped by feature validity, not Python objects."""
    Y = np.asarray(Y, float)
    if Y.ndim == 1:
        Y = Y[:, None]
    beta = np.full((X.shape[1], Y.shape[1]), np.nan)
    nobs = np.isfinite(Y).sum(axis=0)
    patterns, inverse = np.unique(np.isfinite(Y).T, axis=0, return_inverse=True)
    for number, valid in enumerate(patterns):
        columns = np.flatnonzero(inverse == number)
        rank = np.linalg.matrix_rank(X[valid]) if valid.any() else 0
        if valid.sum() <= rank:
            continue
        beta[:, columns] = np.linalg.lstsq(X[valid], Y[np.ix_(valid, columns)], rcond=None)[0]
    return beta, nobs


def apply_adjustment(Y, X, beta, reference, mode="mean"):
    prediction = X @ beta
    if mode == "residuals":
        return Y - prediction
    if mode == "intercept":
        # Historical intercept at uncentred C=0.
        intercept = beta[0] - np.asarray(reference) @ beta[1:]
        return Y - prediction + intercept
    return Y - prediction + beta[0]


def adjusted_features(train_df, test_df, train_y, test_y, transform):
    """Fit solely on training participants, then apply fixed parameters to held-out rows."""
    design = Design.fit(train_df, transform)
    Xtrain, Xtest = design.matrix(train_df), design.matrix(test_df)
    beta, nobs = fit_coefficients(Xtrain, train_y)
    return (apply_adjustment(train_y, Xtrain, beta, design.reference, transform.adjustment),
            apply_adjustment(test_y, Xtest, beta, design.reference, transform.adjustment),
            design, beta, nobs)


def cn_standardise(values, reference):
    mu = np.mean(reference, axis=0)
    sd = np.std(reference, axis=0, ddof=0)
    valid = np.isfinite(sd) & (sd > 0) & np.isfinite(mu)
    z = np.full_like(values, np.nan, dtype=float)
    np.divide(values - mu, sd, out=z, where=valid)
    return z, mu, sd, valid


def hierarchy(values):
    a = np.asarray(values)
    if not np.isfinite(a[:, [0, 15]]).all():
        return np.full(len(a), np.nan)
    return stats.rankdata(a[:, 15], method="average") - stats.rankdata(a[:, 0], method="average")


def correlation(x, y):
    x, y = np.asarray(x), np.asarray(y)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 3 or np.ptp(x[valid]) == 0 or np.ptp(y[valid]) == 0:
        return np.nan
    return float(np.corrcoef(x[valid], y[valid])[0, 1])


def fdr(pvalues, alpha=.05, correction="fdr_bh"):
    p = np.asarray(pvalues, float)
    q = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    if valid.any():
        pv = p[valid]
        if correction == "none":
            q[valid] = pv
        else:
            order = np.argsort(pv, kind="stable")
            adjusted = np.minimum.accumulate((pv[order] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
            out = np.empty(len(pv))
            out[order] = np.minimum(adjusted, 1)
            q[valid] = out
    return q, np.isfinite(q) & (q <= alpha)


def compare(a, b, welch=True):
    """All signed quantities are B minus A. Degenerate effects are never zero-filled."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    result = {"n_a": len(a), "n_b": len(b), "mean_a": np.mean(a) if len(a) else np.nan,
              "mean_b": np.mean(b) if len(b) else np.nan, "t": np.nan, "p": np.nan,
              "cohens_d": np.nan, "hedges_g": np.nan, "status": "ok"}
    result["difference"] = result["mean_b"] - result["mean_a"]
    if min(len(a), len(b)) < 2:
        result["status"] = "insufficient_samples"
        return result
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    if va + vb == 0:
        result["status"] = "zero_variance"
        return result
    test = stats.ttest_ind(b, a, equal_var=not welch)
    result["t"], result["p"] = float(test.statistic), float(test.pvalue)
    degrees = len(a) + len(b) - 2
    pooled = np.sqrt(((len(a) - 1) * va + (len(b) - 1) * vb) / degrees)
    if pooled > 0:
        result["cohens_d"] = result["difference"] / pooled
        result["hedges_g"] = result["cohens_d"] * (1 - 3 / (4 * degrees - 1))
    return result


def omnibus(groups, welch=True):
    samples = [np.asarray(g, float)[np.isfinite(g)] for g in groups]
    ns = np.array([len(g) for g in samples], float)
    result = {"F": np.nan, "p": np.nan, "df1": len(samples) - 1, "df2": np.nan, "eta_squared": np.nan, "status": "ok"}
    if len(samples) < 2 or np.min(ns) < 2:
        result["status"] = "insufficient_samples"
        return result
    means = np.array([g.mean() for g in samples])
    variances = np.array([g.var(ddof=1) for g in samples])
    if np.any(variances <= 0):
        result["status"] = "zero_variance"
        return result
    grand = np.average(means, weights=ns)
    between = np.sum(ns * (means - grand)**2)
    total = sum(np.sum((g - grand)**2) for g in samples)
    result["eta_squared"] = between / total if total else np.nan
    if welch:
        k = len(samples)
        weights = ns / variances
        weighted_mean = np.average(means, weights=weights)
        term = np.sum((1 - weights / weights.sum())**2 / (ns - 1))
        F = np.sum(weights * (means - weighted_mean)**2) / (k - 1) / (1 + 2 * (k - 2) * term / (k*k - 1))
        df2 = (k*k - 1) / (3 * term)
        result.update(F=float(F), df2=float(df2), p=float(stats.f.sf(F, k - 1, df2)))
    else:
        test = stats.f_oneway(*samples)
        result.update(F=float(test.statistic), p=float(test.pvalue), df2=float(ns.sum() - len(samples)))
    return result


def interval(samples):
    a = np.asarray(samples)
    a = a[np.isfinite(a)]
    return np.quantile(a, [.025, .975]) if len(a) >= 2 else np.array([np.nan, np.nan])


def permutation_p(observed, null, centre=0.):
    null = np.asarray(null)
    null = null[np.isfinite(null)]
    if not np.isfinite(observed) or not len(null):
        return np.nan
    return (1 + np.sum(np.abs(null - centre) >= abs(observed - centre))) / (len(null) + 1)
