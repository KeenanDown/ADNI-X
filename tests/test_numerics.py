from dataclasses import replace
from itertools import combinations
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from adni_x.config import Transform
from adni_x.numerics import (regional_sums, network_means, normalise_tdmi, global_tdmi,
    Design, fit_coefficients, apply_adjustment, adjusted_features, fdr, compare, omnibus, cn_standardise)


def test_asymmetric_reductions_and_diagonal():
    a = np.arange(4*4*16, dtype=float).reshape(4, 4, 16)
    original = a.copy()
    expected = a.copy(); expected[np.arange(4), np.arange(4)] = 0
    np.testing.assert_array_equal(regional_sums(a), expected.sum(axis=1))
    assert not np.array_equal(regional_sums(a), expected.sum(axis=0))
    np.testing.assert_allclose(network_means(a, [1, 1, 2, 2])[0], (a[0, 1] + a[1, 0]) / 2)
    assert global_tdmi(a) == a[~np.eye(4, dtype=bool)].sum(axis=-1).mean()
    np.testing.assert_array_equal(a, original)


def test_normalisation_is_distinct_and_invalid_masked():
    a = np.ones((3, 3, 16))
    a[0, 1] = 0
    normalised, count = normalise_tdmi(a)
    assert count == 1
    assert np.isnan(normalised[0, 1]).all()
    assert normalised[1, 0, 0] == 1/16
    _, _, sd, valid = cn_standardise(np.ones((3, 2, 16)), np.ones((3, 2, 16)))
    assert not valid.any() and (sd == 0).all()


def data():
    rng = np.random.default_rng(19)
    n = 50
    return pd.DataFrame({"sex": np.where(np.arange(n) % 2, "Female", "Male"),
        "education": rng.uniform(10, 20, n), "visit_age": rng.uniform(60, 85, n),
        "AGE": rng.uniform(60, 85, n), "ICV": rng.uniform(1e6, 2e6, n)})


def test_known_confound_coefficients_mean_preservation_and_all_toggles():
    df = data()
    names = ["sex", "education", "visit_age", "ICV"]
    for count in range(5):
        for confounds in combinations(names, count):
            t = Transform(confounds=confounds)
            model = Design.fit(df, t)
            X = model.matrix(df)
            beta = np.arange(X.shape[1]*3, dtype=float).reshape(X.shape[1], 3)*1e-4
            beta[0] = [5, 6, 7]
            y = X @ beta
            fitted, _ = fit_coefficients(X, y)
            np.testing.assert_allclose(fitted, beta, atol=1e-9)
            adjusted = apply_adjustment(y, X, fitted, model.reference)
            np.testing.assert_allclose(adjusted.mean(axis=0), y.mean(axis=0), atol=1e-9)
            np.testing.assert_allclose(adjusted, np.broadcast_to(beta[0], adjusted.shape), atol=1e-9)


def test_training_only_design_does_not_change_when_test_changes():
    df = data()
    rng = np.random.default_rng(7)
    y = rng.normal(size=(len(df), 4))
    t = Transform(missing="median")
    first = adjusted_features(df.iloc[:35], df.iloc[35:], y[:35], y[35:], t)
    changed = df.iloc[35:].copy(); changed["ICV"] = 1e9
    second = adjusted_features(df.iloc[:35], changed, y[:35], y[35:]*10, t)
    assert first[2].to_dict() == second[2].to_dict()
    np.testing.assert_array_equal(first[3], second[3])
    np.testing.assert_array_equal(first[0], second[0])


def test_statistics_sign_fdr_degenerate_and_welch():
    a, b = np.array([1., 2, 3, 5]), np.array([5., 6, 9, 12])
    comparison = compare(a, b)
    assert comparison["difference"] > 0 and comparison["t"] > 0 and comparison["hedges_g"] > 0
    np.testing.assert_allclose(comparison["p"], stats.ttest_ind(b, a, equal_var=False).pvalue)
    q, rejected = fdr([.001, .02, .8, np.nan])
    np.testing.assert_allclose(q[:3], [.003, .03, .8])
    assert np.isnan(q[-1]) and not rejected[-1]
    assert compare([1, 1], [2, 2])["status"] == "zero_variance"
    samples = [a, b, np.array([2., 2.5, 6, 9, 17])]
    np.testing.assert_allclose(omnibus(samples, False)["p"], stats.f_oneway(*samples).pvalue)
    # Independent Welch ANOVA reference, implemented by SciPy.
    np.testing.assert_allclose(omnibus(samples, True)["p"], stats.f_oneway(*samples, equal_var=False).pvalue)

