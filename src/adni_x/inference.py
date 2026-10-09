"""Bounded-memory symmetry checks and participant-based coefficient inference."""
import numpy as np
from scipy import stats

from .constants import ATOMS
from .numerics import connection_mask


def validate_symmetric_atoms(data, rows, exclude_diagonal=True, block_features=4096):
    """Check each participant, including mirrored missingness; never symmetrise inputs."""
    p = data.shape[1]
    ii, jj = np.where(connection_mask(p, exclude_diagonal, upper=True))
    flat = data.reshape(len(data), -1)
    report = {"rtol": 1e-7, "atol": 1e-10, "atoms": {}}
    for atom in (0, 15):
        maximum = 0.0
        for start in range(0, len(ii), block_features):
            i, j = ii[start:start + block_features], jj[start:start + block_features]
            a = flat[np.ix_(rows, (i*p+j)*16+atom)]
            b = flat[np.ix_(rows, (j*p+i)*16+atom)]
            valid = np.isfinite(a)
            if not np.array_equal(valid, np.isfinite(b)):
                raise ValueError(f"{ATOMS[atom]} symmetry failed: mirrored finite-value masks differ")
            difference = np.abs(a[valid]-b[valid])
            if len(difference):
                maximum = max(maximum, float(difference.max()))
                tolerance = report["atol"] + report["rtol"]*np.maximum(abs(a[valid]), abs(b[valid]))
                if np.any(difference > tolerance):
                    raise ValueError(f"{ATOMS[atom]} symmetry failed: discrepancy {maximum:g} exceeds tolerance")
        report["atoms"][ATOMS[atom]] = {"max_abs_discrepancy": maximum,
                                       "unique_pairs": len(ii), "participants": len(rows)}
    return report


def hc3_coefficients(X, Y, beta):
    """HC3 diagonal sandwich covariance, with explicit approximate t inference.

    Reuse fitted coefficients. Group features by valid participant rows and do
    not allocate a covariance matrix for every feature. The caller supplies a
    bounded feature block, before confound adjustment.
    """
    X, Y, beta = np.asarray(X, float), np.asarray(Y, float), np.asarray(beta, float)
    if Y.ndim != 2 or beta.shape != (X.shape[1], Y.shape[1]) or len(X) != len(Y):
        raise ValueError("Incompatible coefficient-inference shapes")
    if not np.isfinite(X).all():
        raise ValueError("Coefficient design must be finite")
    output = {name: np.full(beta.shape, np.nan) for name in
              ("standard_error", "t", "p", "ci_lower", "ci_upper")}
    output["status"] = np.full(beta.shape, "insufficient_observations", dtype="U32")
    output["nobs"] = np.isfinite(Y).sum(axis=0)
    output["df_resid"] = np.full(Y.shape[1], np.nan)
    patterns, inverse = np.unique(np.isfinite(Y).T, axis=0, return_inverse=True)
    for number, valid in enumerate(patterns):
        columns = np.flatnonzero(inverse == number)
        design = X[valid]
        rank = np.linalg.matrix_rank(design) if len(design) else 0
        df = int(valid.sum())-rank
        output["df_resid"][columns] = df
        if df <= 0:
            continue
        if rank < X.shape[1]:
            output["status"][:, columns] = "rank_deficient_design"
            continue
        pinv = np.linalg.pinv(design)
        leverage = np.einsum("ij,ji->i", design, pinv)
        if np.any(1-leverage <= 1e-10):
            output["status"][:, columns] = "unit_leverage"
            continue
        residual = Y[np.ix_(valid, columns)]-design@beta[:, columns]
        weights = (residual/(1-leverage[:, None]))**2
        se = np.sqrt((pinv**2)@weights)
        estimable = np.isfinite(se) & (se > 0) & np.isfinite(beta[:, columns])
        status = np.where(np.isfinite(se), "zero_standard_error", "nonfinite_standard_error")
        output["status"][:, columns] = np.where(estimable, "ok", status)
        output["standard_error"][:, columns] = np.where(estimable, se, np.nan)
        t = np.full(se.shape, np.nan)
        np.divide(beta[:, columns], se, out=t, where=estimable)
        output["t"][:, columns] = t
        output["p"][:, columns] = 2*stats.t.sf(abs(t), df)
        width = stats.t.ppf(.975, df)*se
        output["ci_lower"][:, columns] = np.where(estimable, beta[:, columns]-width, np.nan)
        output["ci_upper"][:, columns] = np.where(estimable, beta[:, columns]+width, np.nan)
    return output
