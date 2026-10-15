"""FLZ reference constants (arXiv:0905.2280), one provenance-checked copy.

FLZ Tables 1 and 2, last column 2*lambda_n^(num), at 14 significant
digits.  Transcribed from the published paper, with a digit-for-digit
audit against the fetched PDF in the 2026-07-17 review.

validate.py, test_thooft.py and test_jacobi.py all share this copy.  So
the literature anchor has exactly one place to audit.

The SOLVERS being validated share no code with this file.
"""

from math import factorial

import numpy as np
from scipy.special import digamma, polygamma

FLZ_EVEN = {  # symmetric sector, FLZ even n
    0: 0.73706174629269, 2: 2.7481609123706, 4: 4.7492953810375,
    6: 6.7496294196488, 8: 8.7497715807892, 10: 10.749845089160,
    12: 12.749888008416, 14: 14.749915244446, 16: 16.749933611057,
    18: 18.749946584034, 20: 20.749956088173, 22: 22.749963259761,
    24: 24.749968804883, 26: 26.749973181145, 28: 28.749976695731,
}
FLZ_ODD = {  # antisymmetric sector, FLZ odd n
    1: 1.7537313369175, 3: 3.7510575817054, 5: 5.7504926236487,
    7: 7.7502843971925, 9: 9.7501851352539, 11: 11.750130142515,
    13: 13.750096503972, 15: 15.750074428438, 17: 17.750059159035,
    19: 19.750048157169, 21: 21.750039967130, 23: 23.750033705317,
    25: 25.750028810060, 27: 27.750024910394, 29: 29.750021753287,
}
FLZ = {**FLZ_EVEN, **FLZ_ODD}                 # merged, keyed by FLZ n
FLZ_PADE_LAMBDA0 = 0.737061746292690          # FLZ Eq. (4.37), Pade estimate


def sum_inv_sq_with_tail(evs, offset, n_cut):
    """G2 partial sum sum_n 4/(2 lambda_n)^2 over one CP sector (levels at
    n = offset, offset+2, ... up to n_cut), closed with the analytic tail
    from the leading asymptotic 2 lambda_n -> n + 3/4.  The tail is
    polygamma(1, (n0 + 3/4)/2), with n0 the first level index past the
    summed head."""
    total, n = 0.0, offset
    for two_lam in evs:
        if n > n_cut:
            break
        total += 4.0 / two_lam ** 2
        n += 2
    return total + polygamma(1, (n + 0.75) / 2.0)


def flz_sum_rule(evs, offset, s, n_cut=120, tail="leading",
                 m_tail=2_000_000):
    """FLZ G^(s)+- = sum_n lambda_n^-s over one CP sector (s >= 2), or the
    regularized s = 1 sum G^(1)+- = sum_m [1/lambda_n - 1/(m+1)]
    (FLZ Eq. 4.18; n = 2m + offset), from the computed 2*lambda_n up to
    n_cut and an analytic tail beyond.

    tail = "leading": the polygamma/digamma tail from 2 lambda_n -> n + 3/4
    (the original; limits G2 to ~1e-10 and G1 to ~4e-9).
    tail = "eq16": the transcribed FLZ Eq. (1.6) correction delta_n is
    summed term by term for m_tail levels before the leading tail takes
    over; that brings every s to the ~1e-12 level of FLZ's Table 3."""
    total, n, m = 0.0, offset, 0
    for two_lam in evs:
        if n > n_cut:
            break
        total += (2.0 / two_lam - 1.0 / (m + 1)) if s == 1 else 2.0 ** s / two_lam ** s
        n += 2
        m += 1
    if tail == "eq16":
        from flz_dconstants import flz_delta_asymptotic
        ns = np.arange(n, n + 2 * m_tail, 2, dtype=float)
        two_lam = ns + 0.75 + flz_delta_asymptotic(ns)
        if s == 1:
            ms = np.arange(m, m + m_tail, dtype=float)
            total += float(np.sum(2.0 / two_lam - 1.0 / (ms + 1)))
        else:
            total += float(np.sum(2.0 ** s / two_lam ** s))
        n += 2 * m_tail
        m += m_tail
    elif tail != "leading":
        raise ValueError(f"unknown tail {tail!r}")
    if s == 1:
        a = 0.375 if offset == 0 else 0.875       # lambda_n -> m + a
        return total + float(digamma(m + 1) - digamma(m + a))
    return total + float((-1) ** s * polygamma(s - 1, (n + 0.75) / 2.0)
                         / factorial(s - 1))
