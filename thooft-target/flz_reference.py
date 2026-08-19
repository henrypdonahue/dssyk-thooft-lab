"""FLZ reference constants (arXiv:0905.2280), one provenance-checked copy.

FLZ Tables 1 & 2, last column 2*lambda_n^(num), 14 significant digits,
transcribed from the published paper (digit-for-digit audit against the
fetched PDF, 2026-07-17 review).  Shared by validate.py, test_thooft.py
and test_jacobi.py so the literature anchor has exactly one place to
audit; the SOLVERS being validated share no code with this file.
"""

from scipy.special import polygamma

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
    from the leading asymptotic 2 lambda_n -> n + 3/4:
    tail = polygamma(1, (n0 + 3/4)/2)."""
    total, n = 0.0, offset
    for two_lam in evs:
        if n > n_cut:
            break
        total += 4.0 / two_lam ** 2
        n += 2
    return total + polygamma(1, (n + 0.75) / 2.0)
