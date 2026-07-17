#!/usr/bin/env python3
"""Validation of the massless 't Hooft spectrum against Fateev-Lukyanov-
Zamolodchikov (arXiv:0905.2280): eigenvalue tables, exact sum rules, and the
large-n asymptotic expansion."""

import sys

import numpy as np
from scipy.special import zeta, polygamma

from thooft_spectrum import solve_sector, spectrum
import sine_solver

PI = np.pi

# --- FLZ Tables 1 & 2, last column 2*lambda_n^(num) (14 sig. digits) ----------
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
FLZ_PADE_LAMBDA0 = 0.737061746292690  # FLZ Eq. (4.37), Pade estimate


def compare_tables(num_basis=300):
    ev_sym, _, _ = solve_sector(num_basis, parity=0)
    ev_anti, _, _ = solve_sector(num_basis, parity=1)

    print(f"Basis size per sector K = {num_basis}\n")
    print("Comparison with FLZ Tables 1 & 2 (2*lambda_n):")
    print(f"{'n':>3} {'this work':>20} {'FLZ':>20} {'|diff|':>10}")
    print("-" * 58)
    worst = 0.0
    for n in sorted({**FLZ_EVEN, **FLZ_ODD}):
        mine = ev_sym[n // 2] if n % 2 == 0 else ev_anti[n // 2]
        ref = FLZ_EVEN.get(n, FLZ_ODD.get(n))
        d = abs(mine - ref)
        worst = max(worst, d)
        print(f"{n:>3} {mine:>20.14f} {ref:>20.14f} {d:>10.1e}")
    print("-" * 58)
    print(f"worst absolute deviation over n<=29: {worst:.2e}")
    print(f"\nGround state 2*lambda_0 = {ev_sym[0]:.15f}")
    print(f"FLZ Pade  2*lambda_0 = {FLZ_PADE_LAMBDA0:.15f}  "
          f"(|diff| = {abs(ev_sym[0]-FLZ_PADE_LAMBDA0):.1e})")
    return worst


def check_sum_rules(num_basis=400, n_cut=120):
    """Exact FLZ sum rules (Eq. 1.8):
        G2+ = sum_m 1/lambda_{2m}^2   = 7 zeta(3)
        G2- = sum_m 1/lambda_{2m+1}^2 = 2

    lambda_n = (2 lambda_n)/2, so 1/lambda_n^2 = 4/(2 lambda_n)^2.  The sum is
    dominated by the lowest few levels; we use accurately computed eigenvalues
    for n <= n_cut and close the (small) tail analytically with the leading
    asymptotic 2 lambda_n -> n + 3/4, giving
        tail = sum_{n=n0,n0+2,...} 4/(n+3/4)^2 = polygamma(1, (n0+3/4)/2)."""
    ev_sym, _, _ = solve_sector(num_basis, parity=0)
    ev_anti, _, _ = solve_sector(num_basis, parity=1)

    def sum_inv_sq(evs, offset):
        total = 0.0
        n = offset
        for two_lam in evs:
            if n > n_cut:
                break
            total += 4.0 / two_lam ** 2  # 1/lambda_n^2
            n += 2
        # analytic tail from the next level n (leading asymptotic n+3/4)
        total += polygamma(1, (n + 0.75) / 2.0)
        return total

    G2p = sum_inv_sq(ev_sym, 0)
    G2m = sum_inv_sq(ev_anti, 1)
    print("\nExact sum rules (FLZ Eq. 1.8):")
    print(f"  G2+ = sum 1/lambda_2m^2   = {G2p:.12f}   "
          f"target 7*zeta(3) = {7*zeta(3):.12f}   |diff|={abs(G2p-7*zeta(3)):.1e}")
    print(f"  G2- = sum 1/lambda_2m+1^2 = {G2m:.12f}   "
          f"target 2         = {2.0:.12f}   |diff|={abs(G2m-2.0):.1e}")
    return max(abs(G2p - 7 * zeta(3)), abs(G2m - 2.0))


def crosscheck_independent():
    """Compare the Chebyshev spectrum against the independent sine-basis solver
    at alpha = 0 and two nonzero masses.  This is the check on the quark-mass
    machinery (invisible at alpha = 0)."""
    print("\nIndependent cross-check (sine-basis solver, sine_solver.py):")
    print("  (algebraic convergence -> ~1e-3..1e-4 agreement is expected & confirms")
    print("   the potential term; NOT a precision statement about the main solver)")
    worst = 0.0
    for alpha in (0.0, 0.5, 2.0):
        cheb = [L["two_lambda"]
                for L in spectrum(num_basis=300, num_levels=6, alpha=alpha)]
        sine = sine_solver.solve(num_basis=90, n_quad=700, alpha=alpha, num_levels=6)
        d = max(abs(cheb[n] - sine[n]) for n in range(6))
        worst = max(worst, d if alpha != 0.0 else 0.0)
        tag = "(duality point)" if alpha == 0.0 else ""
        print(f"  alpha={alpha:>4}: max|cheb - sine| over n<6 = {d:.2e}  {tag}")
    return worst


def convergence_study():
    print("\nVariational convergence of ground state 2*lambda_0 (from above):")
    print(f"{'K':>5} {'2*lambda_0':>20} {'change':>12}")
    prev = None
    for K in [2, 4, 8, 16, 32, 64, 128, 256]:
        ev, _, _ = solve_sector(K, parity=0)
        v = ev[0]
        ch = "" if prev is None else f"{v-prev:+.2e}"
        print(f"{K:>5} {v:>20.15f} {ch:>12}")
        prev = v


def check_asymptotics():
    print("\nApproach to leading 't Hooft/WKB asymptotic 2*lambda_n -> n + 3/4:")
    levels = spectrum(num_basis=300, num_levels=30)
    print(f"{'n':>3} {'2*lambda_n':>18} {'(2lam_n) - (n+3/4)':>20}")
    for L in levels[::4]:
        n = L['n']
        v = L['two_lambda']
        print(f"{n:>3} {v:>18.12f} {v - (n + 0.75):>20.3e}")


if __name__ == "__main__":
    worst_tables = compare_tables()
    worst_sumrule = check_sum_rules()
    convergence_study()
    check_asymptotics()
    worst_cross = crosscheck_independent()

    # Pass/fail gate (exit code != 0 breaks CI if the physics regresses).
    checks = {
        "FLZ tables (< 1e-10)": worst_tables < 1e-10,
        "sum rules (< 1e-8)": worst_sumrule < 1e-8,
        "sine cross-check alpha!=0 (< 5e-3)": worst_cross < 5e-3,
    }
    print("\n" + "=" * 40)
    ok = True
    for name, passed in checks.items():
        print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
        ok = ok and passed
    print("=" * 40)
    sys.exit(0 if ok else 1)
