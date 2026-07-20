#!/usr/bin/env python3
r"""
The DSSYK_infinity <-> 't Hooft comparison dictionary.
======================================================

Every convention needed to compare a DSSYK number against the 't Hooft
yardstick, encoded ONCE, with tests -- built (deliberately) before the DSSYK
side of the comparison exists, so that when a number finally appears there is
exactly one place where an index offset, factor of 2, or normalization can be
wrong, and it is unit-tested.

The four load-bearing corrections (top-level README), as code:

  1. Index offset:      n_paper = n_FLZ + 2   (paper adds the non-dynamical
                        photon n=0 and graviton n=1 below the dynamical tower).
  2. Where the match lives:  fixed p, lambda -> 0 ('t Hooft scaling
                        alpha_bar = p^2 fixed), NOT the lambda-fixed
                        double-scaling line.
  3. Normalization:     M_0 = Q exactly;  M_1 = -(p/2) i H  (not the paper's
                        M_1 = H) -- the constant M1_COEFFICIENT below.
  4. Conventions:       the PAPER's lambda is p^2/N; the CHORD literature
                        (Berkooz et al., the machinery any DSSYK solve will
                        use) defines lambda_chord = 2p^2/N and q =
                        exp(-lambda_chord).  A silent factor of 2 here would
                        wreck any quantitative match -- hence the explicit
                        converters.

Mass units.  thooft-target reports the dimensionless eigenvalue
2*lambda_n = M_n^2/(pi g^2) at the duality point alpha = pi m_q^2/g^2 = 0.
Physical masses need a value of g^2 (or of the string tension
tau = g^2 N = p^2 m_q^2), which the duality fixes only through the dictionary;
RATIOS of masses-squared are parameter-free, so the comparison tables below
are built from ratios and from the interleave splittings -- the sharpest
falsifiable fingerprints available today.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

# ---------------------------------------------------------------------------
# couplings and double-scaling conventions
# ---------------------------------------------------------------------------

# correction #3: M_1 = M1_COEFFICIENT * i * H  (the paper says M_1 = H)
M1_COEFFICIENT = lambda p: -p / 2.0


def lambda_paper(N: float, p: float) -> float:
    """The paper's double-scaling parameter lambda = p^2/N (their Sec. 3.1)."""
    return p * p / N


def lambda_chord(N: float, p: float) -> float:
    """The chord-literature convention lambda_chord = 2 p^2/N (Berkooz et al.)."""
    return 2.0 * p * p / N


def chord_from_paper(lam_paper: float) -> float:
    """lambda_chord = 2 * lambda_paper."""
    return 2.0 * lam_paper


def paper_from_chord(lam_chord: float) -> float:
    return 0.5 * lam_chord


def q_parameter(lam_chord: float) -> float:
    """q = exp(-lambda_chord); q -> 1 is the paper's lambda -> 0 corner where
    an exact 't Hooft match is expected."""
    return math.exp(-lam_chord)


def alpha_bar(p: float) -> float:
    """Dimensionless 't Hooft coupling alpha_bar = gbar^2 N = p^2 (paper Sec. 4.1)."""
    return p * p


def gbar_squared(N: float, p: float) -> float:
    """Dimensionless gauge/string coupling gbar^2 = lambda = p^2/N."""
    return lambda_paper(N, p)


def string_tension(m_q_squared: float, p: float) -> float:
    """tau = g^2 N = alpha_bar * m_q^2 = p^2 m_q^2 (paper Sec. 4.2).  Fixed,
    finite tau with p^2 = lambda N forces m_q^2 -> 0 in the double-scaled
    limit: the duality point is the RENORMALIZED-massless 't Hooft model."""
    return alpha_bar(p) * m_q_squared


# ---------------------------------------------------------------------------
# index and quantum-number maps (correction #1)
# ---------------------------------------------------------------------------
def n_paper_from_flz(n_flz: int) -> int:
    """The paper's tower index: its dynamical mesons start at n_paper = 2;
    FLZ's ground state (n_FLZ = 0) is the paper's n = 2 'eta'."""
    if n_flz < 0:
        raise ValueError("FLZ index must be >= 0")
    return n_flz + 2


def n_flz_from_paper(n_paper: int) -> int:
    if n_paper < 2:
        raise ValueError(
            f"n_paper = {n_paper} is non-dynamical (photon n=0 / graviton "
            "n=1); it has no FLZ level")
    return n_paper - 2


def cp_paper(n_paper: int) -> int:
    """CP = (-1)^(n+1) in the paper's indexing (Sec. 4.3).  For n = 0, 1 this
    is the photon (-1) / graviton (+1); those assignments are EXACT operator
    statements on the SYK side (see syk-self-averaging/dirac.py).  For the
    dynamical tower n >= 2 the raw M_n operators are CP mixtures; the
    assignment refers to the CP-resolved correlator channel."""
    return (-1) ** (n_paper + 1)


def flz_parity(n_flz: int) -> str:
    """FLZ wavefunction parity under x -> 1-x: even n symmetric ("sym"), odd n
    antisymmetric ("anti") -- the labels used in reference_spectrum.json.
    Consistency (asserted in tests): symmetric FLZ levels are the paper's
    CP = -1 states, antisymmetric are CP = +1, and the two interleave with
    alternating CP."""
    return "sym" if n_flz % 2 == 0 else "anti"


# ---------------------------------------------------------------------------
# the yardstick: reference spectrum and its parameter-free fingerprints
# ---------------------------------------------------------------------------
_REFERENCE_PATH = Path(__file__).resolve().parent.parent / "thooft-target" / "reference_spectrum.json"


def load_reference(path=None) -> list:
    """The certified massless 't Hooft levels, decorated with both index
    conventions.  Each row: n_flz, n_paper, parity, CP, two_lambda (float;
    = M^2/(pi g^2)), stable_digits."""
    p = Path(path) if path is not None else _REFERENCE_PATH
    with open(p) as fh:
        raw = json.load(fh)
    rows = []
    for r in raw:
        n_flz = r["n"]
        rows.append(dict(
            n_flz=n_flz,
            n_paper=n_paper_from_flz(n_flz),
            parity=r["parity"],
            CP=r["CP"],
            two_lambda=float(r["value"]),
            stable_digits=r["stable_digits"],
        ))
    return rows


def mass_squared_ratios(rows=None) -> list:
    """Parameter-free fingerprint #1: M_n^2 / M_2^2 in PAPER indexing
    (= two_lambda_n / two_lambda_0 in FLZ indexing).  Any DSSYK singlet tower
    must reproduce these numbers as lambda -> 0, whatever the overall units."""
    if rows is None:
        rows = load_reference()
    ground = rows[0]["two_lambda"]
    return [dict(n_paper=r["n_paper"], CP=r["CP"],
                 ratio=r["two_lambda"] / ground) for r in rows]


def interleave_splittings(rows=None) -> list:
    """Parameter-free fingerprint #2: adjacent-level differences
    two_lambda_{n+1} - two_lambda_n (FLZ indexing).  The 't Hooft tower is a
    single interleaved trajectory with alternating CP and near-unit spacing
    (2lam ~ n + 3/4): gaps 1.0167, 0.9944, 1.0029, ... -> 1.  A DSSYK tower
    showing pairwise-degenerate CP doublets instead of this alternation would
    falsify the identification.  (This also settles what the paper's 'two
    degenerate Regge trajectories... on top of one another' must mean: the
    even- and odd-CP trajectories have the same slope and interleave --
    adjacent levels are NOT degenerate.)"""
    if rows is None:
        rows = load_reference()
    return [dict(between=(rows[i]["n_flz"], rows[i + 1]["n_flz"]),
                 CP_pair=(rows[i]["CP"], rows[i + 1]["CP"]),
                 gap=rows[i + 1]["two_lambda"] - rows[i]["two_lambda"])
            for i in range(len(rows) - 1)]


if __name__ == "__main__":
    rows = load_reference()
    print("The 't Hooft yardstick in both index conventions "
          "(two_lambda = M^2/(pi g^2)):\n")
    print(f"{'n_FLZ':>6} {'n_paper':>8} {'CP':>4} {'parity':>7} "
          f"{'two_lambda':>20} {'M^2/M_eta^2':>13}")
    ground = rows[0]["two_lambda"]
    for r in rows[:12]:
        print(f"{r['n_flz']:>6} {r['n_paper']:>8} {r['CP']:>+4d} "
              f"{r['parity']:>7} {r['two_lambda']:>20.15f} "
              f"{r['two_lambda']/ground:>13.9f}")
    print("\nInterleave splittings (fingerprint #2):")
    for s in interleave_splittings(rows)[:8]:
        print(f"  n_FLZ {s['between'][0]}->{s['between'][1]}  "
              f"CP {s['CP_pair'][0]:+d}->{s['CP_pair'][1]:+d}   "
              f"gap = {s['gap']:.10f}")
    print("\nConvention conversions (N=640, p=26):")
    lp = lambda_paper(640, 26)
    lc = lambda_chord(640, 26)
    print(f"  lambda_paper = {lp:.6f}   lambda_chord = {lc:.6f}   "
          f"q = {q_parameter(lc):.6f}")
    print(f"  alpha_bar = p^2 = {alpha_bar(26):.0f}   gbar^2 = {gbar_squared(640, 26):.6f}")
