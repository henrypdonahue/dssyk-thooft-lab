#!/usr/bin/env python3
r"""
ED reference bench for the large-q (lambda -> 0) anchor (road_map rung 17).
===========================================================================

The lambda -> 0 corner is where the DSSYK <-> 't Hooft match is supposed to
become exact, and closed-form large-q SYK correlators are the analytic anchor
there.  Before trusting any transcribed large-q formula, we need
model-generated ground truth: this module measures, by exact diagonalization
with the fast Pauli-string builder,

    G(t)   = < psi_1(t) psi_1(0) >          (T = inf two-point function)
    C_2(t) = < Obar(t) Obar(0) >,  O = -(1/N) sum_a psidot_a psidot_a
                                        (the M_2-type Majorana singlet;
                                         traceless part)

disorder-averaged with errors, along FIXED p with growing N (lambda = p^2/N
falling -- the direction in which large-q formulas, whose corrections are
O(1/q) AND O(lambda), must improve), stored in largeq_bench.json.  The
comparison against the transcribed closed forms lives in
duality/largeq_anchor.py.

(A trap this module hit and now documents: the "obvious" singlet
(i/N) sum psi_a psidot_a is exactly conserved -- the Majorana analog of
M_1 = -(p/2)iH reads  sum_a psi_a [H, psi_a] = -2p H  in our psi^2 = 1
normalization, asserted in the tests -- so it is 100% static and useless as a
dynamical probe.  The first genuinely dynamical member of the Majorana tower
is sum psidot psidot, by the same conservation identity that fixes M_2.)

Conventions: Var(J) = p!/N^{p-1} (papers' Eq. 3.8, script-J = 1), Majoranas
psi_a^2 = 1.  With this normalization <H^2>/dim = C(N,p) p!/N^{p-1} ~ N/p...
times p!-combinatorics -- the bench stores <H^2>/dim per instance so any
formula's energy-scale convention can be matched exactly rather than assumed.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pauli_strings import majorana_hamiltonian_fast
from syk import majorana_operators


def two_point_ed(N: int, p: int, ts, n_inst: int = 8, seed: int = 0):
    """Disorder-averaged G(t) = <psi_1(t) psi_1> and the singlet bilinear
    correlator at T = inf, by full ED.  Returns dict of curves with standard
    errors, plus the measured energy scale <H^2>/dim per instance."""
    rng = np.random.default_rng(seed)
    psis = [P.toarray() for P in majorana_operators(N)]
    dim = 1 << (N // 2)
    ts = np.asarray(ts, float)
    G_rows, C_rows, m2_rows, sf_rows = [], [], [], []
    for _ in range(n_inst):
        H = majorana_hamiltonian_fast(N, p, rng)
        E, V = np.linalg.eigh(H)
        m2_rows.append(float(np.mean(E ** 2)))
        omega = E[:, None] - E[None, :]

        psi1 = V.conj().T @ psis[0] @ V
        # psi1 Hermitian => psi1_ab psi1_ba = |psi1_ab|^2:
        # G(t) = sum_ab e^{i omega_ab t} |psi1_ab|^2 / dim,  G(0) = 1 exactly
        W1 = np.abs(psi1) ** 2 / dim
        G_rows.append([complex(np.sum(np.exp(1j * omega * t) * W1))
                       for t in ts])

        # M_2-type singlet O = -(1/N) sum_a psidot_a psidot_a, traceless part,
        # constructed directly in the eigenbasis (psidot = i[H, psi])
        O = np.zeros((dim, dim), dtype=complex)
        for P in psis:
            Pd = V.conj().T @ P @ V
            Pdot = 1j * (E[:, None] * Pd - Pd * E[None, :])
            O -= Pdot @ Pdot
        O /= N
        O -= np.trace(O) / dim * np.eye(dim)
        WO = np.abs(O) ** 2 / dim
        C_rows.append([float(np.sum(np.cos(omega * t) * WO)) for t in ts])
        sf_rows.append(float(WO[np.abs(omega) < 1e-9].sum() / WO.sum()))
    G = np.array(G_rows)
    C = np.array(C_rows)
    sem = 1.0 / np.sqrt(max(n_inst - 1, 1))
    return dict(
        N=N, p=p, n_inst=n_inst, ts=ts.tolist(),
        G_re=np.mean(G.real, 0).tolist(),
        G_re_err=(np.std(G.real, 0, ddof=1) * sem).tolist(),
        G_im=np.mean(G.imag, 0).tolist(),
        C2=np.mean(C, 0).tolist(),
        C2_err=(np.std(C, 0, ddof=1) * sem).tolist(),
        # static (omega = 0) fraction of the M_2-type singlet correlator: the
        # large-q "Hamiltonian exchange only" structure predicts the dynamical
        # remainder is subleading as lambda -> 0
        static_fraction=float(np.mean(sf_rows)),
        static_fraction_err=float(np.std(sf_rows, ddof=1) * sem),
        m2_mean=float(np.mean(m2_rows)), m2_err=float(np.std(m2_rows, ddof=1) * sem),
    )


if __name__ == "__main__":
    ts = np.linspace(0.0, 6.0, 61)
    rows = []
    print("ED bench for the large-q anchor (T = inf, Var J = p!/N^{p-1}).")
    print("Fixed p, growing N: lambda = p^2/N falls toward the regime where")
    print("the large-q closed forms (corrections O(1/q) AND O(lambda)) apply.\n")
    for (N, p, n_inst) in [(16, 4, 12), (20, 4, 8), (24, 4, 4),
                           (16, 6, 12), (18, 6, 8), (16, 8, 12)]:
        r = two_point_ed(N, p, ts, n_inst=n_inst, seed=40 + N + p)
        rows.append(r)
        print(f"  N={N} p={p} (lam={p*p/N:.2f}): <H^2>/dim = "
              f"{r['m2_mean']:.3f} ± {r['m2_err']:.3f},  G(0) = {r['G_re'][0]:.4f},  "
              f"static frac of C2 = {r['static_fraction']:.3f} ± "
              f"{r['static_fraction_err']:.3f}")
    out = Path(__file__).resolve().parent / "largeq_bench.json"
    with open(out, "w") as fh:
        json.dump(rows, fh)
    print(f"\nwrote {out}")
    print("Ground truth for duality/largeq_anchor.py: any transcribed large-q")
    print("closed form must hit these curves to O(1/p) with no free parameters.")
