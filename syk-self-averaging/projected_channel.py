#!/usr/bin/env python3
"""Finite-N dynamical spectral weight after removing the full commutant of H.

Subtracting Tr(A)/dim removes only the identity. In an energy eigenbasis,
all matrix elements within degenerate energy blocks are stationary. Removing
those blocks is the Hilbert-Schmidt projection off ker(ad_H), independent
of the chosen basis within each degenerate block. This is an ED diagnostic,
not a large-N proof or a mass extraction. Near-degeneracies are classified
using a recorded absolute energy tolerance; inspect tolerance sensitivity.

For A_2 the large-N statement now exists (duality/singlet_tower.py). At
N = infinity the stationary part is exactly the projection onto H^2 - 1
and the dynamical part is a matter chord of dimension 2, with width
2(1 - q^2). The split computed here agrees with it at ED sizes: the
disorder-averaged stationary weight matches the exact hypergeometric H^2
projection, and at p = 6, N = 16-20 the dynamical width is within 2 % of
the N = infinity value (singlet_tower.ed_commutant_split, tested).
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np


def spectral_split(H, A, energy_tol=None):
    """Return identity-subtracted and dynamical moments in H's energy units.

    Symmetrized positive weights are (|A_ij|^2+|A_ji|^2)/(2 dim).
    Undefined dynamical widths are None, never a physical zero width.
    """
    H, A = np.asarray(H), np.asarray(A)
    if H.ndim != 2 or H.shape[0] != H.shape[1] or A.shape != H.shape:
        raise ValueError('H and A must be square matrices of the same shape')
    if not np.all(np.isfinite(H)) or not np.all(np.isfinite(A)):
        raise ValueError('matrices must be finite')
    if not np.allclose(H, H.conj().T, rtol=1e-12, atol=1e-12):
        raise ValueError('H must be Hermitian')
    E, V = np.linalg.eigh(H)
    dim = len(E)
    scale = max(float(np.max(np.abs(E))), np.finfo(float).tiny)
    tol = float(energy_tol) if energy_tol is not None else 64 * dim * np.finfo(float).eps * scale
    if not np.isfinite(tol) or tol < 0:
        raise ValueError('energy_tol must be finite and nonnegative')
    B = V.conj().T @ (A - np.trace(A) / dim * np.eye(dim)) @ V
    W = (np.abs(B) ** 2 + np.abs(B.T) ** 2) / (2 * dim)
    omega = E[:, None] - E[None, :]
    stationary = np.abs(omega) <= tol
    raw0 = float(W.sum())
    zero0 = float(W[stationary].sum())
    dynW = np.where(stationary, 0.0, W)
    dyn0 = float(dynW.sum())
    mu2 = float(np.sum(dynW * omega ** 2))
    mu4 = float(np.sum(dynW * omega ** 4))
    # Eigensolver rotation can leave roundoff off the exact commutant.
    resolved = dyn0 > 1e-24 * max(float(np.sum(np.abs(A) ** 2) / dim), np.finfo(float).tiny)
    return dict(energy_tol=tol, mu0_traceless=raw0, mu0_stationary=zero0,
                mu0_dynamic=dyn0, mu2_dynamic=mu2, mu4_dynamic=mu4,
                stationary_fraction=zero0 / raw0 if raw0 > 0 else None,
                dynamic_resolved=bool(resolved),
                w2_dynamic=mu2 / dyn0 if resolved else None,
                kurt_dynamic=mu4 * dyn0 / mu2 ** 2 if resolved and mu2 > 0 else None)


def run():
    from conserved_points import a_n_tower, dense_majoranas
    from pauli_strings import majorana_hamiltonian_fast
    rows = []
    for N, p in ((8, 4), (10, 4), (12, 4)):
        for seed in (17, 29, 43):
            H = majorana_hamiltonian_fast(N, p, np.random.default_rng(seed))
            tower = a_n_tower(H, dense_majoranas(N), 3)
            for n in (1, 2, 3):
                row = spectral_split(H, tower[n])
                lo = spectral_split(H, tower[n], row['energy_tol'] / 10)
                hi = spectral_split(H, tower[n], row['energy_tol'] * 10)
                row.update(N=N, p=p, n=n, seed=seed,
                           tolerance_drift=max(abs(lo['mu0_dynamic']-row['mu0_dynamic']),
                                               abs(hi['mu0_dynamic']-row['mu0_dynamic'])))
                rows.append(row)
    payload = dict(provenance='projected_channel.py; per-instance normalized trace; full degenerate-block projection; 3 seeds per size; diagnostic, no scaling fit', rows=rows)
    out = Path(__file__).with_suffix('.json')
    out.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    print(f'wrote {out.name}: {len(rows)} rows')
    return payload


if __name__ == '__main__':
    run()
