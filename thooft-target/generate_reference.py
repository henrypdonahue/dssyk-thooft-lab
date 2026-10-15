"""Generate the certified high-precision reference tables.

Default: the original route.  A full mpmath eigensolve at two basis
sizes (K = 160, 220), reporting the agreeing digits per level, written
to reference_spectrum.json (24 levels, 12-15 digits; the table
duality/ consumes).

--polish: the cheaper and sharper route.  Double-precision Ritz vectors
from thooft_spectrum.solve_sector, each polished by an exact-rational
Rayleigh quotient in mpmath (thooft_highprec.rayleigh_polish, second
order in the vector error), on a three-rung basis ladder whose measured
K^-6 truncation law is extrapolated away.  Written to
reference_spectrum_polished.json: 40 levels, each with its measured
rate, its extrapolation correction and that correction's magnitude as
the error bar.  About seven minutes (395 s recorded in the artifact).
"""
from mpmath import mp, nstr, fabs
from thooft_highprec import solve_sector_mp
import json, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    saved_dps = mp.dps
    mp.dps = 50
    try:
        _main_body()
    finally:
        mp.dps = saved_dps


def _main_body():
    NLEV = 24
    results = {}
    for K in (160, 220):
        t = time.time()
        sym = solve_sector_mp(K, 0)
        anti = solve_sector_mp(K, 1)
        results[K] = (sym, anti)
        print(f"K={K} done ({time.time()-t:.0f}s)", flush=True)

    sym1, anti1 = results[160]
    sym2, anti2 = results[220]

    def agree_digits(a, b):
        if a == b:
            return 40
        d = fabs(a - b) / fabs(b)
        return int(-mp.log10(d)) if d > 0 else 40

    rows = []
    print(f"\n{'n':>3} {'par':>5} {'CP':>3}  2*lambda_n (K=220)                     stable")
    print("-"*78)
    for n in range(NLEV):
        if n % 2 == 0:
            v1, v2, par = sym1[n//2], sym2[n//2], "sym"
        else:
            v1, v2, par = anti1[n//2], anti2[n//2], "anti"
        dig = agree_digits(v1, v2)
        cp = (-1)**(n+1)
        rows.append(dict(n=n, parity=par, CP=cp, value=nstr(v2, 24), stable_digits=dig))
        print(f"{n:>3} {par:>5} {cp:>+3d}  {nstr(v2,22):>24}   ~{dig} digits")

    with open(HERE / "reference_spectrum.json", "w") as f:
        json.dump(rows, f, indent=2)
    print("\nwrote reference_spectrum.json")


def polish_main(ladder=(600, 900, 1350), num_levels=40, dps=40,
                out_path=None):
    """The polished table: see the module docstring."""
    from thooft_highprec import extrapolate_ladder, polished_levels
    t0 = time.time()
    vals = {}
    for K in ladder:
        t = time.time()
        vals[K] = polished_levels(K, num_levels, dps)
        print(f"K={K} polished ({time.time() - t:.0f}s)", flush=True)
    saved = mp.dps
    mp.dps = dps
    try:
        rows = []
        for n in range(num_levels):
            v, r, corr = extrapolate_ladder(*(vals[K][n] for K in ladder),
                                            *ladder)
            err = fabs(corr)
            digits = int(mp.floor(-mp.log10(err / fabs(v))))
            rows.append(dict(
                n=n, parity="sym" if n % 2 == 0 else "anti",
                CP=(-1) ** (n + 1), value=nstr(v, 26),
                value_at_top_rung=nstr(vals[ladder[-1]][n], 26),
                truncation_rate=float(r), richardson_correction=float(corr),
                error_estimate=float(err), stable_digits=digits))
            print(f"{n:>3}  {nstr(v, 24):>28}  rate {float(r):.3f}  "
                  f"+-{float(err):.1e}  (~{digits} digits)")
    finally:
        mp.dps = saved
    payload = dict(
        provenance="generate_reference.py --polish: double-precision Ritz "
                   "vectors (thooft_spectrum.solve_sector, reciprocal "
                   "solve), exact-rational Rayleigh quotient in mpmath at "
                   f"{dps} digits, basis ladder K = {list(ladder)}, "
                   "truncation extrapolated with the per-level measured "
                   "rate (K^-6); error_estimate = |extrapolation "
                   "correction|, a conservative bar",
        ladder=list(ladder), dps=dps, runtime_s=round(time.time() - t0, 1),
        rows=rows)
    out = Path(out_path) if out_path else HERE / "reference_spectrum_polished.json"
    with open(out, "w") as f:
        json.dump(payload, f, indent=1)
    print(f"\nwrote {out.name} ({payload['runtime_s']} s)")
    return payload


if __name__ == "__main__":
    if "--polish" in sys.argv:
        polish_main()
    else:
        main()
