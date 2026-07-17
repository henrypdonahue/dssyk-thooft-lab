"""Generate a certified high-precision reference table.
Runs two basis sizes and reports agreeing digits per level."""
from mpmath import mp, nstr, fabs
from thooft_highprec import solve_sector_mp
import json, time

mp.dps = 50
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

with open("reference_spectrum.json", "w") as f:
    json.dump(rows, f, indent=2)
print("\nwrote reference_spectrum.json")
