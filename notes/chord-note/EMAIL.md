# Draft email to the authors of arXiv:2607.05678

*Keep it under a screen. The repo carries the detail; the email earns the
click. `TODO(author)`: make the repo public before sending, and add your
signature.*

---

**Subject:** The Mₙ tower at strict double scaling — a result you may want
before your promised meson computation

Dear Professors Miyashita, Sekino, and Susskind,

I've been stress-testing the DSSYK∞ ↔ 't Hooft proposal numerically, and one
result seems worth sending before your promised computation of the meson
masses appears.

Using the exact chord solution of double-scaled SYK (validated against a
brute-force diagram enumerator, the published closed forms, and finite-N
diagonalization to ~1%), I find that the singlet bilinear channels
Aₙ = Σᵢψᵢ(ad_H)ⁿψᵢ carry **exactly zero propagating spectral weight at
N = ∞, for every λ** — an algebraic cancellation of which A₁ = −2pH is the
n = 1 case, holding to machine precision for n = 2, 3 as well. Random
operators O_Δ propagate normally; the cancellation is specific to bilinears
built from H's own couplings — i.e., precisely your Mₙ tower. Finite-N
diagonalization shows the channels turning on along fixed p, λ→0.

The consequence: the meson tower is a 1/N effect, invisible at every order
of the double-scaled expansion, so the spectrum match must be extracted in
the fixed-p corner — where singlet self-averaging is only polynomial
(violation amplitude N^{−(p−1)/2}, exact). The exact chord machinery cannot
shortcut it. This seems directly relevant to how you set up the promised
computation.

The repository below contains this result and the supporting program: a few
verified corrections to the operator dictionary (e.g. M₁ = −(p/2)iH, forced
by anti-Hermiticity), the exact λ→0 spectral scans, and — taking the static
patch motivation seriously — a quantitative specification of what the
tomperature map must accomplish at OTOC level (it must flip the scrambling
sign; the Euclidean-fold mechanism does not do this even at N = ∞).
Everything is asserted by a test suite (`make test`, 296 tests, ~3 min).

https://github.com/henrypdonahue/dssyk-thooft-lab

I'd welcome any correction or comment — the point of the exercise is to make
the proposal's best test runnable.

Best regards,
Henry Donahue
