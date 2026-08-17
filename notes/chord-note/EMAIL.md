# Draft email to the authors of arXiv:2607.05678

*Keep it under one screen. The repo carries the detail; the email earns
the click. `TODO(author)`: make the repo public before sending, and add
your signature.*

---

**Subject:** The Mₙ tower at strict double scaling — a result relevant to
your promised meson computation

Dear Professors Miyashita, Sekino, and Susskind,

I have been stress-testing the DSSYK ↔ 't Hooft proposal numerically. One
result seems worth sending before your meson-mass computation appears.

Using the exact chord solution of double-scaled SYK, I find that the
singlet bilinear channels Aₙ = Σᵢ ψᵢ (ad_H)ⁿ ψᵢ carry **zero propagating
spectral weight at N = ∞, at every λ**. For n = 1 this is the exact
identity A₁ = −2pH. For n = 2, 3 the chord amplitudes cancel to machine
precision, at every coupling and matter weight tested. Random operators
propagate normally; the cancellation is specific to bilinears built from
H's own couplings — that is, your Mₙ tower. Finite-N diagonalization
shows the channels turning on along fixed p, λ → 0.

The consequence: the meson tower is a 1/N effect, invisible at every
order of the double-scaled expansion. The spectrum match must be
extracted in the fixed-p corner, where singlet self-averaging is only
polynomial (violation amplitude N^−(p−1)/2, exact). The chord machinery
cannot shortcut it. This seems directly relevant to how you set up the
promised computation.

The repository below contains this result and the supporting program:
verified corrections to the operator dictionary (for example
M₁ = −(p/2)iH, forced by anti-Hermiticity), the exact λ → 0 spectral
scans, and — taking the static-patch motivation seriously — a
quantitative statement of what the tomperature map must do at OTOC level
(flip the scrambling sign; the Euclidean fold does not do this, even at
N = ∞). The engine is validated four independent ways, and everything is
asserted by a test suite (`make test`, 353 tests, about 3 minutes).

https://github.com/henrypdonahue/dssyk-thooft-lab

I would welcome any correction or comment. The point of the exercise is
to make the proposal's best test runnable.

Best regards,
Henry Donahue
