# Reproduction entry points for the DSSYK <-> 't Hooft stress-test repo.
#
#   make test            fast test suites of all three modules   (~15 s)
#   make test-all        full suites including slow ED sweeps    (~minutes)
#   make reproduce       fast deterministic checks + figures     (~1 min)
#   make reproduce-slow  ED sweeps behind the saved JSON         (~hours)
#
# Scripts are run from inside their module directory (module READMEs assume
# this).  Dependencies: pip install -e ".[dev]"  (pinned; python 3.9).

PY ?= python3

.PHONY: test test-all reproduce reproduce-slow

test:
	cd thooft-target && $(PY) -m pytest -q -m 'not slow'
	cd syk-self-averaging && $(PY) -m pytest -q -m 'not slow'
	cd duality && $(PY) -m pytest -q -m 'not slow'

test-all:
	cd thooft-target && $(PY) -m pytest -q
	cd syk-self-averaging && $(PY) -m pytest -q
	cd duality && $(PY) -m pytest -q

# Fast, deterministic: the physics gate, the closed-form/identity printouts,
# the single-realization spectroscopy (Nc = 10), and figures from saved JSON.
reproduce:
	cd thooft-target && $(PY) validate.py
	cd syk-self-averaging && $(PY) exact_wick.py
	cd syk-self-averaging && $(PY) dirac.py
	cd syk-self-averaging && $(PY) mn_spectroscopy.py
	cd syk-self-averaging && $(PY) plot_self_averaging.py
	cd syk-self-averaging && $(PY) plot_freeness.py
	cd duality && $(PY) dictionary.py
	cd duality && $(PY) boost_mass.py
	cd duality && $(PY) largeq_anchor.py
	cd duality && $(PY) antiscrambling.py

# The disorder-ensemble sweeps and high-precision reference generation that
# produced the committed JSON/CSV artifacts.  Hours, not minutes.
reproduce-slow:
	cd syk-self-averaging && $(PY) self_averaging.py
	cd syk-self-averaging && $(PY) validate_syk.py
	cd syk-self-averaging && $(PY) largeq_bench.py
	cd syk-self-averaging && $(PY) moments_pipeline.py
	cd syk-self-averaging && $(PY) qed_campaign.py
	cd syk-self-averaging && $(PY) baryons.py
	cd syk-self-averaging && $(PY) mn_spectroscopy.py --big
	cd syk-self-averaging && $(PY) fold_bench.py
	cd syk-self-averaging && $(PY) fold_bench.py --scan
	cd syk-self-averaging && $(PY) freeness.py
	cd thooft-target && $(PY) generate_reference.py
	cd thooft-target && $(PY) jacobi_solver.py
	cd thooft-target && $(PY) thooft_spectrum.py --csv
