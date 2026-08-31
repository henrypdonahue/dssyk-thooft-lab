# Reproduction entry points for the DSSYK <-> 't Hooft stress-test repo.
#
#   make test            fast test suites of all three modules   (~1.5 min)
#   make test-all        full suites including slow ED sweeps    (~4 min)
#   make reproduce       fast deterministic checks + figures     (~50 s)
#   make reproduce-slow  ED sweeps behind the saved JSON         (~hours)
#
# Run the scripts from inside their module directory (module READMEs assume
# this).  Dependencies: pip install -e ".[dev]"  (pinned; python 3.9-3.12,
# reference tables certified on 3.9).

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
	cd thooft-target && $(PY) condensate.py
	cd thooft-target && $(PY) eigenfunctions.py
	cd thooft-target && $(PY) plot_spectrum.py
	cd syk-self-averaging && $(PY) exact_wick.py
	cd syk-self-averaging && $(PY) dirac.py
	cd syk-self-averaging && $(PY) mn_spectroscopy.py
	cd syk-self-averaging && $(PY) plot_self_averaging.py
	cd syk-self-averaging && $(PY) plot_freeness.py
	cd syk-self-averaging && $(PY) plot_pair_moments.py
	cd duality && $(PY) dictionary.py
	cd duality && $(PY) boost_mass.py
	cd duality && $(PY) largeq_anchor.py
	cd duality && $(PY) antiscrambling.py
	cd duality && $(PY) chord_moments.py --run
	cd duality && $(PY) chord_charges.py
	cd duality && $(PY) dirac_chord.py
	cd duality && $(PY) jt_ds_vacuum.py
	cd thooft-target && $(PY) flz_dconstants.py
	cd duality && $(PY) plot_chord.py
	cd duality && $(PY) plot_signflip.py

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
	cd syk-self-averaging && $(PY) mn_majorana_bench.py
	cd syk-self-averaging && $(PY) fold_edge_probe.py
	cd syk-self-averaging && $(PY) moments_fast.py
	cd syk-self-averaging && $(PY) exact_wick_dirac.py
	cd syk-self-averaging && $(PY) exact_pair_moments.py
	cd syk-self-averaging && $(PY) pair_moment_laws.py --deep
	cd duality && $(PY) chord.py --spectra
	cd duality && $(PY) chord_fold.py
	cd duality && $(PY) signflip_probe.py
	cd duality && $(PY) signflip_probe.py --deep
	cd duality && $(PY) antipodal_modes.py
	cd duality && $(PY) diagonal_chord.py --run
	cd duality && $(PY) boost_forward.py
	cd duality && $(PY) boost_inverse.py
	cd duality && $(PY) largeq_kernel.py
	cd syk-self-averaging && $(PY) pair_moment_ladder.py
	cd syk-self-averaging && $(PY) pair_moment_ladder.py --refit 1,6:2,4,6,8,10,12,14,16:18
	cd syk-self-averaging && $(PY) pair_moment_ladder.py --refit 2,4:2,4,6,8,10,12,14,16:18
	cd syk-self-averaging && $(PY) pair_moment_ladder.py --ds-extend 32,36,40,44,48
	cd thooft-target && $(PY) generate_reference.py
	cd thooft-target && $(PY) jacobi_solver.py
	cd thooft-target && $(PY) thooft_spectrum.py --csv
