# Reproduction and release checks

The Makefile is the command inventory. Use its four targets: `test`,
`test-all`, `reproduce`, and `reproduce-slow`. Slow reproduction includes
large ED ensembles, exact rational interpolation, multiprocessing campaigns
and high-precision reference generation. Allow many hours and inspect logs;
the short unit suite does not regenerate those artifacts.

## Preserve the working tree and record the run

After installing `.[dev]` in Python 3.9–3.12, use:

```bash
python scripts/validate.py test-all
python scripts/validate.py reproduce
python scripts/validate.py reproduce-slow
```

The runner reads the actual Makefile recipes, copies tracked and nonignored
files into an isolated snapshot, and executes commands there. It never
copies ignored third-party papers or the Git database. Existing uncommitted
source changes are included, recorded by SHA-256, and left untouched.
The output directory must be new. `--list` prints the commands, and
`--output PATH` chooses another output location.

Each run saves a machine-readable manifest with the source commit, dirty
status, source hashes, Python/platform/package versions, thread limits,
command lines, elapsed times and exit codes. Every command has its own log.
Failures do not prevent later commands from running; any failure makes the
runner exit nonzero. Keep the complete manifest and logs with a release.
Interruptions are recorded; a partially completed run is not a pass.

Generated files stay in the snapshot. At the end of a run the runner
compares every committed artifact the run rewrote with the working-tree
copy it started from, and records the result in the manifest under
`artifact_comparison`. JSON and CSV files are compared value by value.
Numbers are compared by relative difference; a numeric array is measured
against its largest entry, and a value below 1e-14 of the largest
magnitude in its own record counts as a roundoff zero. Runtime fields are
ignored. Text and structural changes are listed. Images are flagged as
changed bytes and not compared.

Review those differences before copying selected outputs back. Runtime
fields, platform rounding and plot metadata need not be identical.
Compare exact rational strings exactly; compare floating-point results
with the producer's justified tolerances. On 2026-10-14 (macOS, Python
3.9) the fast `reproduce` target regenerated every committed artifact to
within 1e-11 relative. A passing
script exit code is not sufficient if the script merely prints a failed
condition: inspect the log and run the corresponding asserting tests after
regeneration.

## Environment

`pyproject.toml` pins direct scientific, test and plotting dependencies.
It does not lock every build/transitive dependency or the BLAS implementation.
The runner records installed versions rather than claiming bitwise identity.
Use `make ... PY=python` to ensure Make uses the activated interpreter.
If the platform's `make` is unavailable, the Python runner runs the same
recipes without a shell or a Make installation.

CI runs the fast suites and the 't Hooft physics gate on the configured
Python versions. The slow campaigns must be separately recorded before a
reproducibility claim. CI passing is not a substitute for them.

## Scientific comparison contracts

- Majorana: ψᵢ² = 1, normalized trace, even interaction order p.
- Paper λ = p²/N; chord λ = 2p²/N. Always label the convention.
- Distinguish fixed-p N→∞, fixed-λ double scaling, and a subsequent λ→0.
- State whether moments use physical J units or E[Tr H²]/dim = 1.
- State ensemble subtraction versus per-realization trace subtraction.
- Distinguish ratios of ensemble moments from ensemble means of ratios.
- Pair cross-correlations need not have positive spectral weights.
- Identity subtraction does not remove the entire conserved sector.
- Exact interpolation in N at a fixed p requires a justified degree bound;
  extrapolating a coefficient formula across finitely many p is a separate
  claim and does not prove a general-p identity.
- An extrapolation spread is a sensitivity diagnostic, not a statistical
  confidence interval. Disorder means require sample counts, seeds, and
  uncertainties with correlations propagated where relevant.

## Before contacting authors

Run the full suites and both reproduction targets, then rerun the suites on
the regenerated snapshot. Reconcile any generated differences; retain the
logs and source manifest. Check the concise note against the exact release
commit. Verify primary-source versions and equation references. Ensure all
new source/data files are tracked, `git diff --check` passes, and the release
has a clean working tree. Do not include ignored local email drafts or
third-party papers accidentally.
