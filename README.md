# Finite Phase Invariants validation package v3

This package matches the **15 numbered results** in the attached
`Finite_Phase_Invariants_WhitePaper.docx`. It extends v2 with Lemmas 5, 14 and 15,
broader Lemma 7 tests, and a coverage map using the attachment's numbering.
The white paper itself has not been modified.

## Run

Requires Python 3.10 or later. Standard library only; no pip installation,
network access, external data or NumPy dependency is required.

```sh
python run_all.py --output validation_results.json
python demo.py
```

Run from the extracted `validation_v3` directory, without `-O`. The runner
refuses execution with disabled assertions. Failed assertions or exceptions
produce a nonzero exit status. `validation_results.json` records the completed
run on Python 3.12.14. Fixed seeds and rational arithmetic make principal
results reproducible; the independent floating-point phase comparison may
vary slightly by platform.

## Files

- `COVERAGE.md`: map for Definitions 1–4 and every numbered result 1–15,
  including test scope and limitations.
- `coverage_map.csv`: compact machine-readable version of the result map.
- `test_lemma_5.py`: exact half-integer rounding and boundary checks.
- `test_lemma_7.py`: weighted quantiles, integrated ES and composed scaling.
- `test_lemma_14.py`: symbolic bare-action metric derivatives and torus fixture.
- `test_lemma_15.py`: exact finite-dimensional local channel invariance.
- `exact_matrix.py`: rational complex matrix operations for Lemma 15.
- `test_algebra.py`, `test_phase.py`, `test_certificate.py`, `test_gating.py`,
  `test_risk.py`: retained finite-core tests, with updated metadata.
- `detector.py`: normalization, full-coordinate embedding, exact winding and
  perturbation certificate.
- `gating.py`, `pipeline.py`, `demo.py`: state machine, configuration binding,
  and a synthetic end-to-end demonstration.
- `run_all.py`, `validation_results.json`: nine-module runner and recorded results.
- `CHANGELOG.md`: changes from v2.

## Core input and certificate contract

Use integers, Fraction objects, or decimal strings for weights, times, gains
and error bounds. Binary floats are rejected by the detector. Graphs must be
nonnegative, square, have zero diagonal and positive total exposure. Stable
node IDs allow matrix reorderings. Without IDs, row order is the fixed registry.

`epsilon` bounds the TOTAL entrywise absolute measurement error. It is supplied
externally, not inferred by the detector. Include measurement, missing-data,
imputation and quantization errors. The default 2x2 lattice has capacity
N(N-1)<=16; larger registries require a prespecified larger rational baseline.

`accepted=True` is conditional certification of a finite model label, subject
to the externally supplied error bound and fixed embedding. It is not a
certificate of market stress. Failure produces an unknown label, not proof
of instability. Exact arithmetic does not make an incorrect premise true.

The gate requires confirmed baseline, activation and recovery runs. Invalid
evidence emits `None`, resets counters, and preserves latent ARMED state.
Consumers must not convert `None` to zero. Freshness, gap and configuration
version checks are enforced. The pipeline hashes its fixed configuration;
instantiate a new pipeline after intentional model changes instead of mutating
configuration during a run.

## Validation limits

All 15 numbered results have mapped checks, but the suite is not a formal proof
assistant. Lemma 14 has local symbolic and one torus check, not arbitrary
continuum integration. Lemma 15 has exact finite-dimensional examples, not an
infinite-dimensional channel proof. General probability statements are checked
on selected finite laws. No market data, continuum-topology equivalence,
literature attribution, or predictive performance is validated.
