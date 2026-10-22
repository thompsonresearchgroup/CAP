# Validation and reproduction

Run commands from the project root. Numerical CAP convergence, resonance model
validation, and molecular physical accuracy are separate checks.

The October 6 repository cleanup was checked with **40 passing regression
tests** under the runtime-checked build, followed by restoration of the optimized
binary. A clean source export built successfully, passed all **14 Python-only
tests**, and ran the bundled automatic N2 example with exactly its named matrix
and log outputs. Relocated resonance fixtures retained their recorded input
hashes and analytically verified kinetic integrals.

## Regression tests

```bash
make test-python
```

This runs grid, potential, resonance-solver and radial-propagation unit tests
using only Python dependencies and the included Lebedev tables.

The full suite additionally builds/runs Fortran, reads matrix records back,
checks Gaussian AO overlaps through I shells, tests thread reproducibility,
and exercises automatic convergence, output filenames and failure handling:

```bash
make test
# Or test a checked build and restore the optimized executable:
bash tests/run_validation.sh
```

The full suite needs the `.mat` fixtures in `tests/gaussian/`. They remain
available locally but are excluded from Git because some contain hundreds of
megabytes of unused two-electron data. To regenerate on a new machine:

```bash
bash tests/gaussian/run_gaussian_references.sh
```

Gaussian and `formchk` must be installed. The script runs jobs serially.
Two formatted checkpoints and five explicitly historical convergence logs
are retained as small reference fixtures. Their tests supplement live tests;
they do not validate a changed executable by themselves.

## Independent multicenter CAP check

With the Gaussian fixtures available:

```bash
make tests/read_cap
OPENBLAS_NUM_THREADS=1 python3 tests/validate_multicenter_cap.py --output tests/results/multicenter
```

Independent primitive-pair Gauss-Hermite integration checks two contracted S
orbitals, including their cross-center CAP element, against the live Becke grid.
The saved October 6 reference passed with maximum S-block error 1.14e-6 bohr²
against a 2e-5 bohr² threshold. [Reference report](validation/multicenter.json).

## Exact-pole resonance benchmark

The benchmark uses a manufactured radial well/barrier, not a molecular
resonance. Its Gaussian basis input files are in `tests/resonance/`. The large
`.mat` files are available locally and excluded from Git. Regenerate if absent:

```bash
cd tests/resonance
# Set GAUSS_SCRDIR to your writable Gaussian scratch directory first.
g16 < he_s48.gjf > he_s48.log
g16 < he_s96.gjf > he_s96.log
cd ../..
OPENBLAS_NUM_THREADS=1 python3 tests/validate_resonances.py --output tests/results/resonance
```

The included `s48_kinetic.txt` and `s96_kinetic.txt` contain kinetic integrals
exported from these same Gaussian inputs. The script independently verifies
their analytic normalization/order before use. Use a fresh output directory
when changing the binary or inputs; stale cached integrals are rejected.

The saved October 6 benchmark covers 48/96 Gaussian functions, onsets 4/5/6 bohr,
two quadrature grids/extents, two eta resolutions and two overlap cutoffs.
Maximum raw position error was **1.96e-4 hartree** and maximum raw width error
**1.35%**, within the declared 3e-4 hartree / 3% limits. Corrected trajectories
are retained as diagnostics; correction is not assumed to improve every result.
[Reference report](validation/resonance.json) · [Figure](validation/resonance.png).

These saved reports describe the executable hashes recorded in them. They are
historical evidence, not a claim that each later change reran the full benchmark.
Re-run the commands above to establish a new executable's scientific results.

## Physical limits

The base automatic driver certifies its numerical checks, not a molecular
resonance. Optional 1D reflection calibration and finite-basis absorption
capacity are screening tests. Compact N2/H2 bases failed the tested necessary
99%-absorption-in-100-a.u. condition; that condition is not a resonance criterion.
A target molecule still needs its electronic Hamiltonian, state identification,
and onset/grid/basis stability in both resonance position and width.
