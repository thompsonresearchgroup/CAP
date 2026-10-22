# Smooth Voronoi CAP

Generate a smooth Voronoi complex absorbing potential (CAP) in a Gaussian
atomic-orbital basis. The program writes the **real, unscaled matrix W**;
a downstream calculation applies `H_eff = H - i eta W`.

The automatic driver selects the onset and converges the integration grid.
A successful calculation produces **one matrix file and one log file**.

## Requirements

- GNU Fortran with OpenMP and GNU Make.
- MQCPack **26.9.3**, built with compatible 8-byte default real/integer flags.
- Intel MKL (`libmkl_rt`) for the current build configuration.
- Python **3.10 or newer**, NumPy and SciPy.
- Gaussian is needed to generate new molecular matrix inputs, but not to run
  the supplied example.

```bash
python3 -m pip install -r requirements.txt
make mqcroot=/path/to/mqcPack MKLLIB=/path/to/mkl/lib
```

On the original workstation, `make` alone uses the installed MQCPack and MKL
paths. Third-party libraries are not bundled.

## Run a calculation

From the project root:

```bash
python3 cap_auto.py examples/inputs/n2_pure.mat --output runs/n2
```

The output directory must be new or empty. This example selects an onset of
9.75 bohr and produces:

```text
runs/n2/
  CAP_n2_pure_rcut_9.75.mat
  CAP_n2_pure_rcut_9.75.log
```

For your own calculation, replace the input with your Gaussian `.mat` file:

```bash
python3 cap_auto.py /path/to/n2_sto-3g.mat --output runs/n2_sto-3g
```

Output names contain the input basename and selected cutoff in bohr. Cutoffs
use up to six decimal places without trailing zeros. A failed calculation
writes `CAP_<input-stem>_failed.log` without publishing a final matrix.

The driver automatically selects `r_cut`, keeps `r_min=0`, estimates `r_max`,
and refines radial/angular resolution and extent. It starts with 200 radial
and 590 Lebedev points. To change the starting resolution:

```bash
python3 cap_auto.py molecule.mat --output runs/molecule --radial 300 --angular 1202
```

The log lists the **final** grid settings, geometry, overlap and CAP matrices,
traces, and convergence checks. Look for `Quadrature converged: True` and
`Exit status: 0`. The basic calculation does not choose eta or certify a
molecular resonance. Add `--diagnostics` only when you need intermediate
matrices and a JSON report.

## Documentation

- [User guide](docs/user-guide.md): input preparation, manual options, log fields,
  matrix conventions, and the CAP formula.
- [Automatic parameters](docs/automation.md): convergence criteria, tolerances,
  optional reflection calibration, and failure statuses.
- [Examples](examples/README.md): included Gaussian input and executed output.
- [Resonance analysis](docs/resonances.md): H/S/W inputs, eta trajectories,
  widths, and lifetimes.
- [Validation and tests](docs/validation.md): what has been checked and how to
  reproduce the tests.
- [Repository contents](docs/repository.md): files to publish and local archives.

## Tests

```bash
make test-python                # no Gaussian matrix fixtures or Fortran build needed
make test                       # full suite; requires the local Gaussian fixtures
bash tests/run_validation.sh    # checked build + full suite; restores optimized binary
```

Large generated test matrices are excluded from Git. On a fresh clone,
regenerate them before the full suite:

```bash
bash tests/gaussian/run_gaussian_references.sh
```

That command requires `g16` and `formchk` on PATH; both commands can be overridden
with `GAUSSIAN_EXE` and `FORMCHK_EXE`. The small N2 example input is bundled.

## Scope and reference

The implementation has been tested against analytical/Gaussian overlaps,
independent multicenter integration, and an exact-pole radial resonance model.
A molecular resonance still needs an appropriate electronic Hamiltonian and
state-specific onset/grid/basis convergence.

The potential follows Sommerfeld and Ehara, *J. Chem. Theory Comput.* **11**,
4627–4633 (2015), [DOI: 10.1021/acs.jctc.5b00465](https://doi.org/10.1021/acs.jctc.5b00465).
