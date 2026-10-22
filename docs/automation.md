# Automatic CAP parameters and convergence

For resonance energies/lifetimes, use the [resonance workflow](resonances.md).
The optional `--energies` calibration in this guide concerns time-dependent
absorption and a 1D reflection model; it does not select a resonance-optimal eta.

`cap_auto.py` now implements density-based onset selection, automatic extent,
independent radial/angular/extent refinement, a final joint check, eta/onset
calibration in a 1D absorber model, and a finite-basis rejection test.
`cap_reflection.py` supplies stationary scattering and adaptively refined
wavepacket tests. Becke remains the default integration scheme.

Requires the normal CAP build, Python, NumPy and SciPy. Matplotlib is needed
only for plots. Build with `make`; run all tests with
`bash tests/run_validation.sh`.

## A calculation from a Gaussian matrix file

For density-initialized, converged AO integrals:

```bash
python3 cap_auto.py tests/gaussian/n2_pure.mat --output my_n2_cap
```

This selects r_min=0, estimates the first r_max from the basis, selects r_cut
from the spin-summed density, and refines the grid. On success it creates
exactly two files: `my_n2_cap/CAP_n2_pure_rcut_9.75.mat` and
`my_n2_cap/CAP_n2_pure_rcut_9.75.log` for this example.
The log contains the final calculation matrices/traces and an automatic
parameter/convergence summary, following the traditional CAP log layout.
Intermediate files are temporary. Add `--diagnostics` to retain `report.json`,
full-precision `.npz` matrices and per-run logs instead.
The output directory must be new or empty. This command does not choose eta.

To also calibrate eta/onset under the approved benchmark assumptions:

```bash
python3 cap_auto.py tests/gaussian/n2_pure.mat --output my_n2_calibration \
  --energies 0.05 0.2 0.5 --duration 100 --model-end 40
```

Energies are hartree, times atomic units, distances bohr. Here the incident
range is 0.05–0.5 hartree, approximately 1.36–13.61 eV; 100 a.u. is approximately
2.42 fs. The specified energies plus a 31-point grid across their interval
are scanned; a denser 61-point/spatially refined check validates the winner.
The scan searches a finite set, not a mathematically global optimum.

**For this supplied N₂ basis, this command intentionally returns exit 3:**
the integrals converge and the 1D model passes, but the molecular basis has
insufficient absorption capacity at the selected eta/time. The report explains
why. It does not silently recommend these settings for molecular dynamics.

## Meaning of the parameters

| Option | Default | Meaning |
|---|---:|---|
| `--density-tail` | 1e-6 | Maximum initial electrons beyond the selected onset |
| `--cap-tol` | 1e-4 | Absolute spectral-norm change in orthonormalized W, bohr² |
| `--overlap-tol` | 1e-6 | Spectral-norm numerical overlap error in the same space |
| `--electron-tol` | 1e-6 | Density integral error; tightened to at most 0.1*density-tail |
| `--overlap-cutoff` | 1e-9 | Relative analytical-overlap eigenvalue cutoff |
| `--radial`, `--angular` | 200, 590 | Initial quadrature sizes |
| `--max-radial`, `--max-runs` | 1600, 50 | Hard limits, never interpreted as convergence |
| `--run-timeout` | 300 | Per-executable timeout, seconds |
| `--threads` | 2 | OpenMP threads; MKL threads fixed at 1 |
| `--energies` | absent | Enables optional model eta/onset calibration |
| `--etas` | 25 logarithmic values, 1e-5…1e-2 | Candidate eta values, hartree/bohr² |
| `--model-end` | 30 | Hard-wall position in the 1D model |
| `--duration` | 100 | Time for initial-loss/capacity screening |
| `--loss-budget` | 1e-4 | Maximum initial-density loss estimate |
| `--reflection-tol` | 1e-3 | Maximum sampled model reflection probability |
| `--required-absorption` | 0.99 | Necessary single-electron finite-basis capacity target |

`--model-end` is **not r_max**. r_max truncates integrals of a fixed Gaussian
basis; it does not put a hard wall into that basis or add outgoing-wave
functions. The model boundary is an explicit independent calibration choice.

## Numerical algorithm

1. Read real restricted/unrestricted Gaussian density records through MQCPack.
   Check the analytical electron count. Restricted one-spin versus spin-summed
   conventions are accepted only when they reproduce the expected count;
   unrestricted alpha and beta densities are summed. Complex/GHF densities
   are rejected. Density matrices/samples must be finite and nonnegative
   within roundoff tolerance.
2. Integrate density on the same Becke grid as W. Choose a density-tail quantile,
   round upward to 0.25 bohr and add a 0.25-bohr safety margin. Onset selection
   protects the supplied density only, not unknown excited/ionic states.
3. Keep a fixed canonical orthogonalizer X from analytical S for every
   comparison. Record discarded overlap eigenvectors. Compare
   `||X.T @ (W_new-W_old) @ X||_2`, numerical S, electron count, and density tail.
4. Independently increase radial count by 1.5, angular order to the next rule,
   and extent by 1.2 while increasing radial count by 1.2. Failed directions
   are refined further. A violated density-tail criterion restarts with a
   larger onset. Resource exhaustion fails explicitly.
5. Require **two consecutive** sets of successful independent probes and
   joint refinements. Density-tail changes must be <=10% of their budget.
   These are reproducible convergence estimates, not rigorous quadrature
   error bounds or proof of basis convergence.
6. Optionally scan the selected onset and offsets +0.5/+1 bohr with eta. Reject
   initial loss estimates `2 eta T Tr(P W)` above budget, model reflection above
   tolerance, or insufficient model-grid resolution. Re-converge the molecular
   integrals when a different onset wins.
7. Test the necessary single-electron capacity bound in the retained basis:
   `absorbed_fraction(T) <= 1-exp(-2 eta lambda_max(W_orth) T)`.
   A failed bound rules out the requested absorption in that space/time.
   A passed bound is **not sufficient** to establish useful absorption;
   actual molecular propagation still requires validation. This bound assumes
   Hermitian physical H and positive W, and applies to a single-electron norm.

Candidate selection now retries alternatives when a coarse reflection winner
fails the refined check or when a later onset's exact matrix fails an optimistic
capacity bound. The report records each attempted candidate's rejection reason.

The static initial-loss estimate is not a bound on the subsequently evolving
state. Nor does a density-tail threshold prove an optimal physical onset.
If no candidate meets all constraints, change the physical/basis setup or the
explicit tolerances, rather than treating failure as success.

## Output and exit codes

- **0:** requested numerical/model checks passed. If model calibration was
  requested, molecular propagation is still not certified; inspect
  `molecular_status` and `molecular_dynamics_validated`.
- **2:** convergence/calibration failed or a calculation limit was reached.
  No final CAP matrix is published; `CAP_<input-stem>_failed.log` explains the failure.
  With `--diagnostics`, individual run files and the JSON report remain.
- **3:** numerical/model checks passed and a converged matrix was written,
  but the finite molecular basis failed the necessary absorption-capacity test.

In diagnostic mode (`--diagnostics`), `report.json` separates `quadrature_converged`, `converged` (requested
numerical/model workflow), `physical_calibration`, `basis_capacity`, and
`molecular_status`. Always inspect the molecular status or process exit code;
`converged=true` alone is not a molecular-dynamics certificate. SHA-256 hashes
identify input, executable and driver. Each `.npz` contains analytical S,
spin-summed P, numerical Snum and W in the original Gaussian AO ordering.

The default `CAP_<input-stem>_rcut_<cutoff>.log` includes the selected run's detailed matrix output,
convergence status, error estimates and tolerances. Physical screening status
is included when calibration is requested. Full provenance and per-run history
are available with `--diagnostics`.
The default successful output directory contains only the matching `.mat` and `.log` files.
The cutoff is in bohr, formatted with up to six decimal places and no trailing
zeros. Diagnostic-mode `report.json` records the final filename in `output_matrix`.

For the Python driver, the Fortran `--automation_data` option writes native
64-bit stream data in the current directory: `cap_data.bin` begins with two
int64 values (AO count, electron count), followed by S, P, Snum, W as four
column-major float64 matrices. `density_samples.bin` contains float64 pairs
(r_WA, quadrature-weighted density). These are local same-machine interchange
files, not a portable archive format. The driver converts to `.npz` and removes
its temporary working directory; normal matrix-file records remain unchanged.

## Reflection and timestep tests

```bash
python3 cap_reflection.py --output my_reflection_test \
  --energies 0.05 0.2 0.5 --onsets 9.75 \
  --etas 0.000237137370566 --end 40 --propagate
```

The stationary problem uses a transparent left lead and hard wall on the
right. The wavepacket test uses Crank–Nicolson propagation with the same
quadratic profile and compares against the exact infinite-line free Gaussian
inside the unabsorbed region. Timestep and spatial spacing are refined
independently/jointly until two consecutive changes in reported observables
are <=`--propagation-tol` (default 1e-5). Limit failures return exit 2.

Wavepacket runs extend long enough for a nominal boundary round trip; they are
stress tests, often longer than the 100-a.u. loss-screening interval. Gaussian
packets have energy tails outside their central-energy interval, so their
reflection/error can exceed the stationary energy-band criterion. The code
reports that outcome instead of labeling it reflection free. The 1D model
cannot certify curved molecular surfaces, oblique incidence, Coulomb tails,
many-electron dynamics or resonance trajectories.

See the [validation summary](validation.md).
