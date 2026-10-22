# CAP resonance energies and lifetimes

The CAP has now passed an end-to-end resonance benchmark against an independent
outgoing-wave pole, including onset, Gaussian-basis, quadrature, eta-scan, and
overlap-cutoff checks. This establishes the tested implementation and model
accuracy; it does not certify a particular molecular resonance without its
electronic Hamiltonian and convergence study.

**Do not use the 99%-absorption-in-100-a.u. screen to accept or reject a
resonance calculation.** That is a time-dependent wavepacket criterion.
Likewise, `cap_auto.py --energies ...` calibrates a 1D reflection model, not
resonance energies or lifetimes. The October 2–3 absorption-capacity failures
do not, by themselves, show that this CAP is unsuitable for resonances.

## Resonance workflow

The integral executable writes real, unscaled W. Form `H(eta)=H-i*eta*W`.
For an isolated resonance eigenvalue `E=E_R-i*Gamma/2`, the width is
`Gamma=-2*Im(E)` and the lifetime is `1/Gamma` in atomic units.
The added `cap_resonance.py` handles dense **real symmetric H, S and W** in
the same basis, canonical overlap orthogonalization, complex eigenvalues,
state tracking, and eta-trajectory analysis.

It computes the complex-symmetric derivative with the c-product:

```
dE/deta = -i (v.T W v)/(v.T v)
U = E - eta*dE/deta
```

The transpose is intentional; replacing it with a conjugate transpose gives
the wrong derivative. Interior minima of `|eta*dE/deta|` identify candidate
stationary points. Endpoints are rejected, neighboring-state overlap is
checked, and the derivative is compared with a finite difference. Corrected
values U and their sampled logarithmic speeds are also reported. The tool
does **not** assume that the correction improves every finite-basis result.
These conventions follow the [CAP trajectory literature](https://iopenshell.usc.edu/pubs/pdf/jpcl-5-3078.pdf)
and the [smooth Voronoi CAP paper](https://pubs.acs.org/doi/10.1021/acs.jctc.5b00465).

1. Choose a molecule, resonance symmetry/state, electronic-structure method,
   basis, and energy zero. A neutral RHF Fock matrix is not a correlated
   electron-attachment Hamiltonian.
2. Generate and converge W at several onsets with the CAP executable. Keep H,
   S and W in exactly the same AO, orbital, or configuration basis. For a
   many-electron CI Hamiltonian, transform and lift the one-body CAP into the
   CI space before using this solver.
3. Save a NumPy archive containing arrays named `H`, `S`, `W`. Use a broad
   eta range and a target complex energy in that Hamiltonian's energy zero.
4. Inspect **all** stationary candidates and branch continuity. Recompute
   after changing the eta range, scan resolution and seed if needed.
5. Compare resonance positions **and widths** across onsets, basis enrichment,
   quadrature refinement and overlap cutoffs. Only then assign an uncertainty
   and report a lifetime. A small eta derivative alone does not prove a
   physical resonance or a complete basis.

This solver does not implement electronic correlation or non-Hermitian
EOM-CC left/right eigenvectors. It supports a one-electron Hamiltonian or a
real symmetric CI-type Hamiltonian whose matrix elements have already been
computed. The Gaussian CAP generator remains usable with other downstream
solvers through its normal matrix records.

## Run the saved, validated example

From the CAP directory, with a new output filename:

```bash
python3 cap_resonance.py \
  examples/inputs/model_matrices.npz \
  --target 0.12 -0.01 --seed-eta 0.001 \
  --eta-min 1e-6 --eta-max 0.1 --eta-points 401 \
  --overlap-cutoff 1e-10 --output my_resonance.json
```

This is a **manufactured radial well/barrier**, not helium's physical resonance.
Helium Gaussian input files supply normalized basis metadata and kinetic
integrals; the physical potential is replaced explicitly by the benchmark
potential. The target is an approximate branch seed, not an eta-selection
reference.

Exit 0 means an interior stationary candidate was found, not that a molecular
resonance has been certified. The output retains
`resonance_validated=false` and
`status=requires_onset_grid_and_basis_convergence`. Exit 2 means no interior
stationary point was found. Widths and lifetimes are reported in hartree/eV
and atomic units/fs, with the factor of two in the width applied explicitly.

## Executed accuracy

The exact model pole is:

| Quantity | Independent outgoing-wave reference |
|---|---:|
| Position | 0.116675732189 hartree = 3.174908414 eV |
| Width | 0.017692458919 hartree = 0.481436333 eV |
| Lifetime | 56.521256 a.u. = 1.367183803 fs |

Six live-CAP calculations (48/96 S functions, onsets 4/5/6 bohr) pass the
declared model tolerances of 3e-4 hartree in position and 3% in width.
Observed maximum errors are **1.96e-4 hartree** and **1.35%**, respectively.
The exact pole is used only for validation; eta is selected from interior
trajectory stationarity. Both raw and corrected results are retained.

See [validation](validation.md), the [saved reference report](validation/resonance.json),
and the [trajectory figure](validation/resonance.png).

## Repeat the implementation and model checks

The full benchmarks require the generated Gaussian fixtures described in
[validation](validation.md). The compact saved-model CLI example above runs
without Gaussian. NumPy and SciPy must be installed.
From the CAP directory, choose fresh output directories:

```bash
bash tests/run_validation.sh
OPENBLAS_NUM_THREADS=1 python3 tests/validate_multicenter_cap.py --output tests/my_validation/multicenter
OPENBLAS_NUM_THREADS=1 python3 tests/validate_resonances.py --output tests/my_validation/resonance
```

The first command tests a runtime-checked build and restores the optimized
executable. The other commands exercise that executable and save independent
integral and exact-pole comparisons. Both must return zero. The resonance
report includes raw and corrected results; its declared accuracy limits apply
to raw stationary energies and widths. Corrections remain diagnostics, with no
assumption that they improve the result. Historical reports are preserved when
`--output` points to a fresh directory.
