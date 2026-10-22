# Smooth Voronoi CAP guide

Start with [the short tutorial](../README.md). Run the commands below from
the project root. This guide describes the maintained MQCPack implementation.

## Recommended automatic run

Build once, then provide a Gaussian binary matrix file containing the geometry,
basis, analytical overlap and real density records:

```bash
make
python3 cap_auto.py molecule.mat --output my_cap
```

The output directory must be new or empty. A successful run writes exactly:

- `my_cap/CAP_molecule_rcut_<cutoff>.mat`: the real, unscaled AO CAP matrix and supporting records.
- `my_cap/CAP_molecule_rcut_<cutoff>.log`: settings, geometry, overlap/CAP matrices, traces and a short
  convergence summary.

For example, `n2_sto-3g.mat` at 9.75 bohr produces
`CAP_n2_sto-3g_rcut_9.75.mat` and `CAP_n2_sto-3g_rcut_9.75.log`.

The driver selects `r_cut` from the supplied density, keeps `r_min=0`, estimates
`r_max`, then checks radial, angular and extent convergence. It does not select
eta unless optional absorption-model calibration is requested. For your own
input, see the [Gaussian input example](../examples/README.md). Run the included
example with `bash examples/run.sh fresh_run` from the project root.

## Reading the log

| Log entry | Meaning |
|---|---|
| `r_cut` | Absorber onset distance, bohr |
| `r_min`, `r_max` | Integration limits, bohr |
| `n_radial`, `n_lebedev` | Final radial and angular resolutions |
| `Analytical Overlap`, `Numerical Overlap` | Reference and integrated AO overlaps |
| `CAP matrix` | Unscaled W; apply `H_eff = H - i eta W` downstream |
| `Numerical CAP Trace` | Trace of W; a summary, not an accuracy test |
| `Quadrature converged` | Whether the automated integration checks passed |
| `overlap_error` | Spectral-norm overlap error in the retained orthonormal space |
| `electron_error` | Error in integrated electron count |
| `tail` | Electrons outside the selected onset |
| `density_expectation` | Initial density expectation `Tr(P W)` |
| `Overlap directions retained` | Basis directions kept after overlap thresholding |
| `Exit status` | 0: requested checks pass; 2: failure; 3: absorption-capacity rejection |

The normal log omits the APO vector dump, its sum/dimensions, internal hashes
and per-run progress. The matrix file retains its existing records for downstream
compatibility. Add `--diagnostics` to retain full convergence/provenance artifacts.
On convergence failure (exit 2), the normal run leaves `CAP_<input-stem>_failed.log` without
publishing a final matrix.

Converged integrals do not establish a molecular resonance or successful
propagation. Resonance analysis is a separate workflow in the
[resonance guide](resonances.md); full automation options are in the
[automation guide](automation.md).

## Definition

Following Sommerfeld and Ehara, JCTC 11, 4627–4633 (2015), equations 6–8
([paper](https://doi.org/10.1021/acs.jctc.5b00465)):

```text
d_A² = |r - R_A|²
w_A  = [d_A² - min_B(d_B²) + a²]^-2
r_WA = sqrt(sum_A(w_A d_A²) / sum_A(w_A))
W(r) = max(0, r_WA - r_cut)²
W_mu,nu = integral chi_mu(r) W(r) chi_nu(r) d³r
H_eff = H - i eta W
```

The paper's smoothing parameter is `a=1 bohr`. The implementation rescales
weights before taking their ratio to avoid large inverse powers. It evaluates
contracted Gaussian AOs through MQCPack, including SP, Cartesian, and real-pure
shells. Gaussian S–I conventions are covered by live overlap tests. It does
not apply eta or perform resonance searches or AO-to-MO transformations.

## Options

| Option | Meaning | Default |
|---|---|---|
| `-f FILE` | List containing `1`, then one Gaussian matrix path | required |
| `--output FILE` | Override standalone matrix filename (e.g. `CAP_Mat.mat` for older scripts) | derived from input and cutoff |
| `--lebedev_grid FILE` | Supplied angular table | required |
| `--n_lebedev N` | Must match the table's supported point count | 1454 |
| `--n_radial N` | Gauss–Legendre radial points | 500 |
| `--r_min X` | Radial lower limit, bohr | 0 |
| `--r_max X` | Radial upper limit, bohr | 20 |
| `--r_cut X` | Smooth Voronoi onset, bohr | 3 |
| `--smoothing_length X` | Positive smoothing length a, bohr | 1 |
| `--becke_grid` | Atom-centered MQCPack Becke integration | **on** |
| `--centroid_grid` | Older shared molecular-centroid grid | off |
| `--auto_rmax` | Basis/geometry initial extent estimate; exclusive with `--r_max` | off |
| `--validate` | Overlap, symmetry, eigenvalue and closure diagnostics | off |
| `--unit_potential` | Integrate W=1 for validation | off |
| `--equidistant_radial` | Composite Simpson comparison rule | off |
| `--help` | Built-in usage | — |

Paths in the list are relative to the process working directory. Angular table
rows are `azimuth_degrees polar_degrees normalized_weight`, with weights summing
to one. The program applies `4*pi` once. Only positive angular weights are
accepted; Lebedev rules with negative weights are rejected. Lengths must be
finite. Invalid input fails instead of writing a new result.

## Grid choice and convergence

Use `r_min=0`. By default, identical finite radial/angular rules are centered
on each atom and partitioned with MQCPack's original Becke radii and weights.
`--becke_grid` explicitly selects this default. `--centroid_grid` instead
centers a common spherical grid at the arithmetic mean of nuclear positions. The radial bounds then refer to each atom. All AOs share
the resulting quadrature: `W = sum_g weight_g W(r_g) chi_g chi_g^T` is positive
semidefinite up to roundoff. There is no pair screening or eigenvalue clipping.

The potential is still smooth Voronoi in both modes; Becke weights partition
the integration domain only. This differs from the paper's primitive-product
centered quadrature. The paper's reported accuracy cannot be assumed for these
grids. The default Becke grid resolves tight cores much better in the tested
molecules. Its radii lookup rejects ghost centers or unsupported elements;
coincident centers are also invalid. Neither mode is certified here for ECPs,
ghost-center CAP definitions, or arbitrary large/extended systems.

For the supplied small molecules, start with Becke, 1202 angular points,
300 radial points, and `r_max=24`. Verify at 3470 angular and 500 radial points
and larger extent. Tight-core overlap convergence does not establish CAP-tail
convergence; compare the entire W matrix. The automatic driver performs independent radial/angular/extent refinements
and joint checks; see [automation](automation.md).

`--auto_rmax` estimates an initial radius for each shell from
`sqrt((50+2*l)/(2*alpha_min))`. It adds the shell-center distance to the centroid
in centroid mode, or the largest distance to any atom-centered grid in Becke
mode, then takes the maximum over shells. This conservatively encloses these
estimated primitive-tail regions, but is **not a rigorous integral error
bound**: contraction coefficients, CAP growth, and radial resolution still
need convergence tests. The chosen radius is printed; `r_cut` is unchanged.
Explicit `--r_max` and `--auto_rmax` together are rejected.

Diffuse functions need a larger extent. A useful initial tail estimate is
`sqrt(50/(2*alpha_min))` bohr; for the centroid grid also include the largest
atom–centroid distance. This is a starting estimate, not an error bound.
Increase radial resolution when increasing the extent. Atom-centered finite
spheres also omit some partitioned tail contributions: converge their radius.

Negative eigenvalues at machine precision are normal numerical roundoff.
Significant negatives fail validation. A coarse positive grid can still give
an inaccurate positive matrix. The reported symmetry check is a consistency
check; use Gaussian overlap and convergence as independent accuracy checks.

## Output and performance

The standalone executable writes `CAP_<input-stem>_rcut_<cutoff>.mat` in
the working directory. The input basename is used without its extension;
cutoffs use up to six decimal places without trailing zeros. Redirect standard
output to your chosen log filename. Preserve files before repeating a run
with the same input name and cutoff. Records are analytical overlap, numerical overlap, the unscaled
CAP under generic/alpha/beta names, and the legacy `APO Vector` linear integral
`integral chi_p W d³r`. That vector is not `C^T W C` or a channel probability.
MQCPack's matrix reader exposes a vector record as a diagonal matrix.

AO evaluation and the potential are computed once per point. Fixed-size blocks
form weighted Gram matrices, with deterministic block reductions and parallel
point evaluation. Matrix storage scales as O(nBasis²), rather than one dense
matrix per thread. AO scratch is O(256*nBasis + threads*nBasis). The centroid
grid streams coordinates; the default MQCPack Becke grid retains coordinates,
weights and partition metadata, using O(nAtoms*nRadial*nAngular) memory.

`make` uses `-O3` without fast-math. Start with 2–4 OpenMP threads and
`MKL_NUM_THREADS=1`; more threads can increase overhead for small systems.
`make strict` enables runtime checks and replaces the production executable;
run `make -B all` afterward to restore the optimized build.

## Validation and dependencies

`bash tests/run_validation.sh` builds with bounds checks, executes current
Gaussian/analytic/PSD/record/thread/input tests, and restores the optimized
binary. Python tests also retain some explicitly historical archived-log
checks. Gaussian fixtures can be regenerated with
`bash tests/gaussian/run_gaussian_references.sh` (serial jobs).

Uses MQCPack 26.9.3 `loadGTOBasisSet`, `basisSetValuesList1`,
`setup_quadrature_gauss_legendre1d`, `becke1988_radii_from_atomic_numbers`,
`gridBeckeMolecular`, and the molecular grid block API, plus MQCPack Gaussian
I/O, matrices and eigensolvers. See the
[MQCPack release notes](https://github.com/MQCPack/mqcPack/blob/master/releaseNotes.txt).
The local library's Gaussian pure transformations and Cartesian ordering were
validated through I; adding library calls without supported semantics is not
necessary. No changes to the shared MQCPack installation are required.

See [automatic parameter selection](automation.md) for
reflection tests, r_cut selection, and the distinction between integration
extent and the spatial reach of the basis.

## Automatic workflow

Use `cap_auto.py` for density-based onset initialization and independent plus
joint numerical refinement. `cap_reflection.py` supplies model scattering and
adaptively refined propagation checks. See [usage](automation.md) and
[validation summary](validation.md), including deliberate rejection of
insufficient molecular basis capacity.
