# Voronoi-CAP Integral Generator (`noCAPs2withCAPs`)

A standalone Fortran program that computes **smooth Voronoi complex absorbing
potential (CAP) integrals** over a Gaussian basis, plus the per-orbital **APO
weight vector** `W_p`, and writes them into a Gaussian-format matrix file for
downstream use (e.g. by SNOWMAN / TDCIS).

The CAP follows the smooth Voronoi-isosurface construction of
**Sommerfeld & Ehara, *J. Chem. Theory Comput.* **11**, 4627 (2015)** (eqs. 6–8):

```
r_WA(r) = sqrt( Σ_a w_a r_a² / Σ_a w_a )
w_a     = 1 / (r_a² − r_nearest² + 1 a0²)²
W(r)    = 0                  if r_WA ≤ r_cut
W(r)    = (r_WA − r_cut)²    if r_WA  > r_cut
```

Each primitive-pair integral is evaluated on a single large grid centered at the
Gaussian product center, using a **Lebedev angular grid** × an **equidistant
radial grid** (leading trapezoidal term of the Euler–Maclaurin formula, as is
standard for quantum-chemistry grids).

> **Note on η (CAP strength).** This generator produces **η-independent**
> integrals — it computes `W(r) = (r_WA − r_cut)²` only. The `−iη` prefactor is
> applied *downstream* (SNOWMAN multiplies by `−iη` when assembling `H − iηW`).

Authors: L. M. Thompson, A. M. Kinyua, P. D. Simpson — Louisville, KY, 2023–2026.

---

## Requirements

Same toolchain as SNOWMAN:

- `gfortran` (Fortran 2008)
- [**mqcPack**](https://github.com/MQCPack) installed at `$HOME/mqcPack/GNU/`
  (provides `lib/libmqc.a` and the module files in `mod/`)
- **Intel MKL** on the linker path (`$MKLROOT` set; links `-lmkl_rt`)
- OpenMP (the integral loop is threaded)

## Build

```bash
make
```

This produces the `noCAPs2withCAPs` executable. The build flags
(`-std=f2008 -fdefault-real-8 -fdefault-integer-8 -fopenmp`) and the mqcPack /
MKL paths live in the `makefile`; edit `mqcroot` there if your mqcPack lives
elsewhere. `make clean` removes the executable and any `*.o` / `*.mod`.

## Input

A small text file (passed with `-f`) listing the Gaussian matrix files to
process. **Line 1** is the number of files; the following lines are the matrix
file paths. See `input.txt.example`:

```
1
H2-_6-31++g.mat
```

The matrix file is a standard Gaussian unformatted matrix file (`.mat`)
containing the basis-set / shell data (shell-to-atom map, shell types,
primitive exponents, contraction coefficients, coordinates, overlap).

## Running

```bash
export OMP_NUM_THREADS=32

./noCAPs2withCAPs \
    -f input.txt \
    --lebedev_grid lebedev/lebedev_3470.txt \
    --r_min 0.0 \
    --r_max 100.0 \
    --n_lebedev 3470 \
    --n_radial 500 \
    --r_cut 6.22
```

This writes **`CAP_Mat.mat`** in the current directory, containing the analytic
overlap, the numerical overlap (grid-quality check), the numerical CAP matrix,
and the APO weight vector `W_p`. Rename it per geometry as you like
(e.g. `CAP_r6.22.mat`).

### Command-line options

| Flag | Default | Meaning |
|---|---|---|
| `-f <file>` | *(required)* | Input file listing the `.mat` file(s) to process |
| `--lebedev_grid <file>` | *(required)* | Path to a `lebedev/lebedev_<N>.txt` angular-grid table |
| `--n_lebedev <N>` | `1454` | Number of Lebedev angular points — **must match** the grid file passed |
| `--n_radial <N>` | `500` | Number of equidistant radial points (≥ 2) |
| `--r_min <X>` | `0.0` | Inner radial bound (a₀, ≥ 0) |
| `--r_max <X>` | `20.0` | Outer radial bound (a₀, > `r_min`) — make it generous for diffuse bases so the CAP tail is captured |
| `--r_cut <X>` | `3.0` | Voronoi isosurface onset (a₀, ≥ 0); `W = 0` inside `r_cut` |
| `--validate` | off | Also report numerical-vs-analytical overlap + symmetry checks (grid quality) |

> `--n_lebedev` must equal the point count of the `--lebedev_grid` file you
> supply. The bundled `lebedev/` folder provides tables from 6 up to 5810
> points; the Sommerfeld–Ehara reference values are 3470 angular × 500 radial.

## Batch / SLURM

`run_CAP.example` is a SLURM submission script showing a full production run
(grid choice, threading, output renaming). Adjust the paths, partition, and the
`R_CUT` / grid parameters for your system, then `sbatch run_CAP.example`.

## Repository contents

```
noCAPs2withCAPs.f03   the generator (single source file)
makefile              build (gfortran + mqcPack + MKL)
lebedev/              Lebedev angular-grid tables (lebedev_<N>.txt)
input.txt.example     example input list
run_CAP.example       example SLURM submission script
```

Data files (`.mat`, `.chk`, `.log`), per-molecule working directories, and the
compiled binary are intentionally **not** tracked — see `.gitignore`.
