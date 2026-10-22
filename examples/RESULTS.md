# CAP example results

The automatic N2 calculation completes with exit 0 and writes exactly
[CAP_n2_pure_rcut_9.75.mat](results/CAP_n2_pure_rcut_9.75.mat) and [CAP_n2_pure_rcut_9.75.log](results/CAP_n2_pure_rcut_9.75.log).

| Setting | Selected value |
|---|---:|
| AO basis size | 36 |
| r_min | 0 bohr |
| r_cut | 9.75 bohr |
| r_max | 32.04015578 bohr |
| Radial points | 450 |
| Angular points | 3470 |

The log contains the overlap and CAP matrices, traces,
and convergence summary. The matrix is real and unscaled; eta is applied
by the downstream Hamiltonian calculation.

Two consecutive independent plus joint grid-refinement checks pass.
The log states `Integral convergence only; eta not selected.` These are
converged CAP integrals, not a molecular resonance or propagation result.
