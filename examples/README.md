# CAP input/output example

This example runs automatic CAP parameter selection for N2 and produces only:

- [results/CAP_n2_pure_rcut_9.75.mat](results/CAP_n2_pure_rcut_9.75.mat): the real, unscaled AO CAP matrix
  and the existing overlap/APO records in one binary matrix file.
- [results/CAP_n2_pure_rcut_9.75.log](results/CAP_n2_pure_rcut_9.75.log): the calculation log, including geometry,
  grid settings, analytical/numerical overlap matrices, CAP matrix,
  traces, and automatic convergence results.

The log uses the traditional matrix/trace layout with current Becke integration
and a concise convergence summary.

## Input

[inputs/n2_pure.mat](inputs/n2_pure.mat) contains the Gaussian geometry, basis,
overlap, and density records. The readable Gaussian input is
[inputs/n2_pure.gjf](inputs/n2_pure.gjf):

```text
%chk=n2_pure.chk
#p rhf/6-31+g(d) 5d 7f scf=conventional output=matrix int=noraff nosymm

N2 diffuse pure reference

0 1
N 0.0 0.0 -0.55
N 0.0 0.0  0.55

n2_pure.mat

```

The Gaussian geometry uses angstroms; CAP distances are in bohr.
Gaussian is not needed to run the supplied example because the `.mat` input
is included. To regenerate it, run `g16 < n2_pure.gjf > n2_pure.log` in a
separate directory containing a copy of the Gaussian input.

## Run

From the CAP directory:

```bash
make
bash examples/run.sh rerun
```

This creates exactly `examples/rerun/CAP_n2_pure_rcut_9.75.mat` and
`examples/rerun/CAP_n2_pure_rcut_9.75.log`. Choose a fresh directory name for another run.
The supplied executed example is in `results/`.

The equivalent direct command is:

```bash
python3 cap_auto.py examples/inputs/n2_pure.mat --output examples/another_run
```

Automatic onset/extent selection and grid convergence are still enabled.
Intermediate matrices and logs are temporary and removed after the run.
A failed convergence run leaves only `CAP_n2_pure_failed.log`, without publishing a final
matrix. Use `--diagnostics` explicitly if you need all intermediate files
and the JSON report for development or a convergence study.

Apply the matrix downstream as `H_eff = H - 1j*eta*W`, with H and W in the
same basis. This example converges the CAP integrals; it does not select eta
or predict an N2 resonance. See [RESULTS.md](RESULTS.md) for selected parameters.

Resonance analysis is a separate operation described in the
[resonance guide](../docs/resonances.md); it is not run by this CAP example script.
