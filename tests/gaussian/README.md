# Gaussian 16 reference jobs

Run serially from this directory with a writable scratch directory. The pure
route uses Gaussian's `5D 7F` convention (pure d and all higher shells); the
Cartesian route uses `6D 10F` (Cartesian d and all higher shells). These route
semantics are verified by the resulting shell-type arrays, not assumed.

```bash
export GAUSS_SCRDIR=/tmp/g16-cap-$USER
mkdir -p "$GAUSS_SCRDIR"
g16 < one_center_spdfghi_pure.gjf > one_center_spdfghi_pure.log
formchk one_center_spdfghi_pure.chk one_center_spdfghi_pure.fchk
g16 < one_center_spdfghi_cartesian.gjf > one_center_spdfghi_cartesian.log
formchk one_center_spdfghi_cartesian.chk one_center_spdfghi_cartesian.fchk
```

Do not run the two Gaussian jobs concurrently. Formatted checkpoints are local
validation artifacts and may be regenerated with the commands above.
Each input also requests a Gaussian unformatted matrix file, which is the
native input consumed by `noCAPs2withCAPs` through MQCPack.

`run_gaussian_references.sh` runs all references serially, including H2 and N2
pure/Cartesian cases and a translated/rotated N2 geometry.

On 2026-10-02, `h2_rotated` and `he2_spdfghi_pure` were generated afresh.
The latter has off-axis centers and all pure shells through I, so cross-center
overlap elements test phases and ordering beyond isolated-shell normalization.
Both are exercised by `tests/test_executable.py`.

The subsequent `he2_spdfghi_cartesian` reference uses Gaussian `guess=only`
to obtain S–I basis/overlap records without storing two-electron integrals.
Its off-axis two-center overlap agrees with the CAP evaluator to 3.32e-13.

## Files included in the repository

Gaussian `.mat` and `.chk` outputs are generated locally and excluded from Git.
The `.gjf` sources and serial generation script are included. Two all-shell
`.fchk` files and five historical CAP overlap logs are retained for reference
regressions. The current executable tests use freshly computed integrals.
The full test procedure is in [validation](../../docs/validation.md).
