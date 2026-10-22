# Radial resonance fixtures

`he_s48.gjf` and `he_s96.gjf` generate Gaussian basis metadata and overlap/kinetic
integrals for the manufactured radial benchmark. They do not represent a
physical helium resonance calculation.

`s48_kinetic.txt` and `s96_kinetic.txt` are exported kinetic matrices with a
matrix-dimensions header. `validate_resonances.py` verifies them against analytic
Gaussian integrals before constructing the well/barrier Hamiltonian.

The corresponding large `.mat` files are local generated fixtures, excluded
from Git. See [benchmark reproduction](../../docs/validation.md).
