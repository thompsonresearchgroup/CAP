# Repository contents and publishing

The maintained project contains:

```text
README.md                    build and quick start
makefile                     optimized/checked builds and tests
requirements.txt             Python dependencies
noCAPs2withCAPs.f03           CAP integral generator
cap_auto.py                  automatic parameters and grid convergence
cap_resonance.py              resonance trajectory analysis
cap_reflection.py            optional model-reflection screening
cap_wavepacket.py            independent radial propagation utilities
lebedev/                     required angular quadrature tables
docs/                        user, automation, resonance and validation guides
examples/                    small runnable inputs and two-file CAP output
tests/                       regression tests and independent benchmarks
```

`.gitignore` excludes executables, caches, newly generated matrices/logs,
large local Gaussian fixtures, calculation directories under `runs/`, and
`.local/`. Explicit exceptions retain the small runnable N2 example, saved
model H/S/W archive, example output, and seven small regression references.
Generate the large test fixtures on a new machine as described in
[validation](validation.md). They are not required for `make test-python`
or the bundled N2 calculation.

Historical scan directories, old reports, cluster-specific scripts, checkpoints,
backup source files, journal PDFs and obsolete diagnostics were moved into:

```text
.local/archive-2026-10-06/
```

`MANIFEST.json` in that directory records each old and new path, including the
maintained files moved into `docs/` and `tests/resonance/`. No historical data
was deleted. The archive is local only and is excluded from Git. Large matrix
fixtures still needed by the tests remain in their fixture directories but
are also excluded from Git. Build products may remain locally for convenience.

For GitHub, use Git's normal add/status workflow so `.gitignore` is respected;
do not upload the entire directory through a file-manager drag-and-drop.
This cleanup does not create a remote repository, push files, or assign a new
software license. No license grant has been added; the authors can supply their
chosen license before public distribution.

A clean source archive is available locally as `.local/CAP-source.tar.gz`.
It contains the files allowed by `.gitignore`, without the local archive,
executables, or large generated fixtures. `.local/source-files.txt` lists its
contents. The package was checked by building and running the bundled example
from a separate clean export.

The scientific references are linked by DOI. Copies of publisher PDFs and
third-party MQCPack/MKL/Gaussian installations are not part of the source package.
