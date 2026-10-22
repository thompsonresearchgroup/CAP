#!/usr/bin/env python3
"""Brute-force Cartesian n^3 check of <s|W|s>, <s|W|p>, <p|W|d>, ... (Sommerfeld's test).

Fixture tests/gaussian/n2_spd_bruteforce.mat: N2 at z=+-1 bohr, one uncontracted
s(0.5), p(0.4), pure d(0.3) per atom. AOs are rebuilt here independently in
Gaussian order (s; px py pz; d0 d+1 d-1 d+2 d-2) and checked against Gaussian's
analytic overlap before any W comparison. The trapezoid grid is refined in n and
box edge; W from the live Becke grid must agree to the brute-force convergence.
"""
import argparse
import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT/'tests/gaussian/n2_spd_bruteforce.mat'
ATOMS = np.array([[0., 0., -1.], [0., 0., 1.]])
RCUT, SMOOTH = 1.0, 1.0


def aos(x, y, z):
    out = []
    for c in ATOMS:
        dx, dy, dz = x-c[0], y-c[1], z-c[2]
        r2 = dx*dx+dy*dy+dz*dz
        es, ep, ed = np.exp(-.5*r2), np.exp(-.4*r2), np.exp(-.3*r2)
        out += [es, dx*ep, dy*ep, dz*ep,
                (2*dz*dz-dx*dx-dy*dy)*ed, dx*dz*ed, dy*dz*ed, (dx*dx-dy*dy)*ed, dx*dy*ed]
    return out


def potential(x, y, z):
    d2 = np.array([(x-a[0])**2+(y-a[1])**2+(z-a[2])**2 for a in ATOMS])
    w = 1/(d2-d2.min(axis=0)+SMOOTH**2)**2
    rwa = np.sqrt((w*d2).sum(axis=0)/w.sum(axis=0))
    return np.maximum(rwa-RCUT, 0)**2


def brute_force(n, edge):
    """Trapezoid (= midpoint-free uniform) rule on [-edge, edge]^3, summed in x-slabs."""
    g = np.linspace(-edge, edge, n)
    h = g[1]-g[0]
    y, z = np.meshgrid(g, g, indexing='ij')
    s = np.zeros((18, 18)); w = np.zeros((18, 18))
    for x in g:  # integrands vanish at the box faces, so end weights are irrelevant
        phi = np.array(aos(x, y, z)).reshape(18, -1)
        v = potential(x, y, z).ravel()
        s += phi@phi.T; w += (phi*v)@phi.T
    s *= h**3; w *= h**3
    norm = 1/np.sqrt(np.diag(s))
    return s*np.outer(norm, norm), w*np.outer(norm, norm)


def read_record(path, name, workdir):
    subprocess.run([str(ROOT/'tests/read_cap'), str(path), name], cwd=workdir,
                   check=True, stdout=subprocess.DEVNULL)
    return np.loadtxt(Path(workdir)/'record.txt', skiprows=1)


def live_w(radial, angular, extent, workdir):
    (Path(workdir)/'in.list').write_text(f'1\n{FIXTURE}\n')
    subprocess.run([str(ROOT/'noCAPs2withCAPs'), '-f', 'in.list', '--lebedev_grid',
                    str(ROOT/f'lebedev/lebedev_{angular}.txt'), '--n_lebedev', str(angular),
                    '--n_radial', str(radial), '--r_max', str(extent), '--r_cut', str(RCUT),
                    '--output', 'live'], cwd=workdir, check=True, stdout=subprocess.DEVNULL)
    return read_record(Path(workdir)/'live.mat', 'Complex Absorbing Potential', workdir)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', type=Path, default=ROOT/'tests/results/cartesian/report.json')
    out = cli.parse_args().output
    with tempfile.TemporaryDirectory() as td:
        s_ref = read_record(FIXTURE, 'OVERLAP', td)
        w_live = live_w(1000, 5810, 20., td)
    rows = []
    for n, edge in [(241, 10.), (361, 10.), (481, 10.), (481, 12.), (601, 12.)]:
        s, w = brute_force(n, edge)
        rows.append(dict(n=n, edge=edge, overlap_error=float(abs(s-s_ref).max()), w=w))
        print(n, edge, rows[-1]['overlap_error'], flush=True)
    w_bf = rows[-1]['w']
    scale = abs(w_bf).max()
    refinement = float(max(abs(rows[-1]['w']-rows[-2]['w']).max(), abs(rows[-2]['w']-rows[-3]['w']).max())/scale)
    error = float(abs(w_live-w_bf).max()/scale)
    # Labelled elements Sommerfeld suggested, atom A = first N, B = second N.
    labels = {'<sA|W|sA>': (0, 0), '<sA|W|sB>': (0, 9), '<sA|W|pzA>': (0, 3), '<sA|W|pzB>': (0, 12),
              '<pxA|W|pxB>': (1, 10), '<pzA|W|d0B>': (3, 13), '<d+1A|W|d+1B>': (5, 14), '<d0A|W|d0A>': (4, 4)}
    elements = {k: dict(brute_force=float(w_bf[i, j]), live=float(w_live[i, j]),
                        relative=float(abs(w_live[i, j]-w_bf[i, j])/abs(w_bf[i, j]))) for k, (i, j) in labels.items()}
    for k, v in elements.items():
        print(f"{k:16s} {v['brute_force']: .10e} {v['live']: .10e} {v['relative']:.1e}")
    passed = bool(rows[-1]['overlap_error'] < 1e-10 and refinement < 1e-6 and error < 1e-6)
    report = dict(passed=passed, max_relative_live_error=error, brute_force_refinement=refinement,
                  elements=elements, grids=[{k: v for k, v in r.items() if k != 'w'} for r in rows])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2)+'\n')
    print('relative max |W_live - W_bf| =', error, ' brute-force refinement =', refinement, ' passed =', passed)
    return 0 if passed else 2


if __name__ == '__main__':
    raise SystemExit(main())
