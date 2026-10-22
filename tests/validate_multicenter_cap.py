#!/usr/bin/env python3
"""Independent primitive-pair Gauss–Hermite check of live multicenter W.

The production code uses atom-centered Becke/Lebedev quadrature. This test
instead integrates Gaussian products about each primitive-product center.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
from scipy.special import roots_hermite

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cap_auto import Runner, Grid, parser, sha256_file


def independent_block(atoms, cutoff, smoothing, order):
    nodes, weights = roots_hermite(order)
    weights /= np.sqrt(np.pi)
    alpha = np.array([1.2, .3])
    coeff = np.array([.7, .4])
    primitive_s = (2*np.sqrt(alpha[:, None]*alpha)/(alpha[:, None]+alpha))**1.5
    coeff /= np.sqrt(coeff@primitive_s@coeff)
    block = np.zeros((2, 2))
    for i in range(2):
        for j in range(i, 2):
            for k, ak in enumerate(alpha):
                for l, al in enumerate(alpha):
                    p = ak+al
                    center = (ak*atoms[i]+al*atoms[j])/p
                    overlap = (2*np.sqrt(ak*al)/p)**1.5*np.exp(-ak*al/p*np.sum((atoms[i]-atoms[j])**2))
                    expectation = 0.
                    # Slice the tensor quadrature to keep memory bounded.
                    for node, weight in zip(nodes, weights):
                        gx = center[0]+node/np.sqrt(p)
                        gy = center[1]+nodes[:, None]/np.sqrt(p)
                        gz = center[2]+nodes[None, :]/np.sqrt(p)
                        d2 = np.array([(gx-a[0])**2+(gy-a[1])**2+(gz-a[2])**2 for a in atoms])
                        raw_weights = 1/(d2-d2.min(axis=0)+smoothing**2)**2
                        rwa = np.sqrt((raw_weights*d2).sum(axis=0)/raw_weights.sum(axis=0))
                        potential = np.maximum(rwa-cutoff, 0)**2
                        expectation += weight*float(weights@potential@weights)
                    block[i, j] += coeff[k]*coeff[l]*overlap*expectation
            block[j, i] = block[i, j]
    return block


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output', type=Path, default=ROOT/'tests/results/multicenter')
    out = cli.parse_args().output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    source = ROOT/'tests/gaussian/he2_spdfghi_pure.mat'
    with (out/'coordinates.log').open('w') as log:
        subprocess.run([str(ROOT/'tests/read_cap'), str(source), 'COORDINATES OF EACH SHELL'],
                       cwd=out, stdout=log, stderr=subprocess.STDOUT, check=True)
    coords = np.loadtxt(out/'record.txt', skiprows=1)
    # MQCPack converts a vector record to a diagonal mqc_matrix in read_cap.
    if coords.shape == (42,42) and np.count_nonzero(coords-np.diag(np.diag(coords))) == 0:
        coords = np.diag(coords).reshape(14,3)
    if coords.shape == (3,14):
        atoms = coords[:, [0,7]].T
    elif coords.shape == (14,3):
        atoms = coords[[0,7], :]
    else:
        raise ValueError(f'Unexpected fixture shell coordinates: {coords.shape}')
    args = parser().parse_args([str(source), '--output', str(out)])
    args.matrix = source
    runner = Runner(args)
    coarse = runner.run(Grid(400,1202,24.), 1.5)
    fine = runner.run(Grid(700,3470,28.), 1.5)
    block = fine['w'][np.ix_([0,50],[0,50])]
    rows = []
    for order in (48,80,120,160):
        independent = independent_block(atoms, 1.5, 1., order)
        row = dict(order=order, matrix=independent.tolist(),
                   maximum_error=float(np.max(abs(independent-block))))
        if rows:
            row['refinement_change'] = float(np.max(abs(independent-np.array(rows[-1]['matrix']))))
        rows.append(row)
        print(row, flush=True)
    live_change = float(np.max(abs(fine['w']-coarse['w'])))
    passed = bool(rows[-1]['maximum_error'] < 2e-5 and rows[-1]['refinement_change'] < 2e-5
                  and live_change < 2e-5)
    report = dict(passed=passed, live_s_block=block.tolist(),
                  maximum_live_matrix_refinement_change=live_change, independent=rows,
                  atoms_bohr=atoms.tolist(), onset=1.5, smoothing_length=1.,
                  input_sha256=sha256_file(source), executable_sha256=sha256_file(ROOT/'noCAPs2withCAPs'),
                  description='Two contracted S AOs, including their off-diagonal multicenter CAP integral')
    (out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    return 0 if passed else 2


if __name__ == '__main__':
    raise SystemExit(main())
