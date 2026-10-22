#!/usr/bin/env python3
"""N2- 2Pi_g shape resonance, CAP/static-exchange (Koopmans) level, Sommerfeld workflow.

H = diag(eps_virt) in the PIG virtual block of neutral RHF N2 (one x/y component),
W = C^T W_ao C. eta = 1e-5..1e-2 in log steps; every trajectory is tracked and
ranked by distance from the real axis, logarithmic speed and orbital compactness <W>.
usage (in a dir holding BASIS.mat/.log from dzp.gjf or big.gjf): n2_anion.py BASIS RCUT [RADIAL ANGULAR]
"""
import re, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.linalg import eig
from scipy.sparse.csgraph import connected_components

CAP = Path(__file__).resolve().parents[2]; EV = 27.211386245988
basis, rcut = sys.argv[1], float(sys.argv[2])
radial, angular = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (500, 1202)


def record(path, name):
    subprocess.run([str(CAP/'tests/read_cap'), str(path), name], check=True, stdout=subprocess.DEVNULL)
    return np.loadtxt('record.txt', skiprows=1)


tag = f'{basis}_r{rcut}_{radial}_{angular}'
if not Path(f'W_{tag}.npy').exists():
    Path('in.list').write_text(f'1\n{Path(basis + ".mat").resolve()}\n')
    subprocess.run([str(CAP/'noCAPs2withCAPs'), '-f', 'in.list', '--lebedev_grid', str(CAP/f'lebedev/lebedev_{angular}.txt'),
                    '--n_lebedev', str(angular), '--n_radial', str(radial), '--auto_rmax', '--r_cut', str(rcut),
                    '--validate', '--output', f'cap_{tag}'], check=True, stdout=open(f'cap_{tag}.log', 'w'))
    np.save(f'W_{tag}.npy', record(f'cap_{tag}.mat', 'Complex Absorbing Potential'))
    Path(f'cap_{tag}.mat').unlink()
w_ao = np.load(f'W_{tag}.npy')
s = record(basis + '.mat', 'OVERLAP')
c = record(basis + '.mat', 'ALPHA MO COEFFICIENTS')
eps = np.diag(record(basis + '.mat', 'ALPHA ORBITAL ENERGIES'))
if np.abs(c.T@s@c - np.eye(c.shape[1])).max() > 1e-8:
    c = c.T
assert np.abs(c.T@s@c - np.eye(c.shape[1])).max() < 1e-8, 'MO orientation'
log = Path(basis + '.log').read_text()
sym = re.findall(r'\((\w+)\)', log.split('Orbital symmetries:')[1].split('The electronic state')[0])
assert len(sym) == len(eps), (len(sym), len(eps))
nocc = 7
pig = [i for i in range(nocc, len(eps)) if sym[i] == 'PIG']
wmo = c[:, pig].T@w_ao@c[:, pig]
h = eps[pig]
# W is D2h-invariant: the degenerate x/y components form two uncoupled identical blocks.
ncomp, comp = connected_components(np.abs(wmo) > 1e-9*np.abs(wmo).max())
blk = np.flatnonzero(comp == comp[0])
other = np.flatnonzero(comp != comp[0])
assert ncomp == 2 and np.allclose(np.sort(h[blk]), np.sort(h[other]), atol=1e-8), 'x/y split'
h, wmo = h[blk], wmo[np.ix_(blk, blk)]
# Input for cap_resonance.py (orthonormal MO basis, so S = I).
np.savez(f'HSW_{tag}.npz', H=np.diag(h), S=np.eye(len(h)), W=wmo)

etas = np.geomspace(1e-5, 1e-2, 151)
traj = []  # traj[k][m] = energy of branch m at eta k, followed by c-product overlap
prev = None
for eta in etas:
    e, v = eig(np.diag(h) - 1j*eta*wmo)
    v /= np.sqrt(np.sum(v*v, axis=0))
    if prev is None:
        order = np.argsort(e.real)
    else:
        ov = np.abs(prev.T@v)
        order = np.argmax(ov, axis=1)
        assert len(set(order)) == len(order), f'branch crossing at eta={eta:g}'
    prev = v[:, order]; traj.append(e[order])
traj = np.array(traj)
speed = np.abs(np.gradient(traj, np.log(etas), axis=0))  # |dE/dln eta|
print(f'basis={basis} rcut={rcut} grid={radial}x{angular} PIG-x block dim={len(h)}')
# Sommerfeld: the resonance stays closest to the real axis and moves slowest.
low = [m for m in range(len(h)) if h[m] < 0.4]
res = min(low, key=lambda m: -traj[:, m].imag.min())
for m in low:
    print(f'  branch {m:2d} eps={h[m]*EV:7.3f} eV  <W>={wmo[m, m]:9.3e}  max Gamma={-2*traj[:, m].imag.min()*EV:8.3f} eV'
          + ('   <- resonance' if m == res else ''))
v = speed[:, res]
mins = [k for k in range(2, len(etas)-2) if v[k] < v[k-1] and v[k] < v[k+1]]
deriv = np.gradient(traj[:, res], etas)
for k in sorted(mins, key=lambda k: v[k]):
    E, U = traj[k, res], traj[k, res] - etas[k]*deriv[k]
    print(f'  eta*={etas[k]:.3e} |dE/dln eta|={v[k]*EV:.3f} eV  raw Er={E.real*EV:.3f} G={-2*E.imag*EV:.3f}'
          f'  corrected Er={U.real*EV:.3f} G={-2*U.imag*EV:.3f} eV')
np.savez(f'traj_{tag}.npz', etas=etas, traj=traj, h=h, wdiag=np.diag(wmo))
