#!/usr/bin/env python3
"""Complex-symmetric CAP eta trajectories for resonance energies and lifetimes.

Input NPZ: real symmetric H, S, W in the same basis, with H in hartree and W
the real unscaled CAP. This is a one-electron/CI matrix solver, not an electronic
structure method. A stationary trajectory alone is not a resonance certificate.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import physical_constants
from scipy.linalg import eig
from scipy.optimize import minimize_scalar

from cap_auto import orthogonalizer, sha256_file

HARTREE_EV = physical_constants['Hartree energy in eV'][0]
TIME_FS = physical_constants['atomic unit of time'][0]*1e15


def resonance_parameters(energy):
    width = -2*float(energy.imag)
    return dict(energy_hartree=float(energy.real), imaginary_energy_hartree=float(energy.imag),
                energy_ev=float(energy.real)*HARTREE_EV, width_hartree=width,
                width_ev=width*HARTREE_EV,
                lifetime_au=1/width if width>0 else None,
                lifetime_fs=TIME_FS/width if width>0 else None)


def prepare(h, s, w, cutoff=1e-9):
    for name, a in [('H',h),('S',s),('W',w)]:
        if a.ndim!=2 or a.shape[0]==0 or a.shape[0]!=a.shape[1] or a.shape!=s.shape:
            raise ValueError(f'{name} must have the same square shape as S')
        if np.iscomplexobj(a) or not np.isfinite(a).all():
            raise ValueError(f'{name} must be finite and real')
        if np.max(abs(a-a.T))>1e-10*max(1.,np.linalg.norm(a,2)):
            raise ValueError(f'{name} must be symmetric')
    if not 0<cutoff<1:
        raise ValueError('Relative overlap cutoff must be between zero and one')
    x,rank=orthogonalizer(s,cutoff)
    ho=x.T@h@x;wo=x.T@w@x
    ho=(ho+ho.T)/2;wo=(wo+wo.T)/2
    if np.linalg.eigvalsh(wo).min() < -1e-9*max(1.,np.linalg.norm(wo,2)):
        raise ValueError('CAP must be positive semidefinite')
    return ho,wo,rank


def eigenstate(h,w,eta,target=None,previous=None):
    energies,vectors=eig(h-1j*eta*w)
    if previous is None:
        j=int(np.argmin(abs(energies-target)))
        overlap=1.
    else:
        overlaps=abs(previous.conj()@vectors)
        j=int(np.argmax(overlaps));overlap=float(overlaps[j])
    energy=energies[j];u=vectors[:,j]
    # Complex-symmetric Hellmann–Feynman derivative uses transpose, not adjoint.
    c_norm=u.T@u
    if abs(c_norm)<1e-10:
        raise ValueError('Nearly self-orthogonal state: eta trajectory is ill-conditioned')
    derivative=-1j*(u.T@w@u)/c_norm
    residual=float(np.linalg.norm((h-1j*eta*w)@u-energy*u))
    return energy,derivative,u,overlap,residual


def trajectory(h,s,w,etas,target,cutoff=1e-9,seed_eta=None):
    """Track from a seed eta in both directions using normalized vector overlap.

    Return raw and first-order-corrected E-eta*dE/deta trajectories, with all
    interior raw-speed minima. Endpoints are never accepted as stationary.
    """
    etas=np.asarray(etas,dtype=float)
    if etas.ndim!=1 or etas.size<7 or not np.isfinite(etas).all() or np.min(etas)<=0 or np.any(np.diff(etas)<=0):
        raise ValueError('Need at least seven finite, positive, increasing eta values')
    if seed_eta is not None and (not np.isfinite(seed_eta) or seed_eta<=0):
        raise ValueError('Seed eta must be finite and positive')
    if not np.isfinite([target.real,target.imag]).all():
        raise ValueError('Target energy must be finite')
    ho,wo,rank=prepare(h,s,w,cutoff)
    log_seed=np.log(seed_eta) if seed_eta is not None else (np.log(etas[0])+np.log(etas[-1]))/2
    anchor=int(np.argmin(abs(np.log(etas)-log_seed)))
    states=[None]*len(etas)
    states[anchor]=eigenstate(ho,wo,etas[anchor],target=target)
    for sequence in (range(anchor+1,len(etas)),range(anchor-1,-1,-1)):
        previous=states[anchor][2]
        for i in sequence:
            states[i]=eigenstate(ho,wo,etas[i],previous=previous)
            previous=states[i][2]
    energy=np.array([x[0] for x in states]);derivative=np.array([x[1] for x in states])
    corrected=energy-etas*derivative
    speed=abs(etas*derivative)
    corrected_speed=abs(etas**2*np.gradient(derivative,etas,edge_order=2))
    candidates=[]
    for i in range(2,len(etas)-2):
        if not (speed[i]<speed[i-1] and speed[i]<speed[i+1]):
            continue
        if min(states[j][3] for j in (i-1,i,i+1))<.9:
            continue
        reference=states[i][2]
        def objective(log_eta):
            strength=np.exp(log_eta)
            e,d,u,overlap,residual=eigenstate(ho,wo,strength,previous=reference)
            return abs(strength*d) if overlap>=.9 else np.inf
        optimized=minimize_scalar(objective,bounds=np.log([etas[i-1],etas[i+1]]),
                                  method='bounded',options={'xatol':1e-8})
        strength=float(np.exp(optimized.x))
        e,d,u,overlap,residual=eigenstate(ho,wo,strength,previous=reference)
        if not optimized.success or overlap<.9 or e.imag>=0:
            continue
        # Independently check the analytic derivative on the same branch.
        delta=strength*1e-3
        ep=eigenstate(ho,wo,strength+delta,previous=u)[0]
        em=eigenstate(ho,wo,strength-delta,previous=u)[0]
        derivative_error=float(abs((ep-em)/(2*delta)-d))
        candidates.append(dict(eta=strength,raw=resonance_parameters(e),
                               corrected=resonance_parameters(e-strength*d),
                               logarithmic_speed=float(abs(strength*d)),
                               derivative_finite_difference_error=derivative_error,
                               eigenpair_residual=residual,coarse_index=i))
    candidates.sort(key=lambda item:item['logarithmic_speed'])
    rows=[dict(eta=float(eta),raw=resonance_parameters(e),corrected=resonance_parameters(c),
               raw_speed=float(v),corrected_speed=float(cv),tracking_overlap=float(state[3]),
               eigenpair_residual=float(state[4]))
          for eta,e,c,v,cv,state in zip(etas,energy,corrected,speed,corrected_speed,states)]
    return dict(rank=rank,seed_eta=float(etas[anchor]),target=[target.real,target.imag],
                trajectory=rows,stationary_candidates=candidates,
                selected_candidate=candidates[0] if candidates else None,
                status='requires_onset_grid_and_basis_convergence' if candidates else 'no_interior_stationary_point',
                resonance_validated=False,
                note='Corrected values are evaluated at raw-speed stationary points; '
                     'correction is not assumed to improve every finite-basis result.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('matrices',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--target',type=float,nargs=2,required=True,metavar=('REAL','IMAG'))
    p.add_argument('--eta-min',type=float,default=1e-6)
    p.add_argument('--eta-max',type=float,default=.1)
    p.add_argument('--eta-points',type=int,default=201)
    p.add_argument('--seed-eta',type=float,default=.001)
    p.add_argument('--overlap-cutoff',type=float,default=1e-9)
    # Continuum-remover CAP (Sajeev et al., JCP 131, 211102 (2009)): H + lambda*W - i*eta*W.
    p.add_argument('--real-cap',type=float,default=0.,metavar='LAMBDA',
                   help='eta-independent real CAP strength lambda (default 0)')
    a=p.parse_args()
    if not np.isfinite(a.real_cap) or a.real_cap<0:
        raise ValueError('Real CAP strength must be finite and >= 0')
    if a.output.exists():
        raise ValueError('Use a new output filename')
    if not np.isfinite([a.eta_min,a.eta_max,a.seed_eta]).all() or not 0<a.eta_min<a.eta_max or a.eta_points<7 or a.seed_eta<=0:
        raise ValueError('Invalid eta range, count, or seed')
    with np.load(a.matrices) as data:
        report=trajectory(data['H']+a.real_cap*data['W'],data['S'],data['W'],np.geomspace(a.eta_min,a.eta_max,a.eta_points),
                          complex(*a.target),a.overlap_cutoff,a.seed_eta)
    report['real_cap_strength']=a.real_cap
    report['provenance']=dict(input_sha256=sha256_file(a.matrices),driver_sha256=sha256_file(__file__))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(report['status'])
    if report['selected_candidate']:
        print(json.dumps(report['selected_candidate'],indent=2))
    return 0 if report['stationary_candidates'] else 2


if __name__=='__main__':
    raise SystemExit(main())
