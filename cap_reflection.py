#!/usr/bin/env python3
"""One-dimensional quadratic-CAP calibration, atomic units.

This is a one-center normal-incidence model, NOT a molecular propagation.
The right boundary is a hard wall; the scattering left lead is transparent.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_banded
from scipy.sparse import diags
from scipy.sparse.linalg import splu


def reflection(energy, eta, onset, end, dx=0.1):
    if not np.isfinite([energy,eta,onset,end,dx]).all():raise ValueError('Nonfinite scattering parameter')
    if not (0 < energy < 2/dx**2 and eta >= 0 and 0 <= onset < end and dx > 0):
        raise ValueError('Require positive resolvable energy, dx, and 0 <= onset < end')
    n=max(4,int(np.ceil(end/dx))); dx=end/n
    x=np.arange(n)*dx; hopping=0.5/dx**2
    kh=np.arccos(1-energy*dx**2)
    diagonal=(2*hopping-energy-1j*eta*np.maximum(x-onset,0)**2).astype(complex)
    diagonal[0]-=hopping*np.exp(1j*kh)
    band=np.zeros((3,n),complex);band[0,1:]=-hopping;band[2,:-1]=-hopping;band[1]=diagonal
    source=np.zeros(n,complex);source[0]=-2j*hopping*np.sin(kh)
    psi=solve_banded((1,1),band,source,check_finite=False)
    return float(abs(psi[0]-1)**2)


def scan(energies,etas,onsets,end,dx=0.1):
    rows=[]
    for onset in onsets:
        for eta in etas:
            r=[reflection(e,eta,onset,end,dx) for e in energies]
            fine=[reflection(e,eta,onset,end,dx/2) for e in energies]
            rows.append(dict(onset=float(onset),eta=float(eta),reflection=fine,
                             worst_reflection=max(fine),discretization_change=max(abs(np.array(r)-fine))))
    return rows


def wavepacket(energy,eta,onset,end,dx=0.15,dt=0.1):
    """CN propagation versus exact infinite-line free Gaussian in x < onset.

    Duration includes a nominal wall round trip. Norm/error use dx measure.
    Spatial/temporal refinement is reported by the caller, never presumed.
    """
    if not np.isfinite([energy,eta,onset,end,dx,dt]).all():raise ValueError('Nonfinite propagation parameter')
    if min(energy,dx,dt) <= 0 or eta < 0 or not 0 <= onset < end:
        raise ValueError('Invalid propagation parameters')
    k=np.sqrt(2*energy);sigma=max(5.0,3/k);x0=-6*sigma
    left=x0-8*sigma;duration=2*(end-x0)/k+4*sigma/k
    n=int(np.ceil((end-left)/dx));dx=(end-left)/n
    x=left+np.arange(1,n)*dx;steps=int(np.ceil(duration/dt));dt=duration/steps
    potential=-1j*eta*np.maximum(x-onset,0)**2
    h=diags([-np.ones(n-2)/(2*dx*dx),1/dx**2+potential,-np.ones(n-2)/(2*dx*dx)],[-1,0,1],format='csc')
    identity=diags(np.ones(n-1),format='csc')
    lhs=splu(identity+0.5j*dt*h);rhs=identity-0.5j*dt*h
    norm=(2*np.pi*sigma*sigma)**(-0.25)
    psi=norm*np.exp(-(x-x0)**2/(4*sigma*sigma)+1j*k*(x-x0))
    initial_norm=float(np.sum(abs(psi)**2)*dx);inner=x<onset
    largest=0.0;max_norm_increase=0.0;oldnorm=initial_norm
    for step in range(1,steps+1):
        psi=lhs.solve(rhs@psi)
        newnorm=float(np.sum(abs(psi)**2)*dx)
        max_norm_increase=max(max_norm_increase,newnorm-oldnorm);oldnorm=newnorm
        if step%max(1,steps//150)==0 or step==steps:
            t=step*dt;z=1+1j*t/(2*sigma*sigma)
            exact=norm/np.sqrt(z)*np.exp(-(x-x0-k*t)**2/(4*sigma*sigma*z)+1j*k*(x-x0)-0.5j*k*k*t)
            error=float(np.sum(abs(psi[inner]-exact[inner])**2)*dx)
            largest=max(largest,error)
    return dict(energy=energy,eta=eta,onset=onset,end=end,dx=dx,dt=dt,duration=duration,
                initial_norm=initial_norm,final_norm=oldnorm,absorbed=initial_norm-oldnorm,
                final_inner_error=error,max_inner_error=largest,max_norm_increase=max_norm_increase)


def converge_wavepacket(energy,eta,onset,end,tolerance=1e-5,max_levels=4):
    dx=.15;dt=.1;previous=wavepacket(energy,eta,onset,end,dx,dt)
    history=[];passes=0
    for level in range(max_levels):
        time=wavepacket(energy,eta,onset,end,dx,dt/2)
        joint=wavepacket(energy,eta,onset,end,dx/2,dt/2)
        metrics=('max_inner_error','final_inner_error','final_norm')
        time_change=max(abs(time[k]-previous[k]) for k in metrics)
        space_change=max(abs(joint[k]-time[k]) for k in metrics)
        passes=passes+1 if max(time_change,space_change)<=tolerance else 0
        history.append(dict(level=level,time_change=time_change,space_change=space_change,
                            consecutive_passes=passes,time_refined=time,joint_refined=joint))
        if passes>=2:return dict(converged=True,tolerance=tolerance,history=history,result=joint)
        previous=joint;dx/=2;dt/=2
    return dict(converged=False,tolerance=tolerance,history=history,result=previous,
                stopping_reason='Propagation refinement limit reached')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--energies',type=float,nargs='+',default=[.05,.2,.5])
    p.add_argument('--etas',type=float,nargs='+',default=[1e-5,3e-5,1e-4,3e-4,.001,.003,.01])
    p.add_argument('--onsets',type=float,nargs='+',default=[3,5,7])
    p.add_argument('--end',type=float,default=30)
    p.add_argument('--propagate',action='store_true')
    p.add_argument('--propagation-tol',type=float,default=1e-5)
    p.add_argument('--max-propagation-levels',type=int,default=4)
    a=p.parse_args()
    if not np.isfinite(a.propagation_tol) or a.propagation_tol<=0 or a.max_propagation_levels<1:p.error('Positive propagation tolerance and limits required')
    if a.output.exists() and any(a.output.iterdir()):p.error('Output directory must be new or empty')
    a.output.mkdir(parents=True,exist_ok=True)
    rows=scan(a.energies,a.etas,a.onsets,a.end);best=min(rows,key=lambda x:x['worst_reflection'])
    result=dict(model='1D one-center quadratic CAP; not a molecular certificate',energies=a.energies,
                end=a.end,scan=rows,best=best,no_cap=[reflection(e,0,best['onset'],a.end) for e in a.energies])
    if a.propagate:
        result['propagation']=[]
        for e in a.energies:
            params=dict(energy=e,eta=best['eta'],onset=best['onset'],end=a.end)
            result['propagation'].append(dict(adaptive=converge_wavepacket(**params,
                tolerance=a.propagation_tol,max_levels=a.max_propagation_levels),
                no_cap=wavepacket(**dict(params,eta=0))))
            print('Propagated',e,flush=True)
    (a.output/'reflection.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(best,indent=2))
    if a.propagate and not all(x['adaptive']['converged'] for x in result['propagation']):
        raise SystemExit(2)

if __name__=='__main__':main()
