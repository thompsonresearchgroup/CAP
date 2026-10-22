#!/usr/bin/env python3
"""End-to-end CAP resonance test against an exact outgoing-wave pole.

The manufactured radial potential is -0.5 for r<2, +0.5 for 2<=r<4,
and zero outside. Gaussian S/T and live executable W are used, but the
Hamiltonian is this model, NOT helium's electronic Hamiltonian.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import root
from scipy.special import gammainc

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from cap_auto import Runner,Grid,parser,sha256_file
from cap_resonance import trajectory,resonance_parameters

BASE=ROOT/'tests/resonance'
OUT=ROOT/'tests/results/resonance'


def outgoing_residual(z):
    energy=complex(*z)
    k=np.sqrt(2*energy);q=np.sqrt(2*(energy+.5));p=np.sqrt(2*(energy-.5))
    # u(0)=0, u'(0)=1; exact propagation through the well and barrier.
    u=np.sin(2*q)/q;d=np.cos(2*q)
    ub=np.cos(2*p)*u+np.sin(2*p)/p*d
    db=-p*np.sin(2*p)*u+np.cos(2*p)*d
    f=db-1j*k*ub
    return [f.real,f.imag]


def reference_pole():
    roots=[]
    for guess in (.05,.1,.2):
        result=root(outgoing_residual,[guess,-.01],tol=1e-11)
        if not result.success or np.linalg.norm(result.fun)>1e-10:
            raise RuntimeError('Exact outgoing matching failed')
        roots.append(complex(*result.x))
    if max(abs(e-roots[0]) for e in roots)>1e-10 or roots[0].imag>=0:
        raise RuntimeError('Inconsistent outgoing poles')
    return roots[0]


def matrices(count,onset,output=OUT):
    source=BASE/f'he_s{count}.mat'
    out=output/f's{count}_r{onset:g}'
    out.mkdir(parents=True,exist_ok=True)
    provenance=dict(input_sha256=sha256_file(source),executable_sha256=sha256_file(ROOT/'noCAPs2withCAPs'),
                    onset=onset,grids=[[1500,302,720.],[1800,590,864.]])
    stamp=out/'integral_provenance.json'
    if stamp.exists():
        if json.loads(stamp.read_text())!=provenance:
            raise ValueError(f'Stale resonance matrices in {out}')
    else:
        if list(out.glob('run_*')):
            raise ValueError(f'Incomplete resonance matrices in {out}')
        args=parser().parse_args([str(source),'--output',str(out),
                                  '--max-radial','2400','--overlap-cutoff','1e-10'])
        args.matrix=source
        runner=Runner(args)
        for nr,na,extent in provenance['grids']:
            runner.run(Grid(nr,na,extent),onset)
        stamp.write_text(json.dumps(provenance,indent=2)+'\n')
    with np.load(out/'run_000.npz') as d:
        s=d['S'];w=d['W']
    with np.load(out/'run_001.npz') as d:
        wf=d['W']
    # Kinetic record was exported from the same Gaussian input by the packet suite.
    kinetic=np.loadtxt(BASE/f's{count}_kinetic.txt',skiprows=1)
    a=np.geomspace(.0001,100,count);p=a[:,None]+a
    expected_s=(2*np.sqrt(a[:,None]*a)/p)**1.5
    expected_t=3*a[:,None]*a/p*expected_s
    if np.max(abs(s-expected_s))>1e-10 or np.max(abs(kinetic-expected_t))>1e-9:
        raise ValueError('Gaussian AO normalization or ordering mismatch')
    # Analytic finite-radius radial Gaussian moments for the physical potential.
    potential=s*(.5*gammainc(1.5,16*p)-gammainc(1.5,4*p))
    h=kinetic+potential
    np.savez_compressed(out/'model_matrices.npz',H=h,S=s,W=wf)
    return h,s,w,wf,out


def run_trajectory(h,s,w,cutoff=1e-10,points=201):
    result=trajectory(h,s,w,np.geomspace(1e-6,.1,points),.12-.01j,cutoff,.001)
    if not result['selected_candidate']:
        raise RuntimeError('No interior stationary point found')
    return result


def difference(first,second):
    a=first['selected_candidate']['raw'];b=second['selected_candidate']['raw']
    return dict(energy_change=abs(a['energy_hartree']-b['energy_hartree']),
                width_change=abs(a['width_hartree']-b['width_hartree']))


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--output',type=Path,default=OUT)
    output=cli.parse_args().output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    exact=reference_pole();reference=resonance_parameters(exact)
    rows=[]
    for count in (48,96):
        for onset in (4.,5.,6.):
            h,s,w,wf,out=matrices(count,onset,output)
            nominal=run_trajectory(h,s,w)
            fine=run_trajectory(h,s,wf)
            fine_eta=run_trajectory(h,s,wf,points=401)
            cutoff=run_trajectory(h,s,wf,cutoff=1e-11)
            for label,result in [('nominal',nominal),('fine_grid',fine),('fine_eta',fine_eta),('cutoff',cutoff)]:
                (out/f'{label}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
            best=fine_eta['selected_candidate']
            raw=best['raw'];corrected=best['corrected']
            checks=dict(grid=difference(nominal,fine),eta=difference(fine,fine_eta),
                        overlap_cutoff=difference(fine,cutoff))
            row=dict(basis_count=count,onset=onset,rank=fine['rank'],selected=best,checks=checks,
                     raw_energy_error=abs(raw['energy_hartree']-exact.real),
                     raw_width_relative_error=abs(raw['width_hartree']/reference['width_hartree']-1),
                     raw_lifetime_relative_error=abs(raw['lifetime_au']/reference['lifetime_au']-1),
                     corrected_energy_error=abs(corrected['energy_hartree']-exact.real),
                     corrected_width_relative_error=abs(corrected['width_hartree']/reference['width_hartree']-1))
            row['passed']=bool(row['raw_energy_error']<=3e-4 and row['raw_width_relative_error']<=.03
                and checks['grid']['energy_change']<=1e-6 and checks['grid']['width_change']<=1e-6
                and checks['eta']['energy_change']<=1e-6 and checks['eta']['width_change']<=1e-6
                and checks['overlap_cutoff']['energy_change']<=3e-4
                and checks['overlap_cutoff']['width_change']<=.03*reference['width_hartree']
                and best['eigenpair_residual']<=1e-7)
            rows.append(row)
            print(row,flush=True)
            (output/'progress.json').write_text(json.dumps(rows,indent=2)+'\n')
    energy_spread=float(np.ptp([r['selected']['raw']['energy_hartree'] for r in rows]))
    width_spread=float(np.ptp([r['selected']['raw']['width_hartree'] for r in rows]))
    passed=all(r['passed'] for r in rows) and energy_spread<=3e-4 and width_spread<=.03*reference['width_hartree']
    report=dict(benchmark_passed=bool(passed),reference=reference,results=rows,
                energy_spread=energy_spread,width_spread=width_spread,
                tolerances=dict(energy_hartree=3e-4,width_relative=.03,grid_and_eta_hartree=1e-6),
                physical_model='Single radial s-wave: V=-0.5 (r<2), +0.5 (2<=r<4), 0 (r>=4)',
                reference_method='Exact piecewise propagation matched to uprime=i*sqrt(2E)*u at r=4',
                selection='Interior minimum of |eta*dE/deta|; no use of exact pole in eta selection',
                molecular_resonance_validated=False,
                provenance=dict(driver_sha256=sha256_file(__file__),
                                trajectory_sha256=sha256_file(ROOT/'cap_resonance.py'),
                                executable_sha256=sha256_file(ROOT/'noCAPs2withCAPs')))
    (output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return 0 if passed else 3


if __name__=='__main__':
    raise SystemExit(main())
