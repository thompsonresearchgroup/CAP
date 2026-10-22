#!/usr/bin/env python3
"""Density-initialized, error-controlled Smooth Voronoi CAP workflow.

All lengths/energies/times are atomic units. Model reflection calibration is
explicitly separate from molecular quadrature convergence.
"""
import argparse
from contextlib import redirect_stdout
from dataclasses import dataclass, asdict
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import numpy as np
from cap_reflection import scan, reflection

ROOT=Path(__file__).resolve().parent
ANGULAR=[302,590,1202,3470,5810]


class NotConverged(RuntimeError): pass


@dataclass(frozen=True)
class Grid:
    radial:int
    angular:int
    extent:float|None


def orthogonalizer(s,cutoff):
    if not np.isfinite(s).all() or np.max(abs(s-s.T))>1e-10:
        raise ValueError('Invalid analytical overlap')
    eig,u=np.linalg.eigh(s);keep=eig>cutoff*eig[-1]
    if eig[0]<-1e-10 or not keep.any(): raise ValueError('Invalid overlap spectrum')
    return u[:,keep]/np.sqrt(eig[keep]),dict(total=len(eig),retained=int(keep.sum()),
        discarded=int((~keep).sum()),relative_cutoff=cutoff,min_eigenvalue=float(eig[0]))


def opnorm(a): return float(np.linalg.norm(a,2))


def density_cutoff(r,mass,tail_budget,step=.25):
    if tail_budget<=0 or step<=0 or not np.isfinite(r).all() or not np.isfinite(mass).all() or np.min(mass)<0:
        raise ValueError('Invalid density samples or tail budget')
    order=np.argsort(r);r=r[order];mass=mass[order]
    tails=np.cumsum(mass[::-1])[::-1]
    indices=np.flatnonzero(tails<=tail_budget)
    quantile=float(r[indices[0]]) if len(indices) else float(r[-1])
    # Round upward and add one step to protect against sample/quantile movement.
    return (math.ceil(quantile/step)+1)*step


class Runner:
    def __init__(self,args):
        self.a=args;self.cache={};self.runs=[];self.x=None;self.rank=None
    def run(self,grid,rcut):
        key=(grid,round(rcut,12))
        if key in self.cache:return self.cache[key]
        if len(self.runs)>=self.a.max_runs:raise NotConverged('Calculation budget reached')
        if grid.radial>self.a.max_radial:raise NotConverged('Radial point limit reached')
        index=len(self.runs); logpath=self.a.output/f'run_{index:03d}.log'
        with tempfile.TemporaryDirectory(prefix='cap-auto-') as td:
            d=Path(td);(d/'input.list').write_text(f'1\n{self.a.matrix}\n')
            command=[str(ROOT/'noCAPs2withCAPs'),'-f','input.list','--becke_grid',
                '--lebedev_grid',str(ROOT/f'lebedev/lebedev_{grid.angular}.txt'),
                '--n_lebedev',str(grid.angular),'--n_radial',str(grid.radial),
                '--r_cut',str(rcut),'--r_min','0','--automation_data']
            command+=['--auto_rmax'] if grid.extent is None else ['--r_max',str(grid.extent)]
            with logpath.open('w') as log:
                subprocess.run(command,cwd=d,env=dict(os.environ,OMP_NUM_THREADS=str(self.a.threads),MKL_NUM_THREADS='1'),
                    stdout=log,stderr=subprocess.STDOUT,check=True,timeout=self.a.run_timeout)
            if grid.extent is None:
                match=re.search(r'Automatic initial r_max \(bohr\) =\s*(\S+)',logpath.read_text())
                grid=Grid(grid.radial,grid.angular,float(match[1]))
            with (d/'cap_data.bin').open('rb') as f:
                nb,ne=np.fromfile(f,dtype=np.int64,count=2)
                arrays=np.fromfile(f,dtype=np.float64)
            if arrays.size!=4*nb*nb:raise ValueError('Invalid automation matrix payload')
            s,p,snum,w=[a.reshape((nb,nb),order='F') for a in np.split(arrays,4)]
            if not all(np.isfinite(a).all() for a in (s,p,snum,w)):raise ValueError('Nonfinite matrices')
            if self.x is None:
                self.x,self.rank=orthogonalizer(s,self.a.overlap_cutoff);self.s=s.copy()
            if np.max(abs(self.s-s))>1e-12:raise ValueError('Analytical overlap changed between runs')
            x=self.x;wo=x.T@w@x
            if np.linalg.eigvalsh(wo).min() < -1e-10*max(1,opnorm(wo)):raise ValueError('Nonpositive CAP')
            if np.linalg.eigvalsh(p).min() < -1e-8:raise ValueError('Nonpositive initial density')
            profile=np.fromfile(d/'density_samples.bin',dtype=np.float64).reshape((-1,2))
            r,mass=profile.T
            if not np.isfinite(profile).all() or np.min(mass)<0:raise ValueError('Invalid density samples')
            electron_error=abs(float(mass.sum())-ne)
            leakage=float(np.sum(p*w.T))
            direct=float(mass@np.maximum(r-rcut,0)**2)
            if abs(direct-leakage)>1e-8*max(1,abs(leakage)):raise ValueError('Density/CAP contraction mismatch')
            onset=density_cutoff(r,mass,self.a.density_tail)
            candidates=[dict(onset=rcut+offset,tail=float(mass[r>rcut+offset].sum()),
                             expectation=float(mass@np.maximum(r-rcut-offset,0)**2)) for offset in (0,.5,1)]
            record=dict(index=index,grid=asdict(grid),rcut=rcut,
                overlap_error=opnorm(x.T@(snum-s)@x),electron_error=electron_error,
                electron_count=float(mass.sum()),expected_electrons=int(ne),density_cutoff=onset,
                tail=float(mass[r>rcut].sum()),density_expectation=leakage,onset_candidates=candidates,
                minimum_cap_eigenvalue=float(np.linalg.eigvalsh(wo).min()))
            artifact=self.a.output/f'run_{index:03d}.mat'
            shutil.copyfile(d/'CAP_Mat.mat',artifact)
            np.savez_compressed(self.a.output/f'run_{index:03d}.npz',S=s,P=p,Snum=snum,W=w)
        self.runs.append(record)
        result=dict(grid=grid,w=w,wo=wo,record=record,artifact=artifact)
        self.cache[key]=result;self.cache[(grid,round(rcut,12))]=result
        print(f'run {index}: rcut={rcut:g}, grid={grid}, overlap={record["overlap_error"]:.2e}',flush=True)
        return result


def converge(runner,grid,rcut,args,history):
    base=runner.run(grid,rcut);passes=0
    while True:
        g=base['grid'];idx=ANGULAR.index(g.angular)
        if idx==len(ANGULAR)-1:raise NotConverged('Angular rule limit prevents another independent refinement')
        nr=int(math.ceil(g.radial*1.5));extent=g.extent*1.2
        probes={
            'radial':Grid(nr,g.angular,g.extent),
            'angular':Grid(g.radial,ANGULAR[idx+1],g.extent),
            'extent':Grid(int(math.ceil(g.radial*1.2)),g.angular,extent)}
        results={axis:runner.run(test,rcut) for axis,test in probes.items()}
        checks={axis:dict(cap_change=opnorm(res['wo']-base['wo']),
             overlap_error=res['record']['overlap_error'],electron_error=res['record']['electron_error'],
             density_tail=res['record']['tail'],density_tail_change=abs(res['record']['tail']-base['record']['tail'])) for axis,res in results.items()}
        def acceptable(res):
            rec=res['record']
            return rec['overlap_error']<=args.overlap_tol and rec['electron_error']<=args.electron_tol and rec['tail']<=args.density_tail
        if any(res['record']['tail']>args.density_tail for res in results.values()):
            new_cut=max([rcut+.25]+[res['record']['density_cutoff'] for res in results.values()])
            history.append(dict(onset_restart_from=rcut,onset_restart_to=new_cut))
            rcut=new_cut;passes=0;base=runner.run(g,rcut)
            continue
        failed=[axis for axis,res in results.items() if checks[axis]['cap_change']>args.cap_tol or checks[axis]['density_tail_change']>args.density_tail*.1 or not acceptable(res)]
        entry=dict(base=base['record']['index'],checks=checks,failed=failed)
        history.append(entry)
        if failed:
            passes=0
            # Advance the failing directions. If base itself is deficient, joint refinement resolves it.
            grid=Grid(nr if 'radial' in failed else probes['extent'].radial if 'extent' in failed else g.radial,
                ANGULAR[idx+1] if 'angular' in failed else g.angular,
                extent if 'extent' in failed else g.extent)
            base=runner.run(grid,rcut)
            continue
        joint=runner.run(Grid(nr,ANGULAR[idx+1],extent),rcut)
        change=opnorm(joint['wo']-base['wo']);entry['joint_change']=change;entry['joint_run']=joint['record']['index']
        if (change<=args.cap_tol and acceptable(base) and acceptable(joint) and
            abs(joint['record']['tail']-base['record']['tail'])<=args.density_tail*.1):passes+=1
        else:passes=0
        entry['consecutive_passes']=passes
        if passes>=2:return joint
        base=joint


def absorption_capacity(w_orth,eta,duration,target):
    """Necessary finite-basis condition, not a sufficient absorption test.

    For Hermitian H and PSD W, N(t)>=exp(-2 eta lambda_max(W) t) N(0).
    """
    maximum=max(0.0,float(np.linalg.eigvalsh(w_orth).max()))
    upper=float(-np.expm1(-2*eta*duration*maximum))
    return dict(maximum_cap_eigenvalue=maximum,maximum_possible_absorbed_fraction=upper,
                required_fraction=target,necessary_condition_passed=upper>=target,
                minimum_eta_for_target=float(-np.log1p(-target)/(2*duration*maximum)) if maximum>0 else None)


def select_calibrated_candidate(admissible, final, runner, initial_grid, args, history, energies):
    """Validate candidates before selection, including their exact onset matrix.

    W at the density-selected onset only bounds W at later onsets from above.
    An optimistic capacity bound must not prevent trying a genuinely capable
    alternative. Likewise a coarse reflection winner may fail the fine check.
    """
    fine_energies=np.unique(np.r_[energies,np.linspace(min(energies),max(energies),61)]).tolist()
    exact={final['record']['rcut']:final}
    for row in admissible:
        row['basis_capacity_upper_bound']=absorption_capacity(
            final['wo'],row['eta'],args.duration,args.required_absorption)
    candidates=sorted(admissible,key=lambda row:(
        not row['basis_capacity_upper_bound']['necessary_condition_passed'],row['worst_reflection']))
    fallback=None
    for row in candidates:
        values=[reflection(e,row['eta'],row['onset'],args.model_end,.025) for e in fine_energies]
        validation=dict(energies=fine_energies,reflection=values,dx=.025,worst_reflection=max(values))
        row['refined_worst_reflection']=max(values)
        if not np.isfinite(values).all() or max(values)>args.reflection_tol:
            row['selection_rejection']='refined_model_reflection'
            continue
        onset=row['onset']
        if onset not in exact:
            exact[onset]=converge(runner,initial_grid,onset,args,history)
        candidate=exact[onset]
        loss=2*row['eta']*args.duration*candidate['record']['density_expectation']
        if candidate['record']['rcut']!=onset or candidate['record']['tail']>args.density_tail:
            row['selection_rejection']='refined_density_onset'
            continue
        if loss>args.loss_budget:
            row['selection_rejection']='refined_initial_loss'
            continue
        capacity=absorption_capacity(candidate['wo'],row['eta'],args.duration,args.required_absorption)
        row['exact_basis_capacity']=capacity
        selected=(row,candidate,validation,loss)
        if capacity['necessary_condition_passed']:
            return selected
        row['selection_rejection']='insufficient_exact_basis_capacity'
        if fallback is None or row['worst_reflection']<fallback[0]['worst_reflection']:
            fallback=selected
        # Once no remaining candidate can pass even the upper bound, only the
        # best refined reflection result is needed for an explicit rejection.
        if not row['basis_capacity_upper_bound']['necessary_condition_passed']:
            break
    if fallback is not None:
        return fallback
    raise NotConverged('No candidate passes refined reflection and final-grid density checks')


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('matrix',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--diagnostics',action='store_true',help='retain per-run matrices, logs and JSON diagnostics')
    p.add_argument('--density-tail',type=float,default=1e-6,help='allowed electrons outside onset')
    p.add_argument('--cap-tol',type=float,default=1e-4,help='absolute spectral norm tolerance in retained orthonormal basis')
    p.add_argument('--overlap-tol',type=float,default=1e-6)
    p.add_argument('--electron-tol',type=float,default=1e-6)
    p.add_argument('--overlap-cutoff',type=float,default=1e-9)
    p.add_argument('--radial',type=int,default=200);p.add_argument('--angular',type=int,choices=ANGULAR,default=590)
    p.add_argument('--max-radial',type=int,default=1600);p.add_argument('--max-runs',type=int,default=50)
    p.add_argument('--threads',type=int,default=2);p.add_argument('--run-timeout',type=int,default=300)
    p.add_argument('--energies',type=float,nargs='+',help='optional 1D model calibration energies, hartree')
    p.add_argument('--etas',type=float,nargs='+',default=np.geomspace(1e-5,.01,25).tolist())
    p.add_argument('--model-end',type=float,default=30,help='explicit 1D model wall position, NOT integration r_max')
    p.add_argument('--duration',type=float,default=100,help='time for initial-density loss screening, atomic units')
    p.add_argument('--loss-budget',type=float,default=1e-4)
    p.add_argument('--required-absorption',type=float,default=.99,help='necessary finite-basis capacity target over duration')
    p.add_argument('--reflection-tol',type=float,default=1e-3)
    return p


def sha256_file(path):
    digest=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def output_stem(matrix,cutoff):
    radius=f'{cutoff:.6f}'.rstrip('0').rstrip('.')
    return f'CAP_{Path(matrix).stem}_rcut_{radius}'


def run_workflow(args):
    args.matrix=args.matrix.resolve();args.output=args.output.resolve()
    if not 0<args.required_absorption<1:raise ValueError('Required absorption must lie between zero and one')
    if not args.matrix.is_file():raise ValueError('Matrix file does not exist')
    positive=[args.density_tail,args.cap_tol,args.overlap_tol,args.electron_tol,args.overlap_cutoff,
              args.radial,args.max_radial,args.max_runs,args.threads,args.run_timeout,args.duration,args.loss_budget,
              args.reflection_tol,args.model_end,*args.etas,*(args.energies or [1])]
    if not all(np.isfinite(positive)) or min(positive)<=0 or args.radial<2:raise ValueError('All tolerances/limits must be finite and positive')
    if args.output.exists() and any(args.output.iterdir()):raise ValueError('Output directory must be new or empty')
    args.electron_tol=min(args.electron_tol,args.density_tail*.1)
    args.output.mkdir(parents=True,exist_ok=True)
    runner=Runner(args);history=[];report=dict(converged=False,quadrature_converged=False,physical_calibration='not_requested',
        settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},history=history,runs=runner.runs)
    report['provenance']=dict(input_sha256=sha256_file(args.matrix),
        executable_sha256=sha256_file(ROOT/'noCAPs2withCAPs'),driver_sha256=sha256_file(__file__),
        reflection_sha256=sha256_file(ROOT/'cap_reflection.py'),numpy_version=np.__version__)
    try:
        pilot=runner.run(Grid(args.radial,args.angular,None),3.0)
        if pilot['record']['electron_error']>max(.01,100*args.electron_tol):
            raise NotConverged('Pilot density integral is insufficiently resolved; increase initial grid')
        rcut=pilot['record']['density_cutoff']
        final=converge(runner,pilot['grid'],rcut,args,history)
        report['quadrature_converged']=True
        rcut=final['record']['rcut']
        report['orthogonalization']=runner.rank
        if args.energies:
            candidates=final['record']['onset_candidates']
            onsets=[x['onset'] for x in candidates if x['onset']<args.model_end]
            if not onsets:raise NotConverged('Density-safe onset is outside specified model domain')
            energies=np.unique(np.r_[args.energies,np.linspace(min(args.energies),max(args.energies),31)]).tolist()
            report['model_energies']=energies
            rows=scan(energies,args.etas,onsets,args.model_end)
            for row in rows:
                candidate=next(x for x in candidates if x['onset']==row['onset'])
                row['initial_loss_estimate']=2*row['eta']*args.duration*candidate['expectation']
                row['admissible']=bool(row['initial_loss_estimate']<=args.loss_budget and
                    candidate['tail']<=args.density_tail and row['worst_reflection']<=args.reflection_tol and
                    row['discretization_change']<=args.reflection_tol*.1)
            report['model_scan']=rows
            admissible=[row for row in rows if row['admissible']]
            if not admissible:
                report['physical_calibration']='failed'
                raise NotConverged('No eta/onset meets model reflection and initial-loss tolerances')
            try:
                best,final,validation,checked_loss=select_calibrated_candidate(
                    admissible,final,runner,pilot['grid'],args,history,energies)
            except NotConverged:
                report['physical_calibration']='failed'
                raise
            rcut=final['record']['rcut']
            report['model_selection']=best
            report['model_validation']=validation
            report['physical_calibration']='1D_model_only'
            report['final_initial_loss_estimate']=checked_loss
        report['molecular_dynamics_validated']=False
        if args.energies:
            report['basis_capacity']=absorption_capacity(final['wo'],report['model_selection']['eta'],args.duration,args.required_absorption)
            report['molecular_status']='requires_propagation_validation' if report['basis_capacity']['necessary_condition_passed'] else 'rejected_insufficient_basis_capacity'
        else:
            report['molecular_status']='not_calibrated'
        matrix_name=output_stem(args.matrix,rcut)+'.mat'
        shutil.copyfile(final['artifact'],args.output/matrix_name)
        report['output_matrix']=matrix_name
        report.update(converged=True,selected=dict(r_min=0,r_max=final['grid'].extent,
            n_radial=final['grid'].radial,n_lebedev=final['grid'].angular,r_cut=rcut,
            eta=report.get('model_selection',{}).get('eta')),final_run=final['record']['index'],
            stopping_reason='Two consecutive independent plus joint grid-refinement checks passed')
    except (NotConverged,ValueError,subprocess.SubprocessError) as error:
        report['stopping_reason']=str(error)
    finally:
        report['orthogonalization']=runner.rank
        (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('converged','physical_calibration','stopping_reason')}),flush=True)
    if 'molecular_status' in report:print('Molecular status:',report['molecular_status'],flush=True)
    if not report['converged']:return 2
    return 3 if report.get('molecular_status')=='rejected_insufficient_basis_capacity' else 0


def main(argv=None):
    args=parser().parse_args(argv)
    if args.diagnostics:
        return run_workflow(args)
    output=args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be new or empty')
    output.mkdir(parents=True,exist_ok=True)
    # Refinement artifacts are scratch data; publish only the selected matrix/log.
    with tempfile.TemporaryDirectory(prefix='cap-convergence-') as td:
        args.output=Path(td)
        # run_workflow requires an empty output directory.
        with tempfile.TemporaryFile(mode='w+') as console:
            try:
                with redirect_stdout(console):
                    status=run_workflow(args)
            except Exception as error:
                console.seek(0)
                (output/f'CAP_{args.matrix.stem}_failed.log').write_text(console.read()+f'\nCAP failed: {error}\n')
                raise
            report=json.loads((args.output/'report.json').read_text())
            index=report.get('final_run')
            logs=sorted(args.output.glob('run_*.log'))
            selected_log=args.output/f'run_{index:03d}.log' if index is not None else None
            matrix_name=report.get('output_matrix')
            log_name=matrix_name[:-4]+'.log' if matrix_name else f'CAP_{args.matrix.stem}_failed.log'
            with (output/log_name).open('w') as log:
                if selected_log is not None:
                    log.write(selected_log.read_text().replace('Written matrix file: CAP_Mat.mat',
                                                               f'Written matrix file: {matrix_name}'))
                else:
                    for path in logs:
                        log.write(path.read_text())
                log.write('\nAutomatic CAP parameter selection\n')
                log.write(f'Input matrix: {args.matrix}\n')
                log.write(f'Quadrature converged: {report["quadrature_converged"]}\n')
                log.write(f'Stopping reason: {report["stopping_reason"]}\n')
                log.write(f'Integration runs: {len(report["runs"])}\n')
                if index is not None:
                    for key in ('overlap_error','electron_error','tail','density_expectation'):
                        log.write(f'  {key} = {report["runs"][index][key]:.12e}\n')
                log.write(f'CAP convergence tolerance: {args.cap_tol:.6e}\n')
                log.write(f'Overlap tolerance: {args.overlap_tol:.6e}\n')
                log.write(f'Density-tail budget: {args.density_tail:.6e}\n')
                rank=report.get('orthogonalization')
                if rank:
                    log.write(f'Overlap directions retained: {rank["retained"]} / {rank["total"]}\n')
                if args.energies:
                    log.write(f'Physical calibration: {report["physical_calibration"]}\n')
                    log.write(f'Molecular status: {report.get("molecular_status", "not_calibrated")}\n')
                    log.write(f'Selected eta: {report.get("selected",{}).get("eta")}\n')
                    log.write('Molecular dynamics validated: False\n')
                else:
                    log.write('Integral convergence only; eta not selected.\n')
                for key in ('basis_capacity',):
                    if key in report:
                        log.write(f'{key.replace("_", " ").capitalize()}:\n')
                        for name,value in report[key].items():
                            log.write(f'  {name} = {value}\n')
                log.write(f'Exit status: {status}\n')
            if matrix_name:
                shutil.copyfile(args.output/matrix_name,output/matrix_name)
                print(f'CAP matrix: {output / matrix_name}')
    print(f'CAP log: {output / log_name}')
    return status

if __name__=='__main__':raise SystemExit(main())
