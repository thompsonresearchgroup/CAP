"""Independent invariants for numerical/physical automation and failure paths."""
from pathlib import Path
import json
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from cap_auto import (Grid, NotConverged, converge, density_cutoff, orthogonalizer, opnorm,
                      absorption_capacity, select_calibrated_candidate)
from cap_reflection import reflection, scan, wavepacket, converge_wavepacket

ROOT=Path(__file__).resolve().parents[1]


class AutomationTests(unittest.TestCase):
    def test_density_onset_protects_requested_tail(self):
        r=np.linspace(0,20,20001);dr=r[1]-r[0]
        # Independent normalized decaying one-electron radial density.
        mass=2*np.exp(-2*r)*dr
        cut=density_cutoff(r,mass,1e-6)
        self.assertLessEqual(mass[r>cut].sum(),1e-6)
        self.assertGreater(cut,np.log(1e6)/2)
        with self.assertRaises(ValueError):density_cutoff(r,-mass,1e-6)

    def test_error_norm_is_invariant_to_ao_rescaling(self):
        s=np.array([[1,.2],[.2,.8]]);w=np.array([[.01,.002],[.002,.03]])
        c=np.diag([10,.02]);x,_=orthogonalizer(s,1e-12);y,_=orthogonalizer(c.T@s@c,1e-12)
        self.assertAlmostEqual(opnorm(x.T@w@x),opnorm(y.T@c.T@w@c@y),places=12)
        _,rank=orthogonalizer(np.diag([1,1e-12]),1e-9)
        self.assertEqual(rank['discarded'],1)

    def test_two_independent_and_joint_refinement_passes(self):
        class ManufacturedRunner:
            def __init__(self):self.count=0
            def run(self,g,cut):
                self.count+=1
                exact_error=1/g.radial**2+1/g.angular**2+np.exp(-g.extent)
                return dict(grid=g,wo=np.array([[2+exact_error]]),record=dict(index=self.count,
                    overlap_error=0,electron_error=0,tail=0,density_cutoff=cut))
        history=[]
        final=converge(ManufacturedRunner(),Grid(200,590,20),3,
            SimpleNamespace(cap_tol=1e-4,overlap_tol=1e-8,electron_tol=1e-8,density_tail=1e-6),history)
        self.assertEqual(history[-1]['consecutive_passes'],2)
        self.assertTrue(all(set(x['checks'])=={'radial','angular','extent'} for x in history))
        self.assertLess(abs(final['wo'][0,0]-2),1e-4)
        with self.assertRaises(NotConverged):
            converge(ManufacturedRunner(),Grid(200,5810,20),3,SimpleNamespace(),[])

    def test_default_output_is_one_matrix_and_one_log(self):
        with tempfile.TemporaryDirectory(prefix='cap-output-') as td:
            for budget,expected_status in ((50,0),(1,2)):
                out=Path(td)/str(budget)
                result=subprocess.run(['python3',str(ROOT/'cap_auto.py'),str(ROOT/'tests/gaussian/h_uhf.mat'),
                    '--output',str(out),'--max-runs',str(budget)],
                    stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=120)
                self.assertEqual(result.returncode,expected_status,result.stdout.decode())
                files={p.name for p in out.iterdir()}
                if expected_status==0:
                    matrices=list(out.glob('CAP_h_uhf_rcut_*.mat'))
                    self.assertEqual(len(matrices),1)
                    log_path=matrices[0].with_suffix('.log')
                    self.assertEqual(files,{matrices[0].name,log_path.name})
                else:
                    log_path=out/'CAP_h_uhf_failed.log'
                    self.assertEqual(files,{log_path.name})
                log=log_path.read_text()
                for heading in ('Smooth Voronoi CAP integral generator','Analytical Overlap',
                                'Numerical Overlap','CAP matrix','Numerical CAP Trace',
                                'Automatic CAP parameter selection'):
                    self.assertIn(heading,log)
                self.assertIn('budget' if expected_status else 'Two consecutive',log)
                for removed in ('Legacy APO Vector','APO Vector W_p','APO Vector Sum',
                                'APO matrix dimensions','Convergence history:','driver_sha256'):
                    self.assertNotIn(removed,log)

    def test_budget_failure_never_publishes_final_matrix(self):
        with tempfile.TemporaryDirectory(prefix='cap-failure-') as td:
            out=Path(td)/'result'
            process=subprocess.run(['python3',str(ROOT/'cap_auto.py'),str(ROOT/'tests/gaussian/n2_pure.mat'),
                '--output',str(out),'--max-runs','1','--diagnostics'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
            self.assertEqual(process.returncode,2,process.stdout.decode())
            report=json.loads((out/'report.json').read_text())
            self.assertFalse(report['converged']);self.assertIn('budget',report['stopping_reason'])
            self.assertFalse(list(out.glob('CAP_*.mat')))
            self.assertAlmostEqual(report['runs'][0]['electron_count'],14,places=6)

    def test_finite_basis_capacity_rejects_ineffective_absorber(self):
        small=absorption_capacity(np.diag([0,.001]),.0003,100,.99)
        self.assertFalse(small['necessary_condition_passed'])
        self.assertAlmostEqual(small['maximum_possible_absorbed_fraction'],1-np.exp(-.00006),places=14)
        large=absorption_capacity(np.diag([0,100]),.001,100,.99)
        self.assertTrue(large['necessary_condition_passed'])

    def test_candidate_selection_uses_exact_onset_capacity(self):
        args=SimpleNamespace(duration=100,required_absorption=.99,model_end=40,
                             reflection_tol=.001,density_tail=1e-6,loss_budget=1e-4)
        base=dict(wo=np.array([[100.]]),record=dict(rcut=3.,tail=0.,density_expectation=0.))
        later=dict(wo=np.array([[.01]]),record=dict(rcut=4.,tail=0.,density_expectation=0.))
        rows=[dict(onset=4.,eta=.001,worst_reflection=1e-5),
              dict(onset=3.,eta=.001,worst_reflection=2e-5)]
        with patch('cap_auto.reflection',return_value=1e-5),patch('cap_auto.converge',return_value=later):
            best,final,_,_=select_calibrated_candidate(rows,base,None,None,args,[],[.05,.5])
        self.assertEqual(best['onset'],3.)
        self.assertIs(final,base)
        self.assertEqual(rows[0]['selection_rejection'],'insufficient_exact_basis_capacity')

    def test_candidate_selection_retries_failed_fine_reflection(self):
        args=SimpleNamespace(duration=100,required_absorption=.99,model_end=40,
                             reflection_tol=.001,density_tail=1e-6,loss_budget=1e-4)
        base=dict(wo=np.array([[100.]]),record=dict(rcut=3.,tail=0.,density_expectation=0.))
        rows=[dict(onset=3.,eta=.001,worst_reflection=1e-5),
              dict(onset=3.,eta=.002,worst_reflection=2e-5)]
        with patch('cap_auto.reflection',side_effect=lambda e,eta,*rest:.01 if eta==.001 else 1e-5):
            best,_,_,_=select_calibrated_candidate(rows,base,None,None,args,[],[.05,.5])
        self.assertEqual(best['eta'],.002)
        self.assertEqual(rows[0]['selection_rejection'],'refined_model_reflection')

    def test_open_shell_density_electron_count(self):
        with tempfile.TemporaryDirectory(prefix='cap-uhf-') as td:
            out=Path(td)/'result'
            result=subprocess.run(['python3',str(ROOT/'cap_auto.py'),str(ROOT/'tests/gaussian/h_uhf.mat'),
                '--output',str(out),'--max-runs','1','--diagnostics'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
            self.assertEqual(result.returncode,2,result.stdout.decode())
            report=json.loads((out/'report.json').read_text())
            self.assertAlmostEqual(report['runs'][0]['electron_count'],1,places=10)
            data=np.load(out/'run_000.npz')
            self.assertAlmostEqual(float(np.trace(data['P']@data['S'])),1,places=10)

    def test_scattering_conservation_and_resolution(self):
        for energy in [.05,.2,.5]:
            self.assertAlmostEqual(reflection(energy,0,5,40),1,places=11)
            r=reflection(energy,.0003,5,40)
            self.assertGreaterEqual(r,0);self.assertLess(r,1)
            self.assertLess(abs(r-reflection(energy,.0003,5,40,.05)),2e-5)
        with self.assertRaises(ValueError):reflection(.1,-1,5,40)
        with self.assertRaises(ValueError):reflection(.1,float('nan'),5,40)

    def test_propagation_budget_does_not_claim_convergence(self):
        result=converge_wavepacket(.5,.001,5,30,tolerance=1e-20,max_levels=1)
        self.assertFalse(result['converged'])
        self.assertIn('limit',result['stopping_reason'])

    def test_time_propagator_unitarity_and_loss(self):
        free=wavepacket(.5,0,5,30,dx=.4,dt=.3)
        cap=wavepacket(.5,.001,5,30,dx=.4,dt=.3)
        self.assertAlmostEqual(free['final_norm'],free['initial_norm'],places=10)
        self.assertLess(cap['max_norm_increase'],1e-12)
        self.assertLess(cap['max_inner_error'],free['max_inner_error']/10)

if __name__=='__main__':unittest.main()
