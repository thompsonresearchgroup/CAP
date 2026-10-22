"""Execute today's binary and read back its Gaussian matrix records.

Run `make test`; NumPy is required. All outputs use disposable directories.
These tests complement (and do not rely on) the archived Gaussian CAP logs.
"""
from pathlib import Path
import os
import math
import re
import subprocess
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ENV = dict(os.environ, OMP_NUM_THREADS='2', MKL_NUM_THREADS='1', OMP_DYNAMIC='FALSE')


def record(directory, name):
    matrix=next(Path(directory).glob('CAP_*.mat'))
    subprocess.run([str(ROOT/'tests/read_cap'), str(matrix), name],
                   cwd=directory, env=ENV, stdout=subprocess.DEVNULL, check=True)
    return np.loadtxt(Path(directory)/'record.txt', skiprows=1, ndmin=2)


def run_case(stem, angular=590, radial=200, extent=24, extra=(), threads=2):
    with tempfile.TemporaryDirectory(prefix='cap-test-') as td:
        d = Path(td)
        (d/'input.list').write_text(f'1\n{ROOT}/tests/gaussian/{stem}.mat\n')
        extent_args=['--r_max',str(extent)] if extent is not None else []
        result = subprocess.run([str(ROOT/'noCAPs2withCAPs'), '-f', 'input.list',
            '--lebedev_grid', str(ROOT/f'lebedev/lebedev_{angular}.txt'),
            '--n_lebedev', str(angular), '--n_radial', str(radial), *extent_args,
            '--validate', *extra], cwd=d, env=dict(ENV, OMP_NUM_THREADS=str(threads)),
            text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        cutoff=float(extra[extra.index('--r_cut')+1]) if '--r_cut' in extra else 3.
        radius=f'{cutoff:.6f}'.rstrip('0').rstrip('.')
        expected=f'CAP_{stem}_rcut_{radius}.mat'
        if '--output' in extra:
            expected=extra[extra.index('--output')+1]
        if result.returncode or {p.name for p in d.glob('CAP_*.mat')}!={expected}:
            raise AssertionError(result.stdout)
        overlap_error = float(re.search(r'max \|S_num - S_analytical\|\s*=\s*(\S+)',result.stdout)[1])
        return {name:record(d,name) for name in ['numerical overlap',
            'Complex Absorbing Potential','Alpha Complex Absorbing Potential',
            'Beta Complex Absorbing Potential','APO Vector']}, overlap_error


class ExecutableTests(unittest.TestCase):
    def test_output_override(self):
        records,error=run_case('h_uhf',extra=['--output','CAP_custom.mat'])
        self.assertLess(error,1e-8)
        self.assertTrue(np.isfinite(records['Complex Absorbing Potential']).all())

    def test_gaussian_all_shells_and_written_records(self):
        for kind in ['pure', 'cartesian']:
            with self.subTest(kind=kind):
                records, error = run_case('one_center_spdfghi_'+kind, extra=['--unit_potential'])
                self.assertLess(error, 3e-12)
                s, w = records['numerical overlap'],records['Complex Absorbing Potential']
                np.testing.assert_array_equal(s,w)
                for name,a in records.items():
                    self.assertTrue(np.isfinite(a).all(), name)
                    if 'Absorbing' in name: np.testing.assert_array_equal(a,w)

    def test_analytic_one_center_r_squared(self):
        records, _ = run_case('one_center_spdfghi_pure', extra=['--r_cut','0'])
        w = records['Complex Absorbing Potential']
        # A normalized, uncontracted l shell has <r^2>=(l+3/2)/(2 alpha).
        offset=5  # contracted S, then the SP shell's four functions
        for l,alpha in [(2,.7),(3,.6),(4,.5),(5,.4),(6,.35)]:
            n=2*l+1
            np.testing.assert_allclose(w[offset:offset+n,offset:offset+n],
                                       np.eye(n)*(l+1.5)/(2*alpha),atol=3e-12,rtol=0)
            offset+=n
        self.assertGreaterEqual(np.linalg.eigvalsh(w).min(),-1e-12)

    def test_analytic_one_center_finite_onset(self):
        records,_=run_case('one_center_spdfghi_pure',radial=1000,
                           extra=['--r_cut','3','--smoothing_length','0.2'])
        w=records['Complex Absorbing Potential']
        # Incomplete Gaussian moments J_n(c)=integral_c^inf r^n exp(-p r^2) dr.
        def moment(n,p,c):
            if n==0: return math.sqrt(math.pi/p)*math.erfc(math.sqrt(p)*c)/2
            if n==1: return math.exp(-p*c*c)/(2*p)
            return c**(n-1)*math.exp(-p*c*c)/(2*p)+(n-1)*moment(n-2,p,c)/(2*p)
        offset=5
        for l,alpha in [(2,.7),(3,.6),(4,.5),(5,.4),(6,.35)]:
            n=2*l+1; power=2*l+2; p=2*alpha; c=3
            expected=(moment(power+2,p,c)-2*c*moment(power+1,p,c)+c*c*moment(power,p,c))/moment(power,p,0)
            np.testing.assert_allclose(w[offset:offset+n,offset:offset+n],
                                       np.eye(n)*expected,atol=3e-7,rtol=0)
            offset+=n

    def test_becke_gaussian_overlap_and_psd(self):
        for stem in ['n2_pure','n2_cartesian','n2_pure_translated_rotated','h2_pure','h2_cartesian','h2_rotated']:
            with self.subTest(stem=stem):
                records,error=run_case(stem,angular=1202,radial=300,extra=['--becke_grid'])
                self.assertLess(error,2e-8)
                w=records['Complex Absorbing Potential']
                np.testing.assert_allclose(w,w.T,atol=2e-14,rtol=0)
                self.assertGreaterEqual(np.linalg.eigvalsh(w).min(),-2e-12)
                self.assertTrue(np.isfinite(records['APO Vector']).all())

    def test_off_axis_multicenter_high_angular_momentum(self):
        for kind in ['pure','cartesian']:
            with self.subTest(kind=kind):
                records,error=run_case('he2_spdfghi_'+kind,angular=1202,radial=300)
                # Cross-center off-diagonal overlaps expose ordering and phases.
                self.assertLess(error,2e-11)
                self.assertGreater(np.linalg.eigvalsh(records['Complex Absorbing Potential']).min(),0)

    def test_default_becke_and_explicit_centroid(self):
        default,_=run_case('n2_pure',radial=60)
        explicit,_=run_case('n2_pure',radial=60,extra=['--becke_grid'])
        centroid,_=run_case('n2_pure',radial=60,extra=['--centroid_grid'])
        for name in default: np.testing.assert_array_equal(default[name],explicit[name])
        self.assertGreater(np.max(np.abs(default['numerical overlap']-centroid['numerical overlap'])),1e-5)

    def test_automatic_extent(self):
        # One-center fixture: most diffuse exponent .22 (SP), l=1.
        extent=math.sqrt(52/(2*.22))
        auto,error=run_case('one_center_spdfghi_pure',extent=None,extra=['--auto_rmax'])
        manual,_=run_case('one_center_spdfghi_pure',extent=extent)
        self.assertLess(error,3e-12)
        for name in auto: np.testing.assert_allclose(auto[name],manual[name],atol=2e-12,rtol=0)

    def test_thread_reproducibility(self):
        one,_=run_case('n2_pure',radial=100,extra=['--becke_grid'],threads=1)
        four,_=run_case('n2_pure',radial=100,extra=['--becke_grid'],threads=4)
        for name in one: np.testing.assert_array_equal(one[name],four[name])

    def test_zero_cap_and_simpson(self):
        records,error=run_case('one_center_spdfghi_pure',radial=801,
                               extra=['--equidistant_radial','--r_cut','100'])
        self.assertLess(error,1e-10)
        self.assertEqual(np.max(np.abs(records['Complex Absorbing Potential'])),0)
        self.assertEqual(np.max(np.abs(records['APO Vector'])),0)

    def test_bad_inputs_fail(self):
        with tempfile.TemporaryDirectory(prefix='cap-bad-') as td:
            d=Path(td); (d/'zero.list').write_text('0\n')
            (d/'many.list').write_text('2\nmissing.mat\nmissing2.mat\n')
            (d/'input.list').write_text(f'1\n{ROOT}/tests/gaussian/n2_pure.mat\n')
            (d/'bad.grid').write_text('nan 0 1\n'*6)
            cases=[['--auto_rmax','--r_max','20'],['--r_max','NaN'],['--smoothing_length','0'],['--r_min','-1'],
                   ['--n_radial','1'],['--n_radial'],['--unknown'],['-f','zero.list'],
                   ['-f','many.list'],['-f','input.list','--lebedev_grid','bad.grid','--n_lebedev','6'],
                   ['-f','input.list','--lebedev_grid',str(ROOT/'lebedev/lebedev_74.txt'),'--n_lebedev','74']]
            for args in cases:
                with self.subTest(args=args):
                    p=subprocess.run([str(ROOT/'noCAPs2withCAPs'),*args],cwd=d,env=ENV,
                                     stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
                    self.assertNotEqual(p.returncode,0,p.stdout.decode(errors='replace'))
                    self.assertFalse(list(d.glob('CAP_*.mat')))

if __name__=='__main__': unittest.main()
