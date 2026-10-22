#!/usr/bin/env python3
"""Independent radial outgoing-packet diagnostics, in atomic units.

For a single CAP center, W=(r-onset)_+^2. The s-wave radial function u(r)
has measure dr, u(0)=0, and free Hamiltonian -1/2 d²/dr². Its odd extension
provides an exact infinite-domain reference. This is not molecular dynamics.
"""
import numpy as np
from scipy.linalg import expm
from scipy.sparse import diags
from scipy.sparse.linalg import splu
from scipy.fft import dst, idst


def radial_packet(r, time, energy, center=2.0, sigma=1.0):
    """Normalized odd-image Gaussian; energy is the carrier energy, not its mean.

    The packet has an energy distribution and can contain initially inward
    components. The exact reference includes their reflection at the origin.
    """
    if energy <= 0 or sigma <= 0 or center <= 0 or time < 0:
        raise ValueError('Require positive energy, width, center and nonnegative time')
    k = np.sqrt(2*energy)
    z = 1 + 1j*time/(2*sigma*sigma)
    def g(x):
        return np.exp(-(x-center-k*time)**2/(4*sigma*sigma*z)
                      + 1j*k*(x-center) - .5j*k*k*time)/np.sqrt(z)
    norm2 = np.sqrt(2*np.pi)*sigma*(-np.expm1(-center**2/(2*sigma*sigma)-2*k*k*sigma*sigma))
    return (g(np.asarray(r))-g(-np.asarray(r)))/np.sqrt(norm2)


def inner_error(r, numerical, exact, onset):
    """Squared radial L2 error including the zero origin and onset endpoints."""
    mask = r <= onset + 1e-12
    return float(np.trapezoid(np.r_[0., abs(numerical[mask]-exact[mask])**2],
                             np.r_[0., r[mask]]))


def radial_grid(energy, eta, onset=6., end=80., dx=.1, dt=.05,
                duration=200., sample_step=1., method='sine'):
    """Sine-spectral split operator or independent Crank–Nicolson reference.

    Both obey Dirichlet conditions at zero and end. Sine kinetic propagation
    has no finite-difference dispersion; split-step time error still requires
    refinement. The caller must check spatial, time and extent convergence.
    """
    if not np.isfinite([energy, eta, onset, end, dx, dt, duration, sample_step]).all():
        raise ValueError('Nonfinite propagation parameter')
    if min(energy, onset, dx, dt, duration, sample_step) <= 0 or eta < 0 or end <= onset:
        raise ValueError('Invalid radial propagation parameter')
    if method not in ('sine', 'cn'):
        raise ValueError('Unknown propagation method')
    intervals = round(end/dx)
    dx = end/intervals
    if abs(onset/dx-round(onset/dx)) > 1e-8:
        raise ValueError('Onset must coincide with a radial grid point')
    sample_count = round(duration/sample_step)
    sample_step = duration/sample_count
    substeps = int(np.ceil(sample_step/dt))
    dt = sample_step/substeps
    r = np.arange(1, intervals)*dx
    u = radial_packet(r, 0, energy)
    if method == 'cn':
        h = diags([-np.ones(len(r)-1)/(2*dx*dx),
                   np.ones(len(r))/dx**2-1j*eta*np.maximum(r-onset, 0)**2,
                   -np.ones(len(r)-1)/(2*dx*dx)], [-1, 0, 1], format='csc')
        identity = diags(np.ones(len(r)), format='csc')
        lhs = splu(identity+.5j*dt*h)
        rhs = identity-.5j*dt*h
        def advance(u):
            return lhs.solve(rhs@u)
    else:
        damping = np.exp(-.5*dt*eta*np.maximum(r-onset, 0)**2)
        phase = np.exp(-.5j*dt*(np.arange(1, intervals)*np.pi/end)**2)
        def advance(u):
            return damping*idst(phase*dst(damping*u, type=1, norm='ortho'), type=1, norm='ortho')
    times = np.arange(sample_count+1)*sample_step
    norms = [float(np.vdot(u, u).real*dx)]
    errors = [0.]
    for t in times[1:]:
        for _ in range(substeps):
            u = advance(u)
        norms.append(float(np.vdot(u, u).real*dx))
        errors.append(inner_error(r, u, radial_packet(r, t, energy), onset))
    return dict(times=times, norms=np.array(norms), errors=np.array(errors),
                dx=dx, dt=dt, end=end, duration=duration, method=method)


def finite_basis_packet(s, kinetic, w, exponents, energy, eta, cutoff=1e-9,
                        onset=6., duration=200., sample_step=1.):
    """Propagate a packet using actual Gaussian S/T and executable CAP W.

    Only normalized, uncontracted, one-center S functions are supported. The
    caller must establish the AO ordering by comparing S and T analytically.
    """
    from cap_auto import orthogonalizer
    x, rank = orthogonalizer(s, cutoff)
    h = x.T@kinetic@x
    wo = x.T@w@x
    # Resolve the initial packet and all tight basis functions independently.
    from scipy.special import gamma
    r = np.arange(.01, 80.01, .01)
    a = np.asarray(exponents)
    basis = (2*(2*a)**1.5/gamma(1.5))**.5*r[:, None]*np.exp(-r[:, None]**2*a)
    q = basis@x
    packet = radial_packet(r, 0, energy)
    c = q.T@packet*.01
    captured = float(np.vdot(c, c).real)
    projection_error = float(np.sum(abs(q@c-packet)**2)*.01)
    if captured <= 0 or not np.isfinite(captured):
        raise ValueError('Packet has no finite basis projection')
    c /= np.sqrt(captured)
    times = np.arange(round(duration/sample_step)+1)*sample_step
    step = expm((-1j*h-eta*wo)*sample_step)
    free = expm(-1j*h*sample_step)
    free_c = c.copy()
    norms, errors, free_errors = [], [], []
    for t in times:
        reference = radial_packet(r, t, energy)
        norms.append(float(np.vdot(c, c).real))
        errors.append(inner_error(r, q@c, reference, onset))
        free_errors.append(inner_error(r, q@free_c, reference, onset))
        c = step@c
        free_c = free@free_c
    return dict(times=times, norms=np.array(norms), errors=np.array(errors),
                free_errors=np.array(free_errors), initial_projection_error=projection_error,
                captured_initial_norm=captured, rank=rank,
                no_cap_final_norm_error=float(abs(np.vdot(free_c, free_c).real-1)))


def summarize(result):
    times, norms, errors = result['times'], result['norms'], result['errors']
    summary = dict(final_norm=float(norms[-1]), absorbed_fraction=float(norms[0]-norms[-1]),
                   maximum_inner_squared_error=float(errors.max()),
                   maximum_sampled_norm_increase=float(max(0, np.diff(norms).max())),
                   duration=float(times[-1]))
    i = np.flatnonzero(np.isclose(times, 100))
    if len(i):
        summary['absorbed_fraction_at_100'] = float(norms[0]-norms[i[0]])
    for key in ('initial_projection_error', 'captured_initial_norm', 'rank', 'no_cap_final_norm_error'):
        if key in result:
            summary[key] = result[key]
    if 'free_errors' in result:
        summary['no_cap_maximum_inner_squared_error'] = float(result['free_errors'].max())
    return summary
