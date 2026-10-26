"""
Python port of the nondimensional firn compaction model of
Kingslake et al. (2022, The Cryosphere 16, 3413-3430), used to generate the
data shown in the explainer video.

Full model: method of lines on the stretched coordinate zh = z / h(t)
(Appendix B), two-point upwind differences, stiff BDF time stepping.
Reduced models: steady-state ODEs (Eqs. 22, 24, 26) and closed-form
solutions of Eqs. 25 and 27 (n = m = 1).
"""
import numpy as np
from scipy.integrate import solve_ivp

# physical constants (Table 2)
G_ACC = 9.81
SPY = 24 * 365 * 3600
R_GAS = 8.31
RHO_I = 918.0
E_C = 60e3
E_G = 42e3
K_C = 9.2e-9
R2_F = 0.01 ** 2
K_G = 1.3e-7 / R2_F
PHI_830 = 1 - 830 / RHO_I


def scales(b0_mpy=0.1, T_s=253.15, z0=100.0, n=1):
    b0 = b0_mpy / SPY
    t0 = z0 / b0
    sigma0 = G_ACC * RHO_I * z0
    r2_0 = (z0 * K_G * R2_F / b0) * np.exp(-E_G / (R_GAS * T_s))
    gamma = (sigma0 / 4) ** (1 - n)
    alpha = r2_0 / (K_C * t0 * sigma0 ** n * np.exp(-E_C / (R_GAS * T_s))) / gamma
    delta = r2_0 / R2_F
    return dict(b0=b0, t0=t0, t0_yr=t0 / SPY, sigma0=sigma0, r2_0=r2_0,
                alpha=alpha, delta=delta, z0=z0, w0_mpy=b0_mpy)


def _cumtrapz(f, dz):
    out = np.empty_like(f)
    out[0] = 0.0
    np.cumsum(0.5 * (f[1:] + f[:-1]) * dz, out=out[1:])
    return out


def _upwind(y, dz):
    g = np.empty_like(y)
    g[1:] = (y[1:] - y[:-1]) / dz
    g[0] = (y[1] - y[0]) / dz
    return g


def z830_from_profile(z, phi):
    """Depth at which porosity first drops to the 830 kg/m^3 value."""
    idx = np.where(phi < PHI_830)[0]
    if len(idx) == 0 or idx[0] == 0:
        return np.nan
    i = idx[0]
    return z[i - 1] + (PHI_830 - phi[i - 1]) * (z[i] - z[i - 1]) / (phi[i] - phi[i - 1])


class FullModel:
    def __init__(self, beta=1.0, rs2=0.0284, alpha=None, delta=None, phi_s=0.5,
                 n=1, m=1, sim_r=True, nu=1.0, N=101, **scale_kw):
        sc = scales(n=n, **scale_kw)
        self.sc = sc
        self.alpha = sc['alpha'] if alpha is None else alpha
        self.delta = sc['delta'] if delta is None else delta
        self.beta, self.rs2, self.phi_s = beta, rs2, phi_s
        self.n, self.m, self.sim_r, self.nu = n, m, sim_r, nu
        self.N = N
        self.zh = np.linspace(0, 1, N)
        self.dz = self.zh[1] - self.zh[0]
        self.w_s = nu * beta / (1 - phi_s)

    def unpack(self, y):
        N = self.N
        return y[:N], y[N:2 * N], y[2 * N:3 * N], y[-1]

    def diagnostics(self, y):
        phi, r2, A, H = self.unpack(y)
        sigma = _cumtrapz(H * (1 - phi), self.dz)        # |sigma|
        v_int = -(H / self.alpha) * sigma ** self.n * phi ** self.m / r2
        w = _cumtrapz(v_int, self.dz) + self.w_s
        return sigma, w

    def rhs(self, t, y):
        N, dz, zh = self.N, self.dz, self.zh
        phi, r2, A, H = self.unpack(y)
        sigma, w = self.diagnostics(y)
        dH = w[-1] - self.beta / (1 - phi[-1])
        dphi = (1 / H) * (_upwind((1 - phi) * w, dz) + dH * zh * _upwind(phi, dz))
        if self.sim_r:
            dr2 = (1 / H) * (dH * zh - w) * _upwind(r2, dz) + (1 - self.delta * r2)
        else:
            dr2 = np.zeros(N)
        dA = 1 + (1 / H) * (dH * zh - w) * _upwind(A, dz)
        dphi[0] = dr2[0] = dA[0] = 0.0
        return np.concatenate([dphi, dr2, dA, [dH]])

    def y0(self):
        zh = self.zh
        phi = (1 - zh) * self.phi_s
        r2 = self.rs2 + (zh if self.sim_r else 0 * zh)
        return np.concatenate([phi, r2, zh.copy(), [1.0]])

    def run(self, t_end=10.0, t_eval=None, steady_tol=1e-5, stop_at_steady=True, rtol=1e-8, atol=1e-8):
        N = self.N

        def steady(t, y):
            return np.max(np.abs(self.rhs(t, y)[:N])) - steady_tol
        steady.terminal = stop_at_steady
        steady.direction = -1
        events = [steady] if stop_at_steady else None
        sol = solve_ivp(self.rhs, (0, t_end), self.y0(), method='BDF', t_eval=t_eval,
                        events=events, rtol=rtol, atol=atol)
        self.sol = sol
        return sol

    def profiles(self, y):
        phi, r2, A, H = self.unpack(y)
        sigma, w = self.diagnostics(y)
        z = H * self.zh
        return dict(z=z, phi=phi, r2=r2, A=A, sigma=-sigma, w=w, H=H,
                    z830=z830_from_profile(z, phi))

    def final(self):
        return self.profiles(self.sol.y[:, -1])


# ---------------------------------------------------------------- ODE models
def ode_model(beta=1.0, rs2=0.0284, alpha=None, delta=None, phi_s=0.5, n=1, m=1,
              z_max=1.0, nz=401, variant='full'):
    """Steady-state ODE models.
    variant: 'full'  -> Eq. 22
             'eq24'  -> sigma=-z, r2=rs2
             'eq26'  -> sigma=-z, delta=0, evolving r2
    """
    sc = scales(n=n)
    alpha = sc['alpha'] if alpha is None else alpha
    delta = sc['delta'] if delta is None else delta
    z = np.linspace(0, z_max, nz)

    def f(zz, y):
        phi, sig, w, r2, A = y
        phi = max(phi, 0.0)
        s = abs(sig) if variant == 'full' else zz
        r2e = r2 if variant in ('full', 'eq26') else rs2
        dphi = -s ** n * phi ** m * (1 - phi) / (alpha * w * r2e)
        dsig = -(1 - phi)
        dw = -s ** n * phi ** m / (alpha * r2e)
        if variant == 'full':
            dr2 = (1 - delta * r2) / w
        elif variant == 'eq26':
            dr2 = 1 / w
        else:
            dr2 = 0.0
        dA = 1 / w
        return [dphi, dsig, dw, dr2, dA]
    y0 = [phi_s, 0.0, beta / (1 - phi_s), rs2, 0.0]
    sol = solve_ivp(f, (0, z_max), y0, t_eval=z, method='RK45', rtol=1e-10, atol=1e-10)
    phi, sig, w, r2, A = sol.y
    return dict(z=z, phi=phi, sigma=sig, w=w, r2=r2, A=A, z830=z830_from_profile(z, phi))


L830 = None


def logit_drop(phi_s=0.5):
    """ln[(phi_s/(1-phi_s)) / (phi_830/(1-phi_830))]"""
    return np.log((phi_s / (1 - phi_s)) / (PHI_830 / (1 - PHI_830)))


def eq25_phi(z, beta, rs2, alpha, phi_s=0.5):
    q = (phi_s / (1 - phi_s)) * np.exp(-z ** 2 / (2 * alpha * beta * rs2))
    return q / (1 + q)


def eq25_z830(beta, rs2, alpha, phi_s=0.5):
    return np.sqrt(2 * alpha * beta * rs2 * logit_drop(phi_s))


def eq27_F(z, beta, rs2):
    """Integral of z/(z + beta rs2) dz: the 'effective compaction' of Eq. 27 (n=1)."""
    br = beta * rs2
    if br == 0:
        return z
    return z - br * np.log1p(z / br)


def eq27_phi(z, beta, rs2, alpha, phi_s=0.5):
    q = (phi_s / (1 - phi_s)) * np.exp(-eq27_F(z, beta, rs2) / alpha)
    return q / (1 + q)


def eq27_z830(beta, rs2, alpha, phi_s=0.5):
    from scipy.optimize import brentq
    target = alpha * logit_drop(phi_s)
    return brentq(lambda z: eq27_F(z, beta, rs2) - target, 1e-9, 50.0)
