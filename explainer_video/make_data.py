"""Generate every dataset plotted in the explainer video (writes data/*.npz)."""
import numpy as np
import scipy.io as sio
from firn_sim import (scales, FullModel, ode_model, eq25_z830, eq27_z830, eq25_phi,
                      eq27_phi, z830_from_profile, PHI_830)

SAVED = '../savedResults/'
sc = scales()
ALPHA = sc['alpha']
RS2_DEF = 2.5e-7 / sc['r2_0']          # (0.5 mm)^2 at the intermediate climate

# --- Fig 1: alpha(Ts) and delta(Ts, b0)
Ts = np.linspace(233, 273, 41)
b0 = np.logspace(-2, 0, 41)
A = np.array([[scales(b, T)['alpha'] for T in Ts] for b in b0])
D = np.array([[scales(b, T)['delta'] for T in Ts] for b in b0])
np.savez('data/fig1.npz', Ts=Ts, b0=b0, alpha=A, delta=D)

# --- Fig 2: time evolution of the full model + steady ODE comparison
t_eval = np.linspace(0, 1.2, 145)
M = FullModel(beta=1, rs2=RS2_DEF)
sol = M.run(t_end=1.2, t_eval=t_eval, stop_at_steady=False)
prof = [M.profiles(sol.y[:, i]) for i in range(len(sol.t))]
keys = ['z', 'phi', 'r2', 'A', 'sigma', 'w']
out = {k: np.array([p[k] for p in prof]) for k in keys}
out['H'] = np.array([p['H'] for p in prof])
out['z830'] = np.array([p['z830'] for p in prof])
out['t'] = sol.t
O = ode_model(beta=1, rs2=RS2_DEF, z_max=1.0, nz=101)
for k in keys:
    out['ode_' + k] = O[k]
out['ode_z830'] = O['z830']
np.savez('data/fig2.npz', **out)
print('fig2 done; steady z830', out['z830'][-1], 'ode', O['z830'])

# --- Fig 3: no grain growth, beta sweep
betas3 = np.linspace(0.5, 10, 20)
f3 = dict(beta=betas3, rs2=RS2_DEF, zgrid=np.linspace(0, 1, 201))
full_phi, full_w, full_z, z830_full = [], [], [], []
ode_phi, ode_w, z830_ode = [], [], []
an_phi, z830_an = [], []
for b in betas3:
    M = FullModel(beta=b, rs2=RS2_DEF, sim_r=False)
    M.run(t_end=20)
    P = M.final()
    full_phi.append(np.interp(f3['zgrid'], P['z'], P['phi'], right=np.nan))
    full_w.append(np.interp(f3['zgrid'], P['z'], P['w'], right=np.nan))
    z830_full.append(P['z830'])
    O = ode_model(beta=b, rs2=RS2_DEF, variant='eq24', z_max=1.0, nz=201)
    ode_phi.append(O['phi']); ode_w.append(O['w']); z830_ode.append(O['z830'])
    an_phi.append(eq25_phi(f3['zgrid'], b, RS2_DEF, ALPHA))
    z830_an.append(eq25_z830(b, RS2_DEF, ALPHA))
f3.update(full_phi=np.array(full_phi), full_w=np.array(full_w), z830_full=np.array(z830_full),
          ode_phi=np.array(ode_phi), ode_w=np.array(ode_w), z830_ode=np.array(z830_ode),
          an_phi=np.array(an_phi), z830_an=np.array(z830_an))
np.savez('data/fig3.npz', **f3)
print('fig3 z830 full', np.round(f3['z830_full'][[0, 9, 19]], 3), 'ode', np.round(f3['z830_ode'][[0, 9, 19]], 3), 'an', np.round(f3['z830_an'][[0, 9, 19]], 3))

# --- Fig 4/5/6a: saved MATLAB results from the paper (21 beta x 21 rs2)
d = sio.loadmat(SAVED + 'f4_results.mat', squeeze_me=True, struct_as_record=False)
beta4, rs4 = d['beta'], d['r_s']
R = d['results_1']
zgrid = np.linspace(0, 1.3, 261)
phi4 = np.full((21, 21, len(zgrid)), np.nan)
w4 = np.full_like(phi4, np.nan)
r24 = np.full_like(phi4, np.nan)
z830_4 = np.zeros((21, 21))
for j in range(21):
    for i in range(21):
        r = R[j, i]
        zr = r.Depth / 100.0          # saved depths are dimensional (m); z0 = 100 m
        phi4[j, i] = np.interp(zgrid, zr, r.Phi, right=np.nan)
        w4[j, i] = np.interp(zgrid, zr, r.W, right=np.nan)
        r24[j, i] = np.interp(zgrid, zr, r.GrainSize, right=np.nan)
        z830_4[j, i] = r.zeta830final
slopes = np.array([np.polyfit(beta4, z830_4[:, i], 1)[0] for i in range(21)])
np.savez('data/fig4.npz', beta=beta4, rs2=rs4, zgrid=zgrid, phi=phi4, w=w4, r2=r24,
         z830=z830_4, slope=slopes)
print('fig4 slopes', slopes[[0, 2, 20]])

# check this port against the paper's own saved result
M = FullModel(beta=beta4[10], rs2=rs4[10]); M.run(t_end=20)
print('port check z830', M.final()['z830'], 'saved', z830_4[10, 10])

# --- Fig 6b, 6c: ODE (Eq 26) and closed form (Eq 27) on the same grid
z26 = np.zeros((21, 21)); z27 = np.zeros((21, 21))
for j, b in enumerate(beta4):
    for i, r in enumerate(rs4):
        z26[j, i] = ode_model(beta=b, rs2=r, variant='eq26', z_max=2.0, nz=801)['z830']
        z27[j, i] = eq27_z830(b, r, ALPHA)
np.savez('data/fig6.npz', beta=beta4, rs2=rs4, full=z830_4, eq26=z26, eq27=z27)
print('fig6 corners full', z830_4[[0, -1]][:, [0, -1]], 'eq26', z26[[0, -1]][:, [0, -1]], 'eq27', z27[[0, -1]][:, [0, -1]])

# --- Fig 7: saved n = 1..4
d7 = sio.loadmat(SAVED + 'f7_results.mat', squeeze_me=True, struct_as_record=False)
np.savez('data/fig7.npz', beta=d7['beta'], rs2=d7['r_s'], z830=d7['zeta_830_full'])

# --- Fig 8: saved thickening / thinning
d8 = sio.loadmat(SAVED + 'f8_results.mat', squeeze_me=True, struct_as_record=False)
Z8 = d8['zeta_830_full']
np.savez('data/fig8.npz', beta=d8['beta'][:15], rs2=d8['r_s_dim'][:15] / d8['p'].r2_0,
         thicken=Z8[:, :, 1], thin=Z8[:, :, 0])

# --- Appendix C: (almost) zero surface grain size -> beta independence
apc = dict(beta=np.array([0.5, 1, 2, 5, 10]), zgrid=np.linspace(0, 1, 201))
for tag, r in [('small', 1e-4), ('large', 0.1)]:
    phis, z8 = [], []
    for b in apc['beta']:
        M = FullModel(beta=b, rs2=r); M.run(t_end=40)
        P = M.final()
        phis.append(np.interp(apc['zgrid'], P['z'], P['phi'], right=np.nan)); z8.append(P['z830'])
    apc['phi_' + tag] = np.array(phis); apc['z830_' + tag] = np.array(z8)
np.savez('data/appC.npz', **apc)
print('appC z830 small', apc['z830_small'], 'large', apc['z830_large'])

# --- grain-growth time evolution for a parcel-tracking animation (beta=1 and 4)
gg = {}
for b in [1.0, 4.0]:
    M = FullModel(beta=b, rs2=0.01); M.run(t_end=20)
    P = M.final()
    for k in ['z', 'phi', 'w', 'r2', 'A']:
        gg[f'{k}_{b:g}'] = P[k]
np.savez('data/gg.npz', **gg)
print('done')
