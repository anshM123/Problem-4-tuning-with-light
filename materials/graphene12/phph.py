"""Three-phonon interaction strengths on the full Nq x Nq grid (no symmetry reduction)
from the MACE cubic force constants, plus calibration of the absolute prefactor against
phono3py's own imaginary self-energy.

Output phph.npz:
  freq (Nq^2, 6) THz;  trip (Nq^2, Ntrip, 2) partner grid indices (q', q'');
  wt (Nq^2, Ntrip) triplet weights;  phi2 (Nq^2, Ntrip, 6, 6, 6) float32 |Phi|^2 (raw);
  C  prefactor so that the linearized decay rate of n_lambda (1/ps) is
       C * sum_t wt sum_{j'j''} phi2 [ (n'+n''+1) g(w-w'-w'') + (n'-n'')(g(w+w'-w'')-g(w-w'+w'')) ]
  with g a normalized Gaussian of width SIG (THz) in ordinary frequency.
"""
import os, sys, time
import numpy as np
import phono3py
from phono3py.phonon3.imag_self_energy import ImagSelfEnergy

HERE = os.path.dirname(os.path.abspath(__file__))
Nq = int(sys.argv[1]) if len(sys.argv) > 1 else 24
SIG = 0.25  # THz (~1 meV)

ph3 = phono3py.load(os.path.join(HERE, "phono3py_params.yaml"),
                    fc2_filename=os.path.join(HERE, "fc2.hdf5"),
                    fc3_filename=os.path.join(HERE, "fc3.hdf5"),
                    produce_fc=False, is_mesh_symmetry=False, log_level=0)
ph3.mesh_numbers = [Nq, Nq, 1]
ph3.init_phph_interaction()
itr = ph3.phph_interaction
grid = ph3.grid
itr.run_phonon_solver()
freqs_bz, _, _ = itr.get_phonons()
grg2bzg = grid.grg2bzg
bzg2grg = grid.bzg2grg
N = Nq * Nq
freq = freqs_bz[grg2bzg[:N]]

trip_all, wt_all, phi_all = [], [], []
t0 = time.time()
for g in range(N):
    itr.set_grid_point(grg2bzg[g])
    itr.run()
    trips, wts, _, _ = itr.get_triplets_at_q()
    trip_all.append(bzg2grg[trips[:, 1:]])
    wt_all.append(wts)
    phi_all.append(itr.interaction_strength.astype(np.float32))
    if g % 96 == 0:
        print(f"gp {g}/{N} ntrip={len(wts)} {time.time()-t0:.0f}s", flush=True)
nt = max(len(w) for w in wt_all)
trip = np.zeros((N, nt, 2), int); wt = np.zeros((N, nt)); phi2 = np.zeros((N, nt, 6, 6, 6), np.float32)
for g in range(N):
    k = len(wt_all[g]); trip[g, :k] = trip_all[g]; wt[g, :k] = wt_all[g]; phi2[g, :k] = phi_all[g]

# ---- calibration against phono3py at T = 300 K for a few grid points
THz_K = 47.99243  # h*1THz/kB
def nB(f, T):
    x = f * THz_K / T
    return np.where(f > 1e-3, 1 / np.expm1(np.maximum(x, 1e-12)), 0.0)

def gauss(x):
    return np.exp(-0.5 * (x / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)

ratios = []
for g in [N // 7 + 3, Nq * (Nq // 3) + Nq // 3, Nq // 2]:
    itr.set_grid_point(grg2bzg[g]); itr.run()
    ise = ImagSelfEnergy(itr)
    ise.set_grid_point(grg2bzg[g]); ise.set_sigma(SIG); ise.set_temperature(300.0)
    ise.run_interaction()
    ise.run_integration_weights()
    ise.run()
    gam = np.array(ise.imag_self_energy)
    w0 = freq[g]
    f1 = freq[trip[g, :len(wt_all[g]), 0]]; f2 = freq[trip[g, :len(wt_all[g]), 1]]
    n1 = nB(f1, 300.0); n2 = nB(f2, 300.0)
    P = phi2[g, :len(wt_all[g])].astype(float)
    mine = np.zeros(6)
    for j in range(6):
        d1 = gauss(w0[j] - f1[:, :, None] - f2[:, None, :])
        d2 = gauss(w0[j] + f1[:, :, None] - f2[:, None, :]) - gauss(w0[j] - f1[:, :, None] + f2[:, None, :])
        term = (n1[:, :, None] + n2[:, None, :] + 1) * d1 + (n1[:, :, None] - n2[:, None, :]) * d2
        mine[j] = np.sum(wt_all[g][:, None, None] * P[:, j] * term)
    sel = (gam > 1e-6) & (mine > 0)
    ratios += list(gam[sel] / mine[sel])
    print("gp", g, "phono3py gamma (THz):", np.round(gam, 5), flush=True)
ratios = np.array(ratios)
print("ratio gamma/mine: mean %.4e  std/mean %.3e" % (ratios.mean(), ratios.std() / ratios.mean()))
# decay rate of n (1/ps) = 2 * (2 pi gamma) = 4 pi gamma  (gamma in THz = half width in ordinary freq)
C = 4 * np.pi * ratios.mean()
np.savez(os.path.join(HERE, "phph.npz"), freq=freq, trip=trip, wt=wt, phi2=phi2, C=C, sig=SIG)
print("saved; C =", C, "total", time.time() - t0)
