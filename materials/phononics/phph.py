"""Three-phonon scattering of LiNbO3 from the MACE cubic force constants (phono3py).

1. Linewidth of the driven A1(TO4) mode at Gamma vs temperature (tetrahedron method,
   meshes 6^3, 8^3, 10^3) -> convergence and the value used.
2. Lifetimes of all 30 branches on the irreducible 8^3 mesh at 300-1100 K: these are the
   relaxation times of the incoherent phonon bath in drive.py.
3. Decay channels of the coherent TO4 mode: every triplet (Gamma, q, -q) and band pair, with
   the Gaussian-smeared golden-rule weights, calibrated against phono3py's own imaginary
   self-energy at the same smearing, then rescaled to the converged tetrahedron linewidth.
   Final states are mapped onto the irreducible mesh classes used for the bath.
4. Real part of the bubble self-energy of TO4 vs temperature (frequency shift from the
   incoherent bath).

Writes phph.npz.  Run time about 1 min.
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
import phono3py
from phono3py.phonon3.imag_self_energy import ImagSelfEnergy

HERE = os.path.dirname(os.path.abspath(__file__))
NQ = 8
SIG = 0.20        # THz Gaussian width for the channel decomposition
TS = np.array([300., 500., 700., 900., 1100.])
THz_K = 47.99243
t0 = time.time()


def load(mesh, sym=True):
    ph3 = phono3py.load(os.path.join(HERE, "phono3py_params.yaml"),
                        fc2_filename=os.path.join(HERE, "fc2.hdf5"),
                        fc3_filename=os.path.join(HERE, "fc3.hdf5"),
                        produce_fc=False, is_mesh_symmetry=sym, log_level=0)
    ph3.mesh_numbers = [mesh] * 3
    return ph3


g = np.load(os.path.join(HERE, "anharmonic.npz"))
iIR = int(g["ia1"][-1])

# ---- 1. TO4 linewidth vs mesh (tetrahedron)
gam_mesh = {}
for mesh in (6, 8, 10):
    ph3 = load(mesh)
    ph3.init_phph_interaction()
    ise = ph3.run_imag_self_energy([0], TS, frequency_points_at_bands=True)
    gg = np.array(ise.gammas).reshape(len(TS), -1)      # one sigma, one grid point: (T, band)
    gam_mesh[mesh] = gg[:, iIR]
    print(f"mesh {mesh}: Gamma_TO4(T) THz = {np.round(gam_mesh[mesh], 4)}  "
          f"FWHM(300K) = {2*gam_mesh[mesh][0]*33.356:.1f} cm-1 ({time.time()-t0:.0f}s)", flush=True)
gam_to4 = gam_mesh[10]
ise_lt = ph3.run_imag_self_energy([0], [10.0], frequency_points_at_bands=True)
gam_to4_10K = float(np.array(ise_lt.gammas).reshape(-1)[iIR])
print(f"10 K: FWHM = {2*gam_to4_10K*33.356:.1f} cm-1", flush=True)

# ---- 2. bath lifetimes on the irreducible 8^3 mesh
ph3 = load(NQ)
ph3.init_phph_interaction()
ph3.run_thermal_conductivity(temperatures=TS, write_kappa=False, log_level=0)
tc = ph3.thermal_conductivity
ir_gp = np.array(tc.grid_points)
wq = np.array(tc.grid_weights)
fq = np.array(tc.frequencies)                    # (nirr, 30)
gq = np.array(tc.gamma)[0]                       # (T, nirr, 30)
assert ir_gp[0] == 0, "first irreducible point must be Gamma"
print(f"bath: {len(ir_gp)} irreducible q of {NQ**3}; kappa_zz(300K) = "
      f"{np.array(tc.kappa)[0,0,2]:.2f} W/mK ({time.time()-t0:.0f}s)", flush=True)

# ---- 3. decay channels of TO4 at Gamma (Gaussian smearing), symmetry-reduced triplets
itr = ph3.phph_interaction
itr.set_grid_point(0)
itr.run()
trips, wts, _, _ = itr.get_triplets_at_q()
freqs_bz, _, _ = itr.get_phonons()
P = itr.interaction_strength[:, iIR].astype(float)      # (ntrip, 30, 30)
f1 = freqs_bz[trips[:, 1]]; f2 = freqs_bz[trips[:, 2]]
w0 = freqs_bz[trips[0, 0], iIR]
from phonopy.phonon.grid import get_ir_grid_points
bz = itr.bz_grid
ir_g, ir_w, ir_map = get_ir_grid_points(bz)
assert np.array_equal(np.array(bz.grg2bzg[ir_g]), ir_gp), "irreducible grids differ"
cls_of = {int(gp): k for k, gp in enumerate(ir_g)}
bz2gr = np.array(bz.bzg2grg)
cls1 = np.array([cls_of[int(ir_map[bz2gr[b]])] for b in trips[:, 1]])
cls2 = np.array([cls_of[int(ir_map[bz2gr[b]])] for b in trips[:, 2]])


def nB(f, T):
    return np.where(f > 1e-3, 1.0 / np.expm1(np.maximum(f * THz_K / T, 1e-12)), 0.0)


def gauss(x):
    return np.exp(-0.5 * (x / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)


Dd = gauss(w0 - f1[:, :, None] - f2[:, None, :])            # decay   w0 -> w1 + w2
Dc = gauss(w0 + f1[:, :, None] - f2[:, None, :])            # w0 + w1 -> w2
# phono3py at the same smearing, 300 K, for the absolute prefactor
ise = ImagSelfEnergy(itr)
ise.set_grid_point(0); ise.set_sigma(SIG); ise.set_temperature(300.0)
ise.run_interaction(); ise.run_integration_weights(); ise.run()
gam_sig = np.array(ise.imag_self_energy)[iIR]
n1 = nB(f1, 300.); n2 = nB(f2, 300.)
Dc21 = gauss(w0 + f2[:, None, :] - f1[:, :, None])         # w0 + w2 -> w1
mine = np.sum(wts[:, None, None] * P * ((1 + n1[:, :, None] + n2[:, None, :]) * Dd
                                         + (n1[:, :, None] - n2[:, None, :]) * (Dc - Dc21)))
Cpref = gam_sig / mine
print(f"TO4 at sigma={SIG} THz: phono3py Gamma = {gam_sig:.4f} THz; prefactor {Cpref:.4e}; "
      f"tetrahedron (10^3) {gam_to4[0]:.4f} THz", flush=True)
# frequency-resolved decay spectrum at 300 K (for the figure): w1 of the lower partner
sel = Dd > 1e-6
print(f"decay channels with weight: {sel.sum()}")

# ---- 4. real self-energy (bubble) of TO4 vs T
shift = None
try:
    _, dlt = ph3.run_real_self_energy([0], TS, frequency_points_at_bands=True, epsilons=[0.1])
    shift = np.array(dlt).reshape(len(TS), -1)[:, iIR]
    print("TO4 bubble shift (THz) vs T:", np.round(shift, 4), flush=True)
except Exception as exc:
    print("real self-energy failed:", exc)

np.savez(os.path.join(HERE, "phph.npz"), TS=TS, gam_to4=gam_to4, gam_to4_10K=gam_to4_10K,
         gam_mesh6=gam_mesh[6], gam_mesh8=gam_mesh[8], gam_mesh10=gam_mesh[10],
         bath_freq=fq, bath_w=wq, bath_gamma=gq, w0=w0, sig=SIG,
         ch_w=wts, ch_P=P.astype(np.float32), ch_f1=f1, ch_f2=f2, ch_c1=cls1, ch_c2=cls2,
         Cpref=Cpref, gam_sig300=gam_sig, shift=shift if shift is not None else np.zeros(len(TS)))
print(f"saved ({time.time()-t0:.0f}s)")
