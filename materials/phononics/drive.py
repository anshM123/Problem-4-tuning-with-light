"""Light-driven A1(TO4) phonon of LiNbO3: Maxwell + coherent modes + incoherent phonon bath.

One set of equations, solved together in time and depth (1D, normal incidence on an x-cut
crystal, pump field E along the polar c axis):

  Maxwell (Yee FDTD, total-field/scattered-field source in vacuum):
      eps0 eps_inf dE/dt = -dH/dx - dP/dt,    mu0 dH/dt = -dE/dx,
      P = (1/Omega) sum_j Zbar_j Q_j                     (DFPT mode polarities, params.py)
  Coherent A1 modes j = TO1..TO4 in every depth cell (MACE energy surface, anharmonic.py):
      Q_j'' = -dV(Q_1..Q_4)/dQ_j - 2 gamma_j Q_j' + Zbar_j E
      V = sum_j V_j(Q_j) + sum_{j<4} sum_{m,n} c^j_mn Q_4^m Q_j^n   (8th/6th-order fits)
      gamma_4(t, z) = 2 pi Gamma_4[n(t, z)]: golden-rule decay of the coherent TO4 mode into
      all phonon pairs of the 8^3 Brillouin-zone mesh (MACE fc3, phph.py), with the
      stimulated factors (1 + n1 + n2) and (n1 - n2) of the CURRENT bath populations.
      Bath-induced frequency shift: w_4^2 -> w_4^2 + 2 w_4 2pi [Delta(T*) - Delta(300 K)]
      (real part of the three-phonon bubble, phph.py).
  Incoherent phonons n_lambda(t, z) (all 30 branches, irreducible 8^3 mesh, 0.2 um depth bins):
      dn_lambda/dt = (pair creation by the coherent decay, channel by channel)
                     - [n_lambda - n0(T*)] / tau_lambda(T*)      (MACE three-phonon lifetimes)
      T*(t, z) fixed by energy conservation of the relaxation term.
  Energy: the power 2 gamma_4 Q_4'^2 removed from the coherent mode enters the bath channel
  by channel; the relaxation term conserves bath energy exactly; the field energy is
  tracked with the Poynting flux.  The budget is printed at the end of each run.

Usage: python drive.py            (all field / fluence scans, writes drive_out.pkl; about 16 min)
       python drive.py one E0 nu  (single run, prints observables)
"""
import os, sys, time
import numpy as np
from params import (E_CH, AMU, EPS0, C0, HBAR, KB, ANG, QSI, ZMODE, EPS_INF, PS0,
                    EXP_A, EXP_B, fluence_to_field)

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 1.0 / (EPS0 * C0 ** 2)
Z0 = MU0 * C0
eV = E_CH

an = np.load(os.path.join(HERE, "anharmonic.npz"))
ph = np.load(os.path.join(HERE, "phph.npz"))
gm = np.load(os.path.join(HERE, "gamma_modes.npz"))
OMEGA = abs(np.linalg.det(gm["cell"])) * ANG ** 3            # primitive cell volume (m^3)
ZBAR = an["zmode"]                                            # e/amu^1/2, signed, TO1..TO4
ZSI = ZBAR * E_CH / np.sqrt(AMU)                              # C/kg^1/2
P8 = an["p8"].copy(); TERMS = an["terms"]; COUP = an["coup"]; SELF = an["selfpot"].copy()
P8[-1] = 0.0; SELF[:, -1] = 0.0          # energy zero at the equilibrium structure
FA1 = an["freqA1"]                                            # THz
W4 = 2 * np.pi * FA1[3] * 1e12

# polynomial helpers: V in eV, Q in Å amu^1/2  ->  force in SI
FSI = eV / QSI      # eV/(Å amu^1/2) -> J/(m kg^1/2)
dP8 = np.polyder(P8)
dSELF = [np.polyder(s) for s in SELF]


CMAT = np.zeros((3, 6, 6))
for j in range(3):
    for (m, n), c in zip(TERMS, COUP[j]):
        CMAT[j, m, n] = c
CMAT_T = CMAT.transpose(0, 2, 1).copy()
MW = np.arange(6)[None, 1:, None]                   # derivative weights m (or n) = 1..5
dSELF_A = np.array([np.polyder(s) for s in SELF])   # (3, deg)


def _horner(c, x):
    y = np.zeros_like(x) + c[..., :1] if c.ndim > 1 else np.full_like(x, c[0])
    for k in range(1, c.shape[-1]):
        y = y * x + (c[..., k:k + 1] if c.ndim > 1 else c[k])
    return y


def forces(Q):
    """-dV/dQ_j for Q (4, N) in Å amu^1/2; returns (4, N) in eV/(Å amu^1/2)."""
    q4 = Q[3]; qj = Q[:3]
    p4 = np.empty((6,) + q4.shape); p4[0] = 1.0
    pj = np.empty((6,) + qj.shape); pj[0] = 1.0
    for k in range(1, 6):
        p4[k] = p4[k - 1] * q4
        pj[k] = pj[k - 1] * qj
    A = np.matmul(CMAT, pj.transpose(1, 0, 2))          # (3, 6, N): A[j,m] = sum_n C[j,m,n] qj^n
    Bn = np.matmul(CMAT_T, p4)                            # (3, 6, N): B[j,n] = sum_m C[j,m,n] q4^m
    F = np.empty_like(Q)
    F[3] = -_horner(dP8, q4) - np.sum(MW * A[:, 1:] * p4[None, :-1], axis=(0, 1))
    F[:3] = -_horner(dSELF_A, qj) - np.sum(MW * Bn[:, 1:] * pj[:-1].transpose(1, 0, 2), axis=1)
    return F


P8_MACE = P8.copy()


def set_ir_potential(kind="mace"):
    """Swap the Q_IR self-potential: 'mace' (default) or 'dft' (published DFT polynomial of
    von Hoegen et al., Methods, with the MACE harmonic frequency kept)."""
    global P8, dP8
    if kind == "mace":
        P8 = P8_MACE.copy()
    else:
        w2 = np.polyfit(an["qs"][np.abs(an["qs"]) <= 0.3], np.polyval(P8_MACE, an["qs"][np.abs(an["qs"]) <= 0.3]), 2)[0] * 2
        P8 = np.array([0.0071 / 5, 0.9008 / 4, 1.56765 / 3, w2 / 2, 0.0, 0.0])
    dP8 = np.polyder(P8)


SELF_MACE = SELF.copy(); CMAT_MACE = CMAT.copy()
EDEP = dict(Ed=0.0, Q0=1.0)          # optional phenomenological depolarizing field (diagnostic only)


def set_surface(kind="mace", Vlit=None, Ed_MVcm=0.0, Q0=1.0):
    """Choose the TO4-TO1 energy surface.
    'mace': the default (MACE frozen-phonon surface).
    'literature': Vlit(QP_dev, QIR_dev) in eV, a function of the deviations of the ferroelectric
    (TO1) and driven (TO4) coordinates from the equilibrium of that surface; it replaces the
    TO1 self-potential, the TO4 self-potential and the TO4-TO1 coupling.  TO2/TO3 stay MACE.
    Ed_MVcm adds the local depolarizing field E = -Ed q1/Q0 of Abalmasov's model."""
    global P8, dP8, SELF, dSELF_A, CMAT, CMAT_T
    SELF = SELF_MACE.copy(); CMAT = CMAT_MACE.copy(); P8 = P8_MACE.copy()
    EDEP["Ed"] = Ed_MVcm * 1e8; EDEP["Q0"] = Q0
    if kind != "mace":
        q1 = np.linspace(-7.0, 3.0, 81); q4 = np.linspace(-2.6, 2.6, 53)
        A, B = np.meshgrid(q1, q4, indexing="ij")
        V = Vlit(A, B)
        s1 = np.polyfit(q1, Vlit(q1, 0 * q1), 8); s1[-1] = 0.0
        s4 = np.polyfit(q4, Vlit(0 * q4, q4), 10); s4[-1] = 0.0
        R = V - np.polyval(s1, A) - np.polyval(s4, B)
        rows = np.array([(B ** m * A ** n).ravel() for m, n in TERMS]).T
        cf, *_ = np.linalg.lstsq(rows, R.ravel(), rcond=None)
        resid = rows @ cf - R.ravel()
        if np.abs(resid).max() > 1e-6 * max(np.abs(V).max(), 1.0):
            print(f"  literature surface fit residual {np.abs(resid).max():.2e} eV")
        SELF[0] = s1; P8 = s4
        CMAT[0] = 0.0
        for (m, n), c in zip(TERMS, cf):
            CMAT[0, m, n] = c
    CMAT_T = CMAT.transpose(0, 2, 1).copy()
    dSELF_A = np.array([np.polyder(x) for x in SELF])
    dP8 = np.polyder(P8)


def potential(Q):
    q4 = Q[3]
    p4 = np.vstack([q4 ** k for k in range(6)])
    V = np.polyval(P8, q4)
    for j in range(3):
        pj = np.vstack([Q[j] ** k for k in range(6)])
        V = V + np.polyval(SELF[j], Q[j]) + np.sum((CMAT[j].T @ p4) * pj, 0)
    return V


# ---- incoherent bath data (phph.py)
TS = ph["TS"]
bf = ph["bath_freq"]; bw = ph["bath_w"].astype(float); bgam = ph["bath_gamma"]   # (T, nirr, 30)
NQTOT = bw.sum()
nirr, nb = bf.shape
okb = bf > 0.3                                    # exclude acoustic modes at Gamma
hw = HBAR * 2 * np.pi * bf * 1e12                 # J
wcell = (bw[:, None] / NQTOT) * np.ones_like(bf)  # weight per cell of each class
tau_tab = np.where(okb[None], 1.0 / (4 * np.pi * np.maximum(bgam, 1e-6) * 1e12), 1.0)
# decay channels of the coherent TO4 mode
chw, chP = ph["ch_w"], ph["ch_P"].astype(float)
f1, f2 = ph["ch_f1"], ph["ch_f2"]; c1, c2 = ph["ch_c1"], ph["ch_c2"]
SIG = float(ph["sig"]); w0 = float(ph["w0"])
gs = lambda x: np.exp(-0.5 * (x / SIG) ** 2) / (np.sqrt(2 * np.pi) * SIG)
nt = len(chw)
I, J1, J2 = np.meshgrid(np.arange(nt), np.arange(nb), np.arange(nb), indexing="ij")
Dd = gs(w0 - f1[:, :, None] - f2[:, None, :])
Dc12 = gs(w0 + f1[:, :, None] - f2[:, None, :])
Dc21 = gs(w0 + f2[:, None, :] - f1[:, :, None])
base = chw[:, None, None] * chP
keep = (base * np.maximum(np.maximum(Dd, Dc12), Dc21)) > 1e-6 * (base * Dd).max()
K = dict(t=I[keep], a=J1[keep], b=J2[keep], Pd=(base * Dd)[keep], P12=(base * Dc12)[keep],
         P21=(base * Dc21)[keep])
K["ca"] = c1[K["t"]]; K["cb"] = c2[K["t"]]
K["fa"] = f1[K["t"], K["a"]]; K["fb"] = f2[K["t"], K["b"]]
# prefactor: phono3py value at sigma, rescaled to the converged tetrahedron linewidth at 300 K
CSCALE = float(ph["Cpref"]) * float(ph["gam_to4"][0]) / float(ph["gam_sig300"])
SHIFT = ph["shift"]                                  # THz, bubble shift of TO4 vs T
# Gamma-point linewidths of TO1..TO3 at 300 K (fixed damping, energy -> bath as heat)
ig = [int(b) for b in an["ia1"]]
gam_lo = np.array([bgam[0, 0, b] for b in ig[:3]])   # irreducible class 0 is Gamma


def nB(f_THz, T):
    x = HBAR * 2 * np.pi * f_THz * 1e12 / (KB * T)
    return 1.0 / np.expm1(np.maximum(x, 1e-9))


# sparse maps from channel flows to bath classes (class index = irreducible q * nb + band)
from scipy import sparse
nK = len(K["Pd"])
ia_ = K["ca"] * nb + K["a"]; ib_ = K["cb"] * nb + K["b"]
sa_ = K["fa"] / (K["fa"] + K["fb"]); ra_ = K["fa"] / w0; rb_ = K["fb"] / w0
_r = np.arange(nK)
def _sp(wa, wb):
    return sparse.csr_matrix((np.r_[wa, wb], (np.r_[_r, _r], np.r_[ia_, ib_])), shape=(nK, nirr * nb))
M_D = _sp(sa_, 1 - sa_)            # decay w0 -> wa + wb (energy split wa : wb)
M_12 = _sp(-ra_, 1 + ra_)          # w0 + wa -> wb
M_21 = _sp(1 + rb_, -rb_)          # w0 + wb -> wa
OK = okb.ravel(); HW = hw.ravel(); WC = wcell.ravel(); BF = np.maximum(bf.ravel(), 1e-3)


def gamma4(n):
    """Gamma_4 (THz) for bath populations n (nbins, ncls) and the channel rates."""
    na = n[:, ia_]; nbb = n[:, ib_]
    rd = K["Pd"] * (1 + na + nbb)
    r12 = K["P12"] * (na - nbb)
    r21 = K["P21"] * (nbb - na)
    tot = rd.sum(1) + r12.sum(1) + r21.sum(1)
    return CSCALE * tot, rd, r12, r21, tot


def tau_T(T):
    """Relaxation times (s) of all classes at bath temperatures T (nbins,)."""
    k = np.clip(np.searchsorted(TS, T) - 1, 0, len(TS) - 2)
    a = np.clip((T - TS[k]) / (TS[k + 1] - TS[k]), 0, 1)[:, None]
    tt = tau_tab.reshape(len(TS), -1)
    return (1 - a) * tt[k] + a * tt[k + 1]


def solve_T(n, tau, T):
    """Temperatures whose Bose populations make the relaxation term energy conserving."""
    target = np.sum(WC * HW * n / tau * OK, 1)
    T = T.copy()
    for _ in range(40):
        n0 = nB(BF[None], T[:, None])
        f = np.sum(WC * HW * n0 / tau * OK, 1) - target
        dn0 = n0 * (n0 + 1) * HW / (KB * T[:, None] ** 2)
        dT = -f / np.sum(WC * HW * dn0 / tau * OK, 1)
        T += dT
        if np.all(np.abs(dT) < 1e-7 * T):
            break
    return T


class Bath:
    """Incoherent phonon populations n (nbins, nirr*nb) per depth bin."""

    def __init__(self, nbins, T0=300.0):
        self.n = np.repeat((nB(BF, T0) * OK)[None], nbins, 0)
        self.T = np.full(nbins, T0)

    def energies(self):
        return np.sum(WC * HW * self.n * OK, 1)          # J per primitive cell

    def step(self, dt, P4, heat_lo):
        """Advance all bins by dt.  P4 (nbins,): power (J/s per cell) removed from the coherent
        TO4 mode, distributed over its decay channels; heat_lo (nbins,): energy (J per cell)
        dissipated by TO1..TO3, added with heat-capacity weights."""
        n = self.n
        E_target = self.energies() + P4 * dt + heat_lo
        G, rd, r12, r21, tot = gamma4(n)
        Eq = np.where(tot > 0, P4 * dt / np.where(tot > 0, tot, 1), 0.0)
        R = (M_D.T @ rd.T + M_12.T @ r12.T + M_21.T @ r21.T).T          # (nbins, ncls)
        n = n + np.where(OK, Eq[:, None] * R / (WC * HW), 0.0)
        tau = tau_T(self.T)
        n0T = nB(BF[None], self.T[:, None])
        cth = WC * HW * HW * n0T * (n0T + 1) * OK
        cth = cth / cth.sum(1, keepdims=True)
        n = n + np.where(OK, heat_lo[:, None] * cth / (WC * HW), 0.0)
        Ts = solve_T(n, tau, self.T)
        n0 = nB(BF[None], Ts[:, None])
        n = (n0 + (n - n0) * np.exp(-dt / tau)) * OK
        resid = E_target - np.sum(WC * HW * n * OK, 1)     # exact energy conservation
        n = n + np.where(OK, resid[:, None] * cth / (WC * HW), 0.0)
        self.n = n
        self.T = Ts
        return G


def run(E0_MVcm, nu_THz, fwhm_fs=150.0, t_end_ps=1.4, record=True, quiet=False, zbar_scale=1.0,
        eps_inf=EPS_INF, bath_on=True, cep=0.0, record_E=False):
    dx = 20e-9
    Lvac, Lcry, Labs = 3e-6, 30e-6, 12e-6
    N = int((Lvac + Lcry) / dx)
    x = np.arange(N) * dx
    ks = 20                                   # TF/SF boundary (E node index)
    kprobe = 8                                # reflected-field probe (scattered region)
    ksurf = int(Lvac / dx)
    cry = np.arange(N) >= ksurf
    nc = int(cry.sum())
    zc = x[cry] - x[ksurf]
    dt = 0.5 * dx / C0
    mur = (C0 * dt - dx) / (C0 * dt + dx)
    t0 = -0.5e-12
    nsteps = int((t_end_ps * 1e-12 - t0) / dt)
    # absorbing conductivity in the deep crystal
    sig = np.zeros(N)
    xa = Lvac + Lcry - Labs
    m = x > xa
    sig[m] = 3e5 * ((x[m] - xa) / Labs) ** 2           # S/m
    epsr = np.where(cry, eps_inf, 1.0)
    ca = (1 - sig * dt / (2 * EPS0 * epsr)) / (1 + sig * dt / (2 * EPS0 * epsr))
    cb = (dt / (EPS0 * epsr * dx)) / (1 + sig * dt / (2 * EPS0 * epsr))
    E = np.zeros(N); H = np.zeros(N)
    Q = np.zeros((4, nc)); V = np.zeros((4, nc))      # Å amu^1/2 and Å amu^1/2 per s
    zs = ZSI * zbar_scale
    E0 = E0_MVcm * 1e8
    om = 2 * np.pi * nu_THz * 1e12
    tau_env = fwhm_fs * 1e-15 / np.sqrt(2 * np.log(2))   # field envelope exp(-t^2/tau^2)... intensity FWHM

    def Einc(xx, t):
        tt = t - (xx - x[ks]) / C0
        return E0 * np.sin(om * tt + cep) * np.exp(-2 * np.log(2) * tt ** 2 / (fwhm_fs * 1e-15) ** 2)

    # bath bins
    dzb = 0.2e-6
    nbins = 30
    bidx = np.minimum((zc / dzb).astype(int), nbins)       # bin nbins = "deep" (not tracked)
    bath = Bath(nbins)
    g4_bins = np.full(nbins + 1, float(ph["gam_to4"][0]))  # THz (deep bin: 300 K value)
    dw2_bins = np.zeros(nbins + 1)
    nsub = 90                                              # bath step = 90 FDTD steps (3 fs)
    acc_P4 = np.zeros(nbins + 1); acc_lo = np.zeros(nbins + 1)
    Pold = (zs[:, None] * (Q * QSI)).sum(0) / OMEGA
    # records
    rec_every = 10
    nrec = nsteps // rec_every + 1
    zr = zc < 2.0e-6
    rec_t = np.zeros(nrec); rec_V4 = np.zeros((nrec, zr.sum())); rec_dP = np.zeros((nrec, zr.sum()))
    rec_Q4s = np.zeros(nrec); rec_Qs = np.zeros((nrec, 4)); rec_Ts = np.zeros((nrec, 4)); rec_G4 = np.zeros((nrec, 3))
    rec_Ecoh = np.zeros((nrec, 3))
    zE = zc < 18e-6
    rec_E = np.zeros((nrec, zE.sum())) if record_E else None
    Uin = 0.0; Uref = 0.0; Uabs = 0.0; Udiss = 0.0; Udiss_tr = 0.0
    irec = 0
    tstart = time.time()
    for it in range(nsteps):
        t = t0 + it * dt
        # --- H update (half step), TF/SF correction
        H[:-1] += -(dt / (MU0 * dx)) * (E[1:] - E[:-1])
        H[ks - 1] += (dt / (MU0 * dx)) * Einc(x[ks], t)
        # --- mechanics: velocity Verlet-like leapfrog with semi-implicit damping
        Ec = E[cry]
        if EDEP["Ed"] != 0.0:
            Ec = Ec - EDEP["Ed"] * Q[0] / EDEP["Q0"]
        F = forces(Q) * FSI + zs[:, None] * Ec[None, :]                  # SI force
        gam = np.empty((4, nc))
        gam[:3] = (2 * np.pi * gam_lo * 1e12)[:, None]
        gam[3] = 2 * np.pi * g4_bins[bidx] * 1e12
        # bath-induced softening/stiffening of TO4
        F[3] += -(dw2_bins[bidx]) * Q[3] * QSI
        VSI = V * QSI
        VSI_new = ((1 - gam * dt) * VSI + dt * F) / (1 + gam * dt)
        Vmid = 0.5 * (VSI + VSI_new)
        pdiss = 2 * gam * Vmid ** 2                                         # W per (m^3/OMEGA) -> J/s per cell
        V = VSI_new / QSI
        Q = Q + dt * V
        # accumulate dissipation per bath bin
        acc_P4 += np.bincount(bidx, weights=pdiss[3], minlength=nbins + 1) * dt
        acc_lo += np.bincount(bidx, weights=pdiss[:3].sum(0), minlength=nbins + 1) * dt
        Udiss += (pdiss.sum() * dt) * dx / OMEGA
        Pnew = (zs[:, None] * (Q * QSI)).sum(0) / OMEGA
        # --- E update
        Jp = np.zeros(N); Jp[cry] = (Pnew - Pold) / dt
        Pold = Pnew
        E1_old, E0_old = E[1], E[0]
        E[1:] = ca[1:] * E[1:] - cb[1:] * (H[1:] - H[:-1]) - cb[1:] * dx * Jp[1:]
        E[ks] += (dt / (EPS0 * dx)) * Einc(x[ks] - 0.5 * dx, t + 0.5 * dt) / Z0
        E[0] = E1_old + mur * (E[1] - E0_old)          # first-order Mur absorbing boundary
        Uabs += np.sum(sig * E ** 2) * dt * dx
        Uin += Einc(x[ks], t) ** 2 / Z0 * dt
        Uref += E[kprobe] ** 2 / Z0 * dt
        # --- bath update every nsub steps
        if (it + 1) % nsub == 0:
            dtb = nsub * dt
            ncell = np.maximum(np.bincount(bidx, minlength=nbins + 1)[:nbins], 1)
            P4c = acc_P4[:nbins] / ncell / dtb               # J/s per cell
            hlo = acc_lo[:nbins] / ncell
            Udiss_tr += (acc_P4[:nbins].sum() + acc_lo[:nbins].sum()) * dx / OMEGA
            if bath_on:
                g4_bins[:nbins] = bath.step(dtb, P4c, hlo)
                dw2_bins[:nbins] = 2 * W4 * 2 * np.pi * 1e12 * (np.interp(bath.T, TS, SHIFT) - SHIFT[0])
            # energy deposited in bins: track in bath (bath energy includes it)
            acc_P4[:] = 0.0; acc_lo[:] = 0.0
        if record and it % rec_every == 0:
            rec_t[irec] = t
            if record_E:
                rec_E[irec] = E[cry][zE]
            rec_V4[irec] = V[3][zr]
            rec_dP[irec] = Pnew[zr]
            rec_Q4s[irec] = Q[3][0]
            rec_Qs[irec] = Q[:, 0]
            rec_Ts[irec] = [bath.T[0], bath.T[3], bath.T[6], bath.T[15]]
            rec_G4[irec] = [g4_bins[0], g4_bins[6], g4_bins[15]]
            Ek = 0.5 * V[:, :zr.sum()] ** 2 * (QSI ** 2) / eV
            Ep = potential(Q[:, :zr.sum()])
            rec_Ecoh[irec] = [(Ek[:, 0].sum() + Ep[0]), (Ek.sum(0) + Ep)[:65].mean(), Ek[3, 0] + np.polyval(P8, Q[3, 0])]
            irec += 1
    # ---- energy budget (J per m^2)
    Ufield = np.sum(0.5 * EPS0 * epsr * E ** 2 + 0.5 * MU0 * H ** 2) * dx
    Umech = np.sum((0.5 * (V * QSI) ** 2).sum(0) + potential(Q) * eV) * dx / OMEGA
    ncell = np.bincount(bidx, minlength=nbins + 1)[:nbins]
    Ebath_ex = np.sum((bath.energies() - Bath(1).energies()[0]) * ncell) * dx / OMEGA
    budget = dict(Uin=Uin, Uref=Uref, Uabs=Uabs, Ufield=Ufield, Umech=Umech, Udiss=Udiss,
                  Ebath_ex=Ebath_ex, Udiss_tr=Udiss_tr)
    if not quiet:
        bal = (Uref + Uabs + Ufield + Umech + Udiss) / Uin
        print(f"E0={E0_MVcm:.1f} MV/cm nu={nu_THz} THz: R={Uref/Uin:.3f} absorbed(deep)={Uabs/Uin:.3f} "
              f"mech={Umech/Uin:.3f} diss={Udiss/Uin:.3f} field={Ufield/Uin:.3f} -> sum {bal:.4f}; "
              f"bath excess / dissipation in tracked bins {Ebath_ex/max(Udiss_tr,1e-30):.6f} "
              f"({time.time()-tstart:.0f}s)", flush=True)
    return dict(t=rec_t[:irec], V4=rec_V4[:irec], dP=rec_dP[:irec], z=zc[zr], Q4s=rec_Q4s[:irec], Qs=rec_Qs[:irec],
                Ts=rec_Ts[:irec], G4=rec_G4[:irec], Ecoh=rec_Ecoh[:irec], budget=budget, dt=dt * rec_every,
                E=rec_E[:irec] if record_E else None, zE=zc[zE])


# ---------------------------------------------------------------- observables
def sh_signal(r, ng=3.03):
    """SH probe signal: integral of dQ4/dt along the 400 nm light line over the first 2 um."""
    t = r["t"]; z = r["z"]; V4 = r["V4"]
    S = np.zeros_like(t)
    for k, zz in enumerate(z):
        S += np.interp(t + zz * ng / C0, t, V4[:, k], left=0, right=0)
    return S * (z[1] - z[0])


def spectrum(t, s, tmin=0.0e-12, tmax=1.2e-12, pad=16):
    m = (t >= tmin) & (t <= tmax)
    y = s[m] * np.hanning(m.sum())
    nf = pad * m.sum()
    F = np.abs(np.fft.rfft(y, nf))
    f = np.fft.rfftfreq(nf, t[1] - t[0]) / 1e12
    return f, F


def peak(f, F, lo, hi):
    m = (f > lo) & (f < hi)
    i = np.argmax(F[m]); ff = f[m]; FF = F[m]
    if 0 < i < len(ff) - 1:
        a, b, c = FF[i - 1], FF[i], FF[i + 1]
        d = 0.5 * (a - c) / (a - 2 * b + c)
        return ff[i] + d * (ff[1] - ff[0]), b
    return ff[i], FF[i]


def polarization_peak(r, lc=1.27e-6, ng=2.30, probe_fs=35.0, sgnP=1.0):
    t = r["t"]; z = r["z"]; dP = r["dP"]
    sel = z < lc
    S = np.zeros_like(t)
    for k in np.where(sel)[0]:
        S += np.interp(t + z[k] * ng / C0, t, dP[:, k], left=0, right=0)
    S /= sel.sum()
    s = probe_fs * 1e-15 / (2 * np.sqrt(2 * np.log(2)))
    ker = np.exp(-0.5 * ((np.arange(-200, 201) * (t[1] - t[0])) / s) ** 2); ker /= ker.sum()
    Ss = np.convolve(S, ker, mode="same")
    ratio = 1 + sgnP * Ss / PS0
    return ratio, t


def rect(r, win_fs=110.0):
    """Largest cycle-averaged (rectified) displacement of TO1..TO4 at the surface."""
    t = r["t"]; Qs = r["Qs"]
    n = max(int(win_fs * 1e-15 / (t[1] - t[0])), 1)
    ker = np.ones(n) / n
    sm = np.array([np.convolve(Qs[:, j], ker, mode="same") for j in range(4)])
    sm = sm[:, n:-n]
    return sm[np.arange(4), np.argmax(np.abs(sm), 1)]          # signed, largest magnitude


def analyse(r, nu):
    S = sh_signal(r)
    f, F = spectrum(r["t"], S)
    f1, A1 = peak(f, F, 10, 24)
    f2, A2 = peak(f, F, 2 * f1 - 4, 2 * f1 + 4)
    f3, A3 = peak(f, F, 3 * f1 - 5, 3 * f1 + 5)
    q4 = r["Q4s"]
    return dict(f1=f1, A1=A1, A2=A2, A3=A3, qmin=q4.min(), qmax=q4.max(), S=S, f=f, F=F)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "one":
        E0 = float(sys.argv[2]); nu = float(sys.argv[3])
        r = run(E0, nu)
        a = analyse(r, nu)
        print({k: v for k, v in a.items() if np.ndim(v) == 0})
        sys.exit()

    t0 = time.time()
    res = {}
    sgnP = float(an["sgnP"])
    # A. von Hoegen et al.: 17.5 THz, field scan
    runsA = []
    for E0 in EXP_A["fields"]:
        r = run(E0, EXP_A["nu"])
        a = analyse(r, EXP_A["nu"])
        ratio, tt = polarization_peak(r, sgnP=sgnP)
        Ecoh = r["Ecoh"]
        runsA.append(dict(E0=E0, f1=a["f1"], A1=a["A1"], A2=a["A2"], A3=a["A3"], qmin=a["qmin"],
                          qmax=a["qmax"], Ecoh_max_surf=Ecoh[:, 0].max(), Ecoh_max_layer=Ecoh[:, 1].max(),
                          U4max=Ecoh[:, 2].max(), Tmax=r["Ts"][:, 0].max(), G4max=r["G4"][:, 0].max(),
                          Pmin=ratio.min(), budget=r["budget"],
                          Q1rect=rect(r)))
        if E0 in (2.0, 20.0):
            res[f"A_{E0:.0f}"] = dict(t=r["t"], S=a["S"], f=a["f"], F=a["F"], Q4s=r["Q4s"], Qs=r["Qs"], Ts=r["Ts"],
                                      G4=r["G4"], Ecoh=r["Ecoh"], ratio=ratio)
        print(f"  A: E0={E0:5.1f}  f1={a['f1']:.3f} THz  Q4 in [{a['qmin']:.3f}, {a['qmax']:.3f}]  "
              f"E_coh(surface) max {Ecoh[:,0].max():.3f} eV  T_bath(surface) max {r['Ts'][:,0].max():.0f} K  "
              f"Gamma4 max {r['G4'][:,0].max():.3f} THz  ({time.time()-t0:.0f}s)", flush=True)
    # sensitivity runs at 20 MV/cm
    sens = {}
    for name, kw in [("Zbar_exp_strength", dict(zbar_scale=np.sqrt(2.55 / 3.07))),
                     ("eps_inf_DFPT", dict(eps_inf=5.8)),
                     ("no_bath", dict(bath_on=False)),
                     ("same_detuning", dict(nu_THz=EXP_A["nu"] * FA1[3] / 19.0)),
                     ("cep_pi2", dict(cep=np.pi / 2))]:
        kw2 = dict(E0_MVcm=20.0, nu_THz=EXP_A["nu"]); kw2.update(kw)
        r = run(**kw2)
        a = analyse(r, kw2["nu_THz"])
        sens[name] = dict(f1=a["f1"], qmin=a["qmin"], qmax=a["qmax"], Ecoh=r["Ecoh"][:, 0].max())
        print(f"  sens {name}: f1={a['f1']:.3f} Q4 [{a['qmin']:.3f},{a['qmax']:.3f}] "
              f"Ecoh {r['Ecoh'][:,0].max():.3f} eV", flush=True)
    # variant: Q_IR self-potential of the published DFT fit (diagnostic, not a fit to data)
    runsA_dft = []
    set_ir_potential("dft")
    for E0 in (2.0, 5.7, 11.6, 17.6, 20.9):
        r = run(E0, EXP_A["nu"], quiet=True)
        a = analyse(r, EXP_A["nu"])
        runsA_dft.append(dict(E0=E0, f1=a["f1"], qmin=a["qmin"], qmax=a["qmax"],
                              Ecoh_max_surf=r["Ecoh"][:, 0].max()))
        print(f"  A(DFT potential): E0={E0:5.1f} f1={a['f1']:.3f} Q4 [{a['qmin']:.3f}, {a['qmax']:.3f}] "
              f"E_coh {r['Ecoh'][:,0].max():.3f} eV ({time.time()-t0:.0f}s)", flush=True)
    runsB_dft = []
    for F in (12.8, 28.7, 51.8, 94.5):
        E0 = fluence_to_field(F)
        r = run(E0, 19.0, quiet=True)
        ratio, tt = polarization_peak(r, sgnP=sgnP)
        runsB_dft.append(dict(nu=19.0, F=F, E0=E0, Pmin=ratio.min()))
        print(f"  B(DFT potential): F={F:6.1f}: (P/P0) min {ratio.min():.3f}", flush=True)
    set_ir_potential("mace")
    # B. Mankowsky et al.: 19 THz fluence scan (+21, 22 THz)
    runsB = []
    for nu, Fl in [(19.0, EXP_B["fluences"]), (21.0, [24.4, 52.4, 90.2]), (22.0, [49.4, 103.7, 173.8])]:
        for F in Fl:
            E0 = fluence_to_field(F)
            r = run(E0, nu, quiet=True)
            ratio, tt = polarization_peak(r, sgnP=sgnP)
            q4 = r["Q4s"]
            runsB.append(dict(nu=nu, F=F, E0=E0, Pmin=ratio.min(), Pmax=ratio.max(),
                              tmin=tt[np.argmin(ratio)], qmin=q4.min(), qmax=q4.max()))
            if nu == 19.0 and F in (12.8, 94.5):
                res[f"B_{F:.0f}"] = dict(t=tt, ratio=ratio)
            print(f"  B: nu={nu} F={F:6.1f} mJ/cm2 (E0={E0:.1f} MV/cm): (P/P0) min {ratio.min():.3f} "
                  f"max {ratio.max():.3f} at {tt[np.argmin(ratio)]*1e12:.2f} ps  ({time.time()-t0:.0f}s)", flush=True)
    import pickle
    with open(os.path.join(HERE, "drive_out.pkl"), "wb") as fh:
        pickle.dump(dict(runsA=runsA, runsB=runsB, sens=sens, traces=res, sgnP=sgnP,
                         runsA_dft=runsA_dft, runsB_dft=runsB_dft), fh)
    print(f"done ({time.time()-t0:.0f}s)")
