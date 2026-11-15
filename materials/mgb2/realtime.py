"""Real-time dynamics of THz-driven MgB2: two-band time-dependent Bogoliubov-de Gennes (G1 level)
with impurity-assisted light coupling, thin-film electrodynamics, and electron-phonon collision
integrals with explicit phonon populations.

Units: energies meV, times ps, fields V/m, vector potential V s/m.  Band 0 = pi, band 1 = sigma.

1. Electrons (one-body Nambu density matrix, per band a)
   Each band is represented in the exact-eigenstate (Anderson) basis of the disordered normal
   metal: N_a levels eps_n uniformly spaced on [-W_a, W_a], paired with their time-reversed
   partners.  The Nambu density matrix R_ij = <psi_j^+ psi_i>, psi_n = (c_n up, c_n down^+), obeys
       i hbar dR/dt = [H(t), R],
       H = [[eps - alpha_a(t) Y_a, -Delta_a], [-Delta_a^*, -eps - alpha_a(t) Y_a]],
       Delta_a = sum_b V_ab N_b delta_b sum_{|eps_n|<wc} R_b[n+N, n]^*.
   The pairing interaction is the phonon-mediated BCS mean field with the same couplings and
   cutoff as twoband.py (Liu-Mazin-Kortus lambda_ab, wc from Tc = 38 K).
   Light: the velocity operator between disorder eigenstates is a random Hermitian matrix
   Y = iK (K real antisymmetric, time-reversal symmetric disorder), <|Y_nm|^2> = L_a(eps_n-eps_m)/N_a,
   L_a(x) = gamma_a^2/(gamma_a^2 + x^2) (Drude width gamma_a = 1/tau_a of SM Table I).  Its strength
   alpha_a = c_a A with c_a^2 = sigma_n,a W_a / (pi hbar nu_a) reproduces the normal-state Drude
   conductivity sigma_n,a (Kubo formula; checked numerically by validate()).  The coupling enters
   as tau_0 (x) Y in Nambu space: it is the impurity-assisted paramagnetic coupling.  At first order
   it gives the Mattis-Bardeen absorption (case-II coherence factors); at second order the dirty-limit
   superfluid stiffness, Abrikosov-Gor'kov depairing by the supercurrent and the A^2 drive of the
   gap amplitude (Higgs mode).  The diamagnetic current -chi_dia,a A is fixed by gauge invariance
   of the normal state (no Meissner response above Tc): chi_dia,a = static paramagnetic
   susceptibility of the same basis in the normal state.
   Current: j = sum_a [ (c_a/V_a) Tr((1 x Y_a) R_a) - chi_dia,a A ].
2. Thin film on a substrate (field inside the 20 nm film, as in pump.py):
       (1 + n_s) E_f(t) = 2 E_in(t) - Z0 d j(t),   dA/dt = -E_f.
3. Collisions (Chang-Scalapino / Eliashberg kinetic equations, case-I coherence factors), evaluated
   in the instantaneous BCS eigenbasis of each level n (qp energy E_n = sqrt(eps_n^2 + |Delta_a|^2),
   occupation f_n = Tr(P_+ R_n)):
     scattering  n(a) -> m(b) + phonon Omega = E_n - E_m:
        J = (2 pi/hbar) d_a d_b nu^0_a a2F_ab(Omega) (1 + x_n x_m - y_n y_m)/2
            [f_n (1-f_m)(n_Omega+1) - f_m (1-f_n) n_Omega]
     recombination / pair breaking  n(a) + m(b) <-> phonon Omega = E_n + E_m:
        J = (2 pi/hbar) d_a d_b nu^0_a a2F_ab(Omega) (1 - x_n x_m + y_n y_m)/2
            [f_n f_m (n_Omega+1) - (1-f_n)(1-f_m) n_Omega]
   with x = eps/E, y = |Delta|/E, nu^0_a the bare DOS, a2F_ab the first-principles two-band
   Eliashberg functions (Golubov et al. 2002; a2f_golubov2002.csv).  df_n = -J/(nu_a d_a) with the
   dressed DOS nu_a = (1 + lambda_a) nu^0_a; each event changes the phonon population of the
   Omega bin by J/(D_ph Omega_bin) (Debye phonon DOS); phonons escape into the substrate,
   dn/dt = -(n - n_T)/tau_esc.  Energy is conserved exactly by the collision step.  The coherence of
   each level (pseudospin component transverse to its quasiparticle axis) decays at half the sum
   of its in- and out-scattering rates.
4. Probe: a weak 10 kV/cm single-cycle probe delayed by tau; sigma(w) is extracted from the
   odd part of the transmitted-field change, (E(+probe) - E(-probe))/2, by the thin-film formula,
   as in the experiment.  The quasi-stationary state is also evaluated with the nonequilibrium
   Mattis-Bardeen formula (computed f(E) and gaps).
5. After the probe window (A -> 0) the levels decouple and the run is continued to 100 ps as
   per-level Anderson pseudospins with the same collision integrals (pseudospin_continue).
Numerics: RK4, dt = 10 fs; collision step every 50 fs; N = 300 (pi, W = 15 meV) and 250 (sigma,
W = 18 meV) levels; single precision for the density matrices (checked against double precision
and against N = 450/375 and a second disorder realization: absorbed energy and sigma at 2 ps
within 2 %).  Level window: converged at 30 kV/cm (W = 30 meV changes the absorbed energy by 1 %);
at strong drive alpha = c A becomes comparable to W and the finite window adds a spurious
screening current, so the absorbed energy is underestimated: W = 30 meV gives +22 % at 47 kV/cm,
+33 % at 68 kV/cm and +82 % at 150 kV/cm (window_check.py).
"""
import os
import numpy as np
import twoband as tb
from twoband import THZ, KB, V as VPAIR, NR, GAM, Z0

HERE = os.path.dirname(os.path.abspath(__file__))
HB = 0.6582119569             # meV ps
MEV = 1.602176634e-22         # J per meV
HBAR_J = 1.054571817e-34
NS, DFILM = 3.1, 20e-9        # MgO index at THz, film thickness (SM)
VCELL = 2.906e-29             # m^3 per MgB2 cell (a = 3.086 A, c = 3.524 A)
# bare band DOS per spin (Golubov et al. 2002: N_pi = 0.410, N_sigma = 0.300 states/(eV cell), both spins)
NU0 = np.array([0.410, 0.300]) / 2 / 1e3 / VCELL          # states / (meV m^3) per spin
# Eliashberg couplings (Golubov et al. 2002), lambda_ab includes final-band DOS
LAM_E = np.array([[0.448, 0.155], [0.213, 1.017]])
ZREN = 1 + LAM_E.sum(1)                                     # 1 + lambda_a (mass renormalization)
NUD = NU0 * ZREN                                            # dressed DOS (gamma_n = 3.1 mJ/mol K^2)
THETA_D = 750.0               # K, Debye temperature (Bud'ko et al., PRL 86, 1877 (2001)), as pump.py
N_AT = 3 / VCELL
TAU_ESC = 20.0                # ps, phonon escape 20 nm MgB2 -> MgO (see PARAMS)
WC = tb.fit_wc()
_a2f = np.loadtxt(os.path.join(HERE, "a2f_golubov2002.csv"), delimiter=",")
A2F_W = _a2f[:, 0]
# A2F[a][b]: a = initial band, b = final band; columns: sigsig, pipi, sigpi, pisig
A2F = [[_a2f[:, 2], _a2f[:, 4]], [_a2f[:, 3], _a2f[:, 1]]]
SIGN = tb.sigma_n_SI()        # S/m per band (film value, twoband.py)
CDT = np.complex64 if os.environ.get("RT_SINGLE", "1") == "1" else np.complex128   # density-matrix precision
RDT = np.float32 if CDT == np.complex64 else np.float64

PARAMS = {
    "lambda_ab (pairing mean field)": "0.28, 0.96, 0.19; N_s/N_p = 0.73 (Liu, Mazin, Kortus PRL 87, 087005 (2001), SM of Katsumi et al.)",
    "wc": f"{WC:.2f} meV from Tc = 38 K",
    "1/tau_pi, 1/tau_sigma": "12.5, 8.33 meV (SM Table I, Fiore et al. PRB 106, 094515 (2022))",
    "sigma_n": "2.3e4 Ohm^-1 cm^-1 at 0.3 THz, 38 K (SM Fig. S1), band split by wp^2 tau",
    "Eliashberg a2F_ab, lambda_ab": "Golubov et al., J. Phys.: Condens. Matter 14, 1353 (2002), Fig. 1 (digitized)",
    "bare DOS": "N_pi = 0.410, N_sigma = 0.300 states/(eV cell) (Golubov et al. 2002)",
    "Debye temperature": "750 K (Bud'ko et al. PRL 86, 1877 (2001))",
    "phonon escape": "tau_esc = 4d/(eta u) = 20 ps for d = 20 nm (u = 7.5 km/s, eta = 0.5); "
                     "MgB2 bolometer fits give 6-42 ps (Chalmers HEB studies)",
    "film": "d = 20 nm, n_MgO = 3.1 (SM)",
    "pump": "single-cycle, exp(-(t/0.55 ps)^2) cos(2 pi 0.6 THz t + 0.3), as pump.py",
    "narrowband pulses": "0.63 THz, Gaussian field envelope FWHM 3 ps (not given in the paper; 1.5 ps tested)",
    "probe": "same single-cycle form, 10 kV/cm peak (SM Sec. IV)",
}


def fermi(E, T):
    return 1.0 / (np.exp(np.clip(E / (KB * T), -700, 700)) + 1.0)


def bose(W, T):
    return 1.0 / np.expm1(np.clip(W / (KB * T), 1e-12, 700))


def d_ph(W):
    """Debye phonon DOS per volume per meV."""
    return 9 * N_AT * W ** 2 / (KB * THETA_D) ** 3


_PN = None


def _pulse_norm():
    global _PN
    if _PN is None:
        t = np.linspace(-3, 3, 60001)
        _PN = np.max(np.abs(np.exp(-(t / 0.55) ** 2) * np.cos(2 * np.pi * 0.6 * t + 0.3)))
    return _PN


def broadband(Epk, t0=0.0):
    pn = _pulse_norm()
    return lambda t: Epk * np.exp(-((t - t0) / 0.55) ** 2) * np.cos(2 * np.pi * 0.6 * (t - t0) + 0.3) / pn


def narrowband(Epk, t0=0.0, f0=0.63, fwhm=3.0):
    """band-pass filtered multicycle pulse: Gaussian field envelope with the given FWHM (ps)."""
    s = fwhm / (2 * np.sqrt(np.log(2)))
    return lambda t: Epk * np.exp(-((t - t0) / s) ** 2) * np.cos(2 * np.pi * f0 * (t - t0))


class Band:
    def __init__(self, a, N, W, seed):
        self.a, self.N, self.W = a, N, W
        self.eps = -W + (np.arange(N) + 0.5) * 2 * W / N
        self.de = 2 * W / N
        rng = np.random.default_rng(seed)
        L = GAM[a] ** 2 / (GAM[a] ** 2 + (self.eps[:, None] - self.eps[None, :]) ** 2)
        # fixed magnitudes sqrt(L/N), random signs: exact Kubo normalization, random-matrix phases
        G = np.triu(rng.choice([-1.0, 1.0], size=(N, N)) * np.sqrt(L / N), 1)
        self.K = np.ascontiguousarray(G - G.T)
        self.Kc = self.K.astype(RDT)
        self.d = np.concatenate([self.eps, -self.eps]).astype(RDT)
        self.win = np.abs(self.eps) < WC
        nu_J = NUD[a] / MEV                                    # per J per m^3
        self.Vm = N / (2 * W * NUD[a])                         # m^3 represented
        cJ = np.sqrt(SIGN[a] * W * MEV / (np.pi * HBAR_J * nu_J))
        self.cJ = cJ
        self.c = cJ / MEV                                      # meV per (V s/m)
        self.chi_dia = self.chi_para_normal(7.0)

    def chi_para_normal(self, T):
        e = self.eps
        f = fermi(e, T)
        de = e[None, :] - e[:, None]
        df = f[:, None] - f[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(np.abs(de) > 1e-12, df / de, 0.0)
        return 2 * self.cJ ** 2 / self.Vm * np.sum(self.K ** 2 * r) / MEV   # A/(V s m)

    def thermal_R(self, D, T):
        N, e = self.N, self.eps
        E = np.sqrt(e ** 2 + abs(D) ** 2)
        f = fermi(E, T)
        R = np.zeros((2 * N, 2 * N), CDT)
        i = np.arange(N)
        R[i, i] = 0.5 + (f - 0.5) * e / E
        R[i + N, i + N] = 0.5 - (f - 0.5) * e / E
        R[i, i + N] = -(f - 0.5) * D / E                   # h_eh = -Delta
        R[i + N, i] = -(f - 0.5) * np.conj(D) / E
        return R

    def anomalous(self, R):
        N = self.N
        return np.sum(np.conj(np.diagonal(R[N:, :N]))[self.win]) * self.de

    def current(self, R):
        N = self.N
        S = R[:N, :N] + R[N:, N:]
        tr = 1j * np.sum(self.Kc * S.T, dtype=np.complex128)
        return (self.cJ / self.Vm) * tr.real

    def rhs(self, R, D, alpha):
        N = self.N
        M = self.d[:, None] * R
        M[:N] += CDT(-D) * R[N:]
        M[N:] += CDT(-np.conj(D)) * R[:N]
        if alpha != 0.0:
            Rv = R.view(RDT)
            KR = np.empty_like(R)
            KR[:N] = (self.Kc @ Rv[:N]).view(CDT)
            KR[N:] = (self.Kc @ Rv[N:]).view(CDT)
            M += CDT(-alpha * 1j) * KR
        M *= CDT(-1j / HB)
        M += M.conj().T
        return M


def gap_eq_grid(bands, T, guess=(2.3, 7.0)):
    D = np.array(guess, float)
    VN = VPAIR * NR[None, :]
    for _ in range(5000):
        F = []
        for b in bands:
            e = b.eps[b.win]
            E = np.sqrt(e ** 2 + D[b.a] ** 2)
            F.append(np.sum(D[b.a] * np.tanh(E / (2 * KB * T)) / (2 * E)) * b.de)
        Dn = VN @ np.array(F)
        if np.max(abs(Dn - D)) < 1e-11:
            return Dn
        D = 0.5 * D + 0.5 * Dn
    return D


class Phonons:
    def __init__(self, T0, dW=0.1, Wmax=50.0):
        self.dW = dW
        self.W = (np.arange(int(Wmax / dW)) + 0.5) * dW
        self.T0 = T0
        self.nT = bose(self.W, T0)
        self.n = self.nT.copy()
        self.D = d_ph(self.W)

    def at(self, Om):
        return np.interp(Om, self.W, self.n)

    def deposit(self, Om, J):
        """distribute events J (per m^3 per ps) at energies Om into bins, conserving number and energy."""
        x = Om / self.dW - 0.5
        i0 = np.clip(np.floor(x).astype(int), 0, len(self.W) - 2)
        w1 = np.clip(x - i0, 0, 1)
        nb = len(self.W)
        out = np.bincount(i0.ravel(), (J * (1 - w1)).ravel(), nb) + np.bincount(i0.ravel() + 1, (J * w1).ravel(), nb)
        return out / (self.D * self.dW)

    def energy(self):
        return np.sum((self.n - self.nT) * self.D * self.W) * self.dW * MEV       # J/m^3


def _a2f(a, b, Om):
    return np.interp(Om, A2F_W, A2F[a][b], left=0.0, right=0.0)


def kinetics(bands, fs, Ds, ph):
    """Chang-Scalapino collision integrals. Returns df/dt per band, (gain, loss) rates, dn/dt."""
    df = [np.zeros(b.N) for b in bands]
    gl = [np.zeros(b.N) for b in bands]
    dn = np.zeros_like(ph.n)
    Es = [np.sqrt(b.eps ** 2 + abs(Ds[b.a]) ** 2) for b in bands]
    for A in bands:
        for B in bands:
            a, b = A.a, B.a
            Ea, Eb = Es[a], Es[b]
            xa, xb = A.eps / Ea, B.eps / Eb
            ya, yb = abs(Ds[a]) / Ea, abs(Ds[b]) / Eb
            fa, fb = fs[a][:, None], fs[b][None, :]
            pre = 2 * np.pi / HB * A.de * B.de * NU0[a]
            # scattering, n (band a) above m (band b)
            Om = Ea[:, None] - Eb[None, :]
            msk = Om > 1e-9
            Omp = np.where(msk, Om, 1.0)
            g = np.where(msk, _a2f(a, b, Omp), 0.0) * 0.5 * (1 + xa[:, None] * xb[None, :] - ya[:, None] * yb[None, :]) * pre
            n = ph.at(Omp)
            out_ = g * fa * (1 - fb) * (n + 1)
            in_ = g * fb * (1 - fa) * n
            J = out_ - in_
            df[a] -= J.sum(1) / (NUD[a] * A.de)
            df[b] += J.sum(0) / (NUD[b] * B.de)
            # rates for coherence decay: level n loses (emission) or gains; level m likewise
            gl[a] += (g * ((1 - fb) * (n + 1) + fb * n)).sum(1) / (NUD[a] * A.de)
            gl[b] += (g * ((1 - fa) * n + fa * (n + 1))).sum(0) / (NUD[b] * B.de)
            dn += ph.deposit(Omp[msk], 2 * J[msk])
            # recombination / pair breaking
            Om2 = Ea[:, None] + Eb[None, :]
            g2 = _a2f(a, b, Om2) * 0.5 * (1 - xa[:, None] * xb[None, :] + ya[:, None] * yb[None, :]) * pre
            n2 = ph.at(Om2)
            J2 = g2 * (fa * fb * (n2 + 1) - (1 - fa) * (1 - fb) * n2)
            df[a] -= J2.sum(1) / (NUD[a] * A.de)
            gl[a] += (g2 * (fb * (n2 + 1) + (1 - fb) * n2)).sum(1) / (NUD[a] * A.de)
            dn += ph.deposit(Om2, J2)
    dn -= (ph.n - ph.nT) / TAU_ESC
    return df, gl, dn


def level_f(b, R, D):
    N = b.N
    i = np.arange(N)
    Ree, Rhh = R[i, i].real, R[i + N, i + N].real
    Reh, Rhe = R[i, i + N], R[i + N, i]
    E = np.sqrt(b.eps ** 2 + abs(D) ** 2)
    trh = b.eps * (Ree - Rhh) - D * Rhe - np.conj(D) * Reh
    return (Ree + Rhh) / 2 + trh.real / (2 * E), E


def apply_collisions(b, R, D, dfn, gam, dt):
    """R_n += df h/E (population change along the quasiparticle axis); the transverse part
    decays at half the sum of the in- and out-scattering rates."""
    N = b.N
    i = np.arange(N)
    E = np.sqrt(b.eps ** 2 + abs(D) ** 2)
    h11, h12, h21 = b.eps / E, -D / E, -np.conj(D) / E
    Ree, Rhh, Reh, Rhe = R[i, i], R[i + N, i + N], R[i, i + N], R[i + N, i]
    tr = Ree + Rhh
    s = (b.eps * (Ree - Rhh) - D * Rhe - np.conj(D) * Reh) / E       # Tr(h R)/E
    par = [tr / 2 + s / 2 * h11, tr / 2 - s / 2 * h11, s / 2 * h12, s / 2 * h21]
    k = np.exp(-0.5 * gam * dt)
    cur = [Ree, Rhh, Reh, Rhe]
    new = [p + k * (c - p) for p, c in zip(par, cur)]
    new[0] = new[0] + dfn * h11
    new[1] = new[1] - dfn * h11
    new[2] = new[2] + dfn * h12
    new[3] = new[3] + dfn * h21
    R[i, i], R[i + N, i + N], R[i, i + N], R[i + N, i] = new
    return R


class Sim:
    def __init__(self, T0=7.0, N=(400, 400), W=(15.0, 20.0), seed=1, collisions=True):
        self.bands = [Band(0, N[0], W[0], seed), Band(1, N[1], W[1], seed + 1)]
        self.T0 = T0
        self.D0 = gap_eq_grid(self.bands, T0)
        self.R = [b.thermal_R(self.D0[b.a], T0) for b in self.bands]
        self.A = 0.0
        self.t = 0.0
        self.ph = Phonons(T0)
        self.collisions = collisions
        self.ncoll = 5            # collision step every ncoll coherent steps
        self._k = 0
        self.Uabs = 0.0
        self.VN = VPAIR * NR[None, :]

    def gaps(self, R):
        F = np.array([b.anomalous(Rb) for b, Rb in zip(self.bands, R)])
        return self.VN @ F

    def current(self, R, A):
        return sum(b.current(Rb) - b.chi_dia * A for b, Rb in zip(self.bands, R))

    def field(self, R, A, Ein):
        j = self.current(R, A)
        return (2 * Ein - Z0 * DFILM * j) / (1 + NS), j

    def deriv(self, R, A, t, Ein_f):
        D = self.gaps(R)
        Ef, j = self.field(R, A, Ein_f(t))
        dR = [b.rhs(Rb, D[b.a], b.c * A) for b, Rb in zip(self.bands, R)]
        return dR, -Ef * 1e-12, Ef, j

    def step(self, dt, Ein_f):
        R, A, t = self.R, self.A, self.t
        k1, a1, E1, j1 = self.deriv(R, A, t, Ein_f)
        h = CDT(0.5 * dt); hf = CDT(dt)
        R2 = [r + h * k for r, k in zip(R, k1)]
        k2, a2, E2, j2 = self.deriv(R2, A + 0.5 * dt * a1, t + 0.5 * dt, Ein_f)
        R3 = [r + h * k for r, k in zip(R, k2)]
        k3, a3, E3, j3 = self.deriv(R3, A + 0.5 * dt * a2, t + 0.5 * dt, Ein_f)
        R4 = [r + hf * k for r, k in zip(R, k3)]
        k4, a4, E4, j4 = self.deriv(R4, A + dt * a3, t + dt, Ein_f)
        self.R = [r + CDT(dt / 6) * (p + 2 * q + 2 * s + u) for r, p, q, s, u in zip(R, k1, k2, k3, k4)]
        self.A = A + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
        # absorbed energy density (J/m^3): int j E dt (Simpson)
        self.Uabs += dt * 1e-12 / 6 * (j1 * E1 + 2 * j2 * E2 + 2 * j3 * E3 + j4 * E4)
        self.t = t + dt
        self._k += 1
        if self.collisions and self._k % self.ncoll == 0:
            self.collide(dt * self.ncoll)
        return E4

    def collide(self, dt):
        D = self.gaps(self.R)
        fs = [level_f(b, Rb, D[b.a])[0] for b, Rb in zip(self.bands, self.R)]
        df, gl, dn = kinetics(self.bands, fs, D, self.ph)
        self.R = [apply_collisions(b, Rb, D[b.a], df[b.a] * dt, gl[b.a], dt) for b, Rb in zip(self.bands, self.R)]
        self.ph.n = self.ph.n + dn * dt

    def snapshot(self):
        D = self.gaps(self.R)
        out = {"t": self.t, "D": D}
        for b, Rb in zip(self.bands, self.R):
            f, E = level_f(b, Rb, D[b.a])
            out[f"f{b.a}"] = f
            out[f"E{b.a}"] = E
        return out

    def qp_density(self):
        """excess quasiparticle density per spin (per m^3) above the initial thermal state."""
        D = self.gaps(self.R)
        out = []
        for b, Rb in zip(self.bands, self.R):
            f, E = level_f(b, Rb, D[b.a])
            out.append(np.sum(f) * b.de * NUD[b.a])
        return np.array(out)

    def copy(self):
        import copy
        s = copy.copy(self)
        s.R = [r.copy() for r in self.R]
        s.ph = copy.copy(self.ph)
        s.ph.n = self.ph.n.copy()
        return s


def run(sim, Ein_f, t1, dt=0.01, rec_every=1, callback=None):
    ts, Efs, Ds, As = [], [], [], []
    k = 0
    while sim.t < t1 - 1e-9:
        Ef = sim.step(dt, Ein_f)
        if k % rec_every == 0:
            ts.append(sim.t); Efs.append(Ef); As.append(sim.A)
            Ds.append(sim.gaps(sim.R))
            if callback:
                callback(sim)
        k += 1
    return np.array(ts), np.array(Efs), np.array(Ds), np.array(As)


# ---------------------------------------------------------------- continuation after the pulse
def level_blocks(sim):
    """per-level 2x2 Nambu blocks (Ree, Reh, Rhe, Rhh) of each band."""
    out = []
    for b, Rb in zip(sim.bands, sim.R):
        N = b.N; i = np.arange(N)
        out.append(np.array([Rb[i, i], Rb[i, i + N], Rb[i + N, i], Rb[i + N, i + N]], dtype=np.complex128))
    return out


def pseudospin_continue(sim, t1, dt=0.01, rec=0.5, ncoll=5):
    """After the pump (A -> 0) the levels decouple: each level is an Anderson pseudospin (its 2x2
    Nambu block) precessing in h_n = eps_n tau_3 - Re(Delta) tau_1 + Im(Delta) tau_2 with the
    self-consistent gaps, plus the same collision integrals and phonon kinetics.  The coherences
    between different levels, which only the light couples, are dropped."""
    bl = level_blocks(sim)
    VN = VPAIR * NR[None, :]
    bands = sim.bands

    def gaps(bl):
        F = np.array([np.sum(np.conj(x[2])[b.win]) * b.de for b, x in zip(bands, bl)])
        return VN @ F

    def rhs(bl):
        D = gaps(bl)
        out = []
        for b, x in zip(bands, bl):
            e = b.eps; d = D[b.a]
            ee, eh, he, hh = x
            # [h, R] for h = [[e, -d], [-d*, -e]]
            c_ee = -d * he + np.conj(d) * eh
            c_eh = 2 * e * eh + d * (ee - hh)
            c_he = -2 * e * he - np.conj(d) * (ee - hh)
            c_hh = -c_ee
            out.append((-1j / HB) * np.array([c_ee, c_eh, c_he, c_hh]))
        return out

    out = []
    t = sim.t; nxt = t; k = 0
    while t < t1 - 1e-9:
        if t >= nxt - 1e-9:
            D = gaps(bl)
            fs = []
            for b, x in zip(bands, bl):
                E = np.sqrt(b.eps ** 2 + abs(D[b.a]) ** 2)
                trh = b.eps * (x[0] - x[3]) - D[b.a] * x[2] - np.conj(D[b.a]) * x[1]
                fs.append(((x[0] + x[3]) / 2 + trh / (2 * E)).real)
            out.append((t, D.copy(), fs, sim.ph.n.copy()))
            nxt += rec
        k1 = rhs(bl); k2 = rhs([x + 0.5 * dt * y for x, y in zip(bl, k1)])
        k3 = rhs([x + 0.5 * dt * y for x, y in zip(bl, k2)]); k4 = rhs([x + dt * y for x, y in zip(bl, k3)])
        bl = [x + dt / 6 * (p + 2 * q + 2 * r + u) for x, p, q, r, u in zip(bl, k1, k2, k3, k4)]
        t += dt; k += 1
        if k % ncoll == 0:
            D = gaps(bl)
            fs = []
            for b, x in zip(bands, bl):
                E = np.sqrt(b.eps ** 2 + abs(D[b.a]) ** 2)
                trh = b.eps * (x[0] - x[3]) - D[b.a] * x[2] - np.conj(D[b.a]) * x[1]
                fs.append(((x[0] + x[3]) / 2 + trh / (2 * E)).real)
            df, gl, dn = kinetics(bands, fs, D, sim.ph)
            nb = []
            for b, x in zip(bands, bl):
                N = b.N; i = np.arange(N)
                R = np.zeros((2 * N, 2 * N), np.complex128)
                R[i, i], R[i, i + N], R[i + N, i], R[i + N, i + N] = x
                R = apply_collisions(b, R, D[b.a], df[b.a] * dt * ncoll, gl[b.a], dt * ncoll)
                nb.append(np.array([R[i, i], R[i, i + N], R[i + N, i], R[i + N, i + N]]))
            bl = nb
            sim.ph.n = sim.ph.n + dn * dt * ncoll
    return out


# ---------------------------------------------------------------- nonequilibrium Mattis-Bardeen
def mb_band_f(fr, D, fE, n=6000):
    """Mattis-Bardeen sigma/sigma_n with a nonequilibrium distribution f(E) (callable)."""
    w = fr * THZ
    if D < 1e-6:
        return 1.0 + 0j
    e = D + np.geomspace(1e-6, 25.0, n); de = np.gradient(e)
    g = (e ** 2 + D ** 2 + w * e) / (np.sqrt(e ** 2 - D ** 2) * np.sqrt((e + w) ** 2 - D ** 2))
    s1 = 2 / w * np.sum((fE(e) - fE(e + w)) * g * de)
    if w > 2 * D:
        e2 = np.linspace(D - w, -D, n + 2)[1:-1]; de2 = e2[1] - e2[0]
        g2 = (e2 ** 2 + D ** 2 + w * e2) / (np.sqrt(e2 ** 2 - D ** 2 + 0j).real * np.sqrt((e2 + w) ** 2 - D ** 2 + 0j).real)
        s1 += 1 / w * np.sum((1 - 2 * fE(e2 + w)) * np.abs(g2) * de2)
    lo = max(D - w, -D)
    e3 = np.linspace(lo, D, n + 2)[1:-1]; de3 = e3[1] - e3[0]
    num = e3 ** 2 + D ** 2 + w * e3
    den = np.sqrt(np.maximum(D ** 2 - e3 ** 2, 1e-30)) * np.sqrt(np.maximum((e3 + w) ** 2 - D ** 2, 1e-30))
    s2 = 1 / w * np.sum((1 - 2 * fE(e3 + w)) * num / den * de3)
    return s1 + 1j * s2


def f_of_E(E, f, D, T):
    """level occupations -> f(E) (average of the two levels +-eps with the same E); thermal tail above the grid."""
    o = np.argsort(E)
    Es, fs = E[o], f[o]
    # bin pairs of equal E
    Eb = np.linspace(abs(D), Es.max(), 200)
    idx = np.digitize(Es, Eb)
    Em, fm = [], []
    for k in np.unique(idx):
        s = idx == k
        Em.append(Es[s].mean()); fm.append(fs[s].mean())
    Em, fm = np.array(Em), np.clip(np.array(fm), 0, 1)
    return lambda e: np.where(e <= Em[-1], np.interp(e, Em, fm), fermi(e, T))


def sigma_mb_state(fr, snap, T0):
    sn = SIGN
    out = 0
    for a in range(2):
        D = abs(snap["D"][a])
        fE = f_of_E(snap[f"E{a}"], snap[f"f{a}"], D, T0)
        mb = mb_band_f(fr, D, fE)
        # finite-purity correction of the superfluid response as in twoband.sigma_total
        out = out + sn[a] * (mb.real + 1j * mb.imag * tb.purity_factor(D, T0, GAM[a]))
    return out


# ---------------------------------------------------------------- probe analysis
def sigma_from_probe(t, dE, Ein_probe, fmin=0.3, fmax=2.2):
    dt = t[1] - t[0]
    n = len(t)
    pad = 8 * n
    Ew = np.fft.rfft(dE, pad)
    Ei = np.fft.rfft(Ein_probe, pad)
    f = np.fft.rfftfreq(pad, dt)
    m = (f > fmin) & (f < fmax)
    sig = (2 * Ei[m] / Ew[m] - (1 + NS)) / (Z0 * DFILM)
    return f[m], sig            # S/m (e^{-i w t} convention handled below)


def probe_sigma(sim, Epump_f, t_probe, Eprobe=10e5, t_end_after=5.0, dt=0.01, callback=None):
    """Weak-probe conductivity of the evolving state.  Three branches from the current state of
    sim: pump only, pump + probe, pump - probe.  The probe response is the odd part
    dE = (E(+p) - E(-p))/2 (removes the even-order pump-probe mixing), tapered to the probe
    window [t_probe - 1.5, t_probe + 3] ps, and inverted with the thin-film formula.
    Returns f (THz), sigma (S/m, physical sign convention), the pump-only sim, and the traces."""
    s_m = sim.copy()
    s_p = sim.copy()
    pr = broadband(Eprobe, t_probe)
    t1 = t_probe + t_end_after
    tp, Ep, _, _ = run(sim, Epump_f, t1, dt, callback=callback)
    _, Eplus, _, _ = run(s_p, lambda t: Epump_f(t) + pr(t), t1, dt)
    _, Eminus, _, _ = run(s_m, lambda t: Epump_f(t) - pr(t), t1, dt)
    dE = 0.5 * (Eplus - Eminus)
    f, sig = sigma_from_probe(tp, dE * probe_window(tp, t_probe), pr(tp))
    return f, np.conj(sig), sim, (tp, Ep, dE)


def probe_window(t, t_probe, before=1.5, after=3.0, edge=0.5):
    x = np.ones_like(t)
    lo, hi = t_probe - before, t_probe + after
    m1 = t < lo + edge
    x[m1] = 0.5 * (1 - np.cos(np.pi * np.clip((t[m1] - lo) / edge, 0, 1)))
    m2 = t > hi - edge
    x[m2] = 0.5 * (1 + np.cos(np.pi * np.clip((t[m2] - hi + edge) / edge, 0, 1)))
    return x
