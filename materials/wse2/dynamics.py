"""Exciton formation and intervalley scattering in monolayer WSe2 (90 K).

Incoherent excitons in three valley classes, resolved in center-of-mass kinetic energy:
  KK   (bright, K-K and K'-K')        g = 2,  mass M_KK
  KK'  (momentum dark, K-K', K'-K)    g = 2,  mass M_KK
  SK   (momentum dark, Sigma-K)       g = 12, mass M_SK
Energies: E_KK from the measured A resonance is the reference (0); the Sigma-K offset
follows from the conduction-band offset and the Rytova-Keldysh binding energies
(excitons.py).  Exciton-phonon scattering: Fermi golden rule with deformation potentials
(intervalley and intravalley optical: dispersionless phonons; intravalley acoustic:
small quanta), detailed balance at the lattice temperature, exact particle conservation.
Coherent excitation: resonant pump injects KK excitons at zero kinetic energy with the
measured 241 fs Gaussian envelope; above-gap pump injects hot excitons at 0.17 eV excess
energy in all valleys in proportion to their density of states.
Losses: exciton-exciton annihilation (bimolecular, coefficient gamma_EEA).
"""
import numpy as np
from scipy.integrate import solve_ivp
from excitons import exciton_table

HBAR = 6.582119569e-16      # eV s
KB = 8.617333262e-5
M0C2 = 0.51099895e6         # eV
HBARC = 1.97326980e-5       # eV cm
RHO = 6.09e-7 * 5.6096e32 / (2.99792458e10) ** 2  # g/cm^2 -> eV s^2/cm^4  (mass density)


def dos2d(M):
    """states per (eV cm^2) for one exciton species of mass M (units of m0)."""
    return M * M0C2 / (2 * np.pi * HBARC ** 2)


def build(params):
    p = dict(params)
    T = p["T"]; kT = KB * T
    de = p["de"]; eg = np.arange(0, p["emax"], de); nE = len(eg)
    sp = ["KK", "KKp", "SK"]
    g = dict(KK=2, KKp=2, SK=12)
    M = dict(KK=p["M_KK"], KKp=p["M_KK"], SK=p["M_SK"])
    E0 = dict(KK=0.0, KKp=p["dE_KKp"], SK=p["dE_SK"])
    # processes: (from, to, z, D [eV/cm], hOmega [eV])
    # processes (from, to, z, D [eV/cm], hOmega [eV]) with first-principles (DFPT)
    # deformation potentials of Jin, Li, Mullen, Kim, PRB 90, 045422 (2014), Tables II/IV:
    # electron K->Q (phonon at Q'): ac 8.2e8 @ 13.0 meV, op 0.8e8 @ 30.1 meV   (3 target valleys)
    # electron K->Q' (phonon at M):  ac 5.7e8 @ 15.8 meV, op 3.2e8 @ 30.0 meV  (3 target valleys)
    # electron K->K' (phonon at K):  ac 1.3e8 @ 16.8 meV, op 0.8e8 @ 29.7 meV
    # electron Q->Q' (Q3/Q5, M3/M4, K' phonons): ac 2.7/1.8/4.2e8 @ ~13-16 meV
    # intravalley optical: K 2.3e8, Q 2.7e8 @ 30.7 meV
    procs = []
    for (D, w) in ((p["D_KQ_ac"], 0.0130), (p["D_KQ_op"], 0.0301)):
        procs += [("KK", "SK", 3, D, w), ("SK", "KK", 0.5, D, w), ("KKp", "SK", 3, D, w), ("SK", "KKp", 0.5, D, w)]
    for (D, w) in ((p["D_KM_ac"], 0.0158), (p["D_KM_op"], 0.0300)):
        procs += [("KK", "SK", 3, D, w), ("SK", "KK", 0.5, D, w), ("KKp", "SK", 3, D, w), ("SK", "KKp", 0.5, D, w)]
    for (D, w) in ((1.3e8, 0.0168), (0.8e8, 0.0297)):
        procs += [("KK", "KKp", 1, D, w), ("KKp", "KK", 1, D, w)]
    for (D, w, z) in ((2.7e8, 0.0130, 2), (1.8e8, 0.0158, 2), (4.2e8, 0.0168, 1)):
        procs += [("SK", "SK", z, D, w)]
    procs += [("KK", "KK", 1, 2.3e8, 0.0307), ("KKp", "KKp", 1, 2.3e8, 0.0307), ("SK", "SK", 1, 2.7e8, 0.0307)]
    # golden-rule prefactor:  W = (2 pi/hbar) * hbar D^2/(2 rho Omega) * dos2d(M_f) * occupation
    W = []
    nB = lambda w: 1 / np.expm1(w / kT)
    for (a, b, z, D, w) in procs:
        Omega = w / HBAR
        pref = (2 * np.pi / HBAR) * HBAR * D ** 2 / (2 * RHO * Omega) * dos2d(M[b]) * z
        for sign, occ in ((+1, nB(w) + 1), (-1, nB(w))):   # emission (+1), absorption (-1)
            shift = E0[a] - E0[b] - sign * w                # eps_f = eps_i + shift
            W.append((sp.index(a), sp.index(b), pref * occ, int(round(shift / de))))
    # acoustic intravalley energy relaxation: quanta of hw_ac with quasi-elastic rate
    wac = p["w_ac"]
    for s in sp:
        # 2D quasi-elastic acoustic momentum-relaxation rate 1/tau = D^2 kT M / (hbar^3 rho v^2)
        Dac = p["Dac_K"] if s != "SK" else p["Dac_S"]
        erg = 1.602176634e-12
        Mg = M[s] * 9.1093837e-28
        rate = (Dac * erg) ** 2 * (kT * erg) * Mg / ((1.054571817e-27) ** 3 * p["rho_g"] * p["v_s"] ** 2)
        i = sp.index(s)
        W.append((i, i, rate * (nB(wac) + 1), int(round(-wac / de))))
        W.append((i, i, rate * nB(wac), int(round(wac / de))))
    return dict(p=p, eg=eg, nE=nE, sp=sp, g=g, M=M, E0=E0, W=W)


def rhs_factory(model, pump):
    nE = model["nE"]; W = model["W"]; p = model["p"]; eg = model["eg"]
    nsp = len(model["sp"])

    def rhs(t, y):
        N = y.reshape(nsp, nE)
        dN = np.zeros_like(N)
        for (a, b, r, sh) in W:
            if sh >= 0:
                src = N[a, :nE - sh] if sh > 0 else N[a]
                flow = r * src
                dN[a, :len(src)] -= flow
                dN[b, sh:sh + len(src)] += flow
            else:
                k = -sh
                src = N[a, k:]
                flow = r * src
                dN[a, k:] -= flow
                dN[b, :nE - k] += flow
        # exciton-exciton annihilation (bimolecular, total density)
        Ntot = N.sum()
        dN -= p["g_EEA"] * Ntot * N
        dN += pump(t)
        return dN.reshape(-1)
    return rhs


def run(params, mode="resonant", tmax=4e-12, nt=181):
    m = build(params)
    p = m["p"]; nE = m["nE"]; eg = m["eg"]; de = p["de"]
    sig = p["pulse_fwhm"] / (2 * np.sqrt(2 * np.log(2)))
    inj = np.zeros((3, nE))
    if mode == "resonant":
        inj[0, 0] = 1.0
    else:
        k = int(round(p["E_hot"] / de))
        w = np.array([m["g"][s] * m["M"][s] * (eg[k] + p["E_hot"] * 0 - m["E0"][s] > -1) for s in m["sp"]])
        for i, s in enumerate(m["sp"]):
            kk = int(round((p["E_hot"] - m["E0"][s]) / de))
            inj[i, max(kk, 0)] = m["g"][s] * m["M"][s]
        inj /= inj.sum()
    n_exc = p["n_exc"] if mode == "resonant" else p["n_exc_hot"]
    pump = lambda t: n_exc * inj * np.exp(-0.5 * (t / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)
    t_eval = np.linspace(-0.5e-12, tmax, nt)
    sol = solve_ivp(rhs_factory(m, pump), (-0.8e-12, tmax), np.zeros(3 * nE), t_eval=t_eval,
                    method="LSODA", rtol=1e-6, atol=1e3)
    N = sol.y.reshape(3, nE, -1)
    tot = N.sum(1)
    return dict(t=sol.t, KK=tot[0], KKp=tot[1], SK=tot[2], N=N, model=m)


PARAMS = dict(
    T=90.0, de=0.0005, emax=0.30,
    M_KK=0.65, M_SK=0.93,
    dE_KKp=0.0, dE_SK=None,       # set below from binding energies + CB offset
    dCB_SK=0.035,                  # Sigma - K conduction band offset (eV), DFT value
    D_KQ_ac=8.2e8, D_KQ_op=0.8e8, D_KM_ac=5.7e8, D_KM_op=3.2e8,   # Jin et al. 2014
    Dac_K=1.1, Dac_S=0.2,          # exciton acoustic D = |D_e - D_h| (eV), Jin et al. 2014
    rho_g=6.2e-7, v_s=3.3e5,       # g/cm^2, cm/s
    w_ac=0.002,                    # energy quantum of intravalley acoustic scattering (eV)
    g_EEA=0.08,                    # cm^2/s, exciton-exciton annihilation (see report)
    n_exc=8e12, n_exc_hot=2e12,    # cm^-2
    pulse_fwhm=241e-15, E_hot=0.17,
)

if __name__ == "__main__":
    t = exciton_table()
    P = dict(PARAMS)
    P["dE_SK"] = P["dCB_SK"] - (t["Eb_SK"] - t["Eb_KK"])
    P["M_KK"], P["M_SK"] = t["M_KK"], t["M_SK"]
    print("E(Sigma-K) - E(K-K) = %.1f meV" % (1e3 * P["dE_SK"]))
    r = run(P)
    K = r["KK"] + r["KKp"]; S = r["SK"]
    for tt in [0.08, 0.4, 0.7, 1.0, 1.3, 1.6, 2.1, 2.6, 3.1, 3.6]:
        i = np.argmin(abs(r["t"] - tt * 1e-12))
        print(f"t={tt:4.2f} ps  K={K[i]:.3e}  S={S[i]:.3e}  K/S={K[i]/max(S[i],1):.2f}")


def run_above_gap(params, tau_form=0.5e-12, Tc0=2000.0, tmax=4e-12, nt=181):
    """Above-gap (2.48 eV) excitation: hot free carriers first.
    Electrons occupy the 2 K and 6 Sigma conduction valleys (masses 0.39, 0.64; Jin et al.
    2014) in quasi-equilibrium at the carrier temperature Tc(t).  Tc relaxes to the
    lattice temperature by optical-phonon emission: the cooling time is the number of
    30.7 meV phonons needed to remove the excess energy times the golden-rule emission
    time for the intravalley optical deformation potential (2.3e8 eV/cm).  Carriers bind
    into excitons in the valley they occupy with formation time tau_form; excitons then
    follow the exciton kinetics of run().  ARPES signal per valley = excitons + carriers."""
    m = build(params)
    p = m["p"]; nE = m["nE"]; de = p["de"]; TL = p["T"]
    dCB = p["dCB_SK"]
    D, w = 2.3e8, 0.0307
    Omega = w / HBAR
    W_em = (2 * np.pi / HBAR) * HBAR * D ** 2 / (2 * RHO * Omega) * dos2d(0.39) * (1 / np.expm1(w / (KB * TL)) + 1)
    tau_c = (KB * (Tc0 - TL) / w) / W_em
    Tc = lambda t: TL + (Tc0 - TL) * np.exp(-max(t, 0.0) / tau_c)
    sig = p["pulse_fwhm"] / (2 * np.sqrt(2 * np.log(2)))
    gen = lambda t: p["n_exc_hot"] * np.exp(-0.5 * (t / sig) ** 2) / (np.sqrt(2 * np.pi) * sig)
    rhs_x = rhs_factory(m, lambda t: 0.0)

    def fK_of(T):
        wK = 2 * 0.39; wS = 6 * 0.64 * np.exp(-dCB / (KB * T))
        return wK / (wK + wS)

    def rhs(t, y):
        nfree = y[0]; N = y[1:]
        T = Tc(t); fK = fK_of(T)
        form = nfree / tau_form
        inj = np.zeros((3, nE))
        k = min(int(round(KB * T / de)), nE - 1)
        inj[0, k] = 0.5 * fK; inj[1, k] = 0.5 * fK; inj[2, k] = 1 - fK
        return np.concatenate([[gen(t) - form], rhs_x(t, N) + (form * inj).reshape(-1)])

    t_eval = np.linspace(-0.5e-12, tmax, nt)
    sol = solve_ivp(rhs, (-0.8e-12, tmax), np.zeros(1 + 3 * nE), t_eval=t_eval, method="LSODA",
                    rtol=1e-6, atol=1e3)
    nf = sol.y[0]
    N = sol.y[1:].reshape(3, nE, -1).sum(1)
    fK = np.array([fK_of(Tc(t)) for t in sol.t])
    return dict(t=sol.t, K=N[0] + N[1] + nf * fK, S=N[2] + nf * (1 - fK), nfree=nf,
                Tc=np.array([Tc(t) for t in sol.t]), tau_c=tau_c)
