"""Owen-Scalapino (mu*) bottleneck state of MgB2 after the strong THz pump.

Pair-breaking phonons and quasiparticles are in mutual equilibrium (Rothwarf-Taylor
bottleneck), so the excitation is stored as an excess quasiparticle density at the bath
temperature T = 7 K, described by a quasiparticle chemical potential mu*:
   f(E) = 1 / (exp((E - mu*)/kT) + 1),  E > 0.
Gaps follow from the two-band gap equations with tanh(E/2T) -> 1 - 2 f(E); the
conductivity from Mattis-Bardeen with the same distribution (finite-purity superfluid
correction as in twoband.py).  mu* is fixed by the absorbed energy density computed in
pump.py (no adjustable parameter); the energy per band is the quasiparticle energy
2 N_a int E f dxi with N_a from the Sommerfeld coefficient split by the DOS ratio.
"""
import numpy as np
from scipy.optimize import brentq
import twoband as tb
from twoband import THZ, KB, V, NR, GAM
import pump

T = 7.0
wc = tb.fit_wc()
KB_J = 1.380649e-23
N_tot = 3 * pump.GAMMA_N / (2 * np.pi ** 2 * KB_J ** 2)     # states / (J mol) per spin
Nw = NR / NR.sum()
meV = 1.602176634e-22


def f_mu(E, mu):
    return 1 / (np.exp(np.clip((E - mu) / (KB * T), -700, 700)) + 1)


def gaps_mu(mu, guess=(2.3, 7.0)):
    D = np.array(guess, float)
    xi = np.linspace(0, wc, 6001)[1:] - wc / 12000
    for _ in range(3000):
        I = []
        for a in range(2):
            E = np.sqrt(xi ** 2 + D[a] ** 2)
            I.append(np.sum((1 - 2 * f_mu(E, mu)) / E) * (wc / 6000))
        Dn = V @ (NR * np.array(I) * D)
        Dn = np.maximum(Dn, 0)
        if np.max(abs(Dn - D)) < 1e-9:
            return Dn
        D = 0.5 * D + 0.5 * Dn
    return D


def qp_energy(mu, D):
    """J/mol stored in quasiparticles above the 7 K thermal state."""
    xi = np.linspace(0, 30, 30001); dx = xi[1] - xi[0]
    out = 0.0
    for a in range(2):
        E = np.sqrt(xi ** 2 + D[a] ** 2)
        out += 2 * Nw[a] * N_tot * np.sum(E * meV * (f_mu(E, mu) - f_mu(E, -1e9 * 0 + 0) * 0)) * dx * meV
    return out


def sigma_mu(f, D, mu):
    old = tb.mb_band
    def mb(fr, Dg, Tt, n=6000):
        w = fr * THZ
        if Dg < 1e-6:
            return 1.0 + 0j
        ff = lambda e: f_mu(np.abs(e), mu) if np.isscalar(e) else f_mu(np.abs(e), mu)
        e = Dg + np.geomspace(1e-6, 60 * max(KB * Tt, Dg), n); de = np.gradient(e)
        g = (e ** 2 + Dg ** 2 + w * e) / (np.sqrt(e ** 2 - Dg ** 2) * np.sqrt((e + w) ** 2 - Dg ** 2))
        s1 = 2 / w * np.sum((ff(e) - ff(e + w)) * g * de)
        if w > 2 * Dg:
            e2 = np.linspace(Dg - w, -Dg, n + 2)[1:-1]; de2 = e2[1] - e2[0]
            g2 = (e2 ** 2 + Dg ** 2 + w * e2) / (np.sqrt(e2 ** 2 - Dg ** 2 + 0j).real * np.sqrt((e2 + w) ** 2 - Dg ** 2 + 0j).real)
            s1 += 1 / w * np.sum((1 - 2 * ff(e2 + w)) * np.abs(g2) * de2)
        lo = max(Dg - w, -Dg)
        e3 = np.linspace(lo, Dg, n + 2)[1:-1]; de3 = e3[1] - e3[0]
        num = e3 ** 2 + Dg ** 2 + w * e3
        den = np.sqrt(np.maximum(Dg ** 2 - e3 ** 2, 1e-30)) * np.sqrt(np.maximum((e3 + w) ** 2 - Dg ** 2, 1e-30))
        s2 = 1 / w * np.sum((1 - 2 * ff(e3 + w)) * num / den * de3)
        return s1 + 1j * s2
    tb.mb_band = mb
    try:
        return tb.sigma_total(f, D, T)
    finally:
        tb.mb_band = old


if __name__ == "__main__":
    U_abs = pump.absorbed(normal=False) * pump.VM          # J/mol, SC-state absorption
    D0 = gaps_mu(0.0); E0 = qp_energy(0.0, D0)
    def resid(mu):
        D = gaps_mu(mu); return qp_energy(mu, D) - E0 - U_abs
    mu = brentq(resid, 0.0, 6.0)
    D = gaps_mu(mu)
    print(f"absorbed {U_abs:.3f} J/mol -> mu* = {mu:.2f} meV; gaps 2Dpi = {2*D[0]/THZ:.2f} THz, 2Dsig = {2*D[1]/THZ:.2f} THz")
    for f in (0.5, 0.8, 1.0):
        s = sigma_mu(f, D, mu) / 100
        print(f"  f={f} THz: sigma1 = {s.real:.2e}, sigma2 = {s.imag:.2e} Ohm^-1 cm^-1")
    print("measured 2 ps after 150 kV/cm: sigma1(0.8) ~2.55e4, sigma2(0.5) ~2.6e4, sigma2(1.0) ~1.5e4")
