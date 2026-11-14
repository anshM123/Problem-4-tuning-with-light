"""Two-band superconductivity of MgB2: equilibrium gaps, linear THz conductivity,
amplitude-mode (Higgs) mediated third-order response, and the strongly pumped state.

Energies in meV, frequencies f in THz (1 THz = 4.135667 meV), temperatures in K.

Gap equations (pi = band 1, sigma = band 2; Liu-Mazin-Kortus couplings as used by
Katsumi et al.):  Delta_a = sum_b V_ab N_b I_b Delta_b,  I_b = int_0^wc dxi tanh(E/2T)/E,
with lambda_ab = V_ab sqrt(N_a N_b): lambda_11 = 0.28, lambda_22 = 0.96, lambda_12 = 0.19,
N_2/N_1 = 0.73.  The cutoff wc is fixed by Tc = 38 K.

Linear conductivity: dirty-limit Mattis-Bardeen per band, sigma_n,a = eps0 wp_a^2 / gamma_a.

Third-order response (amplitude-mode mediated, dirty/London limit): the superfluid weight
of band a is rho_a(Delta_a, T) = pi sigma_n,a Delta_a tanh(Delta_a / 2T) / hbar.  The
diamagnetic energy 1/2 rho_a A^2 is a source 1/2 rho_a' A^2 for the gap amplitudes, which
respond through the two-band amplitude propagator
   D(w) = [ V^-1 - diag(N_a P_a(w)) ]^-1,   P_a(w) = int dxi tanh(E/2T)/E * 4 xi^2 / (4E^2 - (w+i eta)^2),
and the nonlinear current is j3 = - sum_a rho_a' dDelta_a A.  Hence
   chi3(w_D) ~ sum_ab rho_a' D_ab(w_D) rho_b'   (w_D = 2 Omega for TH; 0 and 2 Omega for FH).
"""
import numpy as np
from scipy.optimize import brentq

THZ = 4.135667        # meV per THz
KB = 0.0861733        # meV/K
EPS0 = 8.8541878e-12
Z0 = 376.730313

LAM = np.array([[0.28, 0.19], [0.19, 0.96]])
NR = np.array([1.0, 0.73])
V = LAM / np.sqrt(np.outer(NR, NR))
WP = np.array([5.2, 4.2]) * 1e3          # meV
GAM = np.array([12.5, 8.33])             # meV
HBAR_meVs = 6.582119569e-13


def sigma_n_SI():
    """dc normal-state conductivities (S/m) per band.  The tabulated plasma frequencies
    (5.2, 4.2 eV) with 1/tau = 12.5, 8.33 meV give sigma_dc ~ 6e5 Ohm^-1 cm^-1, an order of
    magnitude above the film's measured 38 K conductivity (~2.3e4 Ohm^-1 cm^-1 at 0.3 THz,
    rho ~ 50 muOhm cm).  We keep the band ratio and scattering rates and rescale the total
    to the measured normal-state value (film input, not fitted to any superconducting data)."""
    wp = WP / 1e3 / 6.582119569e-16         # rad/s
    g = GAM / 1e3 / 6.582119569e-16
    s = EPS0 * wp ** 2 / g
    w = 2 * np.pi * 0.3e12
    re = np.sum(s * g ** 2 / (g ** 2 + w ** 2))
    return s * (2.3e6 / re)


def _I(D, T, wc, n=4000):
    xi = np.linspace(0, wc, n + 1)[1:] - wc / (2 * n)
    E = np.sqrt(xi ** 2 + D ** 2)
    th = np.tanh(E / (2 * KB * T)) if T > 0 else np.ones_like(E)
    return np.sum(th / E) * (wc / n)


def gaps(T, wc, guess=(2.3, 7.0)):
    D = np.array(guess, float)
    for _ in range(4000):
        I = np.array([_I(D[0], T, wc), _I(D[1], T, wc)])
        Dn = V @ (NR * I * D)
        if np.max(abs(Dn - D)) < 1e-9:
            D = Dn
            break
        D = 0.5 * D + 0.5 * Dn
    return np.abs(D)


def Tc_of(wc):
    """largest eigenvalue of V N I(0,T) = 1 at T = Tc."""
    def f(T):
        I0 = _I(1e-9, T, wc)
        M = V * NR[None, :] * I0
        return np.max(np.linalg.eigvals(M).real) - 1
    if f(0.05) < 0:
        return 0.0
    return brentq(f, 0.05, 2000)


def fit_wc(Tc=38.0):
    return brentq(lambda w: Tc_of(w) - Tc, 2, 30)


def mb_band(f, D, T, n=6000):
    """Mattis-Bardeen sigma/sigma_n (complex) at frequency f (THz) for gap D (meV)."""
    w = f * THZ
    if D < 1e-6:
        return 1.0 + 0j
    kT = KB * T
    fermi = lambda e: 1 / (np.exp(np.clip(e / kT, -700, 700)) + 1)
    # sigma1: thermal part + pair breaking
    e = D + np.geomspace(1e-6, 60 * max(kT, D), n) ; de = np.gradient(e)
    g = (e ** 2 + D ** 2 + w * e) / (np.sqrt(e ** 2 - D ** 2) * np.sqrt((e + w) ** 2 - D ** 2))
    s1 = 2 / w * np.sum((fermi(e) - fermi(e + w)) * g * de)
    if w > 2 * D:
        e2 = np.linspace(D - w, -D, n + 2)[1:-1]; de2 = e2[1] - e2[0]
        g2 = (e2 ** 2 + D ** 2 + w * e2) / (np.sqrt(e2 ** 2 - D ** 2 + 0j).real * np.sqrt((e2 + w) ** 2 - D ** 2 + 0j).real)
        s1 += 1 / w * np.sum((1 - 2 * fermi(e2 + w)) * np.abs(g2) * de2)
    # sigma2
    lo = max(D - w, -D)
    e3 = np.linspace(lo, D, n + 2)[1:-1]; de3 = e3[1] - e3[0]
    num = e3 ** 2 + D ** 2 + w * e3
    den = np.sqrt(np.maximum(D ** 2 - e3 ** 2, 1e-30)) * np.sqrt(np.maximum((e3 + w) ** 2 - D ** 2, 1e-30))
    s2 = 1 / w * np.sum((1 - 2 * fermi(e3 + w)) * num / den * de3)
    return s1 + 1j * s2


def purity_factor(D, T, gam):
    """Superfluid weight for arbitrary purity (nonmagnetic impurities, Abrikosov-Gor'kov)
    relative to the dirty-limit (Mattis-Bardeen) value pi*Delta*tau*tanh(Delta/2T):
       rho_s/rho_0 = pi T sum_n Delta^2 / [(w_n^2+Delta^2)(sqrt(w_n^2+Delta^2) + gam/2)]."""
    if D < 1e-9:
        return 1.0
    kT = KB * max(T, 0.3)
    n = np.arange(0, 20000)
    wn = np.pi * kT * (2 * n + 1)
    s = 2 * np.pi * kT * np.sum(D ** 2 / ((wn ** 2 + D ** 2) * (np.sqrt(wn ** 2 + D ** 2) + gam / 2)))
    dirty = np.pi * D / gam * np.tanh(D / (2 * kT))
    return s / dirty


def sigma_total(f, Ds, T):
    """complex sigma (S/m) at frequency f (THz): per-band Mattis-Bardeen with the
    superfluid (inductive) part corrected for finite purity."""
    sn = sigma_n_SI()
    out = 0
    for a in range(2):
        mb = mb_band(f, Ds[a], T)
        out = out + sn[a] * (mb.real + 1j * mb.imag * purity_factor(Ds[a], T, GAM[a]))
    return out


def P_amp(w, D, T, wc, eta=0.02, n=6000):
    xi = np.linspace(0, wc, n + 1)[1:] - wc / (2 * n)
    E = np.sqrt(xi ** 2 + D ** 2)
    th = np.tanh(E / (2 * KB * max(T, 1e-3)))
    return np.sum(th / E * 4 * xi ** 2 / (4 * E ** 2 - (w + 1j * eta) ** 2)) * (wc / n)


def D_amp(w, Ds, T, wc, eta=0.02):
    P = np.array([P_amp(w, Ds[a], T, wc, eta) for a in range(2)])
    K = np.linalg.inv(V) - np.diag(NR * P)
    return np.linalg.inv(K)


def rho_prime(Ds, T):
    """d rho_a / d Delta_a for rho_a = pi sigma_n,a Delta tanh(Delta/2kT) (arbitrary common units)."""
    sn = sigma_n_SI()
    out = []
    for a in range(2):
        D = Ds[a]
        if D < 1e-9:
            out.append(0.0); continue
        x = D / (2 * KB * T)
        out.append(np.pi * sn[a] * (np.tanh(x) + x / np.cosh(x) ** 2))
    return np.array(out)


if __name__ == "__main__":
    wc = fit_wc()
    print(f"cutoff wc = {wc:.1f} meV ({wc/THZ:.1f} THz)")
    for T in (5, 12, 20, 25, 30, 32, 34, 36, 37.5):
        D = gaps(T, wc)
        print(f"T={T:5.1f} K  2Dpi={2*D[0]/THZ:.3f} THz  2Dsig={2*D[1]/THZ:.3f} THz")
