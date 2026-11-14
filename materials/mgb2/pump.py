"""Strong THz pump of MgB2 at 7 K (Supplemental Fig. S3 of Katsumi et al.).

1. Absorbed energy: broadband single-cycle pulse (peak 150 kV/cm, spectrum peaked at
   0.6 THz) incident on 20 nm MgB2 on MgO; field in the film
   E_f(w) = 2 E_in(w) / (1 + n_s + Z0 sigma(w) d), absorbed energy density
   u = (1/pi) int_0^inf Re sigma(w) |E_f(w)|^2 dw.  sigma is evaluated (a) in the 7 K
   superconducting state and (b) in the normal state, bracketing the strongly driven case.
2. Quasi-equilibrium after 2 ps: the energy heats the electrons (two-band BCS heat
   capacity, Sommerfeld gamma_n = 2.7 mJ/mol K^2 split by the band DOS ratio) and, as a
   second bound, electrons + lattice (Debye, Theta_D = 750 K).
3. Complex conductivity at the resulting temperature T* from the two-band model.
"""
import numpy as np
from twoband import fit_wc, gaps, sigma_total, THZ, KB, Z0, NR

NS, DF = 3.1, 20e-9
wc = fit_wc()
VM = 17.46e-6                    # m^3/mol
GAMMA_N = 2.7e-3                 # J/mol K^2
KB_J = 1.380649e-23


def absorbed(Epk=150e5, T=7.0, normal=False):
    dt = 0.005e-12
    t = np.arange(-15e-12, 15e-12, dt)
    E = np.exp(-(t / 0.55e-12) ** 2) * np.cos(2 * np.pi * 0.6e12 * t + 0.3)
    E *= Epk / np.max(abs(E))
    Ew = np.fft.rfft(E) * dt
    f = np.fft.rfftfreq(len(t), dt)
    D = [0, 0] if normal else gaps(T, wc)
    u = 0.0
    for k in range(1, len(f)):
        fT = f[k] / 1e12
        if fT > 5:
            break
        s = sigma_total(fT, D, T if not normal else 40.0)
        Ef = 2 * Ew[k] / (1 + NS + Z0 * s * DF)
        u += s.real * abs(Ef) ** 2 * (2 * np.pi * (f[1] - f[0]))
    return u / np.pi           # J/m^3


def electron_energy(T):
    """two-band BCS electronic energy (J/mol) relative to T=0 (numerical, via entropy)."""
    Ts = np.linspace(0.5, T, 80)
    C = []
    for Tk in Ts:
        C.append(heat_capacity(Tk))
    return np.trapezoid(C, Ts)


def heat_capacity(T):
    """two-band BCS electronic heat capacity (J/mol K):
    C = sum_a gamma_a T * [int_0^inf f(1-f)/kT (E^2 - T/2 dDelta_a^2/dT) dxi] / [(pi^2/6)(kT)^2],
    gamma_a = gamma_n N_a / (N_1 + N_2)."""
    Nw = NR / NR.sum()
    h = 0.05
    D0 = gaps(T, wc)
    D1 = gaps(T + h, wc) if T + h < 38 else np.zeros(2)
    dD2 = (D1 ** 2 - D0 ** 2) / h
    kT = KB * T
    xi = np.linspace(0, 40 * kT, 4000)
    C = 0.0
    for a in range(2):
        E = np.sqrt(xi ** 2 + D0[a] ** 2)
        x = E / kT
        ff = np.exp(-x) / (1 + np.exp(-x)) ** 2
        num = np.trapezoid(ff / kT * (E ** 2 - 0.5 * T * dD2[a]), xi)
        C += Nw[a] * GAMMA_N * T * num / ((np.pi ** 2 / 6) * kT ** 2)
    return C


def lattice_energy(T, T0=7.0, theta=750.0):
    # Debye T^3 regime: U = (3/4) beta T^4 per mole of atoms x 3 atoms
    beta = 12 * np.pi ** 4 * 8.314462618 * 3 / (5 * theta ** 3)
    return 0.25 * beta * (T ** 4 - T0 ** 4)


if __name__ == "__main__":
    us = absorbed(normal=False); un = absorbed(normal=True)
    print(f"absorbed energy density: SC-state sigma {us:.3e} J/m^3, normal-state sigma {un:.3e} J/m^3")
    Ts = np.linspace(7, 40, 34)
    Ue = np.array([electron_energy(T) for T in Ts]); Ue -= Ue[0]
    Ul = np.array([lattice_energy(T) for T in Ts])
    for u, lab in ((us, "SC-state absorption"), (un, "normal-state absorption")):
        um = u * VM
        Te = np.interp(um, Ue, Ts); Tl = np.interp(um, Ue + Ul, Ts)
        print(f"{lab}: {um:.3f} J/mol -> T* = {Te:.1f} K (electrons only), {Tl:.1f} K (electrons+lattice)")
    print("\nconductivity (Ohm^-1 cm^-1) at candidate T*:")
    for T in (7, 20, 25, 28, 30, 32, 34, 36, 38):
        D = gaps(T, wc) if T < 38 else [0, 0]
        row = [sigma_total(f, D, T) / 100 for f in (0.5, 0.8, 1.0)]
        print(f"T={T:4.0f}  s1(0.8)={row[1].real:.2e}  s2(0.5)={row[0].imag:.2e}  s2(1.0)={row[2].imag:.2e}")
    print("measured at 2 ps after 150 kV/cm (7 K): s1 peak ~2.55e4 near 0.8 THz, s2(0.5) ~2.6e4, s2(1.0) ~1.5e4")
