"""Amplitude-mode mediated nonlinear THz response of MgB2 vs temperature.

Narrowband drive at Omega = 0.63 THz:
  chi_TH(T) = rho'^T D(2 Omega) rho'           (third harmonic, 3 Omega)
  chi_FH(T) = rho'^T [2 D(0) + D(2 Omega)] rho' (first harmonic, Omega)
Thin-film radiation and the experiment's normalisation by the transmitted linear field:
  I/E_A^6 ~ |chi|^2 / |1 + n_s + Z0 sigma(w_out) d|^2,   n_s = 3.1 (MgO), d = 20 nm.
Broadband drive (Fig. 2): single-cycle A(t) with spectrum peaked at 0.6 THz; the gap
amplitudes respond to A^2 through D(w); j3(t) = -sum_a rho'_a dDelta_a(t) A(t); the
radiated spectrum is |j3(w)|^2 / |1 + n_s + Z0 sigma(w) d|^2.
Writes nonlinear_out.npz and prints comparisons with the measured curves.
"""
import os
import numpy as np
from twoband import (fit_wc, gaps, D_amp, rho_prime, sigma_total, THZ, Z0)

HERE = os.path.dirname(os.path.abspath(__file__))
NS, DFILM = 3.1, 20e-9
OMEGA = 0.63
wc = fit_wc()
Ts = np.array([7, 9, 12, 15, 17, 20, 22, 24, 25, 27, 28, 30, 32, 34, 36, 37])
FH, TH, gp = [], [], []
for T in Ts:
    D = gaps(T, wc)
    rp = rho_prime(D, T) / 1e8
    D2 = D_amp(2 * OMEGA * THZ, D, T, wc)
    D0 = D_amp(1e-4, D, T, wc)
    chi_th = rp @ D2 @ rp
    chi_fh = rp @ (2 * D0 + D2) @ rp
    s1 = sigma_total(OMEGA, D, T); s3 = sigma_total(3 * OMEGA, D, T)
    FH.append(abs(chi_fh) ** 2 / abs(1 + NS + Z0 * s1 * DFILM) ** 2)
    TH.append(abs(chi_th) ** 2 / abs(1 + NS + Z0 * s3 * DFILM) ** 2)
    gp.append(2 * D[0] / THZ)
FH = np.array(FH) / FH[0]; TH = np.array(TH) / TH[0]
exp_FH = {7.5: 1.0, 12: 0.83, 17: 0.68, 20: 0.65, 23: 0.5, 27.5: 0.32, 31: 0.2, 35: 0.1}
exp_TH = {8: 1.0, 12: 0.9, 17: 0.5, 20: 0.3, 22: 0.2, 25: 0.1}
print(" T    2Dpi   FH(model)  TH(model)")
for i, T in enumerate(Ts):
    print(f"{T:4.0f}  {gp[i]:.3f}   {FH[i]:.3f}      {TH[i]:.3f}")
print("measured FH (12 kV/cm):", exp_FH)
print("measured TH (12 kV/cm):", exp_TH)

# ---- broadband drive -------------------------------------------------------------
dt = 0.01                                   # ps
t = np.arange(-20, 60, dt)
f0 = 0.6
E = np.exp(-(t / 0.55) ** 2) * np.cos(2 * np.pi * f0 * t + 0.3)     # ~single cycle
A = -np.cumsum(E) * dt
freqs = np.fft.rfftfreq(len(t), dt)          # THz
peaks = []
Tb = np.array([7, 12, 17, 22, 27, 30, 32, 36, 37])
spec_store = []
for T in Tb:
    D = gaps(T, wc)
    rp = rho_prime(D, T) / 1e8
    S = np.fft.rfft(A ** 2)
    dD = np.zeros((2, len(freqs)), complex)
    for k, f in enumerate(freqs):
        if f > 4:
            break
        Dm = D_amp(max(f, 1e-4) * THZ, D, T, wc, eta=0.02)
        dD[:, k] = Dm @ rp * S[k]
    dDt = np.fft.irfft(dD, n=len(t))
    j3 = -(rp[:, None] * dDt).sum(0) * A
    J = np.fft.rfft(j3)
    sig = np.array([sigma_total(f, D, T) if 0.2 < f < 2.5 else 0 for f in freqs])
    out = np.abs(J) ** 2 / np.abs(1 + NS + Z0 * sig * DFILM) ** 2
    m = (freqs > 0.3) & (freqs < 2.0)
    pk = freqs[m][np.argmax(out[m])]
    peaks.append(pk)
    spec_store.append(out[m])
    print(f"broadband T={T:4.0f} K: model peak {pk:.2f} THz   (2Dpi = {2*D[0]/THZ:.2f})")
print("measured peaks: 7K 1.05, 12K 1.03, 17K 0.95, 22K 0.9, 27K 0.8, 30K 0.7, 32K 0.65, 36K 0.6, 37K 0.55")
np.savez(os.path.join(HERE, "nonlinear_out.npz"), Ts=Ts, FH=FH, TH=TH, gap_pi=gp, Tb=Tb,
         peaks=np.array(peaks), freqs=freqs[(freqs > 0.3) & (freqs < 2.0)], spectra=np.array(spec_store))
