"""Measured values used for comparison only (see experimental_targets.md for sources)."""
import numpy as np

# von Hoegen et al., Nature 555, 79 (2018), Fig. 5: reconstructed A1 potential (digitized).
# Q in Å amu^1/2, U in eV per unit cell; the authors scaled both axes onto their DFT curve.
POT_Q = np.array([-1.329, -1.178, -1.031, -0.880, -0.826, -0.736, -0.664, -0.582, -0.557, -0.467,
                  -0.377, -0.287, -0.197, -0.061, 0.090, 0.180, 0.223, 0.298, 0.341, 0.377, 0.424,
                  0.449, 0.496, 0.528, 0.600, 0.750])
POT_U = np.array([0.698, 0.524, 0.410, 0.330, 0.291, 0.235, 0.191, 0.177, 0.155, 0.111, 0.069,
                  0.028, 0.022, 0.006, 0.006, 0.042, 0.042, 0.133, 0.111, 0.150, 0.180, 0.249,
                  0.260, 0.294, 0.382, 0.659])
# DFT potential of the same paper (Methods): V = w^2 Q^2/2 + a3 Q^3/3 + a4 Q^4/4 + a5 Q^5/5
DFT_A = dict(nu=19.0, a3=1.56765, a4=0.9008, a5=0.0071)       # eV, Å amu^1/2
QMAX_SOFT = 1.4          # Å amu^1/2 (Fig. 5 caption)
U_MAX = 0.7              # eV per cell at the largest excursion (Methods)
E_DEP = 0.6              # eV per cell deposited (estimate, Methods)
# Ext. Data Fig. 6: fundamental-frequency renormalization vs peak field (binned points)
REN_E = np.array([5.7, 11.6, 17.6, 20.9])
REN_W = np.array([1.004, 1.012, 1.020, 1.036])
REN_ERR = 0.004
NU_TO = 19.0
NU_POL_SH = 16.2         # polariton fundamental seen by the SH probe (THz)
TO4_FWHM_CM = 33.5       # damping of the TO4 oscillator, reflectivity fit (Ext. Data Table 1), congruent
# Stoichiometric LiNbO3, Raman (Schwarz & Maier, PRB 55, 11041 (1997)): 631 cm-1 main line,
# Gamma = Gamma_d [2 n(nu/2) + 1] + Gamma_s [n(nu_t) - n(nu_t + nu)], Gamma_d = 9.3, Gamma_s = 10.4 cm-1,
# nu_t = 251 cm-1  ->  18.9 cm-1 FWHM at 293 K; 9.3 cm-1 at low temperature.
TO4_FWHM_STOICH_293 = 18.9
TO4_FWHM_STOICH_0 = 9.3

# Mankowsky et al., PRL 118, 197601 (2017), Fig. 6a: (P_S/P_S0) at peak vs fluence (digitized)
REV = {19.0: (np.array([1.8, 6.7, 12.8, 20.7, 28.7, 37.8, 51.8, 65.9, 94.5]),
              np.array([0.99, 0.91, 0.79, 0.65, 0.48, 0.34, 0.17, -0.24, -0.35])),
       21.0: (np.array([10.4, 16.5, 24.4, 31.1, 36.0, 52.4, 90.2]),
              np.array([0.95, 0.91, 0.85, 0.81, 0.78, 0.71, 0.60])),
       22.0: (np.array([4.9, 14.6, 24.4, 34.8, 49.4, 67.7, 87.2, 103.7, 124.4, 173.8]),
              np.array([0.99, 0.98, 0.96, 0.96, 0.93, 0.90, 0.88, 0.86, 0.79, 0.70]))}
REV_THRESHOLD = 60.0     # mJ/cm^2 at 19 THz


def dft_potential(Q):
    w2 = (2 * np.pi * DFT_A["nu"] * 1e12) ** 2 * 1e-20 * 1.66053906660e-27 / 1.602176634e-19
    return 0.5 * w2 * Q ** 2 + DFT_A["a3"] * Q ** 3 / 3 + DFT_A["a4"] * Q ** 4 / 4 + DFT_A["a5"] * Q ** 5 / 5


if __name__ == "__main__":
    import os
    np.savetxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), "exp_potential.csv"),
               np.c_[POT_Q, POT_U], delimiter=",", header="Q_IR (A amu^1/2), U (eV); von Hoegen et al. 2018 Fig. 5, digitized",
               fmt="%.3f")
