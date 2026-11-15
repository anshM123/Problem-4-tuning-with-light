"""External inputs for the LiNbO3 nonlinear-phononics model, each with its source.

Everything else (structure relaxation, phonons, anharmonic energy surface, cubic force
constants, three-phonon lifetimes) is computed here with the MACE potential medium-omat-0.
Nothing below is taken from the two nonlinear-phononics experiments being predicted.
"""
import numpy as np

# ---- physical constants (SI)
E_CH = 1.602176634e-19
AMU = 1.66053906660e-27
EPS0 = 8.8541878128e-12
C0 = 2.99792458e8
HBAR = 1.054571817e-34
KB = 1.380649e-23
ANG = 1e-10
THZ_CM = 33.35641          # cm^-1 per THz
QSI = ANG * np.sqrt(AMU)   # 1 Å amu^1/2 in m kg^1/2

# ---- dielectric / polar inputs (cannot be obtained from an interatomic potential)
# Mode polarities p_m,z = sum_k Z*_k,zb u_m(kb) of the four A1 TO modes, DFPT-LDA (ABINIT),
# eigendisplacements normalized to sum_k M_k |u_m(k)|^2 = 1 in atomic units:
# Veithen, Gonze & Ghosez, Phys. Rev. Lett. 93, 187401 (2004) [arXiv:cond-mat/0311240],
# Table III (Born charges and eigenvectors of Veithen & Ghosez, PRB 65, 214302 (2002)).
P3_AU = np.array([3.65, 0.45, 0.67, 3.82]) * 1e-2     # e / sqrt(m_e), TO1..TO4
ZMODE = P3_AU * np.sqrt(AMU / 9.1093837015e-31)        # e / sqrt(amu): 1.56, 0.19, 0.29, 1.63
OMEGA_DFPT_CM = np.array([243.0, 287.0, 355.0, 617.0])  # same work, A1 TO frequencies
# High-frequency permittivity along the polar axis (extraordinary): Barker & Loudon,
# Phys. Rev. 158, 433 (1967).  (DFPT-LDA gives 5.8 [Li et al., Sci. Rep. 6, 34085 (2016)],
# LDA overestimates eps_inf; the measured value is used, sensitivity reported.)
EPS_INF = 4.6
# Spontaneous polarization at 300 K: Wemple, DiDomenico & Camlibel, J. Phys. Chem. Solids
# 29, 1797 (1968).
PS0 = 0.71   # C/m^2
# Nominal ionic charges, used ONLY to fix the sign of each MACE A1 eigenvector relative to
# the DFPT mode polarity (magnitudes come from P3_AU).  For TO1 and TO4, which carry almost
# all of the polar weight, the sign is unambiguous.
ZZZ_SIGN = {"Li": 1.0, "Nb": 5.0, "O": -2.0}

# ---- the two experiments (inputs that define the drive, not results)
# A: von Hoegen et al., Nature 555, 79 (2018): 17.5 THz, 150 fs, CEP stable, up to 20 MV/cm.
EXP_A = dict(nu=17.5, fwhm_fs=150.0, fields=[2.0, 5.7, 8.0, 11.6, 14.8, 17.6, 20.0, 20.9])
# B: Mankowsky et al., PRL 118, 197601 (2017): 19 THz, 150 fs; SH layer 1.27 um.
EXP_B = dict(nu=19.0, fwhm_fs=150.0, lc_um=1.27,
             fluences=[1.8, 6.7, 12.8, 20.7, 28.7, 37.8, 51.8, 65.9, 94.5])


def fluence_to_field(F_mJcm2, fwhm_fs=150.0):
    """Peak field (MV/cm) of a Gaussian pulse with intensity FWHM fwhm_fs and fluence F."""
    I0 = F_mJcm2 * 10.0 / (1.0645 * fwhm_fs * 1e-15)    # W/m^2  (1 mJ/cm^2 = 10 J/m^2)
    return np.sqrt(2 * I0 / (C0 * EPS0)) / 1e8
