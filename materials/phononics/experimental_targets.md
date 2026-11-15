# Nonlinear phononics in LiNbO3: experimental targets

Material: ferroelectric LiNbO3 (R3c, 10 atoms per rhombohedral cell, T_C = 1480 K), room
temperature, x-cut or c-axis-in-plane single crystals, pump polarized along the polar c axis.
The driven mode is the highest-frequency A1 (TO4) phonon, which is both infrared and Raman
active (no inversion centre).

## Why this experiment

- It is the cleanest case of light resonantly driving a single infrared-active phonon far
  into its anharmonic regime (oxygen displacements of about 14 pm, about 0.7 eV per cell).
- The measured quantities are lattice quantities (amplitude, potential shape, frequency
  renormalization, rectified polar displacement), which is what a machine-learned
  interatomic potential plus phono3py can predict. The electronic contribution enters only
  through the linear Born charges and epsilon_inf.
- Two independent experiments from the same group, with the same mode, pulse duration and
  field range, give two independent sets of observables (A and B below).
- Alternatives considered and rejected:
  - SiC, Cartella et al., PNAS 115, 12148 (2018): the parametric gain is controlled by
    the field-induced change of epsilon_inf and of the Born charge (electronic quantities
    alpha and beta), which were fitted to the data there and which MACE cannot provide.
  - SrTiO3, Li et al. (2019) and Kozina et al. (2019): the soft mode is a quantum-
    paraelectric, strongly temperature-renormalized mode. A static MACE potential does
    not describe it without quantum anharmonic (self-consistent) treatment.

## A. von Hoegen, Mankowsky, Fechner, Först & Cavalleri, Nature 555, 79 (2018)

"Probing the interatomic potential of solids with strong-field nonlinear phononics";
arXiv:1708.07659.

Setup:
- Pump: carrier-envelope-phase-stable mid-infrared pulses, 17.5 THz carrier (just red of
  nu_TO = 19 THz), 4 THz bandwidth, 150 fs duration, 3 uJ, about 65 um spot, peak fields up
  to 20 MV/cm. Room temperature.
- Probe: 30 fs, 800 nm. Phase-resolved polarization rotation and second-harmonic (SH)
  intensity. The SH signal is proportional to the lattice velocity dQ/dt integrated over
  the first ~1.3-2 um below the surface along the 400 nm light line (group index 3.03).

Measured (numbers as published; "digitized" means read from the figure, ±3 % of axis range):

| Observable | Value | Source |
|---|---|---|
| Driven mode frequency | nu_TO = 19 THz (632 cm^-1) | text; Ext. Data Table 1 |
| Reflectivity fit, TO4 oscillator | omega_TO = 632 cm^-1, strength 955.9 cm^-1, damping 33.5 cm^-1, eps_inf = 4.4054 | Ext. Data Table 1 (arXiv version) |
| Phonon-polariton fundamentals (phase matching) | 15.3 THz (800 nm probe), 16.2 THz (SH probe) | Fig. 2 caption |
| Harmonics of the polariton in SH/PR spectra | up to n = 5 (75 THz) at 20 MV/cm | text |
| Field dependence of harmonics n = 1, 2, 3 | linear, quadratic, cubic | Fig. 2c |
| Maximum mode excursion at 20 MV/cm | 1.4 Å amu^1/2 (soft side), oxygen displacement ~14 pm | Fig. 5 caption, Methods |
| Potential energy at the largest excursion | 0.7 eV per unit cell | Methods |
| Energy deposited per unit cell (estimate from pulse energy) | 0.6 eV | Methods |
| Reconstructed potential, extreme points | U(-1.33) = 0.70 eV, U(+0.75) = 0.66 eV (digitized, Fig. 5) | Fig. 5 |
| Asymmetry of the potential | soft for Q < 0, stiff for Q > 0; excursion ratio about 1.8 at 0.7 eV | Fig. 5 (digitized) |
| Frequency renormalization of the fundamental vs peak field (binned points) | omega'/omega = 1.004, 1.012, 1.020, 1.036 at 5.7, 11.6, 17.6, 20.9 MV/cm (blue shift); fit sqrt(1 + a E^2) | Ext. Data Fig. 6 (digitized) |

Note on the reconstructed potential: the authors scale the measured U(Q) onto their DFT
potential with one factor B (x axis by B, y axis by B^2). The absolute amplitudes of Fig. 5
therefore inherit the DFT anharmonicity. The scale-free content of the data is the shape.
The DFT fit used in the paper is a3 = 1567.65 meV amu^-3/2 Å^-3, a4 = 900.8 meV amu^-2 Å^-4,
a5 = 7.1 meV amu^-5/2 Å^-5, with V = 1/2 w^2 Q^2 + a3 Q^3/3 + a4 Q^4/4 + a5 Q^5/5 (Methods).

Digitized reconstructed potential (Fig. 5 of the Nature paper; Q in Å amu^1/2, U in eV),
also stored in `exp_potential.csv`:

    Q     U        Q     U
  -1.329 0.698   0.090 0.006
  -1.178 0.524   0.180 0.042
  -1.031 0.410   0.223 0.042
  -0.880 0.330   0.298 0.133
  -0.826 0.291   0.341 0.111
  -0.736 0.235   0.377 0.150
  -0.664 0.191   0.424 0.180
  -0.582 0.177   0.449 0.249
  -0.557 0.155   0.496 0.260
  -0.467 0.111   0.528 0.294
  -0.377 0.069   0.600 0.382
  -0.287 0.028   0.750 0.659
  -0.197 0.022
  -0.061 0.006

## B. Mankowsky, von Hoegen, Först & Cavalleri, Phys. Rev. Lett. 118, 197601 (2017)

"Ultrafast reversal of the ferroelectric polarization"; arXiv:1701.06312.

Setup:
- Pump: 150 fs CEP-stable mid-infrared pulses at 19 THz (5.5 THz FWHM bandwidth), also at
  21 and 22 THz. Fluence up to 95 mJ/cm^2 (19 THz). 1/e field penetration depth 3.2 um.
  5 mm thick x-cut single-domain crystal, room temperature.
- Probe: 35 fs 800 nm. SH generated within the coherence length l_c = 1.27 um below the
  surface. SH intensity proportional to P_S^2; phase-resolved SH gives the sign of P_S.

Peak field from fluence (used here for the comparison; Gaussian pulse with 150 fs intensity
FWHM): E0 = sqrt(2 F / (c eps0 * 1.0645 * 150 fs)), i.e. 95 mJ/cm^2 -> 21.2 MV/cm, consistent
with the 20 MV/cm of experiment A in the same laboratory.

| Observable | Value | Source |
|---|---|---|
| Normalized SH intensity minimum | 0.77 (10 mJ/cm^2), 0.32 (25 mJ/cm^2), 0 (95 mJ/cm^2) | Fig. 2b (digitized) |
| SH recovery | partial recovery within ~0.5 ps, 16 THz polariton oscillations | Fig. 2b |
| (P_S/P_S0) at peak, 19 THz pump | 0.99, 0.91, 0.79, 0.65, 0.48, 0.34, 0.17, -0.24, -0.35 at 1.8, 6.7, 12.8, 20.7, 28.7, 37.8, 51.8, 65.9, 94.5 mJ/cm^2 | Fig. 6a (digitized) |
| (P_S/P_S0) at peak, 21 THz | 0.95, 0.91, 0.85, 0.81, 0.78, 0.71, 0.60 at 10.4, 16.5, 24.4, 31.1, 36.0, 52.4, 90.2 mJ/cm^2 | Fig. 6a (digitized) |
| (P_S/P_S0) at peak, 22 THz | 0.99, 0.98, 0.96, 0.96, 0.93, 0.90, 0.88, 0.86, 0.79, 0.70 at 4.9, 14.6, 24.4, 34.8, 49.4, 67.7, 87.2, 103.7, 124.4, 173.8 mJ/cm^2 | Fig. 6a (digitized) |
| Reversal threshold (19 THz) | 60 mJ/cm^2; transient reversal lasting ~200 fs | text |
| Reversed polarization | about 40 % of P_S0 | text |

## Linear inputs not taken from the nonlinear experiments

| Quantity | Value | Source |
|---|---|---|
| Lattice constants (300 K) | a = 5.14829 Å, c = 13.8631 Å; Li z = 0.2829, O (0.0492, 0.3446, 0.0647) | Abrahams, Reddy & Bernstein, J. Phys. Chem. Solids 27, 997 (1966) |
| A1 TO mode polarities (DFPT, LDA, ABINIT) | p3 = 3.65, 0.45, 0.67, 3.82 (×10^-2 a.u.) for TO1-TO4 | Veithen, Gonze & Ghosez, PRL 93, 187401 (2004), Table III (arXiv:cond-mat/0311240); polarities from Veithen & Ghosez, PRB 65, 214302 (2002) |
| DFPT A1 TO frequencies (LDA) | 243, 287, 355, 617 cm^-1 | same, Table I |
| epsilon_inf (extraordinary, z) | 4.6 | Barker & Loudon, Phys. Rev. 158, 433 (1967) |
| Spontaneous polarization | P_S0 = 0.71 C/m^2 | Wemple, DiDomenico & Camlibel, J. Phys. Chem. Solids 29, 1797 (1968) |
| Measured A1 TO frequencies (Raman, 300 K) | 252, 275, 332, 632 cm^-1 | Barker & Loudon (1967); Schaufele & Weber (1966) |

## C. Lifetime of the driven mode (independent check of the MACE cubic force constants)

| Observable | Value | Source |
|---|---|---|
| A1(TO4) Raman linewidth, stoichiometric LiNbO3 | Gamma(T) = Gamma_d [2n(nu/2)+1] + Gamma_s [n(nu_t) - n(nu_t+nu)], Gamma_d = 9.3 cm^-1, Gamma_s = 10.4 cm^-1, nu_t = 251 cm^-1 -> 9.3 cm^-1 (T -> 0), 18.9 cm^-1 FWHM at 293 K | Schwarz & Maier, Phys. Rev. B 55, 11041 (1997) |
| TO4 oscillator damping, congruent crystal of experiment A | 33.5 cm^-1 (includes defect broadening) | von Hoegen et al., Ext. Data Table 1 |

## D. Literature energy surfaces (used only in the controlled diagnostics, diagnostics.py)

| Source | What it provides for LiNbO3 | Used as |
|---|---|---|
| Abalmasov, arXiv:1910.05020 (2019) | Landau double well for Q_P (nu_P = 7.5 THz, Q_P0 = 3.1 Å amu^1/2 from ab initio displacements); biquadratic coupling c2 Q_P^2 Q_IR^2, c2 = 34 meV amu^-2 Å^-4 from DFT paraelectric/ferroelectric TO4 frequencies; c1 = 10, c3 = 17 meV amu^-2 Å^-4 and E_d = 2.8 MV/cm fitted to experiment B | L1 (first-principles part); L3 (with fitted terms, listed only) |
| von Hoegen et al., Nature 555, 79 (2018), Methods | DFT a3, a4, a5 of Q_IR | L2 = L1 + these |
| Subedi, PRB 92, 214303 (2015) | qualitative statement that LiNbO3 has a large g Q_P Q_IR^2 coupling; no coefficients | not usable |
| Mankowsky et al. (2017), supplement section 3 | dimensionless two-mode model | not usable |
| Hardy et al., arXiv:2312.04732 | DFT-PBE four-A1-mode surface; coefficients not tabulated, fields in arbitrary units | not usable |
