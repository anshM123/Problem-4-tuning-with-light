# Experimental targets for simulation: three ultrafast experiments

Sources read in full (main text, methods, appendices, supplementary material included in the arXiv PDFs):

| Paper | Published | Open-access version used |
|---|---|---|
| A. Madéo et al., dark excitons in monolayer WSe2 | Science 370, 1199 (2020), doi:10.1126/science.aba1029 | arXiv:2005.00241v1 (1 May 2020, 34 pp. incl. SI) |
| B. René de Cotret et al., UEDS phonon populations in graphite | Phys. Rev. B 100, 214115 (2019) | arXiv:1908.02795v3 (9 Nov 2019) |
| B'. (companion, cited as ref. 22 of B) Stern et al., same graphite UEDS data set | Phys. Rev. B 97, 165416 (2018) | arXiv:1708.01251v2 |
| C. Katsumi et al., MgB2 THz 2DCS | Phys. Rev. Lett. 135, 036902 (2025), doi:10.1103/g5rp-6vb1 | arXiv:2411.10852v2 (17 Jul 2025, incl. SM) |

Conventions in this file:
- "Stated" = number printed in the text/caption/table.
- "Graphical" = read off a figure by eye from the arXiv PDF; typical reading error is about +/-5% of full scale on the axis, larger for crowded panels.
- "Derived" = computed here from stated numbers; the arithmetic is shown.
- The arXiv versions can differ from the journal versions. Paper C's paper title is "Amplitude mode in a multi-gap superconductor MgB2 investigated by terahertz two-dimensional coherent spectroscopy". It is a THz 2DCS paper with first-harmonic (FH) and third-harmonic (TH) analysis, not a THG-only paper.

---

## Paper A: Madéo et al., dark excitons in monolayer WSe2 (TR-XUV-µ-ARPES)

### A.1 Sample
- Exfoliated monolayer (1L) WSe2, about 40 x 20 µm, on hBN about 20 nm thick, on n-doped Si. The Si has a few-nm native oxide (Methods; SI §1).
  - The hBN thickness was chosen thick enough to prevent exciton quenching by the substrate, and thin enough to give a conductive path that avoids charging.
  - The monolayer also contacts bulk WSe2 lying directly on Si, which gives a second conduction path.
- Prepared by PDMS dry transfer, then annealed in vacuum at 450 °C for 4 h (SI §1).
- Room-temperature PL of the A exciton: 1.64 eV, with 532 nm CW excitation at 295 K (SI §1, Fig. S1c).
- Measurement temperature: 90 K, in UHV (main text).
- Photoelectrons were collected from a 6 x 40 µm region of the monolayer selected by a field aperture (Methods).
- Lattice constant used in DFT: 3.28 Å (experimental).
  - Derived: |Γ-K| = 4π/(3a) = 1.277 Å⁻¹.

### A.2 Light sources and resolution
- **Laser:** Yb fiber amplifier, 1 MHz, 1030 nm, 230 fs, 100 µJ.
  - 20 µJ pumps an OPA tunable from 320 to 2500 nm, giving 0.1–1 µJ per pulse.
  - OPA bandwidth: 5 nm, which is 12 meV FWHM at 720 nm (Methods; SI §2a).
- **Pump geometry and duration:**
  - p-polarized, 22° angle of incidence. The pump spot is much larger than the sample.
  - Main text says "linearly polarized"; circular polarization is not used.
  - Pulse duration at the sample: about 240 fs. Autocorrelation gives a Gaussian width of 241 fs (Methods; SI §4, Fig. S6a).
  - The K-valley rise under resonant pumping is pulse-width limited (Fig. S6b).
- **Resonant pump:** 1.72 eV (FWHM 5 nm), 14 µJ/cm².
  - Estimated exciton density: 8 x 10¹² cm⁻², assuming 15% absorption.
  - Derived check: 14 µJ/cm² / (1.72 eV x 1.602e-19 J/eV) = 5.1e13 photons/cm². Multiplied by 0.15, this gives 7.6e12 cm⁻².
- **Above-gap pump:** 2.48 eV (FWHM 5 nm; at the C-exciton resonance), 5 µJ/cm².
  - Estimated initial carrier density: 2 x 10¹² cm⁻², assuming 15% absorption.
  - Derived: 1.26e13 photons/cm² x 0.15 = 1.9e12 cm⁻².
- **XUV probe:**
  - Generation: 515 nm (SHG in a 500 µm BBO, about 40% efficiency), focused to 3 x 10¹⁴ W/cm² in a Kr jet at 1 bar backing pressure.
    - The Methods give 70 µJ in → 30 µJ of 515 nm; SI §2b gives 80 µJ → 35 µJ.
  - Selected harmonic: 21.7 eV, using a 200 nm Al filter and a 300 nm Sn bandpass filter.
  - Flux: >10⁴ photons/pulse. Count rate: about 0.1 detected photoelectrons per pulse, i.e. about 10⁵ e⁻/s. This is below the space-charge limit (SI §4, Fig. S5).
  - Delivery: toroidal mirror (f = 500 mm), 2f-2f geometry, 4° grazing incidence.
- **Analyzer:** time-of-flight momentum microscope (SPECS Metis 1000). Real-space imaging resolution is better than 50 nm.
- **Resolution:**
  - Energy: 30 meV, from the 16–84% width of the Au(111) Fermi edge at 90 K (SI §4).
  - Momentum: about 0.01 Å⁻¹ (Methods).
  - Time: instrument response = pump ⊗ probe cross-correlation. Fig. 3c plots it as a dotted black step and the 241 fs Gaussian as the green shaded curve.
- **Energy reference:** E = 0 is the Gaussian-fit peak of the upper spin-split valence band (VB1) at K, set separately at every delay (SI §5, Fig. S7c).
  - The rigid shift corrections are ≤ 60 meV (Fig. S8b). They remove surface photovoltage.
- **Integration times:** 2 min per delay (resonant data set); 3 h per delay (above-gap data set) (SI §6).

### A.3 Exciton energies and momenta
- **Photoemission excitation spectrum (Fig. 2c, top):**
  - Signal integrated from 1 to 3 eV above the VBM versus pump photon energy, scanned from 1.58 to 2.85 eV.
  - Resonances: A at 1.73 eV, B at 2.17 eV, C at 2.48 eV (caption).
  - Graphical: a smaller peak labeled "T" sits just below A, at about 1.70 eV. The "theoretical direct gap" marker is at about 2.08 eV.
- **Bright K-K (A) exciton:** 1.73 ± 0.03 eV (relative to VB1 at K), at the K valley, |k| ≈ 1.26 Å⁻¹.
  - The energy matches the pump and does not change with delay (main text, Fig. 3c inset).
- **Dark Σ-K exciton** (electron at Σ, also called Q or Λ; hole at K): 1.73 eV, at |k| ≈ 0.75 Å⁻¹ (main text).
  - It is degenerate with K-K within the 0.03 eV resolution.
  - A tighter bound from the late-time population ratio at 90 K is |ΔE| < 0.015 eV (main text).
  - SI §7 Boltzmann estimate with DOS ratio ρΣ/ρK = 3 and N_K/N_Σ = 0.5 gives ΔE = −3 meV, i.e. K-K about 3 meV below Σ-K.
  - Derived: kT(90 K) = 7.76 meV, and kT·ln(1.5) = 3.1 meV.
- **Degeneracy counting (SI §7):**
  - K-valley species: 4 (K-K, K'-K', K-K', K'-K).
  - Σ-valley species: 12 (K-Σ1..6, K'-Σ1..6).
  - Exciton mass is about m0 for both. A roughly 30% mass enhancement from binding applies equally to both and cancels in the ratio.
  - The K-K' species are lumped into the "K-valley" signal; they are not resolved separately.
- **Binding energies** (relative to each CB minimum; no thermal correction applied, since the 90 K thermal energy of about 7 meV is below the resolution):
  - K-K: about 390 meV.
  - Σ-K: about 480 meV.
  - Derived CB minima above VB1: CB(K) ≈ 1.73 + 0.39 = 2.12 eV; CB(Σ) ≈ 1.73 + 0.48 = 2.21 eV.
- **Above-gap pump, at zero delay:**
  - Broad distribution centered at about 1.90 eV in both K and Σ valleys (main text, Fig. S10).
  - Positive dispersion with a minimum around 1.75 eV, from Gaussian fits at k = −1.5 to −1.0 Å⁻¹ (SI §8). Possible origins: free carriers, hot excitons, band-gap renormalization.
- **Hole signal:** depletion is visible at the VB maxima at K/K' (Fig. S12a). The valence band broadens after excitation.

### A.4 Population dynamics: resonant 1.72 eV pump (Fig. 3, Fig. S11)
**Definitions.** The exciton density is the ARPES signal integrated over E = 1–2.5 eV and over the 2 K valleys or the 6 Σ valleys of the 1st BZ. It is corrected by computed photoemission matrix-element factors for K versus Σ (GW quasiparticle bands, dipole matrix elements, 0.2 eV Gaussian broadening).

**Stated results**
- At zero delay only K-K excitons are present (Fig. 3a,b; Fig. 2d).
- The K-valley rise coincides with the pump; it is IRF limited.
- The Σ-K population builds up "via scattering from the K-K excitons on a ~400 fs timescale" (main text). SI §9 / Fig. S11 states a delay of about 300 fs between the normalized K and Σ rises.
- K and Σ show "similar recovery time" (no number given).
- Within about 1 ps Σ-K dominates. The K/Σ ratio tends to about 0.5 (inset of Fig. 3c). Equivalently, dark/bright ≈ 2: the dark population is about twice the bright one at late times.
- No rate-equation fit or exciton-phonon scattering time is quoted beyond the ~400 fs (and ~300 fs) scales. The text calls these consistent with literature intervalley exciton-phonon scattering of "few hundred femtoseconds".

**Fig. 3c, graphical.** Axes: x = −0.5 to 4 ps; y = 0–15 a.u.

| Delay (ps) | K-valley (blue) | Σ-valley (red) | K/Σ (from these values) |
|---|---|---|---|
| −0.5 | 0.1 | 0.1 | – |
| −0.2 | 0.4 | 0.25 | – |
| 0.08–0.1 | 7.5 (peak) | 6.55 | 1.15 |
| 0.4 | 6.8 | 9.8 (peak) | 0.69 |
| 0.7 | 5.1 | 9.55 | 0.53 |
| 1.0 | 4.2 | 8.95 | 0.47 |
| 1.3 | 3.5 | 7.8 | 0.45 |
| 1.6 | 3.2 | 6.55 | 0.49 |
| 2.1 | 2.5 | 5.25 | 0.48 |
| 2.6 | 1.75 | 4.45 | 0.39 |
| 3.1 | 1.65 | 3.35 | 0.49 |
| 3.6 | 1.35 | 2.75 | 0.49 |

- Peak Σ / peak K ≈ 1.3.
- Late decay, from an exponential fit of the readings above over 1–3.6 ps: τ ≈ 2.2–2.3 ps for both valleys.
  - K from its peak: τ ≈ 2.0 ps over 0.4–3.6 ps.
  - Σ from its peak: τ ≈ 2.5 ps over 0.4–3.6 ps.
- Fig. 3c inset, ratio K/Σ (graphical): 1.6 (0 ps) → 1.13 (0.3) → 0.68 (0.6) → 0.53 (0.9) → 0.47 (1.2) → 0.44 (1.5) → 0.45–0.49 (1.8–3.8 ps). A dashed reference line sits at about 0.46.
- Fig. 3c inset, center energies (right axis 1.6–1.8+ eV, graphical):
  - K: about 1.75 eV, flat.
  - Σ: about 1.75 eV, drifting to about 1.72–1.73 eV by 3.8 ps.
  - Error bars: about ±0.03 eV.
- Figs. 3d,e: time-resolved spectra, E = 1.0–2.5 eV, delay −1 to 1.5 ps.
  - The K signal is centered at about 1.75 eV and starts at 0.
  - The Σ signal peaks at about 0.4 ps (arrow labeled "K to Σ scattering"; vertical dashed line at about 0.4 ps).
- Fig. S11, normalized K and Σ densities over −0.4 to 1.0 ps (graphical):
  - K reaches 0.5 at about 0 ps.
  - Σ reaches 0.5 at about 0.25–0.3 ps.

### A.5 Population dynamics: above-gap 2.48 eV pump (Fig. 4)
**Definitions.** Density integrated over 1–2 eV (note: a different window from Fig. 3's 1–2.5 eV) over the first BZ.

**Stated results**
- Full relaxation into K-K or Σ-K excitons within 500 fs. The energy-relaxation time from the average photoemission energy is 500 fs.
- Beyond 500 fs the peak energy is constant at about 1.73 eV in both valleys.
- Σ-K appears coincident with K-K (no delayed rise), and dark Σ-K dominates at all delays.
- At long delays, the text describes the quasi-steady state as similar to the resonant case, with a "nearly identical" K/Σ ratio. The graphical values below are lower (about 0.32–0.35 versus about 0.47).

**Fig. 4c upper panel, graphical.** y = 0 to about 7.5 a.u.

| Delay (ps) | K (blue) | Σ (red) |
|---|---|---|
| −0.2 | 0.05 | 0.25 |
| −0.1 | 0.15 | 1.0 |
| 0 | 0.4 | 3.0 |
| 0.1 | 0.7 | 5.6 |
| 0.2 | 0.8 | 5.75 |
| 0.3 | 1.0 | 6.9 |
| 0.4 | 1.15 | 6.9 |
| 0.5 | 1.1 | 6.15 |
| 0.6 | 1.15 | 5.9 |
| 0.7 | 1.1 | 5.6 |
| 1.7 | 0.5 | 1.55 |
| 3.7 | 0.35 | 1.0 |

- Σ/K is about 6–8 at 0.1–0.7 ps and about 3 at 1.7–3.7 ps.
- The 0.7→1.7 ps drop of Σ, by a factor of 3.6, implies an e-folding time of about 0.8 ps over that interval. Late decay over 1.7→3.7 ps: τ ≈ 4.5 ps for Σ and about 5 ps for K. These are coarse estimates from 2–3 points.

**Fig. 4c lower panel, center energy from Gaussian fits, graphical** (± about 0.03 eV)

| Delay (ps) | K center (eV) | Σ center (eV) |
|---|---|---|
| 0 | 1.875 | 1.92 |
| 0.1 | 1.795 | 1.805 |
| 0.2 | 1.795 | 1.80 |
| 0.3–0.4 | 1.77 | 1.775–1.78 |
| 0.5 | 1.765 | 1.765 |
| 0.6–0.7 | 1.752–1.755 | 1.75 |
| 1.7 | 1.753 | 1.742 |
| 3.7 | 1.752 | 1.737 |

**Other panels**
- Fig. 4c inset, ratio K/Σ (graphical): 0.15 (0 ps) → about 0.2 (0.5 ps) → 0.30 (1.6 ps) → 0.32 (3.6 ps). A dashed reference line sits at about 0.45, the resonant-case limit.
- Fig. 4d: time-resolved spectra over −0.2 to 1.3 ps. Arrows mark "exciton formation" from 0 to about 0.3 ps in both valleys. The Σ-valley maximum is near 0.3–0.4 ps at about 1.75 eV.

### A.6 Theory in the preprint
- DFT-LDA (Quantum ESPRESSO) followed by GW (BerkeleyGW):
  - Dielectric cutoff 35 Ry, 18x18x1 k-grid.
  - QP gap converged to 0.05 eV.
  - SOC added perturbatively.
  - Stated GW systematic error: about 100 meV; convergence error < 50 meV.
- These calculations supply the band overlays (dashed lines in Figs. 2–4) and the K-versus-Σ photoemission matrix-element normalization.
- The preprint contains no ab initio exciton-phonon (GW-BSE + phonon) scattering-time calculation. Intervalley exciton-phonon scattering is discussed only by citation (Christiansen et al. PRB 100, 205401 (2019); Brem et al.; Selig et al.). If the journal version adds such calculations, they are not in arXiv v1.
- Video SI: resonant −1.5 to 3.5 ps; above-gap −0.3 to 3.6 ps.

### A.7 Simulation set-up summary (A)
- **Temperature and states:** T = 90 K. States: K-K (energy E0 = 1.73 eV) and Σ-K (E0 + ΔE, with ΔE ≈ +3 meV and |ΔE| < 15 meV).
- **Degeneracies:** 4 vs 12. Combined with the Boltzmann factor, this gives an equilibrium K/Σ ≈ 0.5.
- **Resonant case:** pump 1.72 eV, Gaussian 240 fs FWHM, n0 = 8e12 cm⁻², injected only into K-K. Targets:
  - Σ rise delayed by 300–400 fs.
  - Σ peak at about 0.4 ps.
  - K/Σ: 1.6 → 0.5 within about 1 ps.
  - Common late decay of about 2–2.5 ps.
- **Above-gap case:** 2.48 eV, n0 = 2e12 cm⁻². Hot population at about 1.9 eV relaxes to 1.75 eV within 0.5 ps. Targets:
  - Σ and K appear simultaneously.
  - K/Σ: 0.15 → 0.32.

---

## Paper B: René de Cotret et al., UEDS phonon populations in graphite

### B.1 Sample and excitation
- Freestanding single-crystal natural graphite flake (supplier: Naturally Graphite), mechanically exfoliated to 70 nm (§II A).
  - Samples 10–90 nm thick showed no quantitative difference in Bragg dynamics, so scattering is close to kinematical (§II B).
- Geometry: [001] zone axis. ZA and ZO out-of-plane modes are invisible because q·e = 0.
- **Pump:**
  - 35 fs, 800 nm (1.55 eV), fluence 12 mJ/cm².
  - Pump spot 1 mm x 1 mm FWHM; probed area 500 µm x 500 µm.
  - Stated "absorbed energy density of 8 J m⁻³". This value as printed is not dimensionally plausible for 12 mJ/cm² on 70 nm; J cm⁻³ or larger is likely intended. No excitation (carrier) density is given.
- **Electron probe:**
  - 10⁷ electrons/bunch, 90 keV, 1 kHz.
  - RF-compressed to about 150 fs at the sample (measured with a streak camera).
  - Momentum resolution 0.06 Å⁻¹; real-space resolution < 1 Å.
- **Detector:** Gatan Ultrascan 895 (2.54 cm phosphor, fiber-coupled 2048 x 2048 CCD) at 25 cm.
  - Delay range: −40 to 680 ps. Transient dynamic range: 1:10⁸.
  - Diffuse intensity is about 10⁴–10⁶ times weaker than Bragg. Brightest Bragg ≈ 20,000 counts; average diffuse feature ≈ 0.2 counts.
- **Processing:** patterns six-fold symmetrized (D6h). Initial temperature: 300 K.
- No Bragg peak shift (lattice expansion) up to 680 ps.

### B.2 Computational inputs
- **Structure relaxation:** Quantum ESPRESSO PWscf, 18x18x10 Γ-centered k-mesh.
  - Force threshold 1e-8 Ry/Bohr; energy threshold 1e-14 Ry.
  - Result: a = 2.462 Å, c = 6.837 Å, 4 atoms per cell, P6₃/mmc.
- **Phonons:** DFPT (PHonon) on a 5x3x3 q-grid, self-consistency threshold 1e-18 Ry.
  - Functional: B86bPBE GGA with PAW. Ecut = 100 Ry (wavefunctions), 1200 Ry (charge density).
  - Fermi-Dirac smearing: 0.06 Ry. Dispersion correction: XDM.
- **Branch assignment:** clustering by continuity of (ω, e) along a path from Γ (Appendix A; scikit-ued).
- **Derived BZ geometry:** |Γ-K| = 4π/(3a) = 1.701 Å⁻¹; |Γ-M| = 2π/(√3 a) = 1.473 Å⁻¹.

### B.3 Theory and extraction method
**Scattering intensities**
- Bragg: I0 = Nc Ie |Σs fs e^{−Ws} e^{−iq·Rs}|².
- One-phonon diffuse: I1(q,t) = Nc Ie Σj (n_{j,k} + 1/2)/ω_{j,k} |F1j(q,t)|².
- One-phonon structure factor: |F1j|² = |Σs e^{−Ws} fs/√µs (q·e_{j,s,k})|².
- Debye-Waller factor: Ws = (1/4µs) Σ_{j,k} |a_{j,k}|² |q·e|², with |a|² = (2ħ/(Nc ω)) (n + 1/2).

**Debye-Waller sensitivity test**
- Test: all modes held at 300 K except one, raised to 5000 K (optical) or 1000 K (acoustic).
- Σs Ws changed by at most 1.5% (optical) or 5% (acoustic), uniformly in q.
- Conclusion: |F1j(q,t)| ≈ |F1j(q,t0)|.

**Modes and selection rules**
- Eight in-plane modes: LA, TA, LO1–3, TO1–3.
- Only 37% of q-points have a single mode contributing more than 50% of Σj |F1j|²/ω. Only 1% have one mode above 75% (Fig. 5).

**Inversion**
- At each reduced k, solve ΔI_k(t) = F_k n_k(t) with n_{j,k} expressed as Δn_{j,k}/ω_{j,k}.
- Data: M = 44 Brillouin zones (all in-plane H with |H| ≤ 12 Å⁻¹), N = 8 modes.
- Non-negative least squares (Lawson–Hanson). No free parameters.
- Stable solutions for |k| > 0.45 Å⁻¹. Population uncertainty: about 4e-5 (a.u., from the covariance).
- Γ-point modes, including E2g at Γ, are therefore not directly resolved in Fig. 6.

### B.4 Extracted populations Δn_{j,k}(t) (Fig. 6)
Fig. 6 shows maps for TO2, TA and LA at 500 fs, 1.5 ps, 5 ps and 100 ps, on a common colour scale of 0–1 a.u. The maps are bounded by the |k| = 0.45 Å⁻¹ circle and the BZ edge; the dashed hexagon marks the BZ midpoints. Fig. 2 shows the corresponding raw ΔI(q,τ) at the same four delays, over |q| < 12 Å⁻¹. Only TO2, TA and LA are shown. LO, TO1/TO3 and LO1–3 are solved for but not plotted. ZA and ZO are not measurable.

**Stated**
- Strongly coupled optical phonons: A1' at K and E2g at Γ.
- < 500 fs: TO2 is strongly populated at K, which is the A1' signature.
- By 5 ps: transfer out of TO2 is well underway, and acoustic populations rise along Γ-M.
- < 5 ps: the TA increase along Γ-M agrees with the predicted anharmonic decay of E2g.
- 500 fs: small LA increase at ½K, consistent with E2g decay.
- > 25 ps: TA pools at ⅓M and at M.
- At the mid-BZ (½K, ½M), TA and LA populations are lower than average. This is consistent with Umklapp decay into ZA (Paulatto et al.).
- Allowed three-phonon channels: L → T+T and L → L+T. A purely transverse mode has no three-phonon decay, so TA population builds up.

**Graphical description**

| Mode | 500 fs | 1.5 ps | 5 ps | 100 ps |
|---|---|---|---|---|
| TO2 | Moderate (about 0.4–0.6) near the K corners, weaker mid-zone | Maximum overall: bright (about 0.8–1) spots near K and along Γ-K mid-zone | Reduced, more uniform (about 0.4–0.6) | Similar to 5 ps or slightly lower, nearly uniform |
| TA | Essentially zero | Weak (about 0.2–0.3) along Γ-M mid-zone | Strong bands along Γ-M (about 0.6–0.8) | Strongest (about 1) at M on the BZ edge and near ⅓M just outside the inner circle; depleted at the mid-zone hexagon |
| LA | Weak, about 0.1, diffuse | Ring at small k (0.45 Å⁻¹ up to the mid-zone hexagon), about 0.5 | Brighter small-k ring, about 0.7–0.8 | Brightest small-k ring, about 1, inside the mid-zone hexagon |

### B.5 A1' (K-point TO2) dynamics and the non-thermal lattice model (Fig. 7, Table I)
**Observable**
- Δn_A1'(t): TO2 population integrated in a circular region |k − K| ≤ 0.3 Å⁻¹.
- Error bars: standard error of the pre-t0 data.
- Assumed n ∝ T_ph (high-T limit: n ≈ kBT/ħω − 1/2), with initial T = 300 K.

**Model: three-temperature subset of the NLM** (Waldecker et al. PRX 6, 021003)
- Ce dTe/dt = f(t − t0) − G_{e,A1'}(Te − T_A1') − G_{e,l}(Te − Tl)
- C_A1' dT_A1'/dt = G_{e,A1'}(Te − T_A1') − G_{A1',l}(T_A1' − Tl)
- Cl dTl/dt = G_{e,l}(Te − Tl) + G_{A1',l}(T_A1' − Tl)

**Heat capacities**
- Electronic: Ce(Te) = 13.8 Te + 1.16e-3 Te² + 2.6e-7 Te³, from Nihira & Iwata PRB 68, 134305.
  - Units are not printed. The linear coefficient matches graphite's electronic γ ≈ 13.8 µJ mol⁻¹ K⁻², so molar units are likely. Verify against the source before use.
- Mode heat capacity: Einstein/Debye form, Eq. (16). For A1', D(ω) = δ(ω − ω_A1').
- Effective "lattice" Cl: built from Bonini et al. (PRL 99, 176802) A1' decay probabilities: 9% TA+TA, 36% TA+LA, 55% LA+TA or LO+LA.
- Electronic thermalization: about 100 fs, after which a single Te is assumed.

**Table I (stated), in W m⁻³ K⁻¹**

| Coupling | Value |
|---|---|
| G_{e,A1'} | (6.8 ± 0.3) x 10¹⁷ |
| G_{A1',l} | (8.0 ± 0.5) x 10¹⁷ |
| G_{e,l} | (0.0 ± 6.0) x 10¹⁵ |

**Mode-projected electron-phonon coupling**
- ⟨g²_{e,A1'}⟩ = 0.035 ± 0.001 eV² (Appendix B).
- Appendix B relations:
  - 1/τ_{e,j} = G_{ep,j}/Ce − Σi G_{ep,i}/C_{ph,j}
  - ħ/τ_{e,j} = 2π ⟨g²_{ep,j}⟩ De(ħω_ν − ħω_{j,k}), with ħω_ν = 1.55 eV
  - De(ε) = 2A|ε|/(π(ħvF)²), with vF = 9.06e5 m/s and A the unit-cell area.
- Stated to agree with trARPES (Johannsen 2013, Stange 2015, Rohde 2018) and Na et al. 2019.

**Fig. 7, graphical**
- Axes:
  - Left: Δn_A1' from −0.25 to 1.1 a.u.
  - Right: temperature from about −700 to 5200 K. Δn = 0 corresponds to about 300 K; Δn = 1.0 corresponds to about 4500 K.
  - x: −2 to 10 ps.
- Data:
  - Rise from t0 to the maximum: Δn = 1.0 at about 0.6 ps. Data at 0.35 ps ≈ 0.65; at 0.15 ps ≈ 0.34.
  - Decay: about 0.85–0.95 at 0.8–1.2 ps; about 0.6 at 2 ps; about 0.45 at 3 ps; about 0.3 at 4–5 ps; scattered 0.25–0.4 for 5–10 ps.
- Fit (pink):
  - Peak Δn ≈ 0.96 at about 0.7 ps, i.e. T_A1' ≈ 4300–4400 K.
  - Decays to a plateau of about 0.30 (≈ 1550 K) by about 6–8 ps.
  - Fitted 1/e decay time from the peak: about 1.5–1.7 ps (graphical).
  - The 10–90% rise time is about 0.3–0.4 ps, limited by Te and the 150 fs probe.
- Lattice Tl (orange dotted, the effective temperature of A1' decay products): rises from 300 K starting at about 0.3 ps, reaches its half-rise at about 1.5–2 ps, and plateaus near the A1' value (about 1500 K) by about 6 ps.
- Inset (−200 to 1200 fs, 0 to about 25,000 K):
  - Te peaks at about 25,000 K at about 50 fs, then falls roughly linearly to about 4000 K at about 800 fs, where it meets T_A1'.
  - T_A1' rises to about 4400 K at about 700–800 fs.
  - Tl reaches about 800–1000 K by 1.2 ps.
  - These temperatures follow the model. They are much higher than a simple Ce-based estimate from "8 J m⁻³" would give; see the units caveat in B.1.
- E2g at Γ: not extractable in this paper (|k| > 0.45 Å⁻¹ only). No E2g rise or decay times are quoted. No ZA or ZO data (geometrically invisible).

### B.6 Same data set, time constants from the companion paper B' (Stern et al. PRB 97, 165416, 2018)
Same conditions: 35 fs, 800 nm, 12 mJ/cm², graphite. These are single-exponential fits to raw diffuse intensity ΔI(q,τ), not decomposed populations.

| Feature | Value |
|---|---|
| K-point TO (A1') rise | 280 fs |
| K-point TO (A1') decay | 1.7 ps |
| K-point LO rise (where A1' is symmetry-forbidden) | 730 fs |
| Γ-E2g LO rise (near {200}) | 430 fs (time-resolved Raman at 0.2 mJ/cm² shows a faster rise, matching K-A1') |
| M-point rise | 2.1 ps |
| ½K / ½M TA-LA and TA-TA channel rise | 1.0–1.5 ps |
| Rise near {2̄10} (LO-LA and LA-TA decay of A1') | 1.5 ps |
| Intraband TA/LA thermalization | 34 ps (30–40 ps) |
| Interband LA/TA/ZA relaxation (diffuse-halo decay) | 115 ps |
| EPC into A1'/E2g | 200–400 fs |
| Decay out of A1'/E2g | 1–3 ps |
| Effective T_TO,K | about 1300 K at 1 ps; about 450 K at 10 ps (Bose-Einstein inversion of ΔI) |
| Electron thermalization (literature) | about 50 fs |
| Heat transport out of the probed volume | about 10 µs |

Note: B' estimates about 1300 K at 1 ps. Paper B's NLM curve (Fig. 7) implies about 4000 K at 1 ps. The two papers use different observables and normalizations.

### B.7 Simulation set-up summary (B)
- **Initial state:** 70 nm graphite, T0 = 300 K, 35 fs pump at 1.55 eV, 12 mJ/cm².
- **NLM:** use Table I couplings with C_A1' from Eq. (16), and Ce as above (check units).
- **Targets, Fig. 7:** normalized Δn_A1' peak at about 0.6–0.7 ps; plateau about 0.30 of peak by about 6 ps; decay τ ≈ 1.5–1.7 ps.
- **Targets, B':** A1' rise 280 fs; E2g rise 430 fs; acoustic thermalization 34 ps; interband 115 ps.
- **Targets, Fig. 6 (qualitative ordering):** TO2 at K peaks at about 0.5–1.5 ps. TA at M and ⅓M grows until 100 ps. LA near Γ grows monotonically to 100 ps.

---

## Paper C: Katsumi et al., amplitude mode in MgB2 by THz 2DCS

### C.1 Sample and equilibrium optics (SM §I, Table I, Figs. S1, S2)
- **Sample:**
  - 20 nm MgB2 film on single-crystal MgO(111), grown by HPCVD (Xi group, Temple). Tc = 38 K.
  - Normal-state resistivity: about 50 µΩ cm, giving k_F l ≈ 13.5.
- **Two-component Mattis-Bardeen parameters** (SM Table I, taken from Fiore et al. PRB 106, 094515 (2022)). ωp and τp are held fixed at all T.

| Band | 2Δ(0)/2π (THz) | 2Δ(0) (meV, derived) | ωp (eV) | 1/τp (meV) |
|---|---|---|---|---|
| π | 1.1 | 4.55 | 5.2 | 12.5 |
| σ | 3.4 | 14.1 | 4.2 | 8.33 |

- Derived: Δπ(0) ≈ 2.3 meV; Δσ(0) ≈ 7.0 meV; 2Δπ/kBTc ≈ 1.39; 2Δσ/kBTc ≈ 4.3. Conversion: 1 THz = 4.136 meV.
- **Two-band BCS gap equations** (SM §II):
  - Δα = Σβ Vαβ Δβ F(Δβ), with F(Δβ) = Nβ ∫₀^{ωc} dξ tanh(√(ξ² + Δβ²)/2kBT)/√(ξ² + Δβ²).
  - Parameters (from Liu, Mazin, Kortus PRL 87, 087005):
    - λ11 = V11 N1 = 0.28 (π)
    - λ22 = V22 N2 = 0.96 (σ)
    - λ12 = λ21 = V12 √(N1 N2) = 0.19
    - N2/N1 = 0.73
  - Cut-off "ωc = 1.78 THz" (citing Xi, Rep. Prog. Phys. 71, 116501). This value is as printed. A phonon cutoff of 1.78 THz (7.4 meV) lies below 2Δσ and is likely a misprint, e.g. for about 17.8 THz ≈ 74 meV. Choose ωc so that Tc = 38 K and 2Δ(0) match the table.
- **"Lower-gap temperature curve"** 2Δπ(T)/2π: the main-text Fig. 3c and SM Figs. S2 and S7c show circles from the MB fits of σ(ω) plus the solid BCS curve.
  - Graphical data points: 1.10 (5 K), 1.05 (12 K), 0.87 (20 K), 0.76 (25 K), 0.68 (30 K), 0.57 (32 K), 0.45 (34 K), 0.35 (36 K), 0 (≥ 38 K).
  - The solid BCS curve sits slightly below the points above 20 K: about 0.70 at 25 K, about 0.55 at 30 K.
  - Ω/2π = 0.63 THz crosses 2Δπ at about 27–28 K (BCS curve) or about 31 K (data).
  - 2Ω/2π = 1.26 THz never matches 2Δπ, which is always below 1.1 THz.
- **2Δσ/2π (Fig. S2, BCS curve and STM triangles from Iavarone et al. PRL 89, 187002):** about 3.4 THz at 5 K; about 3.3 at 20 K; about 2.9 at 28 K; about 2.4 at 32 K; about 1.0 at 37 K; 0 at 38 K (graphical).
  - The σ-band resonance Ω = 2Δσ would occur just below Tc, at about 37.5 K.
- **Fig. S1, equilibrium σ(ω), 0.25–3 THz, at T = 5, 12, 20, 25, 30, 32, 34, 36, 38 K (graphical):**
  - σ1 at 5 K: dips to about 0 (slightly negative within noise, about −0.1e4 Ω⁻¹cm⁻¹) near 0.8–1.0 THz, rises to about 0.8e4 at 2–3 THz, with a weak low-frequency upturn below 0.4 THz.
  - σ1 at 38 K: Drude-like, about 2.3e4 at 0.3 THz falling to about 1.0e4 at 3 THz.
  - σ2 at 5 K: diverges ∝ 1/ω (about 8e4 at about 0.3 THz, about 2e4 at 1 THz, about 1e4 at 2–3 THz).
  - σ2 at 38 K: about 1e4, nearly flat.
  - Power transmittance: 5 K about 1e-4 at 0.5 THz rising to about 4.4e-4 at about 1.5 THz; 38 K about 1.2e-4 at 0.5 THz rising to about 2–2.5e-4 at 2–3 THz.

### C.2 THz sources and fields (SM §III)
- Ti:sapphire regenerative amplifier: 800 nm, 100 fs, 9 mJ, 1 kHz.
- Two tilted-pulse-front LiNbO3 THz sources (A and B), combined on a Si beam splitter.
- EO sampling in 1 mm ZnTe(110).
- **Broadband peak fields:** A up to 150 kV/cm, B up to 80 kV/cm. For 2DCS the A pulse is attenuated to equal B.
  - Drive spectra (Fig. 2a top): both peak at about 0.55–0.6 THz and fall to about 0 by 1.6–1.8 THz.
  - Field limit: peak fields must stay below 21 kV/cm to avoid suppressing superconductivity.
- **Narrowband:** Ω/2π = 0.63 THz, using bandpass filters (2 on A, 1 on B). Maximum peak field: 30 kV/cm for both pulses.
  - Fields used: 12, 21 and 30 kV/cm. 25 kV/cm is mentioned for Fig. 3d in the text, but the panels show 12, 21 and 30.
- **Nonlinear signal:** E_NL(t,τ) = E_AB − E_A − E_B, where t is the sampling time and τ the A-B delay.

### C.3 Broadband 2DCS results (Fig. 1, Fig. 2, Figs. S6, S7)
- **Fig. 1d (7 K, 21 kV/cm, t = 0.36 ps):** E_NL(τ) has a single peak of about 1.7 a.u. at τ ≈ 0. It returns to about 0 with no long-lived component.
  - Above 21 kV/cm a long-lived quasiparticle component appears.
- **Fig. 2a bottom, nonlinear power spectrum I_NL at τ = 0 (graphical):**

| T (K) | Peak frequency (THz) | Peak I_NL (a.u.) |
|---|---|---|
| 7 | 1.05 | 2.35 |
| 12 | 1.03 | 1.9 |
| 17 | 0.95 | 1.55 |
| 22 | 0.9 | 0.9 |
| 27 | 0.8 | 0.75 |
| 30 | 0.7 | 0.65 |
| 32 | 0.65 | 0.5 |
| 36 | 0.6 | 0.3 |
| 37 | 0.55 | 0.2 |
| 42 | – | about 0 |

- At 7 K the line spans about 0.6–1.6 THz (FWHM about 0.5 THz). Above about 20 K the spectra broaden greatly.
- Stated: the peak "at 1 THz at 7 K" coincides with 2Δπ from σ(ω).
- **Fig. 2b:** colour map of I_NL(T, f), colour scale 0 to about 2.5, overlaid with the 2Δπ(T) BCS curve.
- **Fig. S7c, Lorentzian peak positions (graphical):**
  - I_NL: 1.05 (7 K), 1.05 (12), 0.97 (17), 0.92 (22), 0.87 (27), 0.80 (30), 0.70 (32), 0.72 (34), 0.72 (36), about 0.62 (37) THz.
  - I_B³ (cube of transmitted B-pulse power): 0.97, 0.93, 0.90, 0.87, 0.85, 0.80, 0.77, 0.70, 0.65, about 0.6, and about 0.6 at 42 K.
  - Below 20 K I_NL follows 2Δπ. Above 20 K it tracks I_B³, which indicates an overdamped mode. I_NL is wider than I_B³.
- **Fig. S7a:** the source spectrum peaks at about 0.6 THz. Multiplied by T(7 K), it peaks at about 1.0 THz, as do the transmitted A and B spectra at 7 K. This is why the S7b,c comparison is needed.
- **Fig. S6 (2D spectrum, 7 K, broadband):** strongest features at (ωt, ωτ)/2π ≈ (1, 0) and along the diagonal near (1, 1). The rephasing signal does not lie exactly on ωτ = −ωt, so homogeneous and inhomogeneous broadening cannot be separated.

### C.4 Narrowband FH and TH results (Fig. 3, Fig. 4, Figs. S4, S5, S8, S9)
**Fig. 3a (7 K, 12 kV/cm, τ = 0)**
- FH peak at 0.63 THz, amplitude about 0.67 a.u. TH peak at 1.9 THz, amplitude about 0.36 a.u. (graphical).
- TH/FH peak ratio: about 0.54.

**Fig. 3b, field dependence at 7 K** (log-log; dashed guides ∝ E⁶)
- Stated: FH follows E⁶ up to 21 kV/cm and TH up to 25 kV/cm, i.e. both are χ⁽³⁾ processes.
- Graphical I_NL (a.u.):

| E (kV/cm) | FH | TH |
|---|---|---|
| 12 | 5e-7 | 2.8e-7 |
| 15 | 2e-6 | 1.2e-6 |
| 18 | 1e-5 | 3.5e-6 |
| 21 | 3.5e-5 | 9e-6 |
| 25 | 2.1e-4 | 1.6e-5 |
| 30 | 4.5e-4 | 2.2e-5 |

- Above about 21–25 kV/cm, FH runs above E⁶ and TH falls below it (saturation).

**Integration windows:** FH 0.3–1.0 THz; TH 1.6–2.2 THz.

**Fig. 3d, raw FH vs T** (curves offset by +1 per field; 12 kV/cm baseline 0, 21 kV/cm baseline 1, 30 kV/cm baseline 2; each roughly normalized to a peak of 1; graphical)

| Field | Low-T value | Peak | Above Tc |
|---|---|---|---|
| 12 kV/cm | 0.28 at 8 K; 0.3 at 17 K; 0.45 at 22 K; 0.7 at 28 K | about 1.0 at about 35–36 K | 0 at ≥ 38 K |
| 21 kV/cm | 0.3 at 8 K | about 1.0 at about 32–33 K | 0 above about 37 K |
| 30 kV/cm | 0.63 at 8 K | about 1.0 at about 25–27 K | 0 above about 36 K |

- Stated: the raw FH peak is not at Ω = 2Δπ, and the peak temperature drops as the field increases.

**Fig. 3e, raw TH vs T** (same offsets; graphical)
- 12 kV/cm: 0.85 at 8 K, 1.0 at 12 K, about 0.65 at 17–20 K, about 0.5 at 25 K, about 0.2–0.4 at 28–35 K, 0 at 38 K.
- 21 and 30 kV/cm: monotonic decrease from about 1.0 at 8 K to 0 at Tc.
- Stated: consistent with earlier THG (Kovalev et al. PRB 104, L140505), where the maximum is toward 2Ω = 2Δπ.

**Fig. S8, transmitted field at 0.63 THz** (12 kV/cm; graphical)
- A pulse: 0.95 (7 K), 1.0 (12), 1.03 (17), 1.1 (20–22), 1.2 (25), 1.3 (28), 1.4 (30), 1.5 (33–35), maximum 1.7 at about 37 K, 1.67 above Tc.
- B pulse: 0.6 (7 K) → 1.05 at Tc (B is about 2x weaker at the detector because of focusing).
- Derived: (E_A(7 K)/E_A(normal))⁶ ≈ (0.95/1.67)⁶ ≈ 0.034. The normalization therefore enhances the 7 K value by about 30x relative to Tc.

**Fig. 4, I_NL/E_A⁶ vs T** (offsets as in Fig. 3; fit curves are solid)
- FH: monotonic increase on cooling at all fields. Graphical values (offset removed) for 12 kV/cm:
  - 1.0 (7.5 K), 0.83 (12), 0.68 (17), 0.65 (20), about 0.5 (22–24), 0.32 (27–28), about 0.2 (30–32), about 0.1 (34–36), 0 (≥ 38 K).
  - Half value at about 23–24 K.
  - 21 and 30 kV/cm have similar shapes, reaching about 0.8–1.0 at 7–12 K.
- TH: also monotonic but steeper. Graphical values for 12 kV/cm:
  - 1.0 (8 K), 0.9 (12), 0.5 (17), 0.3 (20), 0.2 (22), 0.1 (25), about 0 (≥ 28 K).
  - Half value at about 17 K.
  - The 21 and 30 kV/cm TH curves fall to about 0 by about 28–30 K.
- The same trend holds for E_B⁶ normalization (Fig. S9, 12/21/30 kV/cm, each normalized to its maximum).
- NbN contrast: NbN's normalized FH shows a resonant peak at Ω = 2Δ.

**Model fit (Eq. 1), as printed**
- I(ω,T) = I0 Δπ(T)² / [(ω + iδ)² − (2Δπ(T))²]. The modulus is presumably taken for intensity; Ref. [20] (Katsumi et al. PRL 132, 256903) defines the convention.
- Settings: ω = Ω for FH and ω = 2Ω for TH, with Δπ(T) from the two-band BCS curve.
- Fitted damping:
  - 12 kV/cm: δ = 0.55 THz, from a simultaneous FH+TH fit.
  - 21 and 30 kV/cm: same δ. At 30 kV/cm, 2Δπ is scaled by 0.81.
  - NbN comparison: δ = 0.12 THz.
- Interpretation: the π-band amplitude mode is overdamped, attributed to interband coupling and to the higher resonance temperature, about 2x that in NbN.

**Figs. S4/S5 (narrowband, 7 K, 12 kV/cm)**
- E_NL is oscillatory, with essentially no long-lived component.
- The 2D spectrum for τ ≥ 0 shows nonrephasing (Ω, Ω) and rephasing (Ω, −Ω) peaks.

**Leggett mode and theory scope**
- The Leggett mode is discussed only as motivation (refs. 21–25, 32). No Leggett-mode frequency is fitted or reported.
- Theory framework cited: Fiore et al. PRB 106, 094515 (2022). In that work, for MgB2 parameters, the π-band amplitude mode is strongly suppressed and the σ-band one is pronounced at T = 0. Paramagnetic coupling dominates in dirty superconductors (Murotani & Shimano PRB 99, 224510).

### C.5 Supplemental Fig. S3: THz pump / THz probe at 7 K (SM §IV)
**Setup**
- Same apparatus with a weak broadband probe (peak 10 kV/cm).
- δE_probe recorded at a fixed sampling time t = 0.36 ps.

**Fig. S3a, δE_probe(τ)** (graphical; τ = −2 to 5 ps)

| Pump (kV/cm) | Onset | Level |
|---|---|---|
| 150 | τ ≈ 0.4–0.5 ps | Rises to about 6.5 at about 1.7 ps, then stays about 6 to 5 ps (long-lived) |
| 79 | – | About 3.5–4 |
| 68 | – | About 2 |
| 47 | – | About 0.8 |
| 21 (narrowband, with BPF) | – | About 0.2–0.3, nearly undetectable |

**Fig. S3b, maximum δE_probe vs E_pump** (graphical)
- 21 (BPF): 0.2; 47: 1.0; 60: 1.4; 68: 2.1; 79: 3.8; 87: 5.4; 116: 6.6; 150: 6.6.
- Threshold-like growth above about 60 kV/cm; saturation at ≥ about 116 kV/cm.

**Fig. S3c, σ1(ω) at τ = 2 ps, 7 K, 0.5–2.0 THz**, in 10⁴ Ω⁻¹cm⁻¹ (graphical)
- 150 kV/cm: 1.6 (0.5 THz) → maximum about 2.55 at about 0.8 THz → 1.9 (1.1) → 1.5 (1.4) → 1.0 (1.6) → about 0.8 (1.85) → 1.0 (2.0).
  - Stated: "close to that at 40 K in equilibrium".
- 79 kV/cm: 0.7 (0.5) → 1.2 (0.65–0.85) → about 0.9 (1.2) → 0.75 (1.5–1.9) → 1.1 (2.0).
- 68 kV/cm: 0.25 (0.5) → 0.65 (0.65–0.85) → 0.4 (1.0–1.3) → 0.6–0.7 (1.7–2.0).
- 21 kV/cm (BPF): about 0 (0.5) → 0.25 (0.7–0.85) → about 0 (1.0) → 0.6 (1.9) → 0.2 (2.0). This is essentially equal to equilibrium 5 K.
- Equilibrium dashed curves:

| T | σ1 values |
|---|---|
| 40 K | 2.2 (0.5 THz) decreasing smoothly to 1.35 (2.0 THz) |
| 25 K | 1.6 (0.5), about 1.0 (0.9–1.2), about 1.05 (2.0) |
| 20 K | about 1.0 (0.5), 0.55 (1.0), 0.85 (2.0) |
| 12 K | 0.6 (0.5), 0.2 (0.9), 0.75 (2.0) |
| 5 K | 0.2 (0.5), about −0.1 (0.85–1.0), 0.7 (2.0) |

**Fig. S3d, σ2(ω) at τ = 2 ps, 7 K**, in 10⁴ Ω⁻¹cm⁻¹ (graphical)
- 150 kV/cm: 2.6 (0.5 THz) → 1.4 (0.65) → about 1.5–1.6 (0.8–1.2) → 1.0 (1.8) → 0.5 (2.0).
  - Stated: "strongly suppressed but larger than that at 40 K".
- 79 kV/cm: 3.7 (0.5) → 3.0 (0.75) → 2.5 (1.0) → 1.0 (1.8).
- 68 kV/cm: 4.7 (0.5).
- 21 kV/cm (BPF): 5.6 (0.5), equal to 5–12 K equilibrium.
- Equilibrium at 0.5 THz: 5 K 6.0; 12 K about 5.3; 20 K 4.3; 25 K 3.6; 40 K 0.6 (rising to about 1.0 at 1–2 THz).
- All curves converge to about 1.0 near 1.5–2.0 THz.

**Conclusion (stated):** 150 kV/cm strongly depletes the condensate, so the signals at that field come from quasiparticles. The narrowband ≤ 30 kV/cm drive leaves σ1 and σ2 "almost intact" at τ = 2 ps.

### C.6 Simulation set-up summary (C)
- **Film:** 20 nm MgB2 on MgO, Tc = 38 K.
- **Gaps:** two-band BCS with λππ = 0.28, λσσ = 0.96, λπσ = 0.19, Nσ/Nπ = 0.73. Targets: 2Δπ(0) = 1.1 THz (4.5 meV) and 2Δσ(0) = 3.4 THz (14 meV).
- **σ(ω):** two-component Mattis-Bardeen with ωp = 5.2 eV (π) and 4.2 eV (σ); 1/τ = 12.5 meV (π) and 8.33 meV (σ).
- **Drive:** multicycle 0.63 THz at 12, 21 and 30 kV/cm peak (third-order regime ≤ 21 kV/cm). Broadband single-cycle drive peaks at about 0.6 THz.
- **Targets:**
  - Broadband I_NL peak at 1.05 THz at 7 K, following 2Δπ(T) to about 20 K, then washing out.
  - FH and TH ∝ E⁶.
  - Normalized FH and TH monotonic in T, fitted by Eq. 1 with δ = 0.55 THz.
  - Field-dependent gap reduction factor 0.81 at 30 kV/cm.
  - Pump-probe depletion threshold above about 50 kV/cm, with σ1 at 2 ps after 150 kV/cm close to the 40 K normal state.
