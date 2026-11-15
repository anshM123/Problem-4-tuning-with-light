# Black phosphorus: published UED / UEDS reference measurements

This file collects measured values from two publications by the same group, using the same sample and the same pump conditions. It was put together without reading `preregistered_predictions.json` or `results.json`, so the numbers here are independent of those files.

## Sources

- **[Z20]** D. Zahn, P.-N. Hildebrandt, T. Vasileiadis, Y. W. Windsor, Y. Qi, H. Seiler, R. Ernstorfer, "Anisotropic nonequilibrium lattice dynamics of black phosphorus", Nano Lett. 20, 3728 (2020), doi:10.1021/acs.nanolett.0c00734. Values were taken from arXiv:2002.12038v3 (2 Apr 2020). Page numbers refer to the arXiv PDF, which numbers its pages at the top.
  - Data set: Zenodo doi:10.5281/zenodo.3723318. Its `measurement_1_metadata.txt` file describes the 100 K run, and `measurement_2_metadata.txt` describes the 295 K run.
- **[S21]** H. Seiler, D. Zahn, M. Zacharias, P.-N. Hildebrandt, T. Vasileiadis, Y. W. Windsor, Y. Qi, C. Carbogno, C. Draxl, R. Ernstorfer, F. Caruso, "Accessing the anisotropic nonthermal phonon populations in black phosphorus", Nano Lett. 21, 6171 (2021), doi:10.1021/acs.nanolett.1c01786.
  - Values were checked against both arXiv:2006.12873v3 (27 Apr 2021) and the published version.
  - Page numbers are given as `pub p.XXXX` for the journal and `arXiv p.N` for the preprint. The Supporting Information is cited as `SI p.N` and is appended to the arXiv PDF.

## Experimental conditions (common to both papers, 100 K data set)

| Quantity | Value | Source |
|---|---|---|
| Sample | Free-standing exfoliated film, thickness **39 ± 5 nm** (from optical transmission and transfer-matrix calculation; about 80 layers) | Z20 p.2: "we estimate the film thickness to be 39 ± 5 nm"; S21 pub p.6172 / arXiv p.5 |
| Base temperature | **100 K** | Z20 p.1: "All measurements are performed at a base temperature of 100 K"; S21 pub p.6172 |
| Pump | 770 nm (1.61 eV), about 50 fs, linear polarization along armchair, 10 deg off-normal, 4 kHz, spot about 200 µm | Z20 p.1; S21 pub p.6172: "photoexcited with a 50 fs light pulse with energy hν = 1.61 eV"; Zenodo metadata (m1) |
| Incident fluence | **9.8 ± 1.4 mJ/cm²**. S21 states the same value as 98 J/m², without an uncertainty. | Z20 Fig. 2 caption: "incident fluence of (9.8 ± 1.4) mJ/cm2"; S21 SI p.2: "incident fluence on the BP flake is Iinc = 98 J/m2" |
| Probe | Z20: 70 keV electrons, transmission at normal incidence. S21 estimates the electron pulse duration at 200 fs. | Z20 p.1; S21 pub p.6172: "duration of the electron pulse is estimated to be 200 fs" |
| Time resolution / IRF | Z20 models the IRF as a Gaussian with FWHM = 150 fs, and the metadata give a time resolution of about 150 fs. S21 does not state an IRF for its diffuse-scattering fits. | Z20 p.2: "convolving the fit function with a Gaussian with a FWHM of 150 fs" |
| Delay range (Z20) | 104 delays from -8 to +300 ps, scanned 9 times | Zenodo metadata (m1) |

## Main table: measured observables

The two zone-boundary points around the (400) reflection are labelled as follows:
- **X** is reached by moving along the zigzag (h) axis.
- **A** is reached by moving along the armchair (l) axis.

These assignments come from the inset of S21 Fig. 2b.

| Observable | Definition / fit form | Measured value ± uncertainty | Fluence (incident) | Source location (short quote) |
|---|---|---|---|---|
| ΔMSD armchair, fast amplitude A₁ | ΔU₃₃(t) from (00l) Bragg intensities via Eq. 3, I/I₀ = exp[−ΔU₃₃·4π²l²/c²], c = 4.37 Å. Biexponential fit convolved with a 150 fs FWHM Gaussian. The functional form is not written out, but the amplitudes in Å² imply ΔU = A₁(1−e^{−t/τ₁}) + A₂(1−e^{−t/τ₂}). | **0.66 ± 0.03 ×10⁻³ Å²** | 9.8 ± 1.4 mJ/cm² | Z20 Table I, p.2: "biexponential function convolved with a Gaussian" |
| τ₁ armchair (electron-phonon) | Same fit | **0.48 ± 0.05 ps** | 9.8 | Z20 Table I, p.2 |
| ΔMSD armchair, slow amplitude A₂ | Same fit | **1.84 ± 0.03 ×10⁻³ Å²** | 9.8 | Z20 Table I, p.2 |
| τ₂ armchair (phonon-phonon) | Same fit | **22 ± 1 ps** | 9.8 | Z20 Table I, p.2 |
| ΔMSD zigzag, fast amplitude A₁ | ΔU₁₁(t) from (h00) intensities via Eq. 2, a = 3.31 Å. Same fit form and IRF as above. | **0.68 ± 0.02 ×10⁻³ Å²** | 9.8 | Z20 Table I, p.2 |
| τ₁ zigzag | Same fit | **0.58 ± 0.04 ps** | 9.8 | Z20 Table I, p.2 |
| ΔMSD zigzag, slow amplitude A₂ | Same fit | **0.85 ± 0.03 ×10⁻³ Å²** | 9.8 | Z20 Table I, p.2 |
| τ₂ zigzag | Same fit | **20 ± 2 ps** | 9.8 | Z20 Table I, p.2 |
| Fit uncertainty convention | Errors are 68.3 % confidence intervals of the fit. They do not include the fluence or thickness uncertainty. | — | — | Z20 Table I caption, p.2 |
| Asymptotic ΔMSD armchair (t ≫ τ₂) | A₁ + A₂, derived here from Table I | **2.50 ± 0.04 ×10⁻³ Å²**. Fit value at 100 ps: 2.48 ×10⁻³ Å². Figure reading: data plateau about 2.5 ×10⁻³ Å² (±0.05). | 9.8 | Derived from Z20 Table I; consistent with Z20 Fig. 2c, p.3 |
| Asymptotic ΔMSD zigzag (t ≫ τ₂) | A₁ + A₂, derived here | **1.53 ± 0.04 ×10⁻³ Å²**. Fit value at 100 ps: 1.52 ×10⁻³ Å². Figure reading: about 1.5 ×10⁻³ Å² (±0.05). | 9.8 | Derived from Z20 Table I; Fig. 2c, p.3 |
| Late-time anisotropy ΔU₃₃/ΔU₁₁ | Ratio of the asymptotes (derived here) | **1.63 ± 0.05** | 9.8 | Derived. Z20 p.2: "higher MSD increase in the armchair direction" |
| Fast-component anisotropy A₁(ac)/A₁(zz) | Derived | **0.97 ± 0.05** (equal within error) | 9.8 | Z20 p.2: "amplitudes of the initial MSD rise are the same" |
| Slow-component anisotropy A₂(ac)/A₂(zz) | Derived | **2.16 ± 0.08** | 9.8 | Derived from Z20 Table I |
| Fraction of the total ΔMSD in the fast step | A₁/(A₁+A₂), derived | Armchair **0.26**; zigzag **0.44** | 9.8 | Derived from Z20 Table I |
| Non-thermal lattice lifetime | Delay after which MSD normalized to its 100 ps value is the same for both directions | **≈ 60 ps** (stated without uncertainty) | 9.8 | Z20 p.3: "a thermal state is reached after ~60 ps"; Fig. 2d |
| Pre-excitation MSD (100 K) | Not measured. Estimated in Z20 by assuming MSD ∝ T together with the 270 K rise. Read off Fig. 3a. | U₃₃ ≈ 0.92 ×10⁻³ Å², U₁₁ ≈ 0.57 ×10⁻³ Å² (reading ±0.03 ×10⁻³ Å²) | 9.8 | Z20 Fig. 3a, p.4. Model-based estimate, not a measurement. |
| MSD at 100 ps (absolute) | Same estimate | U₃₃ ≈ 3.4 ×10⁻³ Å², U₁₁ ≈ 2.1 ×10⁻³ Å² (reading ±0.05) | 9.8 | Z20 Fig. 3a, p.4 |
| (200) relative intensity, plateau at 100–300 ps | I/I₀ of the Bragg peak, average of the Friedel pair, 2D pseudo-Voigt integrated intensity | **≈ 0.982** (figure reading ±0.003) | 9.8 | Z20 Fig. 2b, p.3 (read off figure) |
| (400) relative intensity, plateau | Same | **≈ 0.911** (±0.003) | 9.8 | Z20 Fig. 2b, p.3 (read off figure) |
| (600) relative intensity, plateau | Same | **≈ 0.834** (±0.003) | 9.8 | Z20 Fig. 2b, p.3 (read off figure) |
| (002) relative intensity, plateau | Same | **≈ 0.976** (±0.003) | 9.8 | Z20 Fig. 2b, p.3 (read off figure) |
| (004) relative intensity, plateau | Same | **≈ 0.910** (±0.003). It reaches the same level as (400) from about 100 ps onward. | 9.8 | Z20 Fig. 2b, p.3 (read off figure) |
| (006) relative intensity, plateau | Same | **≈ 0.844** (±0.003) | 9.8 | Z20 Fig. 2b, p.3 (read off figure) |
| (400) and (004) Bragg dynamics (S21 data) | Biexponential decay. Constants are quoted only as approximate values that point back to Z20. No uncertainties, and no IRF statement. | τ_fast ≈ **500 fs**, τ_slow ≈ **20 ps** | 9.8 (98 J/m²) | S21 pub p.6172 / arXiv p.5: "fast time constant of 500 fs and a slower time constant of 20 ps" |
| (400) / (004) I/I₀ at 2 ps (S21) | Elastic signal, average of the Friedel pair | (400) ≈ **0.961**; (004) ≈ **0.977** (reading ±0.003) | 9.8 | S21 Fig. 2a inset, pub p.6173 (read off figure) |
| (400) / (004) I/I₀ at 10 ps (S21) | Same | (400) ≈ **0.939**; (004) ≈ **0.953** (±0.003) | 9.8 | S21 Fig. 2a, pub p.6173 (read off figure) |
| (400) / (004) I/I₀ at 100 ps (S21) | Same | (400) ≈ **0.921**; (004) ≈ **0.924** (±0.003) | 9.8 | S21 Fig. 2a, pub p.6173 (read off figure) |
| Diffuse intensity at A next to (400): rise | Relative diffuse intensity I(t)/I(t<0) integrated at A. Biexponential fit (rise, then decay). The exact functional form and IRF treatment are not given. | τ_rise = **1.7 ± 0.1 ps** | 9.8 | S21 pub p.6173 / arXiv p.6: "rising time constant of 1.7±0.1 ps, followed by a slower relaxation of 30±2 ps" |
| Diffuse intensity at A next to (400): decay | Same fit | τ_decay = **30 ± 2 ps** | 9.8 | Same sentence. The paper says "similar dynamics are observed at all the investigated A points." |
| Diffuse intensity at A: amplitudes | Read off the fit curve | Maximum ≈ **+1.7 %** (I/I₀ ≈ 1.017) at about 5–10 ps, relaxing to ≈ **+1.4 %** (1.014) by 100 ps. Reading ±0.2 %. | 9.8 | S21 Fig. 2b, pub p.6173 (read off figure) |
| Diffuse intensity at X next to (400): rise | Single-exponential fit (rise only). The exact form, offset and IRF treatment are not given. | τ = **14.3 ± 0.1 ps** | 9.8 | S21 pub p.6173 / arXiv p.6: "exponential function yields a time constant of 14.3 ± 0.1 ps" |
| Diffuse intensity at X: amplitude | Read off the fit curve | Plateau ≈ **+3.9 %** (I/I₀ ≈ 1.039), reached by about 50 ps. Reading ±0.2 %. Below about 15 ps the A and X data overlap. | 9.8 | S21 Fig. 2b, pub p.6173 (read off figure) |
| Elastic vs diffuse timescales | — | The 1.7 ps rise at A does not appear in the Bragg (elastic) analysis | 9.8 | S21 p.6173: "1.7 ps time constant does not appear in an elastic scattering analysis" |
| Diffuse maps | I(Q,t) − I(Q,t<t₀) at 2, 10 and 50 ps. At 2 ps: faint lines along Γ–A. At 50 ps: stronger, diamond-shaped and anisotropic. No numbers are given. | Qualitative only | 9.8 | S21 Fig. 3a–c, pub p.6174 |
| Return to thermal equilibrium (S21 summary) | — | "on time scales of the order of **50–100 ps**" | 9.8 | S21 pub p.6171–6172 / arXiv p.3 |
| Absorbed energy density | From optical constants (Jiang et al. 2018) and film thickness | **380 ± 70 J/cm³** | 9.8 | Z20 Fig. 2 caption, p.3: "absorbed energy density to be (380 ± 70) J/cm3" |
| Absorbed fluence / fraction | Derived here: 380 J/cm³ × 39 nm | ≈ 1.5 mJ/cm² absorbed, about 15 % of the incident fluence | 9.8 | Derived |
| Photoexcited carrier density | n = 10⁻⁶·(1/E_ph)·(I_abs/d)·c·½, using n_arm(800 nm) = 3.19 + 0.29i and transfer matrices | **(7.3 ± 0.9) ×10¹³ cm⁻²** per layer. As a volume density this is 1.42 ×10²¹ cm⁻³, the value used in the S21 theory. | 9.8 | S21 SI p.2; S21 pub p.6172: "ne = 7.3 × 1013 cm-2" |
| Lattice temperature rise (late times) | From the absorbed energy density and the lattice heat capacity (from the vibrational DOS) | **ΔT = 270 ± 50 K**, so T_final ≈ **370 ± 50 K** from a 100 K base | 9.8 | Z20 p.4: "temperature of the system rises by (270 ± 50) K" |
| Final lattice temperature (S21 statement) | Experimental estimate, used to fix the theory's initial T_el | **≈ 300 K** (no uncertainty given) | 9.8 | S21 pub p.6173: "matches the experimental estimate of 300 K" |
| Electron temperature | **Not measured** in either paper. S21 uses T_el⁰ = 7000 K as a model input, chosen so the final lattice temperature comes out at 300 K. | (7000 K: model input, not data) | — | S21 pub p.6173 |
| τ₁ vs fluence, polarization and base temperature | Same biexponential MSD fit | τ₁ ≈ **0.5 ps** in all cases (no table given) | Several fluences; the values are not listed | Z20 p.4: "we obtain values of ~0.5 ps for the fast time constants τ1" |
| τ₂ vs base temperature and fluence | Same fit | τ₂ gets shorter at higher base temperature and at higher excitation density (qualitative) | Not tabulated | Z20 p.4: "acceleration of the phonon thermalization with increasing excitation density" |
| τ₂ at 295 K, armchair | Same fit | **11.5 ± 0.7 ps** | About 2 mJ/cm². Per the Zenodo metadata: about 1.9 mJ/cm², 800 nm, polarization unknown, film 30–60 nm. | Z20 p.4: "τ2 values of (11.5±0.7) ps and (8±1) ps" |
| τ₂ at 295 K, zigzag | Same fit | **8 ± 1 ps** | Same as above | Z20 p.4 (same sentence) |
| Pump photon energy dependence (S21) | Diffuse maps with a 0.59 eV pump at the same assumed carrier density | Same qualitative thermalization. The 2 ps signal is slightly weaker than with the 1.61 eV pump. No fitted numbers. | Not stated | S21 SI p.3 |

## Theory values in S21 (computed, not measured; listed for context only)

| Quantity | Value | Source |
|---|---|---|
| Peak of the vibrational temperature T̃_q near Γ | 1.7 ps | S21 pub p.6174 / arXiv p.11: "reaches a maximum at 1.7 and 2.3 ps" |
| Peak of T̃_q near A | 2.3 ps | Same |
| T̃_q near X | No maximum (rises monotonically) | Same: "no maximum is observed around X" |
| Thermal equilibrium reached in the simulation | By about 40 ps, at T = 300 K. The text summary says "within 50 ps". | S21 pub p.6174–6176, Fig. 4f |
| Dominant electron-phonon channel | High-energy optical branches (ν = 10–12); phonons emitted along Γ–A (zigzag q) | S21 Fig. 4h–k |

## Derived reference curve (Z20 Table I fit evaluated without convolution; valid for t ≫ 0.15 ps)

The fit form is ΔU(t) = A₁(1−e^{−t/τ₁}) + A₂(1−e^{−t/τ₂}). The I/I₀ columns are computed from that ΔU using Eqs. 2 and 3 (a = 3.31 Å, c = 4.37 Å). Comparing them with the measured plateaus in the main table shows how far the reflections are from ideal Q² scaling.

| t (ps) | ΔU₃₃ armchair (10⁻³ Å²) | ΔU₁₁ zigzag (10⁻³ Å²) | ratio ac/zz | I/I₀ (400) | I/I₀ (004) | I/I₀ (600) | I/I₀ (006) |
|---|---|---|---|---|---|---|---|
| 0.5 | 0.47 | 0.41 | 1.13 | 0.976 | 0.985 | 0.948 | 0.966 |
| 1 | 0.66 | 0.60 | 1.10 | 0.966 | 0.978 | 0.925 | 0.952 |
| 2 | 0.81 | 0.74 | 1.10 | 0.958 | 0.974 | 0.909 | 0.942 |
| 5 | 1.03 | 0.87 | 1.19 | 0.951 | 0.966 | 0.894 | 0.926 |
| 10 | 1.33 | 1.01 | 1.31 | 0.943 | 0.957 | 0.877 | 0.906 |
| 20 | 1.76 | 1.22 | 1.44 | 0.932 | 0.944 | 0.854 | 0.877 |
| 50 | 2.31 | 1.46 | 1.58 | 0.919 | 0.926 | 0.827 | 0.842 |
| 100 | 2.48 | 1.52 | 1.63 | 0.916 | 0.921 | 0.821 | 0.831 |

The table below inverts the measured Z20 Fig. 2b plateaus one reflection at a time, using Eqs. 2 and 3. These are derived values, and the uncertainty comes from the ±0.003 figure reading.

| Reflection | ΔU (10⁻³ Å²) |
|---|---|
| (200) | 1.26 ± 0.21 |
| (400) | 1.62 ± 0.06 |
| (600) | 1.40 ± 0.03 |
| (002) | 2.9 ± 0.4 |
| (004) | 2.82 ± 0.10 |
| (006) | 2.28 ± 0.05 |

The published ΔMSD values are weighted averages over the Friedel pairs and reflections of each family: about 2.5 for armchair and 1.53 for zigzag.

## Caveats

1. **Fit definitions differ.**
   - Z20 fits the MSD (ΔU₁₁, ΔU₃₃), averaged over all (h00) or (00l) reflections. The fit is convolved with a 150 fs Gaussian IRF.
   - S21 quotes Bragg-intensity time constants only as approximate values (500 fs, 20 ps) taken from Z20.
   - S21's diffuse fits (biexponential at A, single exponential at X) give no explicit functional form, offset or IRF treatment.
   - A time constant fitted to MSD is not identical to one fitted to −ln I or to I itself. The difference is small here because the Debye–Waller exponent is below 0.2.
2. **Which reflection.**
   - Z20's MSDs combine (200), (400) and (600) for zigzag, and (002), (004) and (006) for armchair.
   - S21's elastic traces are (400) and (004) only, and its diffuse traces are taken in the Brillouin zone around (400) only.
   - The per-reflection plateaus in Z20 Fig. 2b do not follow Q² scaling exactly. For example, (004) ≈ (400) ≈ 0.911, whereas the averaged MSDs predict 0.921 and 0.916.
3. **Late-time Bragg drop differs between the two papers.** For (400) and (004) at 100 ps:
   - Z20 Fig. 2b gives about 0.911.
   - S21 Fig. 2a gives about 0.921–0.924, with the same stated sample, temperature and fluence.
   - S21 does not say whether it uses the same scan or the same intensity-extraction method. Treat about ±0.01 in late-time I/I₀ as the spread between analyses.
4. **Final temperature is inconsistent between the papers.** Z20 gives ΔT = 270 ± 50 K (T_final ≈ 370 K). S21 uses an "experimental estimate of 300 K". Neither paper measures the electron temperature; 7000 K is a model input.
5. **Fluence dependence.**
   - All quantitative 100 K numbers are for 9.8 ± 1.4 mJ/cm² incident (about 1.5 mJ/cm² absorbed, 380 ± 70 J/cm³).
   - The only other quantitative τ₂ values are at 295 K and about 2 mJ/cm² incident, where the pump was 800 nm with unknown polarization (Zenodo metadata).
   - The fluence dependence at 100 K is described only qualitatively: τ₁ stays at about 0.5 ps, and τ₂ gets shorter as the fluence increases.
6. **Figure readings.** All "read off figure" values come from the vector figures in the arXiv and journal PDFs, with the stated reading uncertainty.
   - Z20 Fig. 2b has a broken time axis: it is linear from 0 to 50 ps and compressed from 50 to 300 ps.
   - The plateau values listed are the fitted-curve levels at 100–300 ps.
7. **Published vs preprint.**
   - The S21 numbers are identical in arXiv v3 and the published version. The published version adds the 50 fs pump duration, the 200 fs electron pulse and the Fig. 2a inset.
   - The Z20 numbers come from arXiv v3, posted the same month as publication. The journal PDF was not accessible for line-by-line comparison.
