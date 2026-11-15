# Black phosphorus: preregistered prediction vs measurement

- Predictions: `preregistered_predictions.json`, committed in 9389d30 before the measured values
  were compiled. The file is unchanged.
- Measurements: `reference_measurements.md`, from Zahn et al. 2020 [Z20] and Seiler et al.
  2021 [S21].
- Numbers: `comparison.json`, written by `compare_with_measurements.py`.
- The "post hoc" column is filled only where the preregistered file contains the point the paper
  measured, but under a different label (see the convention audit below). Its values are read
  from the preregistered file and are not recomputed.
- No physical model, coupling or code was changed.

## Comparison

| Observable (fit form as in the paper) | Preregistered | Post hoc (label convention: points matched by coordinates) | Measured | Ratio pre / meas (post / meas) |
|---|---|---|---|---|
| ΔMSD zigzag τ₁ (biexp ⊗ 150 fs, Z20) | 0.81 ps | — | 0.58 ± 0.04 ps | 1.40 |
| ΔMSD zigzag τ₂ | 20.5 ps | — | 20 ± 2 ps | **1.03** |
| ΔMSD armchair τ₁ | 0.12 ps | — | 0.48 ± 0.05 ps | 0.25 |
| ΔMSD armchair τ₂ | 14.0 ps | — | 22 ± 1 ps | 0.64 |
| ΔMSD zigzag A₁ / A₂ (10⁻³ Å²) | 2.35 / 1.48 | — | 0.68 / 0.85 | 3.45 / 1.74 |
| ΔMSD armchair A₁ / A₂ (10⁻³ Å²) | 1.37 / 5.89 | — | 0.66 / 1.84 | 2.07 / 3.20 |
| ΔMSD at 100 ps, zigzag / armchair (10⁻³ Å²) | 3.87 / 7.38 | — | 1.53 / 2.50 | 2.53 / 2.95 |
| Late anisotropy ΔU_ac / ΔU_zz | 1.91 | — | 1.63 ± 0.05 | **1.17** |
| Fast fraction A₁/(A₁+A₂), zigzag / armchair | 0.61 / 0.19 | — | 0.44 / 0.26 | 1.40 / 0.72 |
| I/I₀ plateau (400) / (004) (Z20, 100–300 ps) | 0.800 / 0.784 | — | 0.911 / 0.910 | −ln ratio 2.4 / 2.6 |
| I/I₀ plateau (200) / (002) | 0.946 / 0.941 | — | 0.982 / 0.976 | −ln ratio 3.1 / 2.5 |
| I/I₀ (400) / (004) at 2 ps (S21) | 0.886 / 0.932 | — | 0.961 / 0.977 | −ln ratio 3.1 / 3.0 |
| I/I₀ (400) / (004) at 10 ps (S21) | 0.842 / 0.865 | — | 0.939 / 0.953 | −ln ratio 2.7 / 3.0 |
| Diffuse at S21's A (zigzag, (400) zone): rise | 6.5 ps (labelled A = corner; negative amplitude, no rise-then-decay shape) | **1.93 ps** (point (4.5, 0, 0), labelled X) | 1.7 ± 0.1 ps | 3.8 (**1.14**) |
| Diffuse at S21's A: relaxation | 16.9 ps | 19.2 ps | 30 ± 2 ps | 0.56 (0.64) |
| Diffuse at S21's X (armchair, (400) zone): single exponential | 1.15 ps (labelled X = zigzag) | **11.0 ps** (point (4, 0, 0.5), labelled Z) | 14.3 ± 0.1 ps | 0.08 (**0.77**) |
| Diffuse amplitudes ΔI/I₀ at A (peak) / X (plateau) | 1.12 / 0.72 | 0.85 / 1.35 | 0.017 / 0.039 | not comparable (see note) |
| Non-thermal lifetime (Z20 definition: normalized MSDs coincide) | 26 ps (2 % threshold; 20 ps at 5 %), derived from the preregistered traces | — | ≈ 60 ps | 0.44 |
| Lattice equilibration (preregistered definition: branch-temperature spread < 10 %) | 30 ps | — | ≈ 60 ps (Z20); 50–100 ps (S21) | 0.5 |
| Final lattice temperature | 305 K | — | 370 ± 50 K (Z20); 300 K (S21) | 0.82 |
| Absorbed energy density (input-derived) | 377 J/cm³ | — | 380 ± 70 J/cm³ | 0.99 |

**Coupling sensitivity** (preregistered, ×0.5 / ×2):
- Zigzag τ₁ = 1.06 / 0.61 ps, τ₂ = 21 / 19 ps.
- Armchair τ₁ = 0.18 / 0.12 ps, τ₂ = 16 / 12 ps.
- The measured armchair τ₁ (0.48 ps) lies outside this bracket. The measured zigzag τ₁ lies at the
  ×2 edge.

**Note on diffuse amplitudes.**
- The simulation gives the one-phonon intensity at one q point only.
- The measured intensity also contains Bragg tails, multi-phonon and inelastic background, and
  this background is unknown here. Measured relative changes are therefore diluted by an unknown
  factor.
- Only the time constants are like-for-like.

## Convention audit

**(a) Axes and point labels.**
- **Axes.** The pipeline's Cartesian axes are x = a = 3.31 Å (zigzag), y = b (out of plane) and
  z = c = 4.37 Å (armchair). So (h00) → U_xx = U₁₁ and (00l) → U_zz = U₃₃, exactly as in Z20.
  There is no axis swap.
- **Pipeline labels.** Special points carry the Setyawan–Curtarolo labels (from ASE) of the
  base-centred orthorhombic zone. In the (h0l) plane around (400):
  - X = (4 ± 0.55, 0, 0) along zigzag (grid point 4 ± 0.5);
  - Z = (4, 0, ±0.5) along armchair;
  - A = (4 ± 0.55, 0, ±0.5) at the corner.

  The preregistered file states these coordinates.
- **S21 labels.** S21 labels its points "according to the convention of Ref. 30", which is a
  different convention. The S21 text is explicit on both points:
  - X is along armchair: "no local minima arise in the conduction band along Γ–X and Z–Q
    (armchair direction)".
  - A is along zigzag: relaxation "entails a predominant emission of phonons with momenta q along
    the Γ–A directions", which the conclusions describe as "along the zigzag direction".
- **Consequence.** S21's A is the pipeline's X point type, and S21's X is the pipeline's Z.
- **Conflict with `reference_measurements.md`.** That file reads the Fig. 2b inset as the opposite
  assignment. The two statements in the S21 text are used here.
- **Effect on the comparison.** With points matched by coordinates, both diffuse time constants
  agree within 15–25 %: rise 1.93 against 1.7 ps, and single exponential 11.0 against 14.3 ps. The
  relaxation is 19 against 30 ps. With the literal labels they miss by 4× and 12×.
- **Verdict.** This is a label-convention mismatch, not a code bug. The pipeline output is
  unambiguous through `example_Q_invA`, so no code change was made.

**(b) MSD definition and equilibrium values.**
- **Definition.** The pipeline computes the per-atom tensor
  ⟨u_α u_β⟩ = ħ/(2 M N) Σ (2n+1)/ω Re(e_α e_β*), projects it on Ĝ, and uses I/I₀ = exp(−G·ΔU·G).
  This is Z20's Eqs. 2–3 with the same U₁₁/U₃₃. It is not ⟨|u|²⟩ and not a per-cell sum.
- **Cross-check against phonopy's thermal displacements** (MACE force constants):

  | | phonopy, same grid | phonopy, 24×24×18 | pipeline |
  |---|---|---|---|
  | U_xx / U_yy / U_zz at 100 K (10⁻³ Å²) | 3.19 / 4.55 / 4.74 | 3.24 / 4.71 / 4.92 | 3.36 / 4.33 / 4.85 |
  | Change, 100 → 304 K (10⁻³ Å²) | 3.83 / 6.50 / 6.89 (thermal) | — | 3.87 / 6.29 / 7.38 (at 100 ps) |

  The differences at 100 ps are the residual non-thermal part. There is no bookkeeping error.
- **Where the amplitude gap comes from.**
  - The harmonic MACE lattice gives ΔU(100 → 370 K) = 5.2 ×10⁻³ (zigzag) and 9.4 ×10⁻³ Å²
    (armchair) at full thermal equilibrium. That is 3.4–3.8× the measured asymptotes, even at the
    experiment's own ΔT = 270 K.
  - In the classical Debye limit, dU/dT = 3ħ²/(M k Θ²). The measured slopes correspond to Θ_eff ≈
    900 K (zigzag) and 710 K (armchair); the MACE phonons correspond to 490 K and 370 K.
- **Low-frequency modes.**
  - Modes below 3 THz carry 45 % (zigzag) and 69 % (armchair) of the MACE thermal ΔU.
  - The MACE interlayer shear modes at the zone boundary along b* are at 1.28 and 1.32 THz.
  - A PBE-level potential without dispersion softens these modes. That can inflate ΔU by tens of
    percent, but not by a factor 3.
- **Not resolved here.** The rest of the gap is not explained. Candidates:
  - the absolute Debye–Waller sensitivity of a 39 nm film (multiple scattering);
  - a smaller true lattice heating;
  - an under-estimate of BP's phonon stiffness by MACE, which the Raman frequencies (5 % soft) do
    not show.

  A short literature search found no tabulated equilibrium U₁₁/U₃₃ for BP to settle this.
- **Verdict.** Amplitudes are not a convention issue, so no code change and no post hoc
  amplitude column. The predicted Bragg drops are 2.4–3.1× too large in −ln(I/I₀).

## Hits and misses

**Hits:**
- Zigzag slow time constant: τ₂ = 20.5 against 20 ± 2 ps.
- Order of magnitude of all four MSD time constants.
- The sign and approximate size of the late-time anisotropy (armchair > zigzag, 1.9 against 1.63).
- The diffuse dynamics once points are matched by coordinates:
  - a fast (≈ 2 ps) rise followed by a slower relaxation along zigzag;
  - a slow (≈ 11–14 ps) single-exponential build-up along armchair.

  This means the momentum anisotropy of the hot-phonon population, with emission along zigzag,
  came out right without any material-specific input.
- Absorbed energy density.

**Misses:**
- **Debye–Waller amplitudes.** They are 2.5–3× too large in every direction and at every time.
  This is the largest miss.
- **Armchair fast time constant.** It is 4× too short (0.12 against 0.48 ps), and the coupling
  bracket does not cover it.
- **Armchair slow time constant.** It is 1.6× too short (14 against 22 ps).
- **Diffuse relaxation at S21's A.** It is too fast (19 against 30 ps).
- **Return to equilibrium.** It is about 2× too fast (26–30 against ≈ 60 ps). This is consistent
  with MACE's known over-estimate of the anharmonic decay rates (about 2× in graphite).
- **Fast fraction of the MSD rise.** The direction ordering is right, but the fraction is too
  large along zigzag (0.61 against 0.44).

## What this test says about the pipeline

- **Momentum-resolved energy flow.** The structure-only route predicted which part of the zone
  heats first and which heats last. It also predicted the 20 ps zigzag phonon–phonon time and the
  approximate time scales at the zone boundary.
- **Systematic error 1: phonon–phonon rates too fast.** The pipeline converges to thermal
  equilibrium about 2× too fast. The 14 ps armchair time and the 30 ps equilibration point the
  same way as the graphite comparison, where MACE lifetimes are about 2× short.
- **Systematic error 2: absolute MSD too large.** The pipeline predicts about 3× larger absolute
  Debye–Waller changes than measured. The audit rules out bookkeeping and conventions: the
  harmonic MACE lattice itself gives that much ΔU for the measured temperature rise.
- **Electron–phonon coupling.** The universal tight-binding coupling makes the carrier-to-optical
  energy transfer too fast along armchair (0.12 against 0.48 ps). Halving the coupling changes this
  by only 50 %, so the carrier phase space of the gap-shifted bands is the likely cause, as stated
  in the preregistered accuracy note.
- **Overall.** Time constants are within a factor 2 for 6 of the 8 coordinate-matched constants;
  only armchair τ₁ (4×) and the 1/e-type equilibration estimate (about 2.3×) fall outside.
  Amplitudes are outside the preregistered ±30 %.
- **Process lesson.** Point labels must be defined by coordinates. The preregistration did this,
  which is what made the label mismatch recoverable without changing any number.
