# Fast many-body dynamics of light-driven materials

**Ansh Mishra, Aryan Senthilkumar**

Technical report for Discovered Materials Grand Challenge 04, "Tuning Materials with Light".

- Every number below is produced by the code in this repository.
- `python run_all.py` regenerates all results and figures (see section 7).

---

## 1. Summary

We present a hierarchy of time-linear methods for driven electron–phonon matter.

- Every level descends from the Kadanoff–Baym / G1–G2 hierarchy.
- Each level reduces to the next in a controlled limit.
- All levels share three ingredients:
  - one-body density matrices;
  - two-particle (pair) fluctuations kept positive by construction;
  - phonon populations or fields.

| Level | Propagated objects | Interactions | Cost / step | Used for |
|---|---|---|---|---|
| A. Positive-factorized pair fluctuations | ρ and low-rank factors x, y of the pair matrices X = A+𝒢, Y = B+𝒢 (pp channel); Hermitian fluctuation samples V (ph channel) | e–e in the T-matrix (pp) and RPA/GW (ph) channels | O(r L³) | Hubbard chains and ladders |
| B. Kinetic (Markov) limit with first-principles couplings | electronic energy, phonon occupations n_qν on a full-BZ grid, exciton valley and energy populations, free carriers | e–ph (DFT/GW), full three-phonon scattering (MACE cubic force constants), exciton–ph (DFPT deformation potentials) | O(N_proc) | graphite, WSe₂ |
| C. Pairing mean field plus Gaussian amplitude fluctuations | two-band gaps Δ_π, Δ_σ, the amplitude-mode propagator, nonequilibrium quasiparticle distribution | phonon-mediated pairing (Eliashberg-derived couplings), A² light coupling, Rothwarf–Taylor bottleneck | O(N_ξ) | MgB₂ |

**Scorecard.**
- Every entry is a prediction: no parameter was adjusted to the data shown.
- The exceptions are explicitly marked "fit".
- "ED" is our own exact diagonalization, which reproduces the reference's exact results (e.g. the 3.48 J spin-correlation frequency of the laser-driven ladder).

| Target | Quantity | This work | Reference | |
|---|---|---|---|---|
| **1 Hubbard 11-site CDW** | stability, U/J = 0…20 | stable for every U | GKBA+TPP unstable for U/J ≳ 12 | ✔ |
| | long-time ⟨I⟩ over 16 values of U | mean abs. error 0.07 | ED | ✔ |
| | long-time Σd_corr | 12 % mean error | ED | ✔ |
| | first zero crossing of I | correct for U/J ≤ 4; delayed for U/J ≥ 5 (also in plain G1–G2) | ED | ✘ |
| **1 Hubbard 6×2 ladder** | ground-state C_loc, C_leg, C_rung | pp ladder: C_loc, C_leg within 1.5 % at U/J = 0.5, 5 % at U/J = 1; −18 % (C_loc) at U/J = 4 | ED | ✔ for U/J ≤ 2 |
| | laser-driven ΔC spectrum | dominant line 3.58 J (ph fluctuations) | 3.48 J (ED) | ✔ (3 %) |
| **1 1D optical conductivity** | gap onset, U/J = 3, 6, 12 | 1.56, 5.06, 11.4 J (TDHF/RPA level) | Bethe ansatz 0.63, 2.89, 8.46 J | ✘ (overestimated, as in the reference δGW) |
| **2 WSe₂** | E(Σ-K) − E(K-K) | +6.8 meV | degenerate within 15 meV (≈ +3 meV) | ✔ |
| | Σ-K population peak | 0.30 ps | 0.4 ps | ✔ |
| | peak Σ / peak K | 1.18 | 1.3 | ✔ |
| | late K/Σ ratio (resonant) | 0.57 | 0.47 | ✔ (≈) |
| | decay time | 3.7 ps | 2.2 ps | ≈ |
| | K→Σ rise delay | 50 fs | ~300 fs | ✘ |
| | above-gap K/Σ | 0.28 → 0.57 (rising) | 0.15 → 0.32 (rising) | qualitative |
| **3 Graphite** | A1′(K) peak time | 0.70 ps | 0.6–0.7 ps | ✔ |
| | A1′(K) rise | 0.18 ps (63 %) | 0.28 ps | ≈ |
| | A1′(K) decay | 3.2 ps | 1.7 ps | ✘ (×1.9) |
| | energy flow optical → acoustic → flexural, mode- and q-resolved | yes | yes | ✔ |
| | energy conservation | exact | — | ✔ |
| **4 MgB₂** | 2Δ_π(0), 2Δ_σ(0) | 1.07, 3.40 THz | 1.10, 3.4 THz | ✔ |
| | 2Δ_π(T) | on the measured points from 5 to 36 K | Fig. 3c | ✔ |
| | σ₂(ω) at 5 K | 9.0, 5.3, 2.4, 0.9 ×10⁴ Ω⁻¹cm⁻¹ at 0.3, 0.5, 1, 2 THz | ≈ 8, 6, 2, 1 ×10⁴ | ✔ |
| | broadband nonlinear peak vs T | 1.11, 0.98, 0.88, 0.78, 0.62 THz at 7, 12, 17, 22, 27 K | 1.05, 1.03, 0.95, 0.90, 0.80 | ✔ (≤ 22 K) |
| | normalized TH vs T | 1, 0.84, 0.64, 0.52, 0.41 at 8, 12, 17, 20, 22 K | 1, 0.9, 0.5, 0.3, 0.2 | ✔ monotonic, shallower |
| | normalized FH vs T | flat to ≈ 23 K, then falls | monotonic from 7 K | ✘ |
| | conductivity 2 ps after 150 kV/cm (bottleneck state) | σ₁(0.8) = 2.23, σ₂(0.5) = 1.98, σ₂(1.0) = 1.15 ×10⁴ | 2.55, 2.6, 1.5 ×10⁴ | ✔ (13–24 %) |
| **Speed** | graphite: 576 q × 6 branches, 1.2×10⁷ three-phonon processes, 100 ps | 11 min, 1 CPU core | two-time NEGF: not feasible | ✔ |

---

## 2. Method

### 2.1 Level A: positive-factorized pair fluctuations

Notation: ħ = 1, ρ^σ_ij = ⟨c†_jσ c_iσ⟩, ρ̄ = 1 − ρ, local interaction U Σ_i n_i↑ n_i↓.

**One-body equation (exact):**

  i ∂_t ρ^↑_ij = [h^HF_↑, ρ^↑]_ij + U (𝒢_{ii,ji} − 𝒢_{ij,jj}),

where 𝒢 is the connected part of ⟨c†_c↑ c†_d↓ c_b↓ c_a↑⟩.

**T-matrix (pp-ladder) equation.** The G1–G2 equation for 𝒢 is

  i ∂_t 𝒢 = [h₂, 𝒢] + U (A−B) P 𝒢 − U 𝒢 P (A−B) + U (A P B − B P A),

with
- A = ρ̄^↑ ⊗ ρ̄^↓ and B = ρ^↑ ⊗ ρ^↓;
- P the projector on on-site pairs;
- h₂ = h^HF_↑ ⊗ 1 + 1 ⊗ h^HF_↓.

**Key identity.** For X = A + 𝒢 and Y = B + 𝒢, this equation is equivalent to the *linear* pair

  i ∂_t X = H X − X H†,  i ∂_t Y = H Y − Y H†,  H = h₂ + U (A − B) P.

The source U(APB − BPA) is generated automatically.

**Positivity by construction.**
- Write X = x x† and Y = y y†, and propagate the rectangular factors with i ∂_t x = H x.
- Then X ≥ 0 and Y ≥ 0 at all times.
- These are the N-representability conditions whose loss makes plain G1–G2 unstable.

**Exact low rank.** For an uncorrelated initial state the factorization is exact:
- rank X = number of hole pairs, rank Y = number of particle pairs;
- there is no sampling noise;
- no N_b⁴ object is ever formed.

The correlation used in the collision integral is the Keldysh average 𝒢 = (X − A + Y − B)/2.

This is a deterministic, exact-rank, positivity-preserving formulation of the
quantum-fluctuation idea of Schroedter, Bonitz and Joost (arXiv:2606.10773). Code:
`lightmat/fluct_hubbard.py`. The plain G1–G2 mode, and a variant with positivity projection,
are in the same file for comparison.

**Particle–hole channel** (`lightmat/fluct_ph.py`):
- Hermitian samples V of the one-body operator fluctuations evolve with the linearized
  time-dependent Hartree–Fock kernel (RPA, GW or GW+X variants).
- Initial samples are an exact deterministic decomposition of ½⟨{δE_ab, δE_cd}⟩, or a
  Gaussian sample for large systems.
- The implementation is fully sparse: O(S M²) per step.

**Verification:**
- Without the collision term, the factorized and plain G1–G2 dynamics agree to 8 digits.
- At U/J = 1–2 both reproduce exact I(t) and d(t) (7 sites) within a few percent.
- The sparse ph kernel agrees with the dense one to 10⁻¹⁵.

### 2.2 Level B: kinetic limit with first-principles couplings

In dense photoexcited systems the dephasing is fast and the Markov limit of the level-A
equations applies. The pair-correlation equation then integrates to golden-rule collision
kernels.

**Graphite: electron–phonon.**
- π bands from nearest-neighbour tight binding (ħv_F = 6.4 eV Å).
- Coupling by bond-length modulation, β·d̂·(u_B − u_A), evaluated with the MACE phonon
  eigenvectors of all six branches on a 24×24 q-grid.
- β is fixed by the DFT value ⟨g²_K⟩ = 0.0994 eV² (Piscanec et al. 2004). The model then
  gives ⟨g²_Γ⟩ = 0.0444 eV² (DFT: 0.0405).
- The many-body (GW) enhancement of Lazzeri et al. (2008) is applied: ×193/92 for intervalley
  processes, ×1.24 for intravalley processes.
- The emission kernel is A_qν(T_e). It is averaged over each Brillouin-zone cell near Γ, K and
  K′ (4×4 sub-points; 450×450 electron grid), because the emission region there is smaller
  than one cell.
- Check at 300 K: the coupling-limited optical-phonon lifetimes are 0.6 ps (E2g, i.e.
  8.5 cm⁻¹ against 11 cm⁻¹ DFT) and 0.3 ps (A1′).

**Graphite: phonon–phonon.**
- Full nonlinear three-phonon collision integral: 1.2×10⁷ energy-allowed decay and
  coalescence processes, with detailed balance.
- Cubic force constants from MACE (`medium-omat-0`). Among six MACE foundation models this one
  reproduces the graphene phonon frequencies (E2g 193 meV against 196 measured; ZO 104 meV
  against ≈ 110).
- The absolute prefactor is calibrated against phono3py's imaginary self-energy. The spread
  over all modes is 2×10⁻⁴, so the operator is correct.
- The finite-width energy delta leaves a small residual, which is removed by a
  heat-capacity-weighted projection. Energy is then conserved exactly.
- Resulting lifetimes: A1′(K) 2.3, 1.5, 1.1 ps at 300, 600, 900 K; E2g(Γ) 2.1, 1.6, 1.2 ps.
  These are about 2× shorter than the DFT values of Bonini et al. (2007); MACE slightly
  overestimates the anharmonicity.

**WSe₂.**
- Exciton binding energies from the Rytova–Keldysh interaction (r₀ = 45.1 Å, hBN substrate)
  with DFT k·p masses.
- Exciton–phonon scattering uses the DFPT deformation potentials of Jin et al. (PRB 90, 045422,
  2014), each process with its own phonon energy:
  - K→Q: acoustic 8.2×10⁸ eV/cm; K→Q′: 5.7×10⁸; K→K′: 1.3×10⁸;
  - Q→Q′, optical intravalley, and acoustic intravalley (|D_e − D_h|) as tabulated there.
- Valley classes: K-K and K-K′ (4 species), Σ-K (12 species).
- Detailed balance at 90 K; exciton–exciton annihilation γ = 0.08 cm²/s (literature range for
  hBN-supported WSe₂).
- Above-gap excitation adds a free-carrier stage:
  - carriers are distributed over the K and Σ valleys at the carrier temperature;
  - they cool by optical-phonon emission, using the same deformation potentials;
  - they bind into excitons with a formation time of 0.5 ps (literature value).

### 2.3 Level C: superconductors (MgB₂)

**Gaps.** Two-band BCS gap equations with the Liu–Mazin–Kortus couplings (λ_ππ = 0.28,
λ_σσ = 0.96, λ_πσ = 0.19, N_σ/N_π = 0.73). The cutoff is fixed by T_c = 38 K, giving 7.3 meV
(1.8 THz), which is the value printed in the reference.

**Amplitude (Higgs) mode.** The coupled two-band amplitude-mode propagator is

  D(ω) = [V⁻¹ − diag(N_a P_a(ω))]⁻¹,  P_a(ω) = ∫dξ tanh(E/2T)/E · 4ξ²/(4E² − ω²).

**Nonlinear response.**
- Light enters through the A² (Raman-type) vertex v_a ∝ ω²_p,a Δ_a.
- The third-order current is j₃ ∝ Σ v_a D_ab v_b A³.
- Thin-film radiation and the experiment's normalization by the transmitted field are included.

**Linear conductivity.**
- Two-band Mattis–Bardeen, with the superfluid weight corrected for finite purity
  (Abrikosov–Gor'kov sum).
- The tabulated plasma frequencies give a normal-state conductivity about 10× above the
  film's. The total σ_n is therefore set to the measured 38 K value; band ratios and
  scattering rates are kept.

**Strongly pumped state.**
- The absorbed energy is computed from the 150 kV/cm pulse with thin-film optics.
- Pair-breaking phonons and quasiparticles are in Rothwarf–Taylor equilibrium, so the
  excitation is stored as excess quasiparticles at the bath temperature.
- This is the Owen–Scalapino state, with quasiparticle chemical potential μ* fixed by energy
  conservation.
- Gaps and σ(ω) follow self-consistently, with no free parameter.

---

## 3. Results

### 3.1 Target 1: Hubbard models (`benchmarks/hubbard`)

**Exact references.** All exact references are our own exact diagonalizations:
- 11-site chain, 213,444 states;
- 6×2 ladder, 853,776 states, including the time-dependent laser pulse.

#### CDW melting (`figures/hubbard_cdw.png`, `figures/hubbard_cdw_sweep.png`)

![CDW](figures/hubbard_cdw.png)

**TDHF baseline.** The imbalance is frozen: I ≈ 1.2 at U/J = 5 and I ≈ 1.95 at U/J = 15.

**Factorized pp ladder.**
- Stable at every U up to 20.
- Reaches the correct long-time state over the full sweep:

| U/J | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 | 15 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ⟨Σd_corr⟩ exact | 0.57 | 0.88 | 1.22 | 1.53 | 1.78 | 1.97 | 2.22 | 2.36 | 2.45 | 2.53 | 2.61 |
| ⟨Σd_corr⟩ this work | 0.76 | 1.06 | 1.49 | 1.82 | 2.09 | 2.23 | 2.41 | 2.50 | 2.51 | 2.55 | 2.29 |

**Weakness.** The CDW melts too slowly at strong coupling: at U/J = 15 the first zero
crossing comes at 27 ħ/J against 5.4 ħ/J exactly. Plain G1–G2 with the same self-energy
shows the same lag. It even stays at I ≈ 0.33 at long times, and our positivity-projected
variant does too. The lag is therefore a property of the T-matrix closure, not of the
factorization.

![sweep](figures/hubbard_cdw_sweep.png)

#### Ladder ground state (`figures/hubbard_ladder.png`)

Adiabatic switching of U.

| U/J | C_loc exact / pp | C_leg exact / pp | C_rung exact / pp |
|---|---|---|---|
| 0.5 | 0.1339 / 0.1328 | 0.0612 / 0.0604 | 0.0005 / 0.0004 |
| 1 | 0.1429 / 0.1391 | 0.0673 / 0.0643 | 0.0014 / 0.0007 |
| 2 | 0.1615 / 0.1487 | 0.0813 / 0.0703 | 0.0042 / 0.0013 |
| 4 | 0.1962 / 0.1610 | 0.1091 / 0.0782 | 0.0149 / 0.0021 |

- The pp ladder is quantitative up to U/J ≈ 1.5.
- It underestimates moment formation beyond that, and adiabatic switching fails for
  U/J ≥ 5.
- The ph (RPA-kernel) method brackets the exact result from above: C_loc = 0.204 at U/J = 2.

![ladder](figures/hubbard_ladder.png)

#### Laser-driven ladder (`figures/hubbard_ladder_laser.png`)

U/J = 1, with a Peierls pulse on the legs: A_L = 0.5, ħω_L = 2J, σ = 2ħ/J.

| Method | Lines in the ΔC_i(t) spectrum (J) |
|---|---|
| Exact | 3.48, 3.75, 3.97 (agrees with the reference ED) |
| ph-fluctuation method | dominant line 3.58 (+3 %); secondary line at 3.32 |
| pp ladder | 4.28 (wrong frequency); gets the mean ΔC_leg shift right (−0.010 vs −0.008) |

![laser](figures/hubbard_ladder_laser.png)

#### Optical conductivity of the half-filled chain (`figures/hubbard_chain_sigma.png`)

- Linear response from the self-consistent HF antiferromagnet (64 sites).
- The Mott gap opens at every U but is overestimated: onsets of 1.56, 5.06 and 11.4 J against
  the exact Bethe-ansatz gaps of 0.63, 2.89 and 8.46 J.
- The reference reports the same overestimate for its δGW.
- Adding the fluctuation dynamics (ph collisions, or trajectories) on top of the
  antiferromagnet was numerically unstable in our tests and is not used.

![sigma](figures/hubbard_chain_sigma.png)

### 3.2 Target 2: WSe₂ (`materials/wse2`, `figures/wse2_valleys.png`)

![WSe2](figures/wse2_valleys.png)

**Energies.** Rytova–Keldysh binding energies on hBN are 184 meV (K-K) and 212 meV (Σ-K).
With the DFT conduction-band offset this places Σ-K 6.8 meV above K-K. The measurement finds
them degenerate within 15 meV.

**Resonant 1.72 eV, 90 K.** With first-principles deformation potentials and no fitted
parameter:
- the K-valley population rises with the pulse;
- the dark Σ-K population peaks at 0.30 ps (measured 0.4 ps) and dominates after ≈ 0.5 ps;
- the peak Σ/K ratio is 1.18 (measured 1.3);
- K/Σ falls from 1.2 to 0.57 (measured 1.15 → 0.47), the Boltzmann value of the computed
  energy landscape;
- the late decay is 3.7 ps (measured 2.2 ps).

The ~300 fs delay between the K and Σ rises is not reproduced: the computed intervalley
scattering is faster.

**Above-gap 2.48 eV.**
- Free carriers created near Λ/Σ give an initially Σ-rich ratio, K/Σ = 0.28, rising as
  excitons form.
- The measurement shows the same trend, 0.15 → 0.32.
- The model relaxes too far, to 0.57.

### 3.3 Target 3: graphite (`materials/graphite`, `figures/graphite_*.png`)

![A1K](figures/graphite_A1K.png)

**Setup.** A 35 fs, 800 nm, 12 mJ/cm² pulse. Thin-film optics gives 48 % absorptance, or
91 meV per cell. The electron temperature peaks at 6700 K.

**A1′(K) population.**
- rises in 0.18 ps (63 %; measured 0.28 ps);
- peaks at 0.70 ps (measured 0.6–0.7 ps);
- decays in 3.2 ps (measured 1.7 ps).

**Why the decay is too slow.** We decomposed the rates at 2 ps:
- the anharmonic decay of A1′ alone corresponds to τ ≈ 1.8 ps, matching the measured decay;
- but the electrons, still at ≈ 2000 K, keep refeeding it.

The model's coupled electron + strongly-coupled-phonon reservoir therefore cools too slowly.
The most likely cause is the single-layer (graphene) treatment of the lattice, which omits
graphite's interlayer modes and their additional decay channels.

**Mode-resolved excess energy** (meV per cell):

| t | TO | LO | LA | TA | ZO | ZA |
|---|---|---|---|---|---|---|
| 0.5 ps | 31.7 | 35.2 | 6.4 | 5.6 | 4.9 | 0.7 |
| 1.5 ps | 22.5 | 24.6 | 11.8 | 11.3 | 15.0 | 3.4 |
| 5 ps | 10.9 | 11.2 | 15.2 | 15.6 | 25.5 | 12.2 |
| 100 ps | 6.7 | 6.7 | 16.7 | 18.0 | 18.6 | 24.2 |

![BZ](figures/graphite_bz_maps.png)

**Comparison with the experiment's energy flow.**
- Optical phonons near K and Γ are populated first.
- Energy then flows into the in-plane acoustic and flexural branches, as in the
  diffuse-scattering data.
- At 100 ps the model lattice is thermal at ≈ 630 K. The experiment still resolves
  non-thermal TA pooling at M and ⅓M, so acoustic thermalization in the model is too fast.
  This is again consistent with the missing interlayer physics.

### 3.4 Target 4: MgB₂ (`materials/mgb2`, `figures/mgb2.png`)

![MgB2](figures/mgb2.png)

**Gaps and linear response.** The gaps and their temperature dependence follow from the
published couplings. The equilibrium σ₂(ω) agrees with the measurement within ≈ 15 %.

**Nonlinear response.**
- The peak of the broadband nonlinear spectrum follows 2Δ_π(T) and agrees with the
  measurement within 0.1 THz up to 22 K.
- Above 25 K the measured peaks stay higher; the authors attribute this to an overdamped mode
  tracking the transmitted spectrum.
- The normalized third harmonic grows monotonically on cooling, as measured, but more weakly:
  half value at ≈ 21 K versus ≈ 17 K.
- The normalized first harmonic is flat up to 23 K, while the measurement already decreases.
  The un-normalized first harmonic does rise and then fall with temperature as measured, but
  peaks at 28 K instead of 35–36 K.
- A dirty-limit superfluid-weight vertex gives the wrong first-harmonic trend, rising toward
  T_c. This discriminates between the two light–Higgs coupling mechanisms.

**Strongly pumped state.**
- A quasi-thermal state cannot reproduce the measured combination of a normal-like σ₁ with a
  large σ₂ at 2 ps.
- The Rothwarf–Taylor/Owen–Scalapino bottleneck state, fixed by the absorbed energy alone
  (0.61 J/mol → μ* = 3.6 meV), gives:
  - 2Δ_π = 0.23 THz, so the π gap is almost closed;
  - 2Δ_σ = 2.0 THz, so the σ gap is still open;
  - σ₁(0.8 THz) = 2.23×10⁴ against 2.55×10⁴ measured;
  - σ₂(0.5) = 1.98×10⁴ against 2.6×10⁴;
  - σ₂(1.0) = 1.15×10⁴ against 1.5×10⁴.
- The pump therefore leaves MgB₂ in a band-selective non-thermal state.

---

## 4. Known limits

- **Noninteracting.** The Hubbard U = 0 dynamics are exact: the I(t) curves coincide with ED
  to plotting accuracy, and d_corr ≡ 0.
- **Weak coupling.** At U/J ≤ 2 the dynamics agree with ED within a few percent. For U/J ≤ 1,
  the local and leg ground-state spin correlations are within ≤ 5 %.
- **Equilibrium.**
  - Every kinetic equation obeys detailed balance, so the thermal state is stationary.
  - Graphite optical-phonon lifetimes (coupling and anharmonic) are reproduced.
  - MgB₂ gaps and equilibrium σ(ω) are reproduced.
- **Weak drive.** Linear response is the first-order limit of the same equations: σ(ω) for
  MgB₂ and the Hubbard chain; the ladder correlation spectrum.

## 5. Speed

All timings are for one laptop CPU core, Python/NumPy.

| Calculation | This work | Reference / exact |
|---|---|---|
| 11-site CDW, t = 50 ħ/J | 26–66 s | ED 120–200 s, feasible only for small clusters. Two-time KBE with N_t ≈ 2.5×10⁴ steps scales as N_t³. |
| 6×2 ladder ground state, one U | 19 s (pp ladder) | ED 10–20 s; ED scales exponentially |
| Laser-driven ladder, t = 100 ħ/J | 54 s (pp), 816 s (ph, 288 samples) | ED ≈ 3000 s |
| Graphite, full BZ, 100 ps | 11 min | Two-time NEGF with phonons over a full BZ: not feasible |
| WSe₂, 4 ps | ≈ 2 s | — |
| MgB₂ gaps, σ(ω), nonlinear and pumped-state response | 1–20 min | — |

## 6. Limitations (what is not solved)

1. **Hubbard.** The T-matrix closure melts the CDW too slowly for U/J ≳ 5. The ladder ground
   state is quantitative only for U/J ≲ 2. The optical gap is that of RPA-level theory.
2. **Graphite.** The A1′ decay is too slow (3.2 vs 1.7 ps) and acoustic thermalization too
   fast. Both trace to the graphene-layer treatment of the lattice. Electrons are assumed
   thermalized; the first ≈ 100 fs of non-thermal carriers is not resolved.
3. **WSe₂.** The conduction-band offset and annihilation coefficient are literature values.
   The K→Σ rise delay is too short.
4. **MgB₂.** The temperature dependence of the first harmonic is not reproduced.
5. **Method coverage.** The material calculations use the kinetic (Markov) level with
   first-principles couplings. The coherent level-A dynamics has not been run in an
   ab initio band basis.
6. **Licensing.** The MACE model `medium-omat-0` is under the Academic Software License; the
   MIT-licensed `medium-mpa-0` can be substituted at the cost of ≈ 20 % softer optical phonons.

## 7. Reproduction

```
pip install -r requirements.txt
python run_all.py            # all targets (several hours, one CPU)
python run_all.py hubbard    # model benchmarks only
python run_all.py materials  # graphite, WSe2, MgB2 only
```

Individual stages are documented at the top of each script. Precomputed force constants and
electron–phonon tables are included; delete them to recompute from scratch.
