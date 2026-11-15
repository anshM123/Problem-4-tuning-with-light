# Submission summary: Grand Challenge 04, "Tuning Materials with Light"

**Ansh Mishra, Aryan Senthilkumar**

Full method, equations, results and limitations: [REPORT.md](REPORT.md). Experimental numbers we
compare against, quoted with sources: [refs/experimental_targets.md](refs/experimental_targets.md).

## Method: one solver, three limits

`lightmat/unified.py` propagates, in one time variable, the Nambu one-body density matrix, the
positive-factorized pair correlations, phonon means, covariances and three-phonon correlations, and
the mixed electron–phonon correlations. Electron–electron, electron–phonon and phonon–phonon
interactions and the light field enter the same equations of motion.

| Limit | Production level | Check |
|---|---|---|
| no phonons | A: positive-factorized pair fluctuations (T-matrix, self-interaction-free and doublon–holon closures) | identical to the Level A code, 6×10⁻¹³ |
| Markov + random phase | B: kinetic equations with first-principles couplings (DFT/GW e–ph, MACE three-phonon, DFPT exciton–ph) | collision operators identical, 10⁻¹⁵ |
| pairing mean field + linearization | C: two-band time-dependent Bogoliubov–de Gennes with collisions; amplitude-mode propagator | gaps to 10 digits, Higgs response within 3 % |

All three interactions at once, against exact diagonalization (Hubbard–Holstein dimer with anharmonic
phonons and a light pulse): exact to second order in the e–ph coupling, energy conserved to 3×10⁻⁸.
Graphite runs end to end through the same code (light → coherent Bloch electrons → e–ph → phonons →
three-phonon decay) and reproduces the Level B time scales at equal absorbed energy.

**Speed.** Against our own two-time Kadanoff–Baym solver on the same 11-site chain, time step and
machine: 300× faster at 400 steps; cost N_t^1.21 vs N_t^3.15; about 10⁴× at t = 50 ħ/J (9 s vs
≈ 29 h and 48 GB). Graphite, 1.2×10⁷ three-phonon processes over 100 ps: 11 min on one core.

**No parameter is adjusted to the data shown.** Inputs not computed here are literature values,
listed with sources in REPORT.md.

## Results

| Target | Quantity | This work | Reference |
|---|---|---|---|
| Hubbard 11-site CDW | stability, U/J = 0–20 | stable for all U | GKBA+TPP unstable for U/J ≳ 12 |
| | melting time, 16 values of U | mean error 0.29 ħ/J (U/J = 15: 3.9 vs 4.0) | exact diagonalization |
| Hubbard 6×2 ladder | spin correlations, U/J ≤ 1; laser-driven line | within 1.5–5 %; 3.58 J | ED; 3.48 J |
| 1D optical conductivity | Mott gap, U/J = 6, 12 (doublon–holon closure) | 2.91, 8.45 J | Bethe ansatz 2.89, 8.46 J |
| | 14-site threshold, U/J = 6, 12 | 3.06, 8.63 J | ED 3.06, 8.63 J |
| WSe₂ | "~400 fs" Σ maximum and "~300 fs" K→Σ delay, as the paper defines them | 0.40, 0.31 ps | 0.40, 0.31 ps |
| | K and Σ populations at all measured delays | rms 6.7 % of the Σ peak | Madéo et al. 2020 |
| | late K/Σ; bright-exciton linewidth at 50 K | 0.57; 6.9 meV | 0.47; 6.0 meV |
| Graphite | A1′(K) peak; K diffuse rise | 0.70 ps; 0.29 ps | 0.6–0.7 ps; 0.28 ps |
| | mid-zone rises K/2, M/2; TA pooling | 0.82, 1.31 ps; > 25 ps | 1.0–1.5 ps; > 25 ps |
| MgB₂ | 2Δ_π(0), 2Δ_σ(0); 2Δ_π(T) | 1.07, 3.40 THz; on the data 5–36 K | 1.10, 3.4 THz |
| | normalized FH / TH vs T | within 0.09 / 0.19 | Katsumi et al. 2025 |
| | pump-probe signal vs pump field (real-time dynamics, beam profile) | within 0.22 at 47–150 kV/cm | Fig. S3b |
| | conductivity 2 ps after 150 kV/cm (real-time dynamics, no assumed state) | σ₁(0.8) = 2.05, σ₂(0.5) = 1.63, σ₂(1.0) = 1.59 ×10⁴ | 2.55, 2.6, 1.5 ×10⁴ |
| LiNbO₃, light driving an IR phonon | TO frequencies; harmonic field scaling n = 1, 2, 3 | within 4 %; 0.87, 1.90, 2.89 | 1, 2, 3 |
| Blind: black phosphorus (`run_material.py`) | absorbed energy; zigzag ΔMSD slow time | 377 J/cm³; 20.5 ps | 380 J/cm³; 20 ps |
| | diffuse times at Seiler's A (rise) and X, matched by coordinates | 1.9 ps; 11.0 ps | 1.7 ps; 14.3 ps |

**Known limits recovered:** exact at U = 0; within a few percent of exact for U/J ≤ 2; detailed
balance in every kinetic equation; linear response as the first-order limit of the same equations.

**Open** (REPORT §6):
- Hubbard: U/J = 15 first zero crossing late; optical gap not controlled at U/J = 3.
- WSe₂: late decay 3.7 vs 2.2–2.5 ps.
- Graphite: K-point decay 2.46 vs 1.7 ps.
- MgB₂: real-time field scale 1.3–1.5× too low, so the residual σ₂ at 150 kV/cm is missed.
- LiNbO₃: blue shift of the driven mode and polarization reversal not reproduced; the controls show
  this is not a MACE error.
- Blind test: Debye–Waller amplitudes 2.5–3× too large; armchair fast time 4× too short.

## Reproduce

```
pip install -r requirements.txt
python -m pytest tests -q        # 44 tests, ~1 min
python run_all.py                # every number and figure (many hours, one CPU)
python run_all.py hubbard        # model benchmarks
python run_all.py materials      # graphite, WSe2, MgB2, LiNbO3
python run_all.py unified        # one solver: limits, ED test, graphite end to end
python run_all.py blind          # run_material.py: graphite self-check, black phosphorus
python run_material.py my_crystal.cif --pump-eV 1.55 --fluence-mJcm2 5 --fwhm-fs 50 --T 300
```
