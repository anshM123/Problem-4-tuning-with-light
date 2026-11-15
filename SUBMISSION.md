# Submission summary: Grand Challenge 04, "Tuning Materials with Light"

**Ansh Mishra, Aryan Senthilkumar**

Full method, equations, results and limitations: [REPORT.md](REPORT.md). Experimental numbers we
compare against, quoted with sources: [refs/experimental_targets.md](refs/experimental_targets.md).

## Method

A hierarchy of methods whose cost grows linearly with simulated time, all derived from the
Kadanoff–Baym / G1–G2 hierarchy. Each level reduces to the next in a controlled limit.

| Level | What is propagated | Interactions | Used for |
|---|---|---|---|
| A | one-body density matrix; pair fluctuations kept positive by factorization | e–e: T-matrix, self-interaction-free pp ladder (exact for two particles), RPA/GW | Hubbard chains and ladders |
| B | carrier/exciton populations; phonon occupations on the full Brillouin zone | e–e, e–ph (DFT/GW), full three-phonon scattering (MACE), exciton–ph (DFPT) | graphite, WSe₂ |
| C | two-band gaps, amplitude-mode propagator, non-equilibrium quasiparticles | phonon-mediated pairing, A² light coupling, Rothwarf–Taylor bottleneck | MgB₂ |

**Speed.** Against our own two-time Kadanoff–Baym solver on the same 11-site chain, time step and
machine: 300× faster at 400 steps; cost N_t^1.21 vs N_t^3.15; about 10⁴× at t = 50 ħ/J (9 s vs
≈ 29 h and 48 GB). Graphite, 576 q × 6 branches, 1.2×10⁷ three-phonon processes, 100 ps: 11 min on
one CPU core.

**No parameter is adjusted to the data shown.** Inputs not computed from first principles: WSe₂
conduction-band offset, annihilation coefficient and above-gap formation time (literature);
MgB₂ pairing couplings (published), cutoff from T_c = 38 K, amplitude damping δ = 0.55 THz
(the experiment's own value).

## Results

| Target | Quantity | This work | Reference |
|---|---|---|---|
| Hubbard 11-site CDW | stability, U/J = 0–20 | stable for all U | GKBA+TPP unstable for U/J ≳ 12 |
| | melting time, 16 values of U | mean error 0.29 ħ/J (U/J = 15: 3.9 vs 4.0) | exact diagonalization |
| | long-time imbalance / double occupancy | mean error 0.06 / 12–17 % | exact diagonalization |
| Hubbard 6×2 ladder | ground-state spin correlations, U/J ≤ 1 | within 1.5–5 % | exact diagonalization |
| | laser-driven correlation line | 3.58 J | 3.48 J |
| 1D optical conductivity | main peak, U/J = 3 | 0.62 J | Bethe ansatz 0.63 J |
| WSe₂ | Σ-K − K-K energy | +6.8 meV | degenerate within 15 meV |
| | K and Σ populations at all measured delays | rms 7 % of the Σ peak | Madéo et al. 2020 |
| | Σ peak; peak Σ/K; late K/Σ | 0.30 ps; 1.18; 0.57 | 0.4 ps; 1.3; 0.47 |
| | bright-exciton phonon linewidth, 50 K | 6.9 meV | 6.0 meV |
| Graphite | A1′(K) peak; K diffuse rise | 0.70 ps; 0.29 ps | 0.6–0.7 ps; 0.28 ps |
| | mid-zone rises K/2, M/2 | 0.82, 1.31 ps | 1.0–1.5 ps |
| | TA pooling; diffuse maps at 0.5, 1.5, 5, 100 ps | > 25 ps; same sequence | > 25 ps; Fig. 2 |
| MgB₂ | 2Δ_π(0), 2Δ_σ(0); 2Δ_π(T) | 1.07, 3.40 THz; on the data 5–36 K | 1.10, 3.4 THz |
| | nonlinear peak at 7 K | 1.11 THz | 1.05 THz |
| | normalized FH / TH vs T | within 0.09 / 0.19 | Katsumi et al. 2025 |
| | conductivity 2 ps after 150 kV/cm | within 13–24 % | Katsumi et al. 2025 |

**Known limits recovered:** exact at U = 0; within a few percent of exact for U/J ≤ 2; detailed
balance in every kinetic equation; linear response as the first-order limit of the same equations.

**Open** (REPORT §6): Hubbard U/J = 15 first zero crossing late; WSe₂ intrinsic K→Σ transfer
≈ 20 fs vs the quoted "~400 fs" and late decay 3.7 vs 2.2 ps; graphite K decay 2.46 vs 1.7 ps;
MgB₂ field dependence not separately reproduced, raw FH maximum at 25 K vs 35–36 K.

## Reproduce

```
pip install -r requirements.txt
python run_all.py              # every number and figure (several hours, one CPU)
python run_all.py hubbard      # model benchmarks only
python run_all.py materials    # graphite, WSe2, MgB2 only
python -m pytest tests -q
```
