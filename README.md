# Tuning materials with light: fast many-body dynamics of driven materials

Authors: Ansh Mishra, Aryan Senthilkumar

This repository is our submission to Discovered Materials Grand Challenge 04, "Tuning
Materials with Light". It contains a hierarchy of time-linear methods for laser-driven
materials that:
- treat electron–electron, electron–phonon and phonon–phonon interactions together;
- are orders of magnitude cheaper than two-time Kadanoff–Baym propagation;
- are verified on Hubbard models against exact results;
- are applied to WSe₂, graphite and MgB₂.

**Where to find things**
- **Report:** [REPORT.md](REPORT.md) — method, results, scorecard, known limits, speed and limitations.
- **Experimental targets:** `refs/experimental_targets.md` (quoted with sources).
- **Figures:** `figures/`.

## Layout
| Path | Content |
|---|---|
| `lightmat/` | library: exact diagonalization (`ed.py`), positive-factorized pair-fluctuation NEGF, T-matrix and self-interaction-free strong-coupling closures (`fluct_hubbard.py`), particle–hole fluctuations (`fluct_ph.py`), phase-space trajectories (`ftwa.py`) |
| `benchmarks/hubbard/` | target 1: exact references and fast-method runs (CDW melting, ladder ground state, laser-driven ladder, optical conductivity) |
| `materials/graphite/` | target 3: MACE phonons, electron–phonon tables, full three-phonon kinetics, dynamics |
| `materials/graphite3d/` | the same pipeline for the 4-atom AB graphite cell (12 branches) |
| `materials/graphene12/` | matched monolayer control on the same 12×12 grid as `graphite3d` |
| `materials/wse2/` | target 2: exciton energetics and valley-resolved exciton/carrier kinetics |
| `materials/mgb2/` | target 4: two-band superconductivity, amplitude-mode nonlinear response, pumped bottleneck state |
| `tests/` | `python -m pytest tests -q` |
| `run_all.py` | reproduces every number and figure |

## Quick start
```
pip install -r requirements.txt
python -m pytest tests -q          # ~5 s
python run_all.py hubbard          # model benchmarks
python run_all.py materials        # graphite, WSe2, MgB2
```

## Applying the method to a new driven material

The material pipelines are written so that only the material inputs change.

### Semiconductor or semimetal with hot carriers and phonons (the graphite pipeline)

1. **Structure and phonons.**
   - Edit the cell in `phonons_mace.py` (or `materials/graphite3d/phonons.py`).
   - The script produces harmonic and cubic force constants with a MACE foundation model.
   - Pick the model by comparing its phonons with experiment or DFPT; `phonons_mace.py` prints frequencies at high-symmetry points.
   - Alternatively, supply `fc2.hdf5` and `fc3.hdf5` from DFT (phono3py format).
2. **Electron–phonon coupling.**
   - Replace the π-band tight-binding model and bond-modulation coupling in `eph_tb.py` with the material's bands and coupling vertex.
   - Calibrate the coupling to a DFPT or GW value of ⟨g²⟩.
   - `eph_rates.py` and `eph_rates_avg.py` then tabulate the electron-temperature-dependent generation coefficient A_qν(T_e) on the phonon grid.
3. **Three-phonon kernel.** `phph.py N [Nz]` builds the full nonlinear collision operator. It calibrates the absolute prefactor against phono3py automatically.
4. **Drive.** Set the absorbed energy (thin-film optics in `dynamics.py`) and the pulse duration.
5. **Run.** `dynamics.py` returns T_e(t) and n_qν(t) on the full Brillouin zone. `plots.py` resolves them by branch and wavevector.

### Excitonic semiconductor (the WSe₂ pipeline)

1. In `excitons.py`, set the screening length r₀, the dielectric environment and the band masses. The script computes binding energies for each valley pair.
2. In `dynamics.py`, set the valley energy offsets and the deformation potentials of each process (with its phonon energy, e.g. from DFPT).
3. Choose the pump: resonant (coherent exciton injection) or above-gap (`run_above_gap`, free carriers followed by exciton formation).

### Superconductor (the MgB₂ pipeline)

1. In `twoband.py`, set the band couplings λ_ab, the density-of-states ratio, plasma frequencies and scattering rates.
2. `nonlinear.py` computes amplitude-mode (Higgs) mediated harmonic and broadband responses.
3. `pump.py` and `bottleneck.py` compute the absorbed energy and the resulting nonequilibrium state with its conductivity.

### Correlated lattice models

`lightmat/fluct_hubbard.py`, `fluct_ph.py` and `ftwa.py` accept any hopping matrix and time-dependent Peierls phases or on-site potentials. See `benchmarks/hubbard/*.py` for examples.

## Licenses of external models
- The MACE model `medium-omat-0` is distributed under the Academic Software License.
- The MIT-licensed `medium-mpa-0` can be substituted in the phonon scripts. It gives about 20 % softer optical phonons for graphene.
