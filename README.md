# Tuning materials with light — a fast many-body method for driven materials

Authors: Ansh Mishra, Aryan Senthilkumar

Work toward Discovered Materials Grand Challenge 04 ("Tuning Materials with Light"):
a method for laser-driven materials that treats electron–electron, electron–phonon and
phonon–phonon interactions together, is orders of magnitude cheaper than two-time
Kadanoff–Baym propagation, is verified on Hubbard models against exact results and is
applied to WSe2, graphite and MgB2.

Status: in progress. The full technical report with results, scorecard and limitations is
[REPORT.md](REPORT.md); experimental targets are collected in `refs/experimental_targets.md`.

## Layout
* `lightmat/` – library (exact diagonalization, time-linear factorized NEGF solvers, …)
* `benchmarks/hubbard/` – model-system verification (exact references + fast method)
* `tests/` – small consistency checks
* `refs/` – extracted experimental targets
