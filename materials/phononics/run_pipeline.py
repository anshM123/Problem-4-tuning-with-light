"""Run the LiNbO3 nonlinear-phononics pipeline in order, skipping cached steps.

  python run_pipeline.py          # reuse fc2/fc3 and npz intermediates when present
  python run_pipeline.py --all    # recompute everything (fc3: 1558 MACE supercells, 20-60 min)
"""
import os, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
env = dict(os.environ, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4",
           PYTHONUTF8="1")
steps = [("phonons.py", ["fc3.hdf5", "gamma_modes.npz"]),
         ("anharmonic.py", ["anharmonic.npz"]),
         ("model_check.py", ["model_check.log"]),
         ("phph.py", ["phph.npz"]),
         ("drive.py", ["drive_out.pkl"]),
         ("isotope_check.py", ["isotope_check.log"]),
         ("diagnostics.py", ["diag_part1.pkl", "diag_part2.pkl", "diag_part3.pkl", "diag_part4.pkl"]),
         ("plots.py", [])]
force = "--all" in sys.argv
for script, outputs in steps:
    if not force and outputs and all(os.path.exists(os.path.join(HERE, o)) for o in outputs):
        print(f"[cached] {script}")
        continue
    t0 = time.time()
    with open(os.path.join(HERE, script.replace(".py", ".out")), "w") as log:
        subprocess.run([sys.executable, script], cwd=HERE, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    print(f"[done] {script} ({time.time()-t0:.0f} s)")
