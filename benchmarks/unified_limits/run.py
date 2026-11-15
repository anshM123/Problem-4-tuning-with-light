"""Run every demonstration of lightmat/unified.py and draw figures/unified_limits.png.

python run.py          (about 1 h on one laptop core; logs in data/*.log)
  full_vs_ed.py 20     (a) full level vs exact diagonalization, Hubbard-Holstein dimer   ~30 min
  levelA_limit.py      (b) Level A limit, 11-site CDW, t = 0..50                          ~9 min
  levelB_limit.py 1 2  (c) kinetic limit vs the graphite Level B code                     ~5 min
  levelB_limit.py 3    (c) full coherent solver vs its kinetic limit                      ~10 min
  levelC_limit.py      (d) two-band gap and amplitude mode vs materials/mgb2              ~25 min
  plot.py
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
for args, log in ((["full_vs_ed.py", "20"], "full_vs_ed.log"), (["levelA_limit.py", "5", "50", "0.002"], "levelA_limit.log"),
                  (["levelB_limit.py", "1", "2"], "levelB_part12.log"), (["levelB_limit.py", "3"], "levelB_part3.log"),
                  (["levelC_limit.py", "100", "5"], "levelC_limit.log"), (["plot.py"], None)):
    print(">>", " ".join(args), flush=True)
    with open(os.path.join(HERE, "data", log), "w") if log else open(os.devnull, "w") as f:
        subprocess.run([sys.executable, "-u", *args], cwd=HERE, check=True, stdout=f if log else None,
                       stderr=subprocess.STDOUT if log else None)
    if log:
        print(open(os.path.join(HERE, "data", log)).read())
