"""Film-averaged absorbed energy from the coherent stage at five depths of the 70 nm film.

The standing-wave field profile |E(z)|^2 of the thin-film optics (n = 2.72 + 1.47i) gives, at the
centres z/d = 0.1, 0.3, 0.5, 0.7, 0.9 of five equal slices, field factors (relative to the
depth-rms field) 1.145, 0.976, 0.923, 0.945, 0.979.  Each slice is run with
    STAGE1_ONLY=1 python run.py 300 rta15 <factor> _depth<factor>
(the 300 x 300 grid; the 450 x 450 grid differs by 2 % at factor 1.0, run_depth1.0.log vs run.log).
This script reads the logs and prints the depth average.
"""
import os, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
fac = [1.145, 0.976, 0.923, 0.945, 0.979]
E, Elin = [], []
for f in fac:
    txt = open(os.path.join(HERE, f"run_depth{f}.log")).read()
    m = re.search(r"at 150 fs\): ([\d.]+) meV/cell; linear optics ([\d.]+)", txt)
    E.append(float(m.group(1))); Elin.append(float(m.group(2)))
ref = re.search(r"at 150 fs\): ([\d.]+)", open(os.path.join(HERE, "run_depth1.0.log")).read())
main = re.search(r"at 150 fs\): ([\d.]+)", open(os.path.join(HERE, "run.log")).read())
corr = float(main.group(1)) / float(ref.group(1))
for f, e, el in zip(fac, E, Elin):
    print(f"field factor {f:5.3f}: absorbed {e:6.1f} meV/cell (linear optics {el:6.1f})")
print(f"film average: {np.mean(E):.1f} meV/cell on 300^2 ({np.mean(E) * corr:.1f} with the 450^2 correction); "
      f"linear optics {np.mean(Elin):.1f}; absorptance {0.483 * np.mean(E) * corr / np.mean(Elin):.3f} (linear 0.483)")
