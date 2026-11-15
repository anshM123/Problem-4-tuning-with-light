"""Reproduce every number and figure of REPORT.md.

python run_all.py            # everything (several hours on one laptop CPU)
python run_all.py hubbard    # only the model benchmarks
python run_all.py materials  # only graphite, WSe2, MgB2
Precomputed force constants (materials/graphite/fc2.hdf5, fc3.hdf5) and electron-phonon
tables are in the repository; delete them to recompute from scratch.
"""
import os, sys, subprocess

ROOT = os.path.dirname(os.path.abspath(__file__))


def run(cwd, *args):
    print(">>", cwd, " ".join(args), flush=True)
    subprocess.run([sys.executable, *args], cwd=os.path.join(ROOT, cwd), check=True)


what = sys.argv[1] if len(sys.argv) > 1 else "all"
if what in ("all", "hubbard"):
    H = "benchmarks/hubbard"
    run(H, "ed_cdw.py", "0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "12", "14", "15", "16", "18", "20")
    for U in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "12", "14", "15", "16", "18", "20"):
        run(H, "fast_cdw.py", U, "factorized", "0.002")
    for U in ("5", "15"):
        run(H, "fast_cdw.py", U, "g1g2p", "0.002")
        run(H, "fast_cdw.py", U, "hf", "0.002")
    run(H, "ladder_gs.py", "ed"); run(H, "ladder_gs.py", "fast", "gwx", "20"); run(H, "ladder_pp.py", "40")
    run(H, "ladder_laser.py", "ed"); run(H, "ladder_laser.py", "fast", "gwx"); run(H, "ladder_laser.py", "pp")
    for U in ("3", "6", "12"):
        run(H, "chain_sigma.py", U, "64", "2", "gwx", "0")
    run(H, "plots.py")
    run("benchmarks/speed", "timing.py")
if what in ("all", "materials"):
    G = "materials/graphite"
    if not os.path.exists(os.path.join(ROOT, G, "fc3.hdf5")):
        run(G, "phonons_mace.py")
    if not os.path.exists(os.path.join(ROOT, G, "eph_rates_point.npz")):
        run(G, "eph_rates.py", "24", "450")
        os.replace(os.path.join(ROOT, G, "eph_rates.npz"), os.path.join(ROOT, G, "eph_rates_point.npz"))
    run(G, "eph_rates_avg.py", "4", "450")
    run(G, "phph.py", "24")
    run(G, "dynamics.py"); run(G, "plots.py"); run(G, "diffuse.py")
    run(G, "eph_kernel.py"); run(G, "dynamics_neq.py")
    G3 = "materials/graphite3d"
    run(G3, "phonons.py"); run(G3, "phph.py", "12", "2"); run(G3, "eph.py", "450"); run(G3, "dynamics.py")
    run("materials/wse2", "plots.py")
    M = "materials/mgb2"
    run(M, "nonlinear.py"); run(M, "pump.py"); run(M, "plots.py")
