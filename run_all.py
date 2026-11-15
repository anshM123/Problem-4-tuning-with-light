"""Reproduce every number and figure of REPORT.md.

python run_all.py            # everything (several hours on one laptop CPU)
python run_all.py hubbard    # only the model benchmarks
python run_all.py materials  # only graphite, WSe2, MgB2
Precomputed force constants (materials/graphite/fc2.hdf5, fc3.hdf5) and electron-phonon
tables are in the repository; delete them to recompute from scratch.
"""
import os, sys, subprocess

ROOT = os.path.dirname(os.path.abspath(__file__))


def run(cwd, *args, env=None):
    print(">>", cwd, " ".join(args), env or "", flush=True)
    subprocess.run([sys.executable, *args], cwd=os.path.join(ROOT, cwd), check=True,
                   env={**os.environ, **(env or {})})


what = sys.argv[1] if len(sys.argv) > 1 else "all"
if what in ("all", "hubbard"):
    H = "benchmarks/hubbard"
    run(H, "ed_cdw.py", "0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "12", "14", "15", "16", "18", "20")
    for U in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "12", "14", "15", "16", "18", "20"):
        run(H, "fast_cdw.py", U, "factorized", "0.002")
    for U in ("0", "1", "2", "3", "4", "5", "6", "7", "8", "10", "12", "14", "15", "16", "18", "20"):
        run(H, "fast_cdw.py", U, "sic", "0.002")
    for U in ("5", "15"):
        run(H, "fast_cdw.py", U, "g1g2p", "0.002")
        run(H, "fast_cdw.py", U, "hf", "0.002")
    run(H, "ladder_gs.py", "ed"); run(H, "ladder_gs.py", "fast", "gwx", "20"); run(H, "ladder_pp.py", "40")
    run(H, "ladder_laser.py", "ed"); run(H, "ladder_laser.py", "fast", "gwx"); run(H, "ladder_laser.py", "pp")
    for U in ("3", "6", "12"):
        run(H, "chain_sigma.py", U, "64", "2", "gwx", "0")
    # optical conductivity: exact 14-site references, mean-field and collision kernels on the
    # same chain, and the doublon-holon pair closure
    run(H, "chain_sigma_ed.py", "pbc", "14", "12", "6", "3")
    run(H, "chain_sigma_ed.py", "obc", "14", "12", "6", "3")
    for U in ("3", "6", "12"):
        run(H, "chain_sigma.py", U, "14", "2", "gwx", "0", "_L14")
        run(H, "chain_sigma.py", U, "14", "exact", "gwx", "1", "_L14")
    for flags in (["sw"], [], ["o1"]):
        run(H, "chain_sigma_dh.py", *flags, "14", "3", "6", "12")
    run(H, "chain_sigma_dh.py", "sw", "ring", "128", "3", "6", "12")
    run(H, "chain_sigma_dh.py", "ring", "128", "3", "6", "12")
    run(H, "chain_sigma_compare.py")
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
    # matched monolayer control on the same 12x12 in-plane grid and the same soft-mode cut
    C = "materials/graphene12"
    run(C, "phph.py", "12"); run(C, "eph_rates.py", "12", "450")
    os.replace(os.path.join(ROOT, C, "eph_rates.npz"), os.path.join(ROOT, C, "eph_rates_point.npz"))
    run(C, "eph_rates_avg.py", "4", "450")
    run(C, "dynamics.py", env={"VALID_MIN_THZ": "1.0", "RTOL": "1e-3", "OUT": "ctrl.npz"})
    run("materials/wse2", "plots.py")
    M = "materials/mgb2"
    run(M, "nonlinear.py"); run(M, "pump.py"); run(M, "plots.py")
