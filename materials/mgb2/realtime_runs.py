"""Drivers for realtime.py (MgB2 at 7 K, film of SM Sec. I).

  python realtime_runs.py validate     equilibrium sigma(w) at 7 K from a weak probe vs Mattis-Bardeen
  python realtime_runs.py pump E       broadband pump of peak E kV/cm (E = 'nb21' for the 21 kV/cm
                                       band-pass pulse): dynamics, probe sigma(w) at tau = 2 ps,
                                       kinetic continuation to 100 ps.  Writes rt_pump_<E>.npz
  python realtime_runs.py harm E       narrowband 0.63 THz drive, nonlinear field E_AB - E_A - E_B
                                       at tau = 0 (A = B = E kV/cm).  Writes rt_harm_<E>.npz
  python realtime_plots.py             figure and tables from the stored runs (rt_out/rt_summary.txt)
Environment: RT_N0, RT_N1 (levels per band), RT_SEED, RT_SINGLE (1 = complex64), RT_FWHM
(narrowband envelope, ps), RT_TKIN (continuation end, ps), RT_SUFFIX (output name suffix).
"""
import os
import sys
import time
import numpy as np
import realtime as rt
import twoband as tb

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "rt_out")
os.makedirs(OUT, exist_ok=True)
NB = dict(N=(int(os.environ.get("RT_N0", 300)), int(os.environ.get("RT_N1", 250))), W=(15.0, 18.0),
          seed=int(os.environ.get("RT_SEED", 1)))
SUF = os.environ.get("RT_SUFFIX", "")     # output suffix for convergence/sensitivity runs
T0 = 7.0
T_PUMP = 1.6             # ps, pump centre
TAU = 2.0                # pump-probe delay (Fig. S3c,d)
DT = 0.01


def new_sim(collisions=True):
    return rt.Sim(T0=T0, N=NB["N"], W=NB["W"], seed=NB["seed"], collisions=collisions)


def pump_field(tag):
    if tag == "nb21":
        return rt.narrowband(21e5, T_PUMP + 1.5)
    return rt.broadband(float(tag) * 1e5, T_PUMP)


def validate():
    t0 = time.time()
    s = new_sim()
    f, sig, _, (tp, Ep, dE) = rt.probe_sigma(s, lambda t: 0.0, T_PUMP + TAU, Eprobe=10e5)
    np.savez(os.path.join(OUT, "rt_validate.npz"), f=f, sig=sig, tp=tp, dE=dE, D0=s.D0)
    print(f"grid gaps at 7 K: 2Dpi = {2*s.D0[0]/rt.THZ:.3f} THz, 2Dsig = {2*s.D0[1]/rt.THZ:.3f} THz  ({time.time()-t0:.0f} s)")
    D = tb.gaps(T0, rt.WC)
    for fr in (0.3, 0.5, 0.8, 1.0, 1.5, 2.0):
        k = np.argmin(abs(f - fr)); ref = tb.sigma_total(fr, D, T0)
        print(f"  {fr:.1f} THz: real-time probe s1 = {sig[k].real/100:.3e}, s2 = {sig[k].imag/100:.3e};"
              f"  Mattis-Bardeen s1 = {ref.real/100:.3e}, s2 = {ref.imag/100:.3e}")


def pump(tag, t_kin=float(os.environ.get("RT_TKIN", 100.0))):
    t0 = time.time()
    s = new_sim()
    Ef = pump_field(tag)
    tpr = (T_PUMP + 1.5 if tag == "nb21" else T_PUMP) + TAU
    rec = {"t": [], "D": [], "nqp": [], "Uabs": [], "Uph": [], "A": []}

    def cb(sim):
        if int(round(sim.t / DT)) % 5:
            return
        rec["t"].append(sim.t); rec["D"].append(sim.gaps(sim.R)); rec["nqp"].append(sim.qp_density())
        rec["Uabs"].append(sim.Uabs); rec["Uph"].append(sim.ph.energy()); rec["A"].append(sim.A)

    rt.run(s, Ef, tpr - 2.0, DT, callback=cb)
    s2 = s.copy()
    rt.run(s2, Ef, tpr, DT)
    snap = s2.snapshot()
    f, sig, s_after, (tp, Ep, dE) = rt.probe_sigma(s, Ef, tpr, Eprobe=10e5, callback=cb)
    print(f"[{tag}] probe done at {time.time()-t0:.0f} s", flush=True)
    fm = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])
    sigmb = np.array([rt.sigma_mb_state(x, snap, T0) for x in fm])
    kin = rt.pseudospin_continue(s_after, t_kin, dt=0.01, rec=0.5)
    np.savez(os.path.join(OUT, f"rt_pump_{tag}{SUF}.npz"),
             t=np.array(rec["t"]), D=np.array(rec["D"]), nqp=np.array(rec["nqp"]), Uabs=np.array(rec["Uabs"]),
             Uph=np.array(rec["Uph"]), A=np.array(rec["A"]), f=f, sig=sig, fmb=fm, sigmb=sigmb,
             tp=tp, Ep=Ep, dE=dE, tpr=tpr,
             snap_E0=snap["E0"], snap_f0=snap["f0"], snap_E1=snap["E1"], snap_f1=snap["f1"],
             snap_D=snap["D"], snap_n=s2.ph.n, phW=s.ph.W, phnT=s.ph.nT, phD=s.ph.D,
             kin_t=np.array([k[0] for k in kin]), kin_D=np.array([k[1] for k in kin]),
             kin_f0=np.array([k[2][0] for k in kin]), kin_f1=np.array([k[2][1] for k in kin]),
             kin_n=np.array([k[3] for k in kin]), eps0=s.bands[0].eps, eps1=s.bands[1].eps, Uabs_end=s_after.Uabs)
    print(f"[{tag}] absorbed {s_after.Uabs:.3e} J/m^3 ({s_after.Uabs*17.46e-6:.3f} J/mol); A after pulse {rec['A'][-1]:.2e} V s/m; "
          f"gaps at tau: 2Dpi = {2*abs(snap['D'][0])/rt.THZ:.3f}, 2Dsig = {2*abs(snap['D'][1])/rt.THZ:.3f} THz; "
          f"total {time.time()-t0:.0f} s", flush=True)
    for fr in (0.5, 0.8, 1.0):
        k = np.argmin(abs(f - fr)); km = np.argmin(abs(fm - fr))
        print(f"   {fr} THz: probe s1 = {sig[k].real/100:.3e} s2 = {sig[k].imag/100:.3e};  "
              f"MB(f, Delta at tau) s1 = {sigmb[km].real/100:.3e} s2 = {sigmb[km].imag/100:.3e}", flush=True)


def harm(E, fwhm=float(os.environ.get("RT_FWHM", 3.0))):
    """E_NL = E_AB - E_A - E_B with identical collinear A, B at tau = 0: run fields E and 2E."""
    t0 = time.time()
    tc = 6.0
    t1 = 14.0
    out = {}
    for k, amp in (("A", E), ("AB", 2 * E)):
        s = new_sim()
        ts, Ef, D, A = rt.run(s, rt.narrowband(amp * 1e5, tc, fwhm=fwhm), t1, DT)
        out[k] = Ef; out["D" + k] = D
    ENL = out["AB"] - 2 * out["A"]
    np.savez(os.path.join(OUT, f"rt_harm_{E:g}{SUF}.npz"), t=ts, EA=out["A"], EAB=out["AB"], ENL=ENL,
             DA=out["DA"], DAB=out["DAB"])
    dt = ts[1] - ts[0]
    F = np.fft.rfft(ENL * np.hanning(len(ENL)), 8 * len(ENL)); fr = np.fft.rfftfreq(8 * len(ENL), dt)
    P = np.abs(F) ** 2
    IFH = P[(fr > 0.3) & (fr < 1.0)].sum(); ITH = P[(fr > 1.6) & (fr < 2.2)].sum()
    print(f"[harm {E:g} kV/cm] I_FH = {IFH:.4e}  I_TH = {ITH:.4e}  min Dpi(AB) = {np.min(abs(out['DAB'][:,0])):.3f} meV  ({time.time()-t0:.0f} s)", flush=True)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "validate":
        validate()
    elif what == "pump":
        pump(sys.argv[2])
    elif what == "harm":
        harm(float(sys.argv[2]))
