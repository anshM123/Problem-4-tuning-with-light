"""(d) Level C limit: two-band superconductor (MgB2 couplings) in the unified solver.

Model: two bands of discrete levels xi = +-(j + 1/2) wc / n_b, j < n_b, with n_b = N_b n (so
the level density is proportional to the band density of states N_b of materials/mgb2), and a
phonon-mediated (Eliashberg-derived, Liu-Mazin-Kortus) reduced pairing interaction
    H_V = sum_{k in a, k' in b} V^lev_ab c^+_ku c^+_kd c_k'd c_k'u,   V^lev_ab = -V_ab wc / n,
with V_ab = twoband.V and the cutoff wc of twoband.fit_wc() (Tc = 38 K).  Units: meV, hbar/meV.

1. Equilibrium: the unified solver's Nambu mean field (UnifiedSolver.equilibrium) against the
   Level C gap equation, evaluated with the Level C functions (twoband._I, V, NR) on the same
   levels, and against twoband.gaps (continuum, n = 4000) to show the discretization error.
2. Amplitude (Higgs) mode: a weak pair-field pulse h_b(t) on band b is applied to the
   equilibrium state; the solver's linear response dDelta_a(t) (time-dependent Bogoliubov-de
   Gennes = linearization of eq. 1) is Fourier transformed with window exp(-eta t) and compared
   with the Level C propagator: dDelta = [D(w) V^-1 - 1] h, D = [V^-1 - diag(N P(w + i eta))]^-1,
   P from twoband.P_amp on the same levels.

Usage: python levelC_limit.py [n] [T_K]     -> data/levelC_limit.npz
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "materials", "mgb2"))
import twoband as TB
from lightmat.unified import Model, UnifiedSolver

ETA = 0.15          # meV, response window (= eta in P_amp)
WMAX = 15.0         # meV, upper end of the comparison window (pulse spectrum exp(-w^2 s^2/2) > 3e-4)


def levels(n, wc):
    nb = [int(round(TB.NR[a] * n)) for a in range(2)]
    assert all(abs(nb[a] - TB.NR[a] * n) < 1e-9 for a in range(2)), "n * N_b must be integer"
    xi, band = [], []
    for a in range(2):
        x = (np.arange(nb[a]) + 0.5) * wc / nb[a]
        xi += list(-x[::-1]) + list(x); band += [a] * (2 * nb[a])
    return np.array(xi), np.array(band), nb


def build(n, wc, pulse=None):
    xi, band, nb = levels(n, wc)
    Vlev = -TB.V[band][:, band] * wc / n
    pf = None
    if pulse is not None:
        b, h0, t0, s = pulse
        sel = (band == b).astype(float)
        pf = lambda t: np.diag(sel * h0 * np.exp(-(t - t0) ** 2 / (2 * s ** 2)))
    m = Model(np.diag(xi), U=0.0, V=Vlev, pair_field=pf)
    return m, xi, band, nb


def solver_gaps(n, wc, T):
    m, xi, band, nb = build(n, wc)
    L = len(xi)
    r = np.diag((xi < 0).astype(float))
    s = UnifiedSolver(m, r, r, kappa=0.1 * np.eye(L), ee=None)
    D = s.equilibrium(TB.KB * T, tol=1e-13)
    return np.abs(np.array([np.real(D[band == a]).mean() for a in range(2)])), np.array(
        [np.real(D[band == a]).std() for a in range(2)]), s, band


def levelC_gaps_same_levels(n, wc, T):
    """Level C gap equation (twoband._I, V, NR) on the solver's levels (n_b midpoints)."""
    nb = [int(round(TB.NR[a] * n)) for a in range(2)]
    D = np.array([2.3, 7.0])
    for _ in range(20000):
        I = np.array([TB._I(D[a], T, wc, n=nb[a]) for a in range(2)])
        Dn = TB.V @ (TB.NR * I * D)
        if np.max(abs(Dn - D)) < 1e-13:
            D = Dn; break
        D = 0.5 * D + 0.5 * Dn
    return np.abs(D)


def levelC_response(w, D, T, wc, n, eta=ETA):
    nb = [int(round(TB.NR[a] * n)) for a in range(2)]
    P = np.array([TB.P_amp(w, D[a], T, wc, eta=eta, n=nb[a]) for a in range(2)])
    K = np.linalg.inv(TB.V) - np.diag(TB.NR * P)
    Dm = np.linalg.inv(K)
    return Dm @ np.linalg.inv(TB.V) - np.eye(2)


def higgs(n, wc, T, b=0, h0=1e-3, t0=1.0, s=0.25, tmax=50.0, dt=0.04):
    Dg, _, s_eq, band = solver_gaps(n, wc, T)
    m, xi, band, nb = build(n, wc, pulse=(b, h0, t0, s))
    L = len(xi)
    sol = UnifiedSolver(m, np.zeros((L, L)), np.zeros((L, L)), ee=None)
    sol.G = s_eq.G.copy()
    nt = int(round(tmax / dt))
    ts = np.arange(nt + 1) * dt
    dD = np.zeros((nt + 1, 2))
    base = np.array([np.real(sol.gap()[band == a]).mean() for a in range(2)])
    for k in range(nt + 1):
        g = sol.gap()
        dD[k] = [np.real(g[band == a]).mean() - base[a] for a in range(2)]
        if k < nt:
            sol.step(dt)
    return ts, dD, Dg


def ft(ts, y, w, eta):
    dt = ts[1] - ts[0]
    wt = np.exp(-eta * ts) * dt
    wt[0] *= 0.5; wt[-1] *= 0.5
    return np.array([np.sum(y * wt * np.exp(1j * ww * ts)) for ww in w])


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    T = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
    wc = TB.fit_wc()
    print(f"wc = {wc:.3f} meV; n = {n} levels per side in band pi, {int(round(TB.NR[1]*n))} in sigma; T = {T} K")
    t0 = time.time()
    Ds, Dstd, _, _ = solver_gaps(n, wc, T)
    Dc = levelC_gaps_same_levels(n, wc, T)
    Dcont = TB.gaps(T, wc)
    print(f"gaps (meV): unified solver {Ds[0]:.10f} {Ds[1]:.10f} (spread over levels {Dstd.max():.1e}); "
          f"Level C same levels {Dc[0]:.10f} {Dc[1]:.10f}; Level C continuum {Dcont[0]:.6f} {Dcont[1]:.6f} "
          f"[{time.time()-t0:.1f}s]")
    print(f"  2Delta (THz): solver {2*Ds[0]/TB.THZ:.4f} {2*Ds[1]/TB.THZ:.4f}; continuum {2*Dcont[0]/TB.THZ:.4f} {2*Dcont[1]/TB.THZ:.4f}")
    w = np.linspace(0.5, 20.0, 976)
    out = dict(n=n, T=T, wc=wc, gap_solver=Ds, gap_levelC=Dc, gap_cont=Dcont, w=w, eta=ETA)
    for b in (0, 1):
        t0 = time.time()
        h0, tp, sp_ = 1e-3, 1.0, 0.25
        ts, dD, _ = higgs(n, wc, T, b=b, h0=h0, t0=tp, s=sp_)
        hw = ft(ts, h0 * np.exp(-(ts - tp) ** 2 / (2 * sp_ ** 2)), w, ETA)
        Rs = np.array([ft(ts, dD[:, a], w, ETA) / hw for a in range(2)])        # (2, nw)
        Rc = np.array([levelC_response(ww, Dc, T, wc, n)[:, b] for ww in w]).T
        out[f"t{b}"] = ts; out[f"dD{b}"] = dD; out[f"Rsolver{b}"] = Rs; out[f"RlevelC{b}"] = Rc
        win = w <= WMAX            # the probe pulse carries spectral weight > 3e-4 of its peak here
        rel = np.abs(Rs - Rc)[:, win].max() / np.abs(Rc)[:, win].max()
        for a in range(2):
            ps = w[np.argmax(np.abs(Rs[a]))]; pc = w[np.argmax(np.abs(Rc[a]))]
            print(f"pulse on band {b}, response of band {a}: |R| peak solver {ps:.3f} meV, Level C {pc:.3f} meV "
                  f"(2Delta_pi = {2*Ds[0]:.3f}, 2Delta_sigma = {2*Ds[1]:.3f})")
        print(f"  max |R_solver - R_levelC| / max |R_levelC| over 0.5..{WMAX:g} meV = {rel:.2e}   [{time.time()-t0:.0f}s]")
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    np.savez(os.path.join(HERE, "data", "levelC_limit.npz"), **out)


def summarize(fn=os.path.join(HERE, "data", "levelC_limit.npz")):
    C = np.load(fn); w = C["w"]; win = w <= WMAX
    for b in (0, 1):
        Rs, Rc = C[f"Rsolver{b}"], C[f"RlevelC{b}"]
        print(f"pulse on band {b}: max |R_solver - R_levelC| / max|R_levelC| over 0.5..{WMAX:g} meV = "
              f"{np.abs(Rs - Rc)[:, win].max() / np.abs(Rc)[:, win].max():.2e}; peaks solver "
              f"{[float(w[win][np.argmax(np.abs(Rs[a][win]))]) for a in range(2)]}, Level C "
              f"{[float(w[win][np.argmax(np.abs(Rc[a][win]))]) for a in range(2)]} meV")
