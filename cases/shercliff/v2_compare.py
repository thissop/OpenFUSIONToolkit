"""
V2 capstone metric — compare the MUG transient-conjugate waveform against the
COMSOL Tier6 transient reference, on the pinned common grid.

Protocol (COMSOL 4e2e255, confirmed OFT ba1707a): Hunt conjugate duct, Ha=100,
c_wall=1, Pm=1, rest IC, step-on constant -dp/dx (G) at t=0, observe tobs = 4*tau_A.
Metric: relL2 of the NORMALIZED waveforms v(t,y) and by(t,y) along the Hartmann
centerline cut x=0, PLUS the first-overshoot Alfven timing t_peak/tau_A.

Why normalize: the two codes run different unit systems and forcing amplitudes
(MUG: nu=eta=1, B0=0.1121 for Ha=100, fy=1e-2, tau_A=2/Ha=0.02; COMSOL: B0=1
nondim, G=1, tau_A=2). The problem is the SAME dimensionless transient (Ha, Pm, c
matched). Re_core ~ 1 => near-linear => the transient SHAPE and the Alfven timing
are amplitude-invariant. So we compare in t/tau_A with each field scaled to its
own asymptote — this isolates the transient physics (the discriminating claim),
not the trivially-cancelling amplitude.

Inputs
  --mug     npz from v2_waveform.py: keys t[nt], z[nz], vely[nt,nz], by[nt,nz]
  --comsol  COMSOL reference. Auto-detected:
              .npz with keys t, y (or z), v (or vely), by
              .csv  long format "t,y,v,by" (header row, comma-sep)
  --tau-mug     tau_A in MUG time units       (default 0.02 = 2/Ha at Ha=100)
  --tau-comsol  tau_A in COMSOL time units     (default 2.0)
  --tobs-taua   observation window in tau_A    (default 4.0)
  --norm    'asymptote' (divide each field by |field(t=tobs, y=0)|; default)
            or 'l2' (divide by sqrt(mean field^2 over the window))
  --out     npz with the gridded normalized fields + scalar metrics

Output: relL2_v, relL2_by, t_overshoot_mug/comsol (in tau_A), timing error.
Self-test: `python3 v2_compare.py --selftest` (no data needed).
"""
import argparse
import os

import numpy as np


# ----------------------------------------------------------------------------- loaders
def load_mug(path):
    d = np.load(path)
    t = np.asarray(d["t"], float).ravel()
    y = np.asarray(d["z"], float).ravel()          # Hartmann coord is mesh r2 (z)
    v = np.asarray(d["vely"], float)               # (nt, ny)
    by = np.asarray(d["by"], float)
    return t, y, v, by


def load_comsol(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npz":
        d = np.load(path)
        keys = set(d.files)
        t = np.asarray(d["t"], float).ravel()
        yk = "y" if "y" in keys else "z"
        y = np.asarray(d[yk], float).ravel()
        vk = "v" if "v" in keys else "vely"
        byk = "by" if "by" in keys else ("bx" if "bx" in keys else "by")
        v = np.asarray(d[vk], float)
        by = np.asarray(d[byk], float)
        # accept either (nt,ny) grid or long columns
        if v.ndim == 1:
            t, y, v, by = _long_to_grid(t, y, v, by)
        return t, y, v, by
    # csv long format: t,y,v,by (skip '#' comments and a non-numeric header row)
    rows = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            try:
                rows.append([float(v) for v in s.split(",")])
            except ValueError:
                continue                        # column-name header
    raw = np.array(rows)
    return _long_to_grid(raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3])


def _long_to_grid(t, y, v, by):
    """Reshape long columns (t,y,v,by) into (nt,ny) grids on the unique axes."""
    tu = np.unique(t)
    yu = np.unique(y)
    V = np.full((tu.size, yu.size), np.nan)
    B = np.full((tu.size, yu.size), np.nan)
    ti = {val: i for i, val in enumerate(tu)}
    yi = {val: i for i, val in enumerate(yu)}
    for k in range(t.size):
        V[ti[t[k]], yi[y[k]]] = v[k]
        B[ti[t[k]], yi[y[k]]] = by[k]
    return tu, yu, V, B


# ----------------------------------------------------------------------------- core
def resample(t_src, y_src, F_src, t_grid, y_grid):
    """Bilinear resample F_src(t_src,y_src) -> (t_grid,y_grid). Assumes sorted axes."""
    # interp along y for each source time, then along t
    tmp = np.empty((t_src.size, y_grid.size))
    for i in range(t_src.size):
        tmp[i] = np.interp(y_grid, y_src, F_src[i])
    out = np.empty((t_grid.size, y_grid.size))
    for j in range(y_grid.size):
        out[:, j] = np.interp(t_grid, t_src, tmp[:, j])
    return out


def normalize(F, y_grid, mode):
    """Scale a field to make amplitude comparable across codes."""
    if mode == "l2":
        s = np.sqrt(np.mean(F ** 2))
    else:  # asymptote: |field at last time, centerline y=0|
        j0 = int(np.argmin(np.abs(y_grid)))
        s = abs(F[-1, j0])
    return F / s if s > 0 else F, s


def first_overshoot(tau_grid, centerline, smooth=5):
    """Time (in tau_A) of the first genuine Alfven overshoot of the centerline
    signal rising from rest: the first local max of the SMOOTHED signal whose
    value exceeds the asymptote (final value) by a prominence margin. Smoothing
    makes it robust to measurement noise; the overshoot-above-asymptote condition
    is the physical definition (an overshoot rises above its own steady state).
    Falls back: if the response is overdamped (no peak above asymptote), returns
    the time it first reaches 99% of the asymptote."""
    c = centerline
    if smooth > 1:                              # centered moving average
        k = np.ones(smooth) / smooth
        c = np.convolve(c, k, mode="same")
    asymp = np.mean(c[-max(3, c.size // 20):])  # tail average = steady value
    scale = np.max(np.abs(c)) or 1.0
    prom = 0.01 * scale                         # ignore sub-1% wiggles
    for i in range(1, c.size - 1):
        if (c[i] > c[i - 1] and c[i] >= c[i + 1]
                and c[i] > asymp + prom):
            return _subgrid_peak(tau_grid, c, i)
    # overdamped: no overshoot -> first crossing of 99% asymptote
    hit = np.where(c >= 0.99 * asymp)[0]
    return tau_grid[hit[0]] if hit.size else tau_grid[int(np.argmax(c))]


def _subgrid_peak(t, c, i):
    """Parabolic sub-grid refinement of a peak at index i (removes the ~1-cell
    bias that smoothing/discretization puts on the overshoot time)."""
    if i <= 0 or i >= c.size - 1:
        return t[i]
    denom = c[i - 1] - 2 * c[i] + c[i + 1]
    if denom == 0:
        return t[i]
    delta = 0.5 * (c[i - 1] - c[i + 1]) / denom   # in grid units, |delta|<=0.5
    return t[i] + delta * (t[i + 1] - t[i - 1]) / 2.0


def compare(tm, ym, vm, bym, tc, yc, vc, byc,
            tau_mug, tau_comsol, tobs_taua, norm):
    # nondimensional time
    taum, tauc = tm / tau_mug, tc / tau_comsol
    # common grid: 201 pts in t/tau_A over [0, tobs_taua]; y on overlap of both
    tg = np.linspace(0.0, tobs_taua, 201)
    ylo, yhi = max(ym.min(), yc.min()), min(ym.max(), yc.max())
    yg = np.linspace(ylo, yhi, 129)

    VM = resample(taum, ym, vm, tg, yg)
    BM = resample(taum, ym, bym, tg, yg)
    VC = resample(tauc, yc, vc, tg, yg)
    BC = resample(tauc, yc, byc, tg, yg)

    VMn, _ = normalize(VM, yg, norm); VCn, _ = normalize(VC, yg, norm)
    BMn, _ = normalize(BM, yg, norm); BCn, _ = normalize(BC, yg, norm)

    def rel(a, b):
        return np.linalg.norm(a - b) / np.linalg.norm(b)

    relL2_v = rel(VMn, VCn)
    relL2_by = rel(BMn, BCn)

    j0 = int(np.argmin(np.abs(yg)))
    to_m = first_overshoot(tg, VMn[:, j0])
    to_c = first_overshoot(tg, VCn[:, j0])

    return dict(relL2_v=relL2_v, relL2_by=relL2_by,
                t_overshoot_mug=to_m, t_overshoot_comsol=to_c,
                timing_err=abs(to_m - to_c),
                tg=tg, yg=yg, VMn=VMn, VCn=VCn, BMn=BMn, BCn=BCn)


# ----------------------------------------------------------------------------- selftest
def selftest():
    """Analytic damped-overshoot surrogate; MUG = COMSOL in different units +
    small noise. relL2 should be tiny, timing error ~0."""
    ny = 65
    y = np.linspace(-1, 1, ny)
    prof = (1 - y ** 2)                       # parabolic-ish spatial shape
    tauA_c, tauA_m = 2.0, 0.02

    def surrogate(tphys, tauA, amp, noise=0.0, rng=None):
        s = tphys / tauA
        # rise to 1 with a damped Alfven overshoot at s~1
        env = 1 - np.exp(-0.7 * s)
        osc = 0.25 * np.exp(-0.4 * s) * np.sin(np.pi * s)
        sig = amp * (env + osc)
        F = np.outer(sig, prof)
        if noise and rng is not None:
            F = F + noise * amp * rng.standard_normal(F.shape)
        return F

    rng = np.random.default_rng(0)
    tc = np.linspace(0, 8, 400)               # COMSOL units, tobs=8=4tau
    tm = np.linspace(0, 0.08, 350)            # MUG units, tobs=0.08=4tau
    vc = surrogate(tc, tauA_c, amp=1.0)
    byc = surrogate(tc, tauA_c, amp=0.3)
    vm = surrogate(tm, tauA_m, amp=137.0, noise=0.01, rng=rng)   # diff amp + noise
    bym = surrogate(tm, tauA_m, amp=41.0, noise=0.01, rng=rng)

    r = compare(tm, y, vm, bym, tc, y, vc, byc, tauA_m, tauA_c, 4.0, "asymptote")
    print("SELFTEST (identical physics, different units+amp, 1% noise):")
    print(f"  relL2_v  = {r['relL2_v']:.4f}   (expect small, ~noise level)")
    print(f"  relL2_by = {r['relL2_by']:.4f}   (expect small)")
    print(f"  overshoot t/tau_A: MUG={r['t_overshoot_mug']:.3f} "
          f"COMSOL={r['t_overshoot_comsol']:.3f}  err={r['timing_err']:.3f}")
    ok = r["relL2_v"] < 0.05 and r["relL2_by"] < 0.05 and r["timing_err"] < 0.1
    # negative control: shift MUG timing by 0.5 tau_A -> relL2 must blow up
    tm2 = tm + 0.5 * tauA_m
    r2 = compare(tm2, y, vm, bym, tc, y, vc, byc, tauA_m, tauA_c, 4.0, "asymptote")
    print(f"  neg-control (0.5 tau_A shift): relL2_v={r2['relL2_v']:.3f} "
          f"(expect >> above)")
    ok = ok and r2["relL2_v"] > 3 * r["relL2_v"]
    print(f"  SELFTEST {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mug")
    ap.add_argument("--comsol")
    ap.add_argument("--tau-mug", type=float, default=0.02)
    ap.add_argument("--tau-comsol", type=float, default=2.0)
    ap.add_argument("--tobs-taua", type=float, default=4.0)
    ap.add_argument("--norm", choices=["asymptote", "l2"], default="asymptote")
    ap.add_argument("--out", default="v2_compare.npz")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        raise SystemExit(0 if selftest() else 1)
    if not (args.mug and args.comsol):
        raise SystemExit("need --mug and --comsol (or --selftest)")

    tm, ym, vm, bym = load_mug(args.mug)
    tc, yc, vc, byc = load_comsol(args.comsol)
    r = compare(tm, ym, vm, bym, tc, yc, vc, byc,
                args["tau_mug"] if isinstance(args, dict) else args.tau_mug,
                args.tau_comsol, args.tobs_taua, args.norm)

    print("V2 transient-conjugate comparison (Ha=100, c=1, Pm=1):")
    print(f"  relL2 v(t,y)  = {r['relL2_v']:.4e}")
    print(f"  relL2 by(t,y) = {r['relL2_by']:.4e}")
    print(f"  first-overshoot Alfven timing (t/tau_A): "
          f"MUG={r['t_overshoot_mug']:.3f}  COMSOL={r['t_overshoot_comsol']:.3f}  "
          f"err={r['timing_err']:.3f}")
    np.savez(args.out, **{k: v for k, v in r.items()
                          if k in ("tg", "yg", "VMn", "VCn", "BMn", "BCn",
                                   "relL2_v", "relL2_by",
                                   "t_overshoot_mug", "t_overshoot_comsol")})
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
