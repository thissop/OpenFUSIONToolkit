#!/usr/bin/env python3
"""P2 headline cross-validation: lay MUG-full over COMSOL-full in a shared nondim.

box 1 (COMSOL) reports a dimensional SI full-solution time series per observable; MUG runs its
own unit system. We reduce BOTH to dimensionless form on their natural scales (time/tau_eta,
velocity/v_A, current/(B0*a/mu0), Joule-rate/(B0^2*a^2/(mu0*tau_eta))) and compare:
  (1) normalized-shape overlay  U(t)/U_peak vs t/tau_eta   -- convention-INDEPENDENT truth check
  (2) absolute nondim peaks + Joule impulse                -- needs the shared ruler (box 1 to confirm)

box 1 side: p1_tauq_c_grid.csv, case tq1ms_c1_full (headline: tau_q=1ms, c=1, Ha=2645.7, tw/a=0.2).
MUG side:   xmhd_2d.moments (t, Fl_net, Fl_abs, EJf_rate, Umax) from a P2 run (fluid observables).

Usage:
  p2_compare.py --comsol p1_tauq_c_grid.csv [--mug xmhd_2d.moments]
Prints box-1 nondim reference immediately; adds MUG comparison if --mug given.
"""
import argparse, csv, math

MU0 = math.pi * 4e-7

def scales(a, B0, rho, lam):
    """Natural nondim scales from (half-width a, Hartmann field B0, density rho, mag. diffusivity lam).
    Ruler confirmed by box 1 (entry m): v_ref = lam/a (resistive velocity unit, = a/tau_eta), NOT v_A.
    All of I_ref/F_ref/E_ref below reproduce box 1's dimensionless refs to the digit."""
    tau_eta = a * a / lam
    return dict(
        a=a, B0=B0, rho=rho, lam=lam,
        tau_eta=tau_eta,
        v_ref=lam / a,                            # velocity unit (box 1: lam/a; MUG's a/tau_eta)
        I_ref=B0 * a / MU0,                       # current per unit depth
        EJ_ref=B0 * B0 * a * a / (MU0 * tau_eta), # Joule power (2D, per depth)
        imp_ref=B0 * B0 * a * a / MU0,            # Joule impulse = E_ref = EJ_ref * tau_eta
    )

# box 1's blind thin-wall (tw/a -> 0) targets, pre-registered before P2 landed (entry m).
# MUG is thin-wall, so score against THESE, not the resolved tw/a=0.2 values.
THIN_WALL_TARGETS = {
    "P2": dict(Umax=1.258071e-1, Iw=2.365607e0, Fl=5.802652e-1, EJf=1.175023e0, EJw=2.699945e-1),
    "P1": dict(Umax=1.043762e-2, Iw=2.876072e-2, Fl=7.204439e-3, EJf=1.971893e-2, EJw=7.802187e-3),
}

# box 1 dimensional SI parameters (from the driver / entry (h))
COMSOL = scales(a=0.1, B0=1.0, rho=9486.0, lam=1.136844)
# MUG unit system (verified: B0=9.0325e-4 -> Ha=2645.706, Pm=9.273e-8)
MUG = scales(a=1.0, B0=9.0324610e-4, rho=1.000228, lam=1.0)


def load_comsol(path, case="tq1ms_c1_full"):
    t, Iw, Fla, EJw, EJf, Um = [], [], [], [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            if not row["case_id"].startswith(case):
                continue
            t.append(float(row["t"]));  Iw.append(float(row["Iw_abs"]))
            Fla.append(float(row["Fl_abs"])); EJw.append(float(row["EJw_rate"]))
            EJf.append(float(row["EJf_rate"])); Um.append(float(row["Umax"]))
    return dict(t=t, Iw=Iw, Fl_abs=Fla, EJw=EJw, EJf=EJf, Umax=Um)


def load_mug(path):
    # columns: t Fl_net Fl_abs EJf_rate Umax [Iw EJw_rate]  (last two present on wall-observable runs)
    t, Fln, Fla, EJf, Um, Iw, EJw = [], [], [], [], [], [], []
    for ln in open(path):
        if ln.startswith("#") or not ln.split():
            continue
        p = ln.split()
        t.append(float(p[0])); Fln.append(float(p[1])); Fla.append(float(p[2]))
        EJf.append(float(p[3])); Um.append(float(p[4]))
        if len(p) >= 7:
            Iw.append(float(p[5])); EJw.append(float(p[6]))
    d = dict(t=t, Fl_net=Fln, Fl_abs=Fla, EJf=EJf, Umax=Um)
    if Iw:
        d["Iw"] = Iw; d["EJw"] = EJw
    return d


def trapz(y, x):
    return sum((y[k] + y[k-1]) * 0.5 * (x[k] - x[k-1]) for k in range(1, len(x)))


def reduce_side(name, d, s, has_wall):
    t = d["t"]
    tnd = [ti / s["tau_eta"] for ti in t]
    Und = [u / s["v_ref"] for u in d["Umax"]]
    EJf_nd = [e / s["EJ_ref"] for e in d["EJf"]]
    peakU = max(Und); tpeakU = tnd[Und.index(peakU)]
    impEJf = trapz(d["EJf"], t) / s["imp_ref"]
    out = dict(name=name, t_final=tnd[-1], peakU=peakU, tpeakU=tpeakU, impEJf=impEJf,
               peakFl=max(f / (s["B0"]*s["B0"]*s["a"]/MU0) for f in d["Fl_abs"]))
    if d.get("Iw"):  # box 1 always; MUG on wall-observable runs
        out["peakIw"] = max(i / s["I_ref"] for i in d["Iw"])
        out["impEJw"] = trapz(d["EJw"], t) / s["imp_ref"]
        out["impEJtot"] = out["impEJf"] + out["impEJw"]
        out["EJf_frac"] = out["impEJf"] / out["impEJtot"]
    return out


def show(o):
    print(f"  [{o['name']}]  t_final/tau_eta = {o['t_final']:.4f}")
    print(f"    peak Umax/v_ref        = {o['peakU']:.6e}  (at t/tau_eta={o['tpeakU']:.4f})")
    print(f"    peak Fl_abs (nondim)   = {o['peakFl']:.6e}")
    print(f"    Joule impulse EJf*     = {o['impEJf']:.6e}")
    if "peakIw" in o:
        print(f"    peak Iw/(B0 a/mu0)     = {o['peakIw']:.6e}")
        print(f"    Joule impulse EJw*     = {o['impEJw']:.6e}   (wall)")
        print(f"    Joule impulse EJtot*   = {o['impEJtot']:.6e}   [EJf frac = {o['EJf_frac']*100:.1f}%]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comsol", required=True)
    ap.add_argument("--mug")
    ap.add_argument("--case", default="P2", choices=["P1", "P2"])
    a = ap.parse_args()

    comsol_case = {"P2": "tq1ms_c1_full", "P1": "tq300ms_c1_full"}[a.case]
    print(f"=== COMSOL-full (box 1, {a.case} resolved tw/a=0.2), reduced to nondim ===")
    c = reduce_side("COMSOL-full", load_comsol(a.comsol, comsol_case), COMSOL, has_wall=True)
    show(c)
    print(f"  NOTE fluid vs wall Joule: EJf is {c['EJf_frac']*100:.1f}% of total "
          f"=> {'fluid dominates but wall non-negligible' if c['EJf_frac']<0.95 else 'fluid-dominated'}")

    if a.mug:
        tgt = THIN_WALL_TARGETS[a.case]
        print(f"\n=== MUG-full ({a.case}), reduced to nondim ===")
        m = reduce_side("MUG-full", load_mug(a.mug), MUG, has_wall=False)
        show(m)
        print(f"\n=== MUG-thin vs box-1 BLIND thin-wall targets ({a.case}); pre-reg: Umax<5%, intEJ few% ===")
        du = abs(m["peakU"] - tgt["Umax"]) / tgt["Umax"] * 100
        de = abs(m["impEJf"] - tgt["EJf"]) / tgt["EJf"] * 100
        df = abs(m["peakFl"] - tgt["Fl"]) / tgt["Fl"] * 100
        print(f"    peak Umax/v_ref : MUG {m['peakU']:.4e} vs target {tgt['Umax']:.4e}  -> {du:.1f}%  {'PASS' if du<5 else 'CHECK'}")
        print(f"    Joule imp EJf*  : MUG {m['impEJf']:.4e} vs target {tgt['EJf']:.4e}  -> {de:.1f}%")
        print(f"    peak Fl (nondim): MUG {m['peakFl']:.4e} vs target {tgt['Fl']:.4e}  -> {df:.1f}% (box 1: don't lean on Fl)")
        if m.get("peakIw") is not None:
            di = abs(m["peakIw"] - tgt["Iw"]) / tgt["Iw"] * 100
            dw = abs(m["impEJw"] - tgt["EJw"]) / tgt["EJw"] * 100
            tot_t = tgt["EJf"] + tgt["EJw"]
            dtot = abs(m["impEJtot"] - tot_t) / tot_t * 100
            print(f"    peak Iw/I_ref   : MUG {m['peakIw']:.4e} vs target {tgt['Iw']:.4e}  -> {di:.1f}%  (box 1: ~15% expected)")
            print(f"    Joule imp EJw*  : MUG {m['impEJw']:.4e} vs target {tgt['EJw']:.4e}  -> {dw:.1f}%")
            print(f"    Joule imp TOTAL : MUG {m['impEJtot']:.4e} vs target {tot_t:.4e}  -> {dtot:.1f}%  {'PASS' if dtot<5 else 'CHECK'}")
        else:
            print(f"    Iw + wall-Joule: not in this run (fluid-only moments).")


if __name__ == "__main__":
    main()
