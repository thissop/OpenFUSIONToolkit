"""
Recompute the Hartmann three-way headline DIRECTLY FROM THE STORED CSVs.

Packaging-session integrity pass (2026-06-30). Does NOT re-extract OpenFOAM VTK
(no meshio dependency) and does NOT copy numbers from any prior markdown. It reads:

  - mhdFoam + analytic:  cases/hartmann/results/profile_Ha{1,5,10,20,50}.csv
        columns: y, u_analytic_norm, u_mhdfoam_norm, u_mug_norm   (already centerline-normalized)
  - resolved MUG:        cases/hartmann/mug/Ha_*_profile.csv and
                         cases/hartmann/mug_ginsburg/Ha_*/Ha_*_profile.csv
                         (raw y, ux; the real MUG runs)

and recomputes, on each source's own grid:
  - L2  = sqrt( int (u_code - u_ana)^2 dy / int u_ana^2 dy )   on centerline-normalized profiles
  - Linf= max|u_code - u_ana|
  - u_mean/u_c and its rel. error vs the exact analytic scalar u_mean_over_uc(Ha)

For MUG it also reports Ha_fit (curve_fit of the Hartmann shape) vs Ha_pred (from the .profile
header) and the L2 against analytic at BOTH integer Ha and fitted Ha, to reconcile the two
numbers that appeared in the prior docs.

Rebuilds: results/error_table.md, results/overlay_Ha*.png, results/flowrate_vs_Ha.png,
and results/report/*.png (smplotlib).  Honest blanks where a source produced no run.
"""
import os, glob, re
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)
import analytic_hartmann as ah

RESULTS = os.path.join(HERE, "results")
MUG = os.path.join(HERE, "mug")
A = 1.0

OF_HA = [1, 5, 10, 20, 50]          # Ha values with a stored mhdFoam results CSV


def detect_mug_ha():
    """Return integer Ha values with stored resolved-MUG profile CSVs."""
    paths = glob.glob(os.path.join(MUG, "Ha_*_profile.csv"))
    paths.extend(glob.glob(os.path.join(HERE, "mug_ginsburg", "Ha_*", "Ha_*_profile.csv")))
    ha_values = set()
    for path in paths:
        match = re.search(r"Ha_(\d+)_profile\.csv$", path)
        if match:
            ha_values.add(int(match.group(1)))
    return sorted(ha_values)


MUG_HA = detect_mug_ha()

# Prior published numbers (for the integrity diff only; NOT used in any output)
PRIOR_OF = {  # from cases/hartmann/results/error_table.md (committed)
    1:  (1.045e-03, 1.080e-03, 0.6783, 1.298e-03),
    5:  (1.951e-04, 2.943e-04, 0.8118, 1.079e-03),
    10: (6.902e-04, 1.194e-03, 0.9012, 1.283e-03),
    20: (5.793e-04, 1.464e-03, 0.9511, 1.156e-03),
    50: (5.971e-04, 2.485e-03, 0.9810, 1.067e-03),
}


def trapz(y, x):
    f = getattr(np, "trapezoid", None) or np.trapz
    return f(y, x)


def norm_uc(y, u):
    uc = np.interp(0.0, y, u)
    return u / uc, uc


def umean_uc(y, u):
    uc = np.interp(0.0, y, u)
    return (trapz(u, y) / (y[-1] - y[0])) / uc


def metrics(y, u_code_norm, Ha):
    """L2, Linf, u_mean/u_c, relerr — all on centerline-normalized profiles."""
    ana = ah.u_profile(y, A, Ha)
    ana_n, _ = norm_uc(y, ana)
    code_n, _ = norm_uc(y, u_code_norm)   # idempotent if already normalized
    d = code_n - ana_n
    l2 = np.sqrt(trapz(d**2, y) / trapz(ana_n**2, y))
    linf = float(np.max(np.abs(d)))
    muc = umean_uc(y, u_code_norm)
    muc_ana = ah.u_mean_over_uc(Ha)
    return l2, linf, muc, abs(muc - muc_ana) / muc_ana


def load_of_csv(ha):
    path = os.path.join(RESULTS, f"profile_Ha{ha}.csv")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.strip().split(",")
            if len(p) < 3 or p[2] == "":
                continue
            rows.append((float(p[0]), float(p[1]), float(p[2])))
    if not rows:
        return None
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2]   # y, u_ana_norm, u_of_norm


def load_mug_csv(ha):
    path = mug_csv_path(ha)
    if not os.path.exists(path):
        return None
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    o = np.argsort(d[:, 0])
    return d[o, 0], d[o, 1]


def mug_csv_path(ha):
    candidates = [os.path.join(MUG, f"Ha_{ha}_profile.csv")]
    candidates.extend(
        sorted(glob.glob(os.path.join(HERE, "mug_ginsburg", f"Ha_{ha}_*", f"Ha_{ha}_profile.csv")))
    )
    for path in candidates:
        if os.path.exists(path):
            return path
    return candidates[0]


def mug_ha_pred(ha):
    """Read Ha_pred from the matching .profile header."""
    csv_path = mug_csv_path(ha)
    csv_dir = os.path.dirname(csv_path)
    cand = []
    if os.path.basename(csv_dir) != "mug":
        cand.append(os.path.join(csv_dir, "hartmann_mug.profile"))
    cand.append(os.path.join(MUG, f"Ha_{ha}", "hartmann_mug.profile"))
    if ha == 5:
        cand.append(os.path.join(MUG, "hartmann_mug.profile"))
    for path in cand:
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    if "Ha_pred" in line:
                        return float(line.split("=")[1])
    return None


def hartmann_shape(y, Ha, C):
    return C * (1.0 - np.cosh(Ha * y / A) / np.cosh(Ha))


def main():
    print("=" * 78)
    print("INTEGRITY RECOMPUTE — Hartmann three-way headline (from CSVs, packaging pass)")
    print("=" * 78)

    of_rows = {}
    print("\n--- mhdFoam (recomputed from results/profile_Ha*.csv) ---")
    print(f"{'Ha':>3} {'L2':>11} {'Linf':>11} {'umean/uc':>9} {'relerr':>10}   prior(L2,relerr) -> delta")
    for ha in OF_HA:
        d = load_of_csv(ha)
        if d is None:
            print(f"{ha:>3}  (no CSV)")
            continue
        y, _ana, uof = d
        l2, linf, muc, relerr = metrics(y, uof, ha)
        of_rows[ha] = (l2, linf, muc, relerr)
        pl2, pli, pmuc, pre = PRIOR_OF[ha]
        print(f"{ha:>3} {l2:11.3e} {linf:11.3e} {muc:9.4f} {relerr:10.3e}   "
              f"prior L2={pl2:.3e} dL2={abs(l2-pl2):.1e}  prior relerr={pre:.3e}")

    mug_rows = {}
    print("\n--- resolved MUG (recomputed from stored MUG profile CSVs) ---")
    print(f"{'Ha':>3} {'Ha_pred':>8} {'Ha_fit':>8} {'fit/pred':>8} "
          f"{'L2(int)':>10} {'L2(fit)':>10} {'Linf':>10} {'umean/uc':>9} {'relerr':>10}")
    for ha in MUG_HA:
        d = load_mug_csv(ha)
        if d is None:
            print(f"{ha:>3}  (no CSV)")
            continue
        y, ux = d
        Ha_pred = mug_ha_pred(ha)
        # fit Hartmann shape (Ha, amplitude) like analyze_mug.py
        uc0 = np.interp(0.0, y, ux)
        try:
            popt, _ = curve_fit(hartmann_shape, y, ux,
                                p0=[Ha_pred or ha, uc0], maxfev=20000)
            Ha_fit = abs(popt[0])
        except Exception as e:
            Ha_fit = float("nan")
        l2_int, linf, muc, relerr = metrics(y, ux, ha)          # vs analytic at integer Ha
        l2_fit, _, _, _ = metrics(y, ux, Ha_fit)                # vs analytic at fitted Ha
        mug_rows[ha] = dict(Ha_pred=Ha_pred, Ha_fit=Ha_fit, l2_int=l2_int,
                            l2_fit=l2_fit, linf=linf, muc=muc, relerr=relerr)
        ratio = Ha_fit / Ha_pred if Ha_pred else float("nan")
        print(f"{ha:>3} {Ha_pred:8.4f} {Ha_fit:8.4f} {ratio:8.4f} "
              f"{l2_int:10.3e} {l2_fit:10.3e} {linf:10.3e} {muc:9.4f} {relerr:10.3e}")

    # ---------- rebuild error_table.md (complete, honest blanks) ----------
    all_ha = sorted(set(OF_HA) | set(MUG_HA))
    lines = []
    lines.append("# Hartmann three-way verification — error table")
    lines.append("")
    lines.append("Profiles normalized by centerline velocity. L2/Linf are on the normalized "
                 "profile vs the analytic referee; u_mean/u_c is the discriminating scalar "
                 "(2/3 at Ha=0, ->1 plug core). MUG and mhdFoam errors are each computed on "
                 "that code's own grid vs analytic at the **integer** Ha.")
    lines.append("")
    def set_text(values):
        return "{" + ",".join(str(v) for v in values) + "}"

    shared_ha = sorted(set(OF_HA) & set(MUG_HA))
    lines.append(
        f"Coverage: mhdFoam was run at Ha in {set_text(OF_HA)}; the resolved-MUG channel "
        f"(drag OFF, body-force driven) has stored runs at Ha in {set_text(MUG_HA)}. "
        f"Both codes are present at Ha in {set_text(shared_ha)}. Blank (—) means no run "
        "was produced for that source at that Ha; no value is ever synthesized."
    )
    lines.append("")
    lines.append("| Ha | u_mean/u_c (analytic) | mhdFoam L2 | mhdFoam Linf | mhdFoam u_mean/u_c | mhdFoam relerr | MUG L2 | MUG Linf | MUG u_mean/u_c | MUG relerr |")
    lines.append("|----|----|----|----|----|----|----|----|----|----|")

    def e(x): return "—" if x is None else f"{x:.3e}"
    def f4(x): return "—" if x is None else f"{x:.4f}"

    for ha in all_ha:
        muc_ana = ah.u_mean_over_uc(ha)
        o = of_rows.get(ha)
        m = mug_rows.get(ha)
        ol2, oli, omuc, ore = (o if o else (None, None, None, None))
        if m:
            ml2, mli, mmuc, mre = m["l2_int"], m["linf"], m["muc"], m["relerr"]
        else:
            ml2 = mli = mmuc = mre = None
        lines.append(f"| {ha} | {muc_ana:.4f} | {e(ol2)} | {e(oli)} | {f4(omuc)} | "
                     f"{e(ore)} | {e(ml2)} | {e(mli)} | {f4(mmuc)} | {e(mre)} |")
    lines.append("")
    lines.append("MUG Ha-mapping detail (fit of the Hartmann shape to the computed profile, "
                 "vs the Ha predicted from inputs Ha = B0 a / sqrt(mu0 eta rho nu)):")
    lines.append("")
    lines.append("| Ha (pred) | Ha (fit) | fit/pred | L2 vs analytic(int Ha) | L2 vs analytic(fit Ha) |")
    lines.append("|----|----|----|----|----|")
    for ha in MUG_HA:
        m = mug_rows.get(ha)
        if not m:
            continue
        lines.append(f"| {m['Ha_pred']:.4f} | {m['Ha_fit']:.4f} | "
                     f"{m['Ha_fit']/m['Ha_pred']:.4f} | {m['l2_int']:.3e} | {m['l2_fit']:.3e} |")
    lines.append("")
    lines.append("_Regenerated by cases/hartmann/recompute_headline.py (packaging integrity pass, "
                 "directly from stored CSVs; mhdFoam from results/profile_Ha*.csv, MUG from "
                 "mug/Ha_*_profile.csv and mug_ginsburg/Ha_*/Ha_*_profile.csv)._")
    with open(os.path.join(RESULTS, "error_table.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nwrote", os.path.join(RESULTS, "error_table.md"))


if __name__ == "__main__":
    main()
