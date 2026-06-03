"""
Presentation-quality report plots for the three-way Hartmann verification.
Styled with smplotlib (classic AAS/matplotlib serif scientific look).

Produces (in cases/hartmann/results/report/):
  1. threeway_overlay_Ha5.png  — analytic vs mhdFoam vs MUG at Ha=5, with a
     residual inset (difference from analytic) so the agreement is quantitative.
  2. umean_uc_vs_Ha.png         — discriminating scalar u_mean/u_c vs Ha:
     analytic curve + mhdFoam points + MUG points (true cross-code scaling).
  3. error_vs_Ha.png            — profile L2 error and |u_mean/u_c| rel. error
     vs Ha for mhdFoam and MUG (quantitative agreement).
  4. flowrate_vs_Ha.png         — analytic Hartmann flow-rate suppression Q* ~ 1/Ha.

MUG points are read from cases/hartmann/mug/Ha_<n>_profile.csv (binned y,ux), which
analyze_mug.py writes; OpenFOAM profiles via foamToVTK+meshio (extract_profile.py).
"""
import os, glob, re
import numpy as np

import smplotlib  # noqa: F401  (applies the style on import)
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)
import analytic_hartmann as ah
from extract_profile import get_profile

OPENFOAM = os.path.join(HERE, "openfoam")
MUG = os.path.join(HERE, "mug")
OUT = os.path.join(HERE, "results", "report")
A = 1.0
OF_HA = [1, 5, 10, 20, 50]


def norm_uc(y, u):
    uc = np.interp(0.0, y, u)
    return u / uc, uc


def umean_uc(y, u):
    uc = np.interp(0.0, y, u)
    return (np.trapezoid(u, y) / (y[-1] - y[0])) / uc


def l2_vs_analytic(y, u, Ha):
    un, _ = norm_uc(y, u)
    ana = ah.u_profile(y, A, Ha); ana, _ = norm_uc(y, ana)
    d = un - ana
    return np.sqrt(np.trapezoid(d**2, y) / np.trapezoid(ana**2, y))


def load_of():
    out = {}
    for ha in OF_HA:
        case = os.path.join(OPENFOAM, f"Ha_{ha}")
        if not os.path.isdir(case):
            continue
        res = get_profile(case, xtarget=10.0)
        if res:
            y, ux = res
            out[ha] = (y, ux)
    return out


def load_mug():
    out = {}
    for f in glob.glob(os.path.join(MUG, "Ha_*_profile.csv")):
        m = re.search(r"Ha_(\d+)_profile", os.path.basename(f))
        if not m:
            continue
        ha = int(m.group(1))
        d = np.loadtxt(f, delimiter=",", skiprows=1)
        y, ux = d[:, 0], d[:, 1]
        o = np.argsort(y)
        out[ha] = (y[o], ux[o])
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    of = load_of()
    mug = load_mug()
    print("OpenFOAM Ha:", sorted(of), " MUG Ha:", sorted(mug))

    # ---- 1. three-way overlay at Ha=5 with residual inset ----
    if 5 in of and 5 in mug:
        yA = np.linspace(-A, A, 400)
        anaA, _ = norm_uc(yA, ah.u_profile(yA, A, 5))
        yo, uo = of[5]; uon, _ = norm_uc(yo, uo)
        ym, um = mug[5]; umn, _ = norm_uc(ym, um)
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        ax.plot(anaA, yA, "-", color="k", lw=1.6, label="analytic", zorder=1)
        ax.plot(uon, yo, "o", ms=4.5, mfc="none", mec="C0", label="mhdFoam", zorder=2)
        ax.plot(umn, ym, "s", ms=4.0, mfc="none", mec="C3", label="MUG (resolved)", zorder=3)
        ax.set_xlabel(r"$u/u_{\rm centerline}$")
        ax.set_ylabel(r"$y/a$")
        ax.set_title(r"Hartmann channel, $\mathrm{Ha}=5$ (three-way)")
        ax.legend(loc="center left")
        # residual inset: difference from analytic on each code's grid
        axin = fig.add_axes([0.62, 0.17, 0.30, 0.28])
        anaO, _ = norm_uc(yo, ah.u_profile(yo, A, 5))
        anaM, _ = norm_uc(ym, ah.u_profile(ym, A, 5))
        axin.plot(yo, (uon - anaO), "o", ms=2.5, mfc="none", mec="C0")
        axin.plot(ym, (umn - anaM), "s", ms=2.0, mfc="none", mec="C3")
        axin.axhline(0, color="k", lw=0.6)
        axin.set_title("residual vs analytic", fontsize=8)
        axin.set_xlabel(r"$y/a$", fontsize=7)
        axin.tick_params(labelsize=6)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, "threeway_overlay_Ha5.png"), dpi=150)
        plt.close(fig)

    # ---- 2. u_mean/u_c vs Ha (cross-code scaling) ----
    Hc = np.logspace(0, np.log10(60), 200)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(Hc, [ah.u_mean_over_uc(h) for h in Hc], "-", color="k", lw=1.4,
            label="analytic")
    if of:
        hs = sorted(of); ax.plot(hs, [umean_uc(*of[h]) for h in hs], "o",
                                 ms=7, mfc="none", mec="C0", label="mhdFoam")
    if mug:
        hs = sorted(mug); ax.plot(hs, [umean_uc(*mug[h]) for h in hs], "s",
                                  ms=6, mfc="none", mec="C3", label="MUG")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\mathrm{Ha}$")
    ax.set_ylabel(r"$u_{\rm mean}/u_{\rm centerline}$")
    ax.set_title(r"Discriminating scalar: $2/3$ (Poiseuille) $\to 1$ (plug core)")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "umean_uc_vs_Ha.png"), dpi=150)
    plt.close(fig)

    # ---- 3. error vs Ha ----
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    if of:
        hs = sorted(of)
        ax.plot(hs, [l2_vs_analytic(*of[h], h) for h in hs], "o-", ms=7, mfc="none",
                mec="C0", color="C0", label="mhdFoam  $L_2$")
    if mug:
        hs = sorted(mug)
        ax.plot(hs, [l2_vs_analytic(*mug[h], h) for h in hs], "s-", ms=6, mfc="none",
                mec="C3", color="C3", label="MUG  $L_2$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\mathrm{Ha}$")
    ax.set_ylabel(r"normalized-profile $L_2$ error vs analytic")
    ax.set_title("Profile error vs analytic (both codes $\\lesssim 10^{-3}$)")
    ax.legend(loc="upper left")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "error_vs_Ha.png"), dpi=150)
    plt.close(fig)

    # ---- 4. analytic flow-rate suppression ----
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    Hf = np.array([1, 2, 5, 10, 20, 50], float)
    ax.loglog(Hf, [ah.flow_rate_star(h) for h in Hf], "o-", color="k", ms=6,
              label=r"analytic $Q^*$")
    ax.loglog(Hf, 1.0 / Hf, "--", color="0.5", label=r"$1/\mathrm{Ha}$ reference")
    ax.set_xlabel(r"$\mathrm{Ha}$")
    ax.set_ylabel(r"$Q^* = \mathrm{Ha}^{-1}[1-\tanh(\mathrm{Ha})/\mathrm{Ha}]$")
    ax.set_title("Hartmann flow-rate suppression")
    ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "flowrate_vs_Ha.png"), dpi=150)
    plt.close(fig)

    print("wrote plots to", OUT)
    for p in sorted(glob.glob(os.path.join(OUT, "*.png"))):
        print("  ", p)


if __name__ == "__main__":
    main()
