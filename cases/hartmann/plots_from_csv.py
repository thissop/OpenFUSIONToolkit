"""
Regenerate all Hartmann result plots DIRECTLY FROM THE STORED CSVs (packaging pass).

No meshio / OpenFOAM dependency: mhdFoam+analytic come from results/profile_Ha*.csv,
resolved MUG from mug/Ha_{2,3,5}_profile.csv. Rebuilds:
  results/overlay_Ha*.png            (analytic vs mhdFoam vs MUG, per Ha that has data)
  results/flowrate_vs_Ha.png         (analytic Q* ~ 1/Ha)
  results/report/threeway_overlay_Ha5.png, umean_uc_vs_Ha.png, error_vs_Ha.png,
  results/report/flowrate_vs_Ha.png  (smplotlib-styled)
"""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
import sys; sys.path.insert(0, HERE)
import analytic_hartmann as ah

RESULTS = os.path.join(HERE, "results")
REPORT = os.path.join(RESULTS, "report")
MUG = os.path.join(HERE, "mug")
A = 1.0
OF_HA = [1, 5, 10, 20, 50]
MUG_HA = [2, 3, 5]


def trapz(y, x):
    f = getattr(np, "trapezoid", None) or np.trapz
    return f(y, x)


def norm_uc(y, u):
    uc = np.interp(0.0, y, u)
    return u / uc, uc


def umean_uc(y, u):
    uc = np.interp(0.0, y, u)
    return (trapz(u, y) / (y[-1] - y[0])) / uc


def l2_vs_analytic(y, u, Ha):
    un, _ = norm_uc(y, u)
    ana, _ = norm_uc(y, ah.u_profile(y, A, Ha))
    d = un - ana
    return np.sqrt(trapz(d**2, y) / trapz(ana**2, y))


def load_of(ha):
    path = os.path.join(RESULTS, f"profile_Ha{ha}.csv")
    if not os.path.exists(path):
        return None
    rows = [l.strip().split(",") for l in open(path).read().splitlines()[1:]]
    yy = [(float(p[0]), float(p[2])) for p in rows if len(p) >= 3 and p[2] != ""]
    if not yy:
        return None
    a = np.array(yy)
    return a[:, 0], a[:, 1]          # y, u_of_norm


def load_mug(ha):
    path = os.path.join(MUG, f"Ha_{ha}_profile.csv")
    if not os.path.exists(path):
        return None
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    o = np.argsort(d[:, 0])
    return d[o, 0], d[o, 1]


# ---------------- basic overlays (plain matplotlib) ----------------
def basic_overlays():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    for ha in sorted(set(OF_HA) | set(MUG_HA)):
        of = load_of(ha)
        mug = load_mug(ha)
        if of is None and mug is None:
            continue
        yA = np.linspace(-A, A, 400)
        anaN, _ = norm_uc(yA, ah.u_profile(yA, A, ha))
        plt.figure(figsize=(6, 4.2))
        plt.plot(anaN, yA, "k-", lw=2, label="analytic")
        if of is not None:
            ofn, _ = norm_uc(of[0], of[1])
            plt.plot(ofn, of[0], "o", ms=3, mfc="none", color="C0", label="mhdFoam")
        if mug is not None:
            mn, _ = norm_uc(mug[0], mug[1])
            plt.plot(mn, mug[0], "s", ms=3, color="C3", label="MUG")
        plt.xlabel("u / u_centerline"); plt.ylabel("y / a")
        plt.title(f"Hartmann profile, Ha = {ha}")
        plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
        plt.savefig(os.path.join(RESULTS, f"overlay_Ha{ha}.png"), dpi=130)
        plt.close()
    # flow-rate
    Ha = np.array([1, 2, 5, 10, 20, 50], float)
    plt.figure(figsize=(6, 4.2))
    plt.loglog(Ha, [ah.flow_rate_star(h) for h in Ha], "k-o", label="analytic Q*")
    plt.loglog(Ha, 1.0 / Ha, "k--", alpha=0.5, label="1/Ha reference")
    plt.xlabel("Ha"); plt.ylabel("Q* = (1/Ha)[1 - tanh(Ha)/Ha]")
    plt.title("Hartmann flow-rate suppression")
    plt.legend(); plt.grid(alpha=0.3, which="both"); plt.tight_layout()
    plt.savefig(os.path.join(RESULTS, "flowrate_vs_Ha.png"), dpi=130)
    plt.close()


# ---------------- report-quality (smplotlib) ----------------
def report_plots():
    os.makedirs(REPORT, exist_ok=True)
    try:
        import smplotlib  # noqa
    except Exception as ex:
        print("smplotlib unavailable, using default style:", ex)
    import matplotlib.pyplot as plt

    of = {h: load_of(h) for h in OF_HA if load_of(h) is not None}
    mug = {h: load_mug(h) for h in MUG_HA if load_mug(h) is not None}

    # 1. three-way overlay Ha=5 + residual inset
    if 5 in of and 5 in mug:
        yA = np.linspace(-A, A, 400)
        anaA, _ = norm_uc(yA, ah.u_profile(yA, A, 5))
        yo, uo = of[5]; uon, _ = norm_uc(yo, uo)
        ym, um = mug[5]; umn, _ = norm_uc(ym, um)
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        ax.plot(anaA, yA, "-", color="k", lw=1.6, label="analytic")
        ax.plot(uon, yo, "o", ms=4.5, mfc="none", mec="C0", label="mhdFoam")
        ax.plot(umn, ym, "s", ms=4.0, mfc="none", mec="C3", label="MUG (resolved)")
        ax.set_xlabel(r"$u/u_{\rm centerline}$"); ax.set_ylabel(r"$y/a$")
        ax.set_title(r"Hartmann channel, $\mathrm{Ha}=5$ (three-way)")
        ax.legend(loc="center left")
        axin = fig.add_axes([0.62, 0.17, 0.30, 0.28])
        anaO, _ = norm_uc(yo, ah.u_profile(yo, A, 5))
        anaM, _ = norm_uc(ym, ah.u_profile(ym, A, 5))
        axin.plot(yo, uon - anaO, "o", ms=2.5, mfc="none", mec="C0")
        axin.plot(ym, umn - anaM, "s", ms=2.0, mfc="none", mec="C3")
        axin.axhline(0, color="k", lw=0.6)
        axin.set_title("residual vs analytic", fontsize=8)
        axin.set_xlabel(r"$y/a$", fontsize=7); axin.tick_params(labelsize=6)
        fig.tight_layout(); fig.savefig(os.path.join(REPORT, "threeway_overlay_Ha5.png"), dpi=150)
        plt.close(fig)

    # 2. u_mean/u_c vs Ha
    Hc = np.logspace(0, np.log10(60), 200)
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(Hc, [ah.u_mean_over_uc(h) for h in Hc], "-", color="k", lw=1.4, label="analytic")
    if of:
        hs = sorted(of); ax.plot(hs, [umean_uc(*of[h]) for h in hs], "o", ms=7,
                                 mfc="none", mec="C0", label="mhdFoam")
    if mug:
        hs = sorted(mug); ax.plot(hs, [umean_uc(*mug[h]) for h in hs], "s", ms=6,
                                  mfc="none", mec="C3", label="MUG")
    ax.set_xscale("log"); ax.set_xlabel(r"$\mathrm{Ha}$")
    ax.set_ylabel(r"$u_{\rm mean}/u_{\rm centerline}$")
    ax.set_title(r"Discriminating scalar: $2/3$ (Poiseuille) $\to 1$ (plug core)")
    ax.legend(loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(REPORT, "umean_uc_vs_Ha.png"), dpi=150)
    plt.close(fig)

    # 3. error vs Ha
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    if of:
        hs = sorted(of); ax.plot(hs, [l2_vs_analytic(*of[h], h) for h in hs], "o-",
                                 ms=7, mfc="none", mec="C0", color="C0", label=r"mhdFoam $L_2$")
    if mug:
        hs = sorted(mug); ax.plot(hs, [l2_vs_analytic(*mug[h], h) for h in hs], "s-",
                                  ms=6, mfc="none", mec="C3", color="C3", label=r"MUG $L_2$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"$\mathrm{Ha}$")
    ax.set_ylabel(r"normalized-profile $L_2$ error vs analytic")
    ax.set_title(r"Profile error vs analytic (both codes $\lesssim 10^{-3}$)")
    ax.legend(loc="upper left"); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(REPORT, "error_vs_Ha.png"), dpi=150)
    plt.close(fig)

    # 4. analytic flow-rate
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    Hf = np.array([1, 2, 5, 10, 20, 50], float)
    ax.loglog(Hf, [ah.flow_rate_star(h) for h in Hf], "o-", color="k", ms=6, label=r"analytic $Q^*$")
    ax.loglog(Hf, 1.0 / Hf, "--", color="0.5", label=r"$1/\mathrm{Ha}$ reference")
    ax.set_xlabel(r"$\mathrm{Ha}$")
    ax.set_ylabel(r"$Q^* = \mathrm{Ha}^{-1}[1-\tanh(\mathrm{Ha})/\mathrm{Ha}]$")
    ax.set_title("Hartmann flow-rate suppression")
    ax.legend(); ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(REPORT, "flowrate_vs_Ha.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    basic_overlays()
    report_plots()
    import glob
    print("Regenerated:")
    for p in sorted(glob.glob(os.path.join(RESULTS, "*.png")) +
                    glob.glob(os.path.join(REPORT, "*.png"))):
        print("  ", os.path.relpath(p, HERE))
