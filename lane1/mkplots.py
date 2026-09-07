"""Figures for REPORT_lane1.md."""
import math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import engine as E
import ahsweep as A
import sweep as W
import seeds as S

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figs")
os.makedirs(FIG, exist_ok=True)


def ah_panel():
    tab = W.seed_table()
    names = list(tab)
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    for ax, name in zip(axes.ravel(), names):
        L = tab[name]["L"]
        phi = W.best_phi(L)
        r = A.evaluate(L, phi, n=320)
        f, s, b, st = r
        run = A.longest_run(st, b)
        ss, bb = s[run], b[run]
        ax.plot(ss, bb, lw=1.1, color="#1f4e79")
        for e in f["extrema"]:
            ax.plot(e["s"], e["b"], "o", ms=5, color="#c0392b")
        if len(f["extrema"]) == 2:
            h = sorted(e["b"] for e in f["extrema"])
            ax.axhspan(h[0], h[1], color="#c0392b", alpha=0.12)
        ax.set_xscale("log")
        ax.set_title(f"{name}: {f['n_extrema']} interior extrema\n"
                     f"3-cycle window {f['wiggle_range']:.2e} in b", fontsize=9)
        ax.set_xlabel("s (ray from the focus)", fontsize=8)
        ax.set_ylabel(r"$\beta^*(s)$", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.suptitle("Andronov-Hopf curves of the nine fat seeds.  A horizontal line "
                 "meets the curve once per limit cycle;\nthe shaded band is the "
                 "rotation window in which three cycles exist.", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(FIG, "ah_curves.png"), dpi=130)
    plt.close(fig)


def row4_panel():
    sd = S.cherkas_seed(4)
    L = E.local10(sd["vec12"], (1.0, -1.0))
    phi = np.pi
    f, s, b, st = A.evaluate(L, phi, n=400)
    run = A.longest_run(st, b)
    ss, bb = s[run], b[run]
    x = 1.0 + np.cos(phi) * ss
    m = (x >= 0.55) & (x <= 0.92)

    coef = [8.89863, 4.39482, -13.5991, 22.9703, -22.4248, 11.9886, -2.72941]
    xx = np.linspace(0.6, 0.9, 2001)
    ah = sum(c * xx ** k for k, c in enumerate(coef))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.6))
    a1.plot(x[m], bb[m], lw=1.3, color="#1f4e79")
    for e in f["extrema"]:
        xe = 1.0 + np.cos(phi) * e["s"]
        if 0.55 <= xe <= 0.92:
            a1.axvline(xe, color="#c0392b", ls="--", lw=0.9)
            a1.text(xe, a1.get_ylim()[0], f" {xe:.4f}", fontsize=8, color="#c0392b")
    a1.set_title(r"this engine: $\beta^*$ (uniform rotation) on $y=-1$", fontsize=10)
    a1.set_xlabel("x"); a1.set_ylabel(r"$\beta^*$")

    a2.plot(xx, ah, lw=1.3, color="#7d3c98")
    d = np.gradient(ah, xx)
    for i in np.where(d[:-1] * d[1:] < 0)[0]:
        a2.axvline(xx[i], color="#c0392b", ls="--", lw=0.9)
        a2.text(xx[i], a2.get_ylim()[0], f" {xx[i]:.4f}", fontsize=8, color="#c0392b")
    a2.set_title("Cherkas-Artes-Llibre 2003, published AH(x) for row 4", fontsize=10)
    a2.set_xlabel("x"); a2.set_ylabel(r"$a_{11}$")
    fig.suptitle("Cherkas row 4: the two interior extrema land at x = 0.6207, 0.7981 "
                 "against the paper's 0.6238, 0.8056.\nThe two curves are different "
                 "rotation parameters of the same field, so only the shape and the "
                 "extremum positions are comparable.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(os.path.join(FIG, "row4_vs_published.png"), dpi=130)
    plt.close(fig)


def displacement_panel():
    sd = S.cherkas_seed(4)
    L = E.local10(sd["vec12"], (1.0, -1.0))
    phi = np.pi
    s = np.geomspace(1e-3, 0.5337 * (1 - 1e-9), 500)
    D, st, noise, T = E.d_curve_noisy(L, phi, s)
    ok = st == 0
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.plot(s[ok], D[ok], lw=1.1, color="#1f4e79", label="D(s, 0)")
    ax.fill_between(s[ok], -noise[ok], noise[ok], color="#c0392b", alpha=0.35,
                    label="two-tolerance noise band")
    ax.axhline(0, color="k", lw=0.6)
    br = E.count_sign_changes(s, D, st, noise)
    for (s1, s2, d1, d2) in br:
        ax.axvline(0.5 * (s1 + s2), color="#27ae60", ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1e-9)
    ax.set_xlabel("s (ray from the focus A=(1,-1) towards -x)")
    ax.set_ylabel("D")
    ax.set_title(f"Cherkas row 4: displacement on the section, "
                 f"{len(br)} sign changes above the noise band", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "row4_displacement.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    ah_panel(); row4_panel(); displacement_panel()
    print("figures written to", FIG)
