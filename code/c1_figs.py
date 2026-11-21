"""H35 1c figures -> results/fig/*.pdf.  Usage: python c1_figs.py <results_dir>"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]; F = os.path.join(R, "fig"); os.makedirs(F, exist_ok=True)
L = lambda f: json.load(open(os.path.join(R, f), encoding="utf-8"))["result"]
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"
C_NEW, C_OLD, C_BASE = "#1f6feb", "#8250df", "#bf8700"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED})

# Fig 1: false alarm (720 h window) per monitor, UCI AQ boot168 vs synthetic AR0.8
e1s, e1a = L("c1_e1_syn.json"), L("c1_e1_aq.json")
mons = list(e1s["AR1_0.8"].keys())
fig, ax = plt.subplots(figsize=(6.4, 3.2))
y = np.arange(len(mons))
for j, (lab, d, col) in enumerate([("AR(1) φ=0.8", e1s["AR1_0.8"], C_OLD), ("UCI CO(GT), bootstrap 168 h", e1a["AQ_CO(GT)|boot168"], C_NEW)]):
    v = np.array([d[m]["rate_wilson95"][0] for m in mons]); lo = np.array([d[m]["rate_wilson95"][1] for m in mons]); hi = np.array([d[m]["rate_wilson95"][2] for m in mons])
    ax.errorbar(v, y + (j - 0.5) * 0.3, xerr=[v - lo, hi - v], fmt="o", ms=3.5, color=col, label=lab, capsize=0, lw=1)
ax.axvline(0.05, color=INK, lw=0.8, ls="--"); ax.text(0.055, len(mons) - 0.4, "α = 0.05", color=INK, fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(mons); ax.set_xscale("symlog", linthresh=0.01); ax.set_xlim(-0.001, 1.05)
ax.set_xlabel("false-alarm probability in a 720-step window (Wilson 95%)"); ax.grid(axis="x", color=GRID, lw=0.5)
ax.legend(frameon=False, fontsize=8, loc="lower right"); fig.tight_layout(); fig.savefig(os.path.join(F, "fig1_false_alarm.pdf")); plt.close(fig)

# Fig 2: delay ratio vs alpha, F2 measured vs predicted overhead, F1
e2 = L("c1_e2_delay.json"); alphas = ["0.05", "0.001", "1e-06"]
fig, ax = plt.subplots(figsize=(4.8, 3.2)); xs = np.array([0.05, 1e-3, 1e-6])
for name, d in e2.items():
    if not name.startswith("AR1") and not name.startswith("ARMA"): continue
    r = [d["cond_delay_nu240"]["F2_bpi_wlc"][a]["ratio_to_Dstar"] for a in alphas]
    ax.plot(xs, r, color=C_NEW, lw=1, alpha=0.55, marker="o", ms=2.5)
    r1 = [d["cond_delay_nu240"]["F1_sn_block"][a]["ratio_to_Dstar"] for a in alphas]
    ax.plot(xs, r1, color=C_OLD, lw=1, alpha=0.55, marker="s", ms=2.5)
p = [list(e2.values())[0]["pred_overhead"][a] for a in alphas]
ax.plot(xs, p, color=INK, lw=1.6, ls="--")
ax.text(1.3e-6, p[2] * 0.82, "predicted overhead", color=INK, fontsize=8)
ax.text(1.3e-6, 2.6, "BPI-WLC (F2)", color=C_NEW, fontsize=8); ax.text(1.3e-6, 9.5, "self-normalized blocks (F1)", color=C_OLD, fontsize=8)
ax.axhline(1, color=MUTED, lw=0.8); ax.set_xscale("log"); ax.invert_xaxis(); ax.set_yscale("log")
ax.set_xlabel("α"); ax.set_ylabel("mean delay / Tartakovsky D*"); ax.grid(color=GRID, lw=0.5)
fig.tight_layout(); fig.savefig(os.path.join(F, "fig2_delay_ratio.pdf")); plt.close(fig)

# Fig 3: Gas drift batch-level mean residual (16 sensors)
g = L("c1_gas.json")["mean_resid"]
fig, ax = plt.subplots(figsize=(5.2, 3.0)); b = np.arange(1, 11)
for k, row in g.items():
    ax.plot(b, [row[str(i)][0] for i in b], color=MUTED, lw=0.8, alpha=0.7)
med = np.median(np.array([[row[str(i)][0] for i in b] for row in g.values()]), 0)
ax.plot(b, med, color=C_NEW, lw=2); ax.text(10.1, med[-1], "median of 16", color=C_NEW, fontsize=8, va="center")
ax.axhline(0, color=INK, lw=0.6); ax.set_xticks(b); ax.set_xlabel("batch (≈ month index, not equally spaced)"); ax.set_ylabel("mean normalized residual")
ax.grid(color=GRID, lw=0.5); fig.tight_layout(); fig.savefig(os.path.join(F, "fig3_gas_batch_drift.pdf")); plt.close(fig)
print("figs ok")
