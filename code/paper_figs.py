"""H35 manuscript figures 1-4 from results/numbers.json only.  Usage: python paper_figs.py <numbers.json> <outdir>"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P = json.load(open(sys.argv[1], encoding="utf-8"))["paper"]; OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"
CH = ["CO(GT)", "C6H6(GT)", "NOx(GT)", "NO2(GT)"]; CHL = {"CO(GT)": "CO", "C6H6(GT)": "C$_6$H$_6$", "NOx(GT)": "NO$_x$", "NO2(GT)": "NO$_2$"}
BL = ["24", "72", "168"]; MK = {"24": "o", "72": "s", "168": "^"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK})

# ---- Figure 1: false alarm per monitor, 4 channels x 3 block lengths
groups = [("Heuristic / fixed-sample", [("M1_cusum", "CUSUM"), ("M2_pagehinkley", "Page-Hinkley"), ("M3_adwin", "ADWIN (per-step bound)"), ("M4_repeated_p", "Repeated $p<0.05$")]),
          ("Nominal guarantee at all times", [("M5_eprocess", "Betting e-process"), ("M6_ctm", "Conformal test martingale")]),
          ("Corrections", [("F1_sn_block", "Self-normalised blocks (F1)"), ("F2_bpi_wlc", "Adaptive block plug-in (F2)"), ("B1c_frozen_lrv_cusum", r"Frozen-LRV CUSUM, $\kappa=0.1$"), ("B2c_ar_residual_cusum", "AR-residual CUSUM")])]
rows = []; labels = []; seps = []
for gname, ms in groups:
    for m, lab in ms: rows.append(m); labels.append(lab)
    seps.append(len(rows))
fig, ax = plt.subplots(figsize=(6.6, 4.4))
cols = {"CO(GT)": "#0969da", "C6H6(GT)": "#bf8700", "NOx(GT)": "#8250df", "NO2(GT)": "#1a7f37"}
for i, m in enumerate(rows):
    for j, c in enumerate(CH):
        for k, b in enumerate(BL):
            v = P["t2"][c][b][m]; y = i + (j - 1.5) * 0.17 + (k - 1) * 0.045
            ax.plot([v[1], v[2]], [y, y], color=cols[c], lw=0.8); ax.plot(v[0], y, MK[b], color=cols[c], ms=3.2, mfc="white" if b == "24" else cols[c])
for s in seps[:-1]: ax.axhline(s - 0.5, color=GRID, lw=0.8)
ax.axvline(0.05, color=INK, lw=0.9, ls="--"); ax.text(0.052, -0.9, r"nominal $\alpha=0.05$", fontsize=8, color=INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels); ax.invert_yaxis(); ax.set_xscale("symlog", linthresh=0.01); ax.set_xlim(-0.0005, 1.05)
ax.set_xlabel("False-alarm probability in a 720-h in-tolerance window (Wilson 95 % interval)")
h = [plt.Line2D([], [], color=cols[c], marker="o", ls="", ms=4, label=CHL[c]) for c in CH] + [plt.Line2D([], [], color=MUTED, marker=MK[b], ls="", ms=4, mfc="white" if b == "24" else MUTED, label=f"block {b} h") for b in BL]
ax.legend(handles=h, frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.42, -0.14), ncol=7, columnspacing=0.8, handletextpad=0.3); ax.grid(axis="x", color=GRID, lw=0.5)
fig.set_size_inches(6.6, 4.8); fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_false_alarm.pdf")); plt.close(fig)

# ---- Figure 2: delay / D* versus alpha (synthetic in-scope nulls)
xs = np.array([0.05, 1e-3, 1e-6]); al = ["0.05", "0.001", "1e-06"]
fig, ax = plt.subplots(figsize=(4.8, 3.3))
for k in P["insc"]:
    ax.plot(xs, [P["delay_ratio_by_alpha"][k]["F2_bpi_wlc"][a] for a in al], color="#0969da", lw=1, alpha=0.6, marker="o", ms=2.5)
    ax.plot(xs, [P["delay_ratio_by_alpha"][k]["F1_sn_block"][a] for a in al], color="#8250df", lw=1, alpha=0.6, marker="s", ms=2.5)
ax.plot(xs, [P["pred_overhead"][a] for a in al], color=INK, lw=1.6, ls="--")
ax.plot([], [], color="#0969da", marker="o", ms=3, label="F2, 7 null models"); ax.plot([], [], color="#8250df", marker="s", ms=3, label="F1, 7 null models")
ax.plot([], [], color=INK, ls="--", lw=1.6, label="predicted overhead of F2"); ax.legend(frameon=False, fontsize=7, loc="lower left")
ax.axhline(1, color=MUTED, lw=0.8); ax.set_xscale("log"); ax.invert_xaxis(); ax.set_yscale("log")
ax.set_xlabel(r"nominal level $\alpha$"); ax.set_ylabel(r"mean delay / information bound $D^*$"); ax.grid(color=GRID, lw=0.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_delay_ratio.pdf")); plt.close(fig)

# ---- Figure 3: burn-in sweep
n0s = ["336", "1000", "3000", "10000"]; x = np.array([336, 1000, 3000, 10000])
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.9), sharex=True)
for k, v in P["e8"].items():
    real = k.startswith("AQ_"); col = "#bf8700" if real else "#0969da"; ls = "--" if real else "-"
    axs[0].plot(x, [v[n]["viol"][0] for n in n0s], color=col, ls=ls, lw=1, marker="o", ms=2.5, alpha=0.8)
    axs[1].plot(x, [v[n]["fa_1overalpha"][0] for n in n0s], color=col, ls=ls, lw=1, marker="o", ms=2.5, alpha=0.8)
axs[1].axhline(0.05, color=INK, lw=0.9, ls="--")
axs[0].set_ylabel("P(variance under-estimated)"); axs[1].set_ylabel(r"false alarm of F2 with threshold $1/\alpha$")
for a in axs: a.set_xscale("log"); a.set_xlabel("in-control burn-in $n_0$ (h)"); a.grid(color=GRID, lw=0.5)
axs[0].plot([], [], color="#0969da", label="synthetic (7)"); axs[0].plot([], [], color="#bf8700", ls="--", label="UCI channels (4)"); axs[0].legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_burnin.pdf")); plt.close(fig)

# ---- Figure 4: false alarm vs delay at alpha = 0.05
fig, ax = plt.subplots(figsize=(5.4, 3.6))
sty = {"F1": ("#8250df", "s", "F1"), "F2": ("#0969da", "o", "F2"), "B1c": ("#cf222e", "v", r"frozen-LRV CUSUM $\kappa=0.1$"),
       "B1c_k0.5": ("#1a7f37", "D", r"frozen-LRV CUSUM $\kappa=0.5$"), "B1c_hcal": ("#bf8700", "^", "frozen-LRV CUSUM, raised threshold")}
keys = P["insc"] + [f"AQ_{c}|boot168" for c in CH]
for m, (col, mk, lab) in sty.items():
    xs_, ys_ = [], []
    for k in keys:
        fa = P["t4_fa"][k][m][0] if not k.startswith("AQ_") else max(P["t4_fa"][k.replace("168", b)][m][0] for b in BL)
        xs_.append(fa); ys_.append(P["t4_delay"][k][m]["mean_delay"] / 60.0)
    ax.plot(xs_, ys_, mk, color=col, ms=4, alpha=0.85, label=lab, ls="")
for k, v in P["oos"].items():
    ax.plot(v["B1c_k0.5"]["FA"][0], v["B1c_k0.5"]["mean_delay"] / 60.0, "D", mfc="white", mec="#1a7f37", ms=6)
    ax.plot(v["F2_bpi_wlc"]["FA"][0], v["F2_bpi_wlc"]["mean_delay"] / 60.0, "o", mfc="white", mec="#0969da", ms=6)
ax.plot([], [], "D", mfc="white", mec=MUTED, ls="", label="unseen null models (open markers)")
ax.axvline(0.05, color=INK, lw=0.9, ls="--"); ax.set_xscale("symlog", linthresh=0.005); ax.set_xlim(-0.0002, 0.1)
ax.set_xlabel("false-alarm probability (UCI: worst of three block lengths)"); ax.set_ylabel(r"mean delay / $D^*$ at $\alpha=0.05$")
ax.legend(frameon=False, fontsize=7, loc="upper right"); ax.grid(color=GRID, lw=0.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_fa_delay.pdf")); plt.close(fig)
print("figs ok")
