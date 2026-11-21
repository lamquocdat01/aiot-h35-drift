"""H35: numbers.json -> manuscript/numbers.tex (macros) + manuscript/tables/*.tex (no colour).
Usage: python make_tex.py <numbers.json> <manuscript_dir>"""
import json, os, sys

N = json.load(open(sys.argv[1], encoding="utf-8")); P = N["paper"]; M = sys.argv[2]
os.makedirs(os.path.join(M, "tables"), exist_ok=True)
CH = ["CO(GT)", "C6H6(GT)", "NOx(GT)", "NO2(GT)"]; CHL = {"CO(GT)": "CO", "C6H6(GT)": r"C$_6$H$_6$", "NOx(GT)": r"NO$_x$", "NO2(GT)": r"NO$_2$"}
BL = ["24", "72", "168"]
NULL_L = {"AR1_0": r"AR(1), $\phi=0$", "AR1_0.5": r"AR(1), $\phi=0.5$", "AR1_0.8": r"AR(1), $\phi=0.8$", "AR1_0.95": r"AR(1), $\phi=0.95$",
          "ARMA11_0.8_+0.4": r"ARMA(1,1), $0.8/{+0.4}$", "ARMA11_0.8_-0.5": r"ARMA(1,1), $0.8/{-0.5}$", "ARMA21_1.2_-0.4_+0.3": r"ARMA(2,1)"}


def p3(x): return f"{x:.3f}"
def pc(x): return f"{100 * x:.0f}\\,\\%" if x >= 0.095 else f"{100 * x:.1f}\\,\\%"
def ci(v): return f"{v[0]:.3f} [{v[1]:.3f}, {v[2]:.3f}]"


mac = {}
r = P["rq1_guaranteed_range"]; mac["RqOneLo"] = pc(r[0]); mac["RqOneHi"] = pc(r[1]); mac["RqOneMinLowerCI"] = p3(P["rq1_guaranteed_min_lowerCI"])
for b in BL:
    mac[f"RqOneLoB{ {'24':'A','72':'B','168':'C'}[b] }".replace(" ", "")] = pc(P["rq1_by_block"][b][0])
    mac[f"RqOneHiB{ {'24':'A','72':'B','168':'C'}[b] }".replace(" ", "")] = pc(P["rq1_by_block"][b][1])
mac["AdwinLo"] = pc(P["rq1_adwin_range"][0]); mac["AdwinHi"] = pc(P["rq1_adwin_range"][1])
t3 = P["t3"]
mac["DeseasonE"] = p3(t3["deseason_block24"]["eprocess"][0]); mac["DeseasonEci"] = ci(t3["deseason_block24"]["eprocess"])
mac["DeseasonC"] = p3(t3["deseason_block24"]["ctm"][0]); mac["IidE"] = p3(t3["iid"]["eprocess"][0]); mac["IidC"] = p3(t3["iid"]["ctm"][0])
mac["IidAdwin"] = p3(t3["adwin_iid"][0]); mac["VarHour"] = pc(t3["var_hour_profile"]); mac["RawE"] = p3(t3["raw_block24"]["eprocess"][0]); mac["RawC"] = p3(t3["raw_block24"]["ctm"][0])
mac["BlockMeanE"] = p3(t3["blockmean24_eproc"][0]); mac["WhitenE"] = p3(t3["whiten_phihat_eproc"][0])
mac["AcfCO"] = f"{P['meta']['CO(GT)']['acf1']:.2f}"
mac["IntelN"] = str(P["intel_n"]); mac["IntelFtwoExceed"] = str(P["intel_F2_n_exceed"])
f1max = max(P["t4_fa"][k]["F1"][0] for k in P["t4_fa"]); f2max = max(P["t4_fa"][k]["F2"][0] for k in P["t4_fa"])
mac["FoneMax"] = p3(f1max); mac["FtwoMax"] = p3(f2max)
c = P["e7_counts"]
for v, tag in (("B1c_k0.1", "Kzero"), ("B1c_k0.25", "Kquarter"), ("B1c_k0.5", "Khalf"), ("B1c_k0.1_hcal", "Hcal"), ("F2_bpi_wlc", "Ftwo")):
    for a, at in (("0.05", "A"), ("0.001", "B")):
        mac[f"Keep{tag}{at}"] = str(c[f"{v}|{a}"]["keep"]); mac[f"KeepFast{tag}{at}"] = str(c[f"{v}|{a}"]["keep_and_faster"])
        mac[f"RatioLo{tag}{at}"] = f"{c[f'{v}|{a}']['delay_ratio_range'][0]:.2f}"; mac[f"RatioHi{tag}{at}"] = f"{c[f'{v}|{a}']['delay_ratio_range'][1]:.2f}"
mac["HcalA"] = f"{P['e7_h']['0.05']:.2f}"; mac["HcalB"] = f"{P['e7_h']['0.001']:.1f}"
mac["FtwoExceedB"] = str(P["f2_fa_1e-3_n_exceed"]); mac["FtwoBmax"] = f"{max(v[0] for v in P['f2_fa_1e-3'].values()):.4f}"
oo = P["oos"]; ks = list(oo)
mac["OosAfa"] = p3(oo[ks[0]]["B1c_k0.5"]["FA"][0]); mac["OosBfa"] = p3(oo[ks[1]]["B1c_k0.5"]["FA"][0])
mac["OosAratio"] = f"{oo[ks[0]]['B1c_k0.5']['mean_delay'] / oo[ks[0]]['F2_bpi_wlc']['mean_delay']:.2f}"
mac["OosBratio"] = f"{oo[ks[1]]['B1c_k0.5']['mean_delay'] / oo[ks[1]]['F2_bpi_wlc']['mean_delay']:.2f}"
mac["OosUciMax"] = p3(P["oos_uci_B1c_FA_max"])
os_ = P["out_of_scope"]
mac["LongFtwo"] = p3(os_["longmem"]["F2_bpi_wlc"][0]); mac["LongFone"] = p3(os_["longmem"]["F1_sn_block"][0]); mac["LongBone"] = p3(os_["longmem"]["B1c_frozen_lrv_cusum"][0])
mac["TfiveFtwo"] = p3(os_["AR0.8_t5"]["F2_bpi_wlc"][0]); mac["SkewFtwo"] = p3(os_["AR0.8_skew"]["F2_bpi_wlc"][0])
e8 = P["e8"]; syn = [k for k in e8 if not k.startswith("AQ_")]
mac["EviolLoInit"] = p3(min(e8[k]["336"]["viol"][0] for k in syn)); mac["EviolHiInit"] = p3(max(e8[k]["336"]["viol"][0] for k in syn))
mac["EviolMaxLong"] = p3(max(e8[k]["10000"]["viol"][0] for k in syn))
mac["FaInvMaxInit"] = p3(max(e8[k]["336"]["fa_1overalpha"][0] for k in syn)); mac["FaInvNexcInit"] = str(sum(e8[k]["336"]["fa_1overalpha"][0] > 0.05 for k in syn))
mac["FaInvMaxThree"] = p3(max(e8[k]["3000"]["fa_1overalpha"][0] for k in e8))
mac["PredA"] = f"{P['pred_overhead']['0.05']:.2f}"; mac["PredC"] = f"{P['pred_overhead']['1e-06']:.2f}"
rA = [P["delay_ratio_by_alpha"][k]["F2_bpi_wlc"]["1e-06"] for k in P["insc"]]; mac["FtwoRatioCLo"] = f"{min(rA):.2f}"; mac["FtwoRatioCHi"] = f"{max(rA):.2f}"
mac["EthreeFtwoMax"] = p3(P["e3_F2_fa_max"])
mac["HthreeFaster"] = str(P["h3_f2_faster"][0]); mac["HthreeTotal"] = str(P["h3_f2_faster"][1]); mac["HiidA"] = f"{P['h_iid_blocks_a05']:.2f}"
mac["NumbersKeys"] = str(len(N))
mac["FixMaxPct"] = f"{100 * max(f1max, f2max):.1f}\,\%"
mac["IidMaxPct"] = f"{100 * max(t3['iid']['eprocess'][0], t3['iid']['ctm'][0]):.1f}\,\%"
mac["FasterLo"] = f"{100 * (1 - c['B1c_k0.5|0.05']['delay_ratio_range'][1]):.0f}"; mac["FasterHi"] = f"{100 * (1 - c['B1c_k0.5|0.05']['delay_ratio_range'][0]):.0f}"
with open(os.path.join(M, "numbers.tex"), "w", encoding="utf-8") as f:
    f.write("% generated by code/make_tex.py from results/numbers.json -- do not edit\n")
    for k, v in mac.items():
        f.write("\\newcommand{\\%s}{%s}\n" % (k, v))


def wt(name, body):
    open(os.path.join(M, "tables", name), "w", encoding="utf-8").write("% generated by code/make_tex.py -- do not edit\n" + body)


# Table 1 data
rows = []
for ch in CH:
    m = P["meta"][ch]
    rows.append(f"{CHL[ch]} & {m['sensor'].replace('PT08.', 'PT08.')} & {m['n_pool_h']} & {m['n_168h_nonoverlap']} & {m['acf1']:.2f} & {m['acf24']:.2f} & {m['sigma2_LR_b72']:.1f} / {m['sigma2_LR_b168']:.1f} & {100 * m['ref_missing_frac_full']:.1f} \\\\")
wt("t1_data.tex", "\\begin{tabular}{llrrrrrr}\n\\hline\nReference & Sensor & Null pool (h) & 168-h blocks & ACF(1 h) & ACF(24 h) & $\\hat\\sigma^2_{\\mathrm{LR}}$ (72 / 168 h) & Missing (\\%) \\\\\n\\hline\n" + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
# Table 2 RQ1
ML = [("M1_cusum", "CUSUM", "heuristic"), ("M2_pagehinkley", "Page--Hinkley", "heuristic"), ("M3_adwin", "ADWIN", "per step"), ("M4_repeated_p", "Repeated $p<0.05$", "fixed sample"),
      ("M5_eprocess", "Betting e-process", "all times"), ("M6_ctm", "Conformal test martingale", "all times"), ("F1_sn_block", "F1 self-normalised blocks", "all times$^{a}$"),
      ("F2_bpi_wlc", "F2 adaptive block plug-in", "window, asympt."), ("B1c_frozen_lrv_cusum", "Frozen-LRV CUSUM ($\\kappa=0.1$)", "asympt. ($m\\to\\infty$)"), ("B2c_ar_residual_cusum", "AR-residual CUSUM", "heuristic")]
rows = []
for m, lab, gu in ML:
    cells = [" / ".join(f"{P['t2'][ch][b][m][0]:.2f}" for b in BL) for ch in CH]
    rows.append(f"{lab} & {gu} & " + " & ".join(cells) + " \\\\")
wt("t2_rq1.tex", "\\begin{tabular}{llcccc}\n\\hline\nMonitor & Nominal guarantee & " + " & ".join(CHL[c] for c in CH) + " \\\\\n\\hline\n" + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
# Table 3 RQ2
t3 = P["t3"]
body = ("\\begin{tabular}{lcc}\n\\hline\nNull construction (CO channel, 720-h window) & Betting e-process & Conformal test martingale \\\\\n\\hline\n"
        f"Block bootstrap 24 h, raw residual & {ci(t3['raw_block24']['eprocess'])} & {ci(t3['raw_block24']['ctm'])} \\\\\n"
        f"Block bootstrap 24 h, daily cycle removed & {ci(t3['deseason_block24']['eprocess'])} & {ci(t3['deseason_block24']['ctm'])} \\\\\n"
        f"Independent resampling (control) & {ci(t3['iid']['eprocess'])} & {ci(t3['iid']['ctm'])} \\\\\n"
        f"Block means (24 h), long-run variance from calibration month & {ci(t3['blockmean24_eproc'])} & -- \\\\\n"
        f"Whitening with estimated AR(1) coefficient & {ci(t3['whiten_phihat_eproc'])} & -- \\\\\n\\hline\n\\end{{tabular}}\n")
wt("t3_rq2.tex", body)
# Table 4a FA corrections
cols = [("F1", "F1"), ("F2", "F2"), ("B1c", "$\\kappa=0.1$"), ("B1c_k0.25", "$\\kappa=0.25$"), ("B1c_k0.5", "$\\kappa=0.5$"), ("B1c_hcal", "raised $h$"), ("B2c", "AR-res.")]
rows = []
for k in P["insc"]:
    rows.append(NULL_L[k] + " & " + " & ".join(f"{P['t4_fa'][k][c][0]:.3f}" for c, _ in cols) + " \\\\")
for ch in CH:
    for b in BL:
        k = f"AQ_{ch}|boot{b}"
        rows.append(f"{CHL[ch]}, block {b} h & " + " & ".join(f"{P['t4_fa'][k][c][0]:.3f}" for c, _ in cols) + " \\\\")
wt("t4a_fa.tex", "\\begin{tabular}{lccccccc}\n\\hline\n & & & \\multicolumn{4}{c}{Frozen-LRV CUSUM} & \\\\\nNull / channel & " + " & ".join(l for _, l in cols) + " \\\\\n\\hline\n" + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
# Table 4b delay
dc = [("F1", "F1"), ("F2", "F2"), ("B1c", "$\\kappa=0.1$"), ("B1c_k0.5", "$\\kappa=0.5$"), ("B1c_hcal", "raised $h$")]
rows = []
for k in P["insc"] + [f"AQ_{c}|boot168" for c in CH]:
    lab = NULL_L[k] if k in NULL_L else CHL[k[3:].split("|")[0]] + ", block 168 h"
    rows.append(lab + " & " + " & ".join(f"{P['t4_delay'][k][c]['mean_delay']:.0f}" for c, _ in dc) + " \\\\")
wt("t4b_delay.tex", "\\begin{tabular}{lccccc}\n\\hline\n & & & \\multicolumn{3}{c}{Frozen-LRV CUSUM} \\\\\nNull / channel & " + " & ".join(l for _, l in dc) + " \\\\\n\\hline\n" + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
# Table 5 burn-in
n0s = ["336", "1000", "3000", "10000"]; rows = []
for k in P["e8"]:
    lab = NULL_L[k] if k in NULL_L else CHL.get(k[3:].split("|")[0], k) + " (bootstrap)"
    rows.append(lab + " & " + " & ".join(f"{P['e8'][k][n]['viol'][0]:.3f} / {P['e8'][k][n]['fa_1overalpha'][0]:.3f}" for n in n0s) + " \\\\")
wt("t5_burnin.tex", "\\begin{tabular}{lcccc}\n\\hline\nNull & $n_0=336$ & 1000 & 3000 & 10000 \\\\\n\\hline\n" + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
# Table 6 E9
if "e9" in P:
    rows = []
    for ch in CH:
        for cval in ("0.0", "0.05", "0.1", "0.2"):
            cells = []
            for m in ("M5_eprocess", "M6_ctm", "F2_bpi_wlc", "B1c_k0.5"):
                cells.append(" / ".join(f"{P['e9'][f'{ch}|c{cval}|boot{b}'][m][0]:.2f}" for b in BL))
            acf = P["e9"][f"{ch}|c{cval}|boot24"]["acf1"]
            rows.append(f"{CHL[ch]} & {cval} & {acf:.2f} & " + " & ".join(cells) + " \\\\")
    wt("t6_e9.tex", "\\begin{tabular}{lcccccc}\n\\hline\nChannel & $c$ & ACF(1 h) & Betting e-process & Conformal t.m. & F2 & Frozen-LRV, $\\kappa=0.5$ \\\\\n\\hline\n" + "\n".join(rows) + "\n\\hline\n\\end{tabular}\n")
    e9 = P["e9"]
    mac2 = {"EnineGuarMin": p3(min(e9[k][m][0] for k in e9 if "|c0.2|" in k for m in ("M5_eprocess", "M6_ctm"))),
            "EnineFtwoMax": p3(max(e9[k]["F2_bpi_wlc"][0] for k in e9)), "EnineBoneMax": p3(max(e9[k]["B1c_k0.5"][0] for k in e9)),
            "EnineBoneNexc": str(sum(e9[k]["B1c_k0.5"][0] > 0.05 for k in e9)), "EnineCells": str(len(e9)),
            "EnineBoneWorst": max(((e9[k]["B1c_k0.5"][0], k) for k in e9))[1].replace("(GT)", "").replace("|c", ", $c=").replace("|boot", "$, block ") + " h"}
    with open(os.path.join(M, "numbers.tex"), "a", encoding="utf-8") as f:
        for k, v in mac2.items(): f.write("\\newcommand{\\%s}{%s}\n" % (k, v))
print("tex ok", len(mac))
