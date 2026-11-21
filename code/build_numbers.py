"""Build results/numbers.json from the S1 result files.  Usage: python -I build_numbers.py <results_dir> <out.json>"""
import json, sys, os
import numpy as np
R = sys.argv[1]
L = lambda f: json.load(open(os.path.join(R, f), encoding="utf-8"))
A = L("pilot_A.json"); B = L("pilot_B.json"); D = L("s1b_a_decompose.json")
if "robustness_boot_block" not in D:
    D.update(L("s1b_a_robust.json"))
a = L("s1b_thm1_check.json"); b = L("s1b_thm2_sim.json"); c = L("s1b_prop0_cex.json")
r = A["rules"]; e = D["eprocess_false_alarm"]; u = a["uci_b1_selfnormalized"]
n = {"_note": "S1 numbers, 21-11-2026; [mean, lo95, hi95]", "alpha": 0.05,
 "A_rule_i_N16_tau_short_lifetimeFDR": r["i"]["tau_lifetimeFDR"]["short"],
 "A_rule_i_N16_tau_first_lifetimeFDR": r["i"]["tau_lifetimeFDR"]["first"],
 "A_rule_i_N1_tau_short_lifetimeFDR": A["rule_i_N1"]["tau_lifetimeFDR"]["short"],
 "A_rule_ii_tau_short_lifetimeFDR": r["ii"]["tau_lifetimeFDR"]["short"],
 "A_rule_ii_tau_short_currentFDR": r["ii"]["tau_currentFDR"]["short"],
 "A_rule_iii_tau_max_lifetimeFDR": max(r["iii"]["tau_lifetimeFDR"].values()),
 "A_rule_iii_delay_j1_j10": [r["iii"]["delay_by_epoch"]["1"], r["iii"]["delay_by_epoch"]["10"]],
 "A_rule_i_delay_j1_j10": [r["i"]["delay_by_epoch"]["1"], r["i"]["delay_by_epoch"]["10"]],
 "A_cexA_lifetimeFDP_tau": A["counterexample_A"][0]["lifetime_FDP_at_tau_exact"], "A_cexB_EZtau": 1.5,
 "B_null_block24_eprocess": B["part1_null_audit"]["block24"]["eprocess"]["rate_wilson95"],
 "B_null_block24_ctm": B["part1_null_audit"]["block24"]["ctm"]["rate_wilson95"],
 "B_null_iid_eprocess": B["part1_null_audit"]["iid"]["eprocess"]["rate_wilson95"],
 "B_null_iid_ctm": B["part1_null_audit"]["iid"]["ctm"]["rate_wilson95"],
 "B_bias_cusum_pre_signed": B["part2_replay"]["cusum_pre"]["signed_bias_seeds"][0],
 "B_dMAE_eproc_pre_vs_sched": B["part2_replay"]["eproc_pre"]["dMAE_vs_sched_seeds"][0],
 "B_dMAE_eproc_post_vs_sched": B["part2_replay"]["eproc_post"]["dMAE_vs_sched_seeds"][0],
 "S1b_a_eproc_deseason": e["1_deseason"]["rate_wilson95"], "S1b_a_eproc_block24": e["2_block24_raw"]["rate_wilson95"],
 "S1b_a_eproc_whiten_phihat": e["3_whiten_phihat_raw"]["rate_wilson95"],
 "S1b_a_block24_bootblk168": D["robustness_boot_block"]["bootblk168_blockmean24"][0],
 "S1b_a_block48_bootblk168": D["robustness_boot_block"]["bootblk168_blockmean48"][0],
 "T1_lemma1_grid_max": a["lemma1_grid_two_point"]["max"],
 "T1_uci_b24g24_boot24_72_168": [u["bootblk24|b24|g24"], u["bootblk72|b24|g24"], u["bootblk168|b24|g24"]],
 "T1_uci_b24g0_boot168": u["bootblk168|b24|g0"],
 "T1_cex_longmem_g24": a["synthetic"]["longmem|b24|g24|d0=0|shift=0"]["rate_wilson95"],
 "T1_cex_skew_b1": a["synthetic"]["skew_iid|b1|g0|d0=0|shift=0"]["rate_wilson95"],
 "T2_N16_iv_EsupFDP": b["N16_iv"]["E_sup_FDP"], "T2_N16_iii_EsupFDP": b["N16_iii"]["E_sup_FDP"],
 "T2_N16_ivb_dynamic_EsupFDP": b["N16_iv_nonpred"]["E_sup_FDP"],
 "T2_N16_delay_j10_iii_iv": [b["N16_iii"]["delay_by_epoch"]["10"], b["N16_iv"]["delay_by_epoch"]["10"]],
 "T2_N4_delay_j10_iii_iv": [b["N4_iii"]["delay_by_epoch"]["10"], b["N4_iv"]["delay_by_epoch"]["10"]],
 "T2_cex_const_gamma_N16": b["N16_iv_const"]["E_sup_FDP"], "T2_cex_raw_N16": b["N16_iv_raw"]["E_sup_FDP"],
 "P0_cex_harmonic_N16": c["N16_iii_harm"]["E_sup_FDP"], "P0_cex_raw_N16": c["N16_iii_raw"]["E_sup_FDP"],
 "branch": "TSP giu: (c),(d) dong; rui ro 'gia tang' cao -> S2 quyet"}
# --- Q1ROUTE additions (optional files)
if os.path.exists(os.path.join(R, "q0_deseason_ctm.json")):
    q0 = L("q0_deseason_ctm.json")
    n["Q0_deseason_eprocess"] = q0["deseason|block24"]["eprocess"]; n["Q0_deseason_ctm"] = q0["deseason|block24"]["ctm"]
    n["Q0_var_explained_hour_profile"] = q0["var_explained_by_hour_profile"]
if os.path.exists(os.path.join(R, "q2_adaptive.json")):
    q2 = L("q2_adaptive.json"); fa = q2["false_alarm_720h"]; dl = q2["delay"]
    n["Q2_bpi_FA_max_over_nulls"] = max(v["bpi_wlc"]["rate_wilson95"][0] for v in fa.values())
    n["Q2_bpi_FA_AR0.95"] = fa["AR1_0.95"]["bpi_wlc"]["rate_wilson95"]
    n["Q2_bpi_FA_B1boot168"] = fa["B1_boot168"]["bpi_wlc"]["rate_wilson95"]
    n["Q2_bpi_delay_ratio_AR0.8_by_alpha"] = {a: v["ratio_to_Dstar"] for a, v in dl["AR1_0.8"]["bpi_wlc"].items()}
    n["Q2_sn_delay_ratio_AR0.8_by_alpha"] = {a: v["ratio_to_Dstar"] for a, v in dl["AR1_0.8"]["sn_block_thm1"].items()}
    n["Q2_predicted_overhead_by_alpha"] = {a: v["ratio_predicted_overhead"] for a, v in dl["AR1_0.8"]["bpi_wlc"].items()}
if os.path.exists(os.path.join(R, "q2_selfcheck.json")):
    sc = L("q2_selfcheck.json")
    n["Q2_cex_noUnion_AR0.95"] = sc["a1_A=1/alpha|AR0.95"]; n["Q2_cex_longmem"] = sc["a4_longmem"]
    n["Q2_cex_nonlagged_delay"] = sc["b1_nonlagged_estimator"]
# --- 1c additions
if os.path.exists(os.path.join(R, "c1_e1_syn.json")):
    C = lambda f: L(f)["result"]
    s1, aq, d2, a0, e4, e5, e3 = C("c1_e1_syn.json"), C("c1_e1_aq.json"), C("c1_e2_delay.json"), C("c1_e2_arl0.json"), C("c1_e4.json"), C("c1_e5.json"), C("c1_e3.json")
    insc = ["AR1_0", "AR1_0.5", "AR1_0.8", "AR1_0.95", "ARMA11_0.8_+0.4", "ARMA11_0.8_-0.5", "ARMA21_1.2_-0.4_+0.3"]
    n["C1_H2_F2_FA_max_inscope"] = max(s1[k]["F2_bpi_wlc"]["rate_wilson95"][0] for k in insc)
    n["C1_H2_F2_FA_max_AQ_boot168"] = max(aq[k]["F2_bpi_wlc"]["rate_wilson95"][0] for k in aq if k.endswith("boot168"))
    n["C1_H1_guaranteed_FA_AQ_boot168_min"] = min(aq[k][m]["rate_wilson95"][1] for k in aq if k.endswith("boot168") for m in ("M3_adwin", "M5_eprocess", "M6_ctm"))
    n["C1_H3_F2_ratio_1e-6_inscope_max"] = max(d2[k]["cond_delay_nu240"]["F2_bpi_wlc"]["1e-06"]["ratio_to_Dstar"] for k in insc)
    n["C1_H3_pred_overhead_1e-6"] = d2["AR1_0"]["pred_overhead"]["1e-06"]
    n["C1_B1c_FA_max_inscope"] = max(s1[k]["B1c_frozen_lrv_cusum"]["rate_wilson95"][0] for k in insc)
    n["C1_B2c_FA_max_AQ_boot168"] = max(aq[k]["B2c_ar_residual_cusum"]["rate_wilson95"][0] for k in aq if k.endswith("boot168"))
    n["C1_delay_ratio_AR0.8_alpha0.05_F2_F1_B1c"] = [d2["AR1_0.8"]["cond_delay_nu240"][m]["0.05"]["ratio_to_Dstar"] for m in ("F2_bpi_wlc", "F1_sn_block", "B1c_frozen_lrv_cusum")]
    n["C1_ARL0_AR0.8_F2_B1c"] = [a0["AR1_0.8"]["F2_bpi_wlc"]["ARL0_lower_bound_mean"], a0["AR1_0.8"]["B1c_frozen_lrv_cusum"]["ARL0_lower_bound_mean"]]
    n["C1_H6_eps_max_p95_AR0.8"] = e4["AR1_0.8"]["eps_max_p95"]
    n["C1_H6_P_under_vb_AR0.8"] = e4["AR1_0.8"]["P_any_under_vb"]
    n["C1_H7_R3a_current_FDR"] = e5["R3a_eBH_reset"]["current_FDR_at_tau_short"]
    n["C1_out_longmem_F2_FA"] = s1["longmem"]["F2_bpi_wlc"]["rate_wilson95"]
    n["C1_E3_F2_FA_max_grid"] = max(v["AR1_0.8"]["FA"][0] for k, v in e3.items() if k.startswith("F2"))
    n["C1_verdict"] = "H1 PASS, H2 PASS, H3 PASS, H6 FAIL, H7 FAIL"
# --- E7/E8/S3-0 additions
if os.path.exists(os.path.join(R, "c1_e7.json")):
    e7 = L("c1_e7.json")["result"]; e8 = L("c1_e8.json")["result"]
    n["E7_hcal"] = {a: e7["calib_" + a]["h"] for a in ("0.05", "0.001")}
    n["E7_F2_FA_max_alpha1e-3"] = max(v["0.001"]["F2_bpi_wlc"][0] for v in e7["false_alarm"].values())
    n["E7_B1c_k0.5_delay_over_F2_AR0.8_a0.05"] = round(e7["delay"]["AR1_0.8"]["cond_delay_nu240"]["B1c_k0.5"]["0.05"]["mean_delay"] / e7["delay"]["AR1_0.8"]["cond_delay_nu240"]["F2_bpi_wlc"]["0.05"]["mean_delay"], 3)
    n["E8_Elow_violation_AR0.8_by_n0"] = {k: v["P_lower_violation"][0] for k, v in e8["AR1_0.8"]["by_n0"].items()}
    n["E8_F2_FA_1overalpha_AR0.8_by_n0"] = {k: v["F2_FA_A=1/alpha"][0] for k, v in e8["AR1_0.8"]["by_n0"].items()}
if os.path.exists(os.path.join(R, "s30_oos.json")):
    s3 = L("s30_oos.json")
    for k in ("AR2_0.5_0.3", "ARMA22_1.0_-0.3_|_0.4_0.2"):
        n["S30_" + k.split("_")[0] + "_B1c_k0.5_FA"] = s3[k]["B1c_k0.5"]["FA"]
        n["S30_" + k.split("_")[0] + "_delay_B1c_F2"] = [s3[k]["B1c_k0.5"]["mean_delay"], s3[k]["F2_bpi_wlc"]["mean_delay"]]
    n["S30_verdict"] = s3["verdict_primary"]
n["branch"] = "Q2 - MST (S2 chot 21-11)"
# --- paper block (F2): every number used in the manuscript tables/text
def _paper():
    P = {}
    if not os.path.exists(os.path.join(R, "c1_e1_aq.json")):
        return None
    C1 = lambda f: L(f)["result"]
    aq, syn, d2, a0, e7, e8, e3 = C1("c1_e1_aq.json"), C1("c1_e1_syn.json"), C1("c1_e2_delay.json"), C1("c1_e2_arl0.json"), C1("c1_e7.json"), C1("c1_e8.json"), C1("c1_e3.json")
    intel = C1("c1_e1_intel.json"); s3 = L("s30_oos.json"); q0 = L("q0_deseason_ctm.json"); sa = L("s1b_a_decompose.json")
    P["meta"] = L("aq_meta.json")
    chans = ["CO(GT)", "C6H6(GT)", "NOx(GT)", "NO2(GT)"]; blocks = ["24", "72", "168"]
    mons = ["M1_cusum", "M2_pagehinkley", "M3_adwin", "M4_repeated_p", "M5_eprocess", "M6_ctm", "F1_sn_block", "F2_bpi_wlc", "B1c_frozen_lrv_cusum", "B2c_ar_residual_cusum"]
    P["t2"] = {c: {b: {m: aq[f"AQ_{c}|boot{b}"][m]["rate_wilson95"] for m in mons} for b in blocks} for c in chans}
    g = [P["t2"][c][b][m] for c in chans for b in blocks for m in ("M5_eprocess", "M6_ctm")]
    P["rq1_guaranteed_range"] = [min(x[0] for x in g), max(x[0] for x in g)]
    P["rq1_guaranteed_min_lowerCI"] = min(x[1] for x in g)
    P["rq1_by_block"] = {b: [min(P["t2"][c][b][m][0] for c in chans for m in ("M5_eprocess", "M6_ctm")),
                             max(P["t2"][c][b][m][0] for c in chans for m in ("M5_eprocess", "M6_ctm"))] for b in blocks}
    P["rq1_adwin_range"] = [min(P["t2"][c][b]["M3_adwin"][0] for c in chans for b in blocks), max(P["t2"][c][b]["M3_adwin"][0] for c in chans for b in blocks)]
    ib = [k for k in intel if k.endswith("boot168")]
    P["intel_guaranteed_frac_lowerCI_gt_alpha"] = round(float(np.mean([min(intel[k]["M5_eprocess"]["rate_wilson95"][1], intel[k]["M6_ctm"]["rate_wilson95"][1]) > 0.05 for k in ib])), 3)
    P["intel_F2_n_exceed"] = int(sum(intel[k]["F2_bpi_wlc"]["rate_wilson95"][1] > 0.05 for k in ib)); P["intel_n"] = len(ib)
    P["t3"] = {"raw_block24": q0["raw|block24"], "deseason_block24": q0["deseason|block24"], "iid": q0["raw|iid"],
               "var_hour_profile": q0["var_explained_by_hour_profile"], "blockmean24_eproc": sa["eprocess_false_alarm"]["2_block24_raw"]["rate_wilson95"],
               "whiten_phihat_eproc": sa["eprocess_false_alarm"]["3_whiten_phihat_raw"]["rate_wilson95"],
               "adwin_iid": L("pilot_B.json")["part1_null_audit"]["iid"]["adwin"]["rate_wilson95"]}
    insc = ["AR1_0", "AR1_0.5", "AR1_0.8", "AR1_0.95", "ARMA11_0.8_+0.4", "ARMA11_0.8_-0.5", "ARMA21_1.2_-0.4_+0.3"]
    P["insc"] = insc
    fa = e7["false_alarm"]

    def row(k, src):
        return {"F1": src[k]["F1_sn_block"]["rate_wilson95"], "F2": fa[k]["0.05"]["F2_bpi_wlc"], "B1c": fa[k]["0.05"]["B1c_k0.1"],
                "B1c_k0.25": fa[k]["0.05"]["B1c_k0.25"], "B1c_k0.5": fa[k]["0.05"]["B1c_k0.5"], "B1c_hcal": fa[k]["0.05"]["B1c_k0.1_hcal"],
                "B2c": src[k]["B2c_ar_residual_cusum"]["rate_wilson95"]}
    P["t4_fa"] = {k: row(k, syn) for k in insc}
    for c in chans:
        for b in blocks:
            P["t4_fa"][f"AQ_{c}|boot{b}"] = row(f"AQ_{c}|boot{b}", aq)
    dl = e7["delay"]
    P["t4_delay"] = {k: {"F1": d2[k]["cond_delay_nu240"]["F1_sn_block"]["0.05"], "F2": dl[k]["cond_delay_nu240"]["F2_bpi_wlc"]["0.05"],
                         "B1c": dl[k]["cond_delay_nu240"]["B1c_k0.1"]["0.05"], "B1c_k0.5": dl[k]["cond_delay_nu240"]["B1c_k0.5"]["0.05"],
                         "B1c_hcal": dl[k]["cond_delay_nu240"]["B1c_k0.1_hcal"]["0.05"]} for k in insc + [f"AQ_{c}|boot168" for c in chans]}
    P["t4_arl0"] = {k: {m: a0[k][m] for m in ("F1_sn_block", "F2_bpi_wlc", "B1c_frozen_lrv_cusum", "B2c_ar_residual_cusum", "M5_eprocess")} for k in ["AR1_0.8", "AQ_CO(GT)|boot168"]}
    P["oos"] = {k: {m: s3[k][m] for m in ("B1c_k0.5", "F2_bpi_wlc")} for k in ("AR2_0.5_0.3", "ARMA22_1.0_-0.3_|_0.4_0.2")}
    P["oos_uci_B1c_FA_max"] = max(s3[k]["B1c_k0.5"]["FA"][0] for k in s3 if isinstance(s3[k], dict) and k.startswith("AQ_"))
    P["e7_h"] = {a: e7["calib_" + a]["h"] for a in ("0.05", "0.001")}

    def keep(v, a, fk):
        return all(fa[x][a][v][0] <= float(a) for x in fk)
    units = [(k, [k], k) for k in insc] + [(f"AQ_{c}", [f"AQ_{c}|boot{b}" for b in blocks], f"AQ_{c}|boot168") for c in chans]
    cnt = {}
    for a in ("0.05", "0.001"):
        for v in ("F2_bpi_wlc", "B1c_k0.1", "B1c_k0.25", "B1c_k0.5", "B1c_k0.1_hcal"):
            kp = sum(keep(v, a, fk) for _, fk, _ in units)
            fs = sum(keep(v, a, fk) and v != "F2_bpi_wlc" and dl[d]["cond_delay_nu240"][v][a]["ci95"][1] < dl[d]["cond_delay_nu240"]["F2_bpi_wlc"][a]["ci95"][0] for _, fk, d in units)
            rr = [dl[d]["cond_delay_nu240"][v][a]["mean_delay"] / dl[d]["cond_delay_nu240"]["F2_bpi_wlc"][a]["mean_delay"] for _, _, d in units]
            cnt[f"{v}|{a}"] = {"keep": kp, "keep_and_faster": fs, "delay_ratio_range": [round(min(rr), 2), round(max(rr), 2)]}
    P["e7_counts"] = cnt
    P["f2_fa_1e-3"] = {k: fa[k]["0.001"]["F2_bpi_wlc"] for k in fa}
    P["f2_fa_1e-3_n_exceed"] = int(sum(v[1] > 0.001 for v in P["f2_fa_1e-3"].values()))
    P["out_of_scope"] = {k: {m: syn[k][m]["rate_wilson95"] for m in ("F1_sn_block", "F2_bpi_wlc", "B1c_frozen_lrv_cusum", "M5_eprocess")} for k in ("AR0.8_t5", "AR0.8_skew", "longmem")}
    P["e8"] = {k: {n0: {"viol": v["P_lower_violation"], "fa_1overalpha": v["F2_FA_A=1/alpha"], "eps95": v["eps_max_p95"]} for n0, v in e8[k]["by_n0"].items()} for k in e8}
    P["delay_ratio_by_alpha"] = {k: {m: {a: d2[k]["cond_delay_nu240"][m][a]["ratio_to_Dstar"] for a in ("0.05", "0.001", "1e-06")} for m in ("F1_sn_block", "F2_bpi_wlc")} for k in insc}
    P["pred_overhead"] = d2["AR1_0"]["pred_overhead"]
    P["e3_F2_fa_max"] = max(v["AR1_0.8"]["FA"][0] for k, v in e3.items() if k.startswith("F2"))
    keys11 = insc + [f"AQ_{c}|boot168" for c in chans]
    cmp = [d2[k]["cond_delay_nu240"]["F2_bpi_wlc"][a]["mean_delay"] < d2[k]["cond_delay_nu240"]["F1_sn_block"][a]["mean_delay"] for k in keys11 for a in ("0.05", "0.001", "1e-06")]
    P["h3_f2_faster"] = [int(sum(cmp)), len(cmp)]
    P["h_iid_blocks_a05"] = round(L("c1_e7.json")["result"]["h_iid"]["0.05"], 2)
    if os.path.exists(os.path.join(R, "e9_ref_noise.json")):
        e9 = L("e9_ref_noise.json")
        P["e9"] = {k: v for k, v in e9.items() if "|" in k}
    return P


_pp = _paper()
if _pp is not None:
    n["paper"] = _pp

json.dump(n, open(sys.argv[2], "w", encoding="utf-8"), indent=1)
print(len(n), "keys")
