"""Build results/numbers.json from the S1 result files.  Usage: python -I build_numbers.py <results_dir> <out.json>"""
import json, sys, os
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
n["branch"] = "Q1-DL co dieu kien (E) + Gauss; S2 quyet" if "Q2_bpi_FA_AR0.95" in n else n["branch"]
json.dump(n, open(sys.argv[2], "w", encoding="utf-8"), indent=1)
print(len(n), "keys")
