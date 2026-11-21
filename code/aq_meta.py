"""H35: metadata of the 4 referenced UCI Air Quality channels for Table 1. Usage: python aq_meta.py <out.json>"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import c1_data as D

aq = D.aq_channels(); df = D.load_aq(); out = {}
for (sx, sy) in D.AQ_PAIRS:
    v = aq[sy]; p = v["pool"]
    acf = lambda k: round(float(np.corrcoef(p[:-k], p[k:])[0, 1]), 3)
    s2 = float(np.mean([p[i:i + 168].sum() ** 2 / 168 for i in range(0, len(p) - 167, 24)]))
    s2_72 = float(np.mean([p[i:i + 72].sum() ** 2 / 72 for i in range(0, len(p) - 71, 24)]))
    out[sy] = {"sensor": sx, "n_pool_h": int(len(p)), "n_168h_nonoverlap": int(len(p) // 168), "acf1": acf(1), "acf24": acf(24),
               "sigma2_LR_b168": round(s2, 2), "sigma2_LR_b72": round(s2_72, 2), "res_sd_units": round(v["sd"], 3),
               "ref_missing_frac_full": round(float(df[sy].isna().mean()), 3), "slope": round(v["coef"][1], 6), "intercept": round(v["coef"][0], 3)}
json.dump(out, open(sys.argv[1], "w", encoding="utf-8"), indent=1); print(out)
