"""H35 Q1ROUTE Q-0: false alarm of e-process and CTM on B1 null, raw vs diurnal cycle removed (no blocking).
Same pool/bootstrap as pilot_B / s1b_a_decompose (block-24 circular bootstrap, H=720, 5 seeds x 400).
Usage: python q0_deseason_ctm.py <AirQualityUCI.csv> <out.json>
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P


def main():
    df = P.load(sys.argv[1])
    s1 = df["PT08.S1(CO)"].to_numpy(float); co = df["CO(GT)"].to_numpy(float); idx = df.index
    march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
    roos = np.full(len(s1), np.nan)
    for d in np.unique(days[march]):
        te = march & (days == d); tr = march & (days != d)
        c = P.fit(s1[tr], co[tr]); roos[te] = P.pred(c, s1[te]) - co[te]
    m = march & np.isfinite(roos)
    z = roos[m]; h = idx.hour.to_numpy()[m]
    z = (z - z.mean()) / z.std()
    prof = np.array([z[h == k].mean() for k in range(24)])
    zds = z - prof[h]; zds = (zds - zds.mean()) / zds.std()
    out = {"config": {"alpha": P.ALPHA, "seeds": P.SEEDS, "H": P.H, "series_per_seed": P.BPS},
           "var_explained_by_hour_profile": round(float(1 - np.var(z - prof[h]) / np.var(z)), 3)}
    for name, pool in (("raw", z), ("deseason", zds)):
        for mode, blk in (("block24", 24), ("iid", 1)):
            he, hc = [], []
            for s in P.SEEDS:
                rng = np.random.default_rng(s)
                X = P.cbb(pool, P.BPS, P.H, blk, rng)
                he.append(P.mon_eprocess(X) >= 0); hc.append(P.mon_ctm(X, rng) >= 0)
            he = np.concatenate(he); hc = np.concatenate(hc)
            out[f"{name}|{mode}"] = {"eprocess": P.wilson(int(he.sum()), len(he)), "ctm": P.wilson(int(hc.sum()), len(hc))}
            print(name, mode, out[f"{name}|{mode}"], flush=True)
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
