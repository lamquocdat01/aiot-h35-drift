"""H35 S1-1b(a) robustness: block-mean / whitened e-process under longer bootstrap blocks (72, 168 h).
(Formerly run inline on 21-11; same code.)  Usage: python -I s1b_a_robust.py <AirQualityUCI.csv> <out.json>"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pilot_B as P, s1b_a_decompose as D
df = P.load(sys.argv[1])
s1 = df["PT08.S1(CO)"].to_numpy(float); co = df["CO(GT)"].to_numpy(float); idx = df.index
march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
roos = np.full(len(s1), np.nan)
for d in np.unique(days[march]):
    te = march & (days == d); tr = march & (days != d); c = P.fit(s1[tr], co[tr]); roos[te] = P.pred(c, s1[te]) - co[te]
z = roos[march & np.isfinite(roos)]; z = (z - z.mean()) / z.std()
out = {}
for bb in (72, 168):
    for b in (24, 48):
        sdb = np.array([z[i:i + b].mean() for i in range(0, len(z) - b + 1, b)]).std(ddof=1)
        fr = []; fw = []
        for s in P.SEEDS:
            rng = np.random.default_rng(s); X = P.cbb(z, P.BPS, P.H, bb, rng)
            fr.append(D.eproc_first(D.block_means(X, b, sdb)))
            if b == 24: fw.append(D.eproc_first(D.whiten_online(X), burn=48))
        out[f"bootblk{bb}_blockmean{b}"] = [D.rate(np.concatenate(fr)), round(float(sdb), 3), len(z) // b]
        if fw: out[f"bootblk{bb}_whiten"] = D.rate(np.concatenate(fw))
json.dump({"robustness_boot_block": out}, open(sys.argv[2], "w", encoding="utf-8"), indent=1)
print(out)
