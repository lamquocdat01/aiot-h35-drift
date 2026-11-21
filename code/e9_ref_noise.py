"""H35 E9 (Amendment 8): sensitivity of the false-alarm audit to reference-analyzer noise.
White Gaussian noise sd = c * sigma_res is added to the reference y during calibration month, residual pool rebuilt
(leave-one-day-out), then false alarm in W = 720 h for M3, M5, M6, F1, F2, B1c(kappa=0.5).
Usage: python e9_ref_noise.py <out.json>"""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P, q2_adaptive as Q, c1_data as D, c1_run as C, c1_e78 as E

SEEDS = C.SEEDS; N0, W, A = C.N0, C.W, 0.05
CS = (0.0, 0.05, 0.1, 0.2)


def pool_with_noise(df, sx, sy, c, rng):
    idx = df.index; march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
    x = df[sx].to_numpy(float); y = df[sy].to_numpy(float)

    def oos(yv):
        r = np.full(len(x), np.nan)
        for d in np.unique(days[march]):
            te = march & (days == d); tr = march & (days != d)
            cf = D._fit(x[tr], yv[tr])
            if cf is not None: r[te] = cf[0] + cf[1] * x[te] - yv[te]
        return r[march & np.isfinite(r)]
    base = oos(y); sres = base.std()
    yn = y + (rng.standard_normal(len(y)) * c * sres if c > 0 else 0.0)
    p = oos(yn); p = (p - p.mean()) / p.std()
    acf1 = float(np.corrcoef(p[:-1], p[1:])[0, 1])
    return p, acf1


def main():
    t0 = time.time(); C.thresholds([A]); df = D.load_aq(); res = {"c": CS}
    h = C.TH[A]["cusum_blocks"]
    for ci_, (sx, sy) in enumerate(D.AQ_PAIRS):
        for k, c in enumerate(CS):
            rng = np.random.default_rng(20261400 + 10 * ci_ + k)
            pool, acf1 = pool_with_noise(df, sx, sy, c, rng)
            for bb in C.BOOTS:
                acc = {m: [] for m in ("M3_adwin", "M5_eprocess", "M6_ctm", "F1_sn_block", "F2_bpi_wlc", "B1c_k0.5")}
                for s in SEEDS:
                    r2 = np.random.default_rng(s); X = P.cbb(pool, 200, N0 + W, bb, r2)
                    o = C.m_monitors(X, r2, which=("M3", "M5", "M6"))
                    o["F1_sn_block"] = Q.sn_block(X, [A])[A]; o["F2_bpi_wlc"] = Q.bpi_wlc(X, [A])[A]
                    M, ends = E.b1c_path(X, 0.5); o["B1c_k0.5"] = E.first_alarm(M, ends, h)
                    for m, f in o.items(): acc[m].append((f >= 0) & (f <= W))
                key = f"{sy}|c{c}|boot{bb}"
                res[key] = {"acf1": round(acf1, 3), **{m: P.wilson(int(np.concatenate(v).sum()), len(np.concatenate(v))) for m, v in acc.items()}}
                print(key, {m: res[key][m][0] for m in acc}, "acf1", res[key]["acf1"], flush=True)
    res["runtime_s"] = round(time.time() - t0, 1)
    json.dump(res, open(sys.argv[1], "w", encoding="utf-8"), indent=1)


if __name__ == "__main__":
    main()
