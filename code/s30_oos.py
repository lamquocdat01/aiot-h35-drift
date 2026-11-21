"""H35 S3-0: out-of-sample check of B1c kappa=0.5 vs F2 (Amendment 7). Usage: python s30_oos.py <out.json>"""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P, q2_adaptive as Q, c1_data as D, c1_run as C, c1_e78 as E

SEEDS = [20261300 + k for k in range(5)]
NEW = {"AR2_0.5_0.3": dict(phi=(0.5, 0.3)), "ARMA22_1.0_-0.3_|_0.4_0.2": dict(phi=(1.0, -0.3), theta=(0.4, 0.2))}
N0, W, A = C.N0, C.W, 0.05


def ci(d):
    return [round(float(d.mean() - 1.96 * d.std(ddof=1) / np.sqrt(len(d))), 1), round(float(d.mean() + 1.96 * d.std(ddof=1) / np.sqrt(len(d))), 1)]


def run_case(gen, sig2):
    h = C.TH[A]["cusum_blocks"]; fa = {"B1c_k0.5": [], "F2_bpi_wlc": []}
    for s in SEEDS:
        rng = np.random.default_rng(s); X = gen(rng, 400, N0 + W)
        M, ends = E.b1c_path(X, 0.5); f = E.first_alarm(M, ends, h); fa["B1c_k0.5"].append((f >= 0) & (f <= W))
        f2 = Q.bpi_wlc(X, [A])[A]; fa["F2_bpi_wlc"].append((f2 >= 0) & (f2 <= W))
    S = np.sqrt(2 * sig2 * np.log(1 / A) / 60.0)
    dl = {"B1c_k0.5": [], "F2_bpi_wlc": []}
    for s in SEEDS:
        rng = np.random.default_rng(s + 1000); X = gen(rng, 200, N0 + 240 + 3000); X[:, N0 + 240:] += S
        M, ends = E.b1c_path(X, 0.5); dl["B1c_k0.5"].append(E.first_alarm(M, ends, h)); dl["F2_bpi_wlc"].append(Q.bpi_wlc(X, [A])[A])
    out = {"S": round(float(S), 4), "sigma2_LR": round(float(sig2), 3)}
    for k in fa:
        f = np.concatenate(dl[k]); d = (f[f >= 240] - 240).astype(float)
        out[k] = {"FA": P.wilson(int(np.concatenate(fa[k]).sum()), len(np.concatenate(fa[k]))),
                  "mean_delay": round(float(d.mean()), 1), "ci95": ci(d), "pre_change": round(float(((f >= 0) & (f < 240)).mean()), 4),
                  "miss": round(float((f < 0).mean()), 4), "ratio_to_Dstar": round(float(d.mean() / 60.0), 3)}
    b, f2 = out["B1c_k0.5"], out["F2_bpi_wlc"]
    out["B1c_keeps_alpha"] = b["FA"][0] <= A
    out["B1c_faster_nonoverlap"] = b["ci95"][1] < f2["ci95"][0]
    return out


def main():
    t0 = time.time(); C.thresholds([A]); res = {"seeds": SEEDS, "h_iid": C.TH[A]["cusum_blocks"]}
    for name, kw in NEW.items():
        sig2 = C.psi_lrv_vb(**kw)[0]
        res[name] = run_case(lambda rng, R, n, kw=kw: Q.gen_arma(rng, R, n, **kw), sig2); print(name, res[name], flush=True)
    aq = D.aq_channels()
    for c, v in aq.items():
        pool = v["pool"]; sig2 = float(np.mean([pool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(pool) - 167, 24)]))
        for bb in C.BOOTS:
            res[f"AQ_{c}|boot{bb}"] = run_case(lambda rng, R, n, bb=bb, pool=pool: P.cbb(pool, R, n, bb, rng), sig2)
            print(c, bb, {k: res[f"AQ_{c}|boot{bb}"][k] for k in ("B1c_keeps_alpha", "B1c_faster_nonoverlap")}, res[f"AQ_{c}|boot{bb}"]["B1c_k0.5"]["FA"], flush=True)
    res["verdict_primary"] = "XAC NHAN NGOAI MAU" if all(res[n]["B1c_keeps_alpha"] and res[n]["B1c_faster_nonoverlap"] for n in NEW) else "VO - kappa can hieu chuan theo du lieu"
    res["runtime_s"] = round(time.time() - t0, 1)
    json.dump(res, open(sys.argv[1], "w", encoding="utf-8"), indent=1); print(res["verdict_primary"])


if __name__ == "__main__":
    main()
