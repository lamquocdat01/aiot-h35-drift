"""H35 E7 (B1c calibrated: inflation kappa, raised threshold) and E8 (burn-in sweep) per Amendment 6.
Usage: python c1_e78.py e7|e8 <out.json>"""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P, q2_adaptive as Q, c1_data as D, c1_run as C

SEEDS = C.SEEDS; CAL_SEEDS = [20261200 + k for k in range(5)]
N0, W = C.N0, C.W
ALPHAS = [0.05, 1e-3]
INSC = list(C.SYN_IN)


def b1c_path(X, kappa, n0=N0, b=Q.B, g=Q.G):
    """Page CUSUM path on frozen-LRV standardized blocks; returns (M (R,L), block-end hours)."""
    R, n = X.shape
    st = Q.s2_tilde(C.autocov0(X[:, :n0]), n0, b) / (1 + Q.KAPPA) * (1 + kappa)
    step = b + g; L = (n - n0 - b) // step + 1
    sp = np.zeros(R); sm = np.zeros(R); M = np.zeros((R, L))
    for l in range(L):
        s = n0 + l * step; z = X[:, s:s + b].sum(1) / np.sqrt(b * st)
        sp = np.maximum(0, sp + z - 0.5); sm = np.maximum(0, sm - z - 0.5); M[:, l] = np.maximum(sp, sm)
    return M, np.arange(L) * step + b


def first_alarm(M, ends, h):
    hit = M >= h
    return np.where(hit.any(1), ends[hit.argmax(1)], -1)


def gen_null(name, rng, R, n, pool=None, blk=168):
    return P.cbb(pool, R, n, blk, rng) if pool is not None else Q.gen_arma(rng, R, n, **C.SYN_IN[name])


def wil(h): return P.wilson(int(h.sum()), len(h))


def e7():
    C.thresholds(ALPHAS); aq = D.aq_channels()
    h_iid = {a: C.TH[a]["cusum_blocks"] for a in ALPHAS}
    out = {"h_iid": h_iid}
    # ---- (b) calibrate raised threshold on CAL seeds, worst in-scope null
    nrep = {0.05: 400, 1e-3: 2000}
    h_cal = {}
    for a in ALPHAS:
        maxes = {}
        for name in INSC:
            ms = []
            for s in CAL_SEEDS:
                rng = np.random.default_rng(s); X = gen_null(name, rng, nrep[a], N0 + W)
                M, _ = b1c_path(X, 0.1); ms.append(M.max(1))
            maxes[name] = np.concatenate(ms)
        h = h_iid[a]
        while max(float((v >= h).mean()) for v in maxes.values()) > a: h = round(h + 0.05, 2)
        h_cal[a] = h
        out[f"calib_{a}"] = {"h": h, "worst_null": max(maxes, key=lambda k: float((maxes[k] >= h_iid[a]).mean())),
                             "cal_FA_at_h": {k: round(float((v >= h).mean()), 5) for k, v in maxes.items()}}
        print("CAL", a, out[f"calib_{a}"], flush=True)
    variants = {"B1c_k0.1": (0.1, h_iid), "B1c_k0.25": (0.25, h_iid), "B1c_k0.5": (0.5, h_iid), "B1c_k0.1_hcal": (0.1, h_cal)}
    out["variants"] = {k: {"kappa": v[0], "h": {str(a): v[1][a] for a in ALPHAS}} for k, v in variants.items()}
    # ---- false alarm (W window)
    fa_rep = {0.05: 400, 1e-3: 2000}
    FA = {}
    cases = [(k, None, None) for k in INSC] + [(f"AQ_{c}|boot{bb}", v["pool"], bb) for c, v in aq.items() for bb in C.BOOTS]
    for name, pool, bb in cases:
        rec = {}
        for a in ALPHAS:
            hs = {k: [] for k in list(variants) + ["F2_bpi_wlc"]}
            for s in SEEDS:
                rng = np.random.default_rng(s)
                X = gen_null(name, rng, fa_rep[a], N0 + W, pool, bb)
                for k, (kap, hh) in variants.items():
                    M, ends = b1c_path(X, kap); f = first_alarm(M, ends, hh[a]); hs[k].append((f >= 0) & (f <= W))
                f2 = Q.bpi_wlc(X, [a])[a]; hs["F2_bpi_wlc"].append((f2 >= 0) & (f2 <= W))
            rec[str(a)] = {k: wil(np.concatenate(v)) for k, v in hs.items()}
        FA[name] = rec
        print("FA", name, {a: {k: v[0] for k, v in r.items()} for a, r in rec.items()}, flush=True)
    out["false_alarm"] = FA
    # ---- delays (same data as c1 E2: seeds +1000, S with D*(0.05)=60)
    e2 = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "c1_e2_delay.json"), encoding="utf-8"))["result"]
    DL = {}
    dcases = [(k, None) for k in INSC] + [(f"AQ_{c}|boot168", v["pool"]) for c, v in aq.items()]
    for name, pool in dcases:
        sig2 = C.psi_lrv_vb(**C.SYN_IN[name])[0] if pool is None else float(np.mean([pool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(pool) - 167, 24)]))
        S = np.sqrt(2 * sig2 * np.log(20) / 60.0); I = S ** 2 / (2 * sig2)
        rec = {}
        for nu, lab in ((240, "cond_delay_nu240"), (0, "ARL1_nu0")):
            acc = {k: {a: [] for a in ALPHAS} for k in variants}
            for s in SEEDS:
                rng = np.random.default_rng(s + 1000)
                X = gen_null(name, rng, 200, N0 + nu + 3000, pool, 168); X[:, N0 + nu:] += S
                for k, (kap, hh) in variants.items():
                    M, ends = b1c_path(X, kap)
                    for a in ALPHAS: acc[k][a].append(first_alarm(M, ends, hh[a]))
            r = {}
            for k in variants:
                r[k] = {}
                for a in ALPHAS:
                    f = np.concatenate(acc[k][a]); post = f >= nu; d = (f[post] - nu).astype(float)
                    r[k][str(a)] = {"mean_delay": round(float(d.mean()), 1), "ci95": [round(float(d.mean() - 1.96 * d.std(ddof=1) / np.sqrt(len(d))), 1), round(float(d.mean() + 1.96 * d.std(ddof=1) / np.sqrt(len(d))), 1)],
                                    "pre_change": round(float(((f >= 0) & (f < nu)).mean()), 4), "miss": round(float((f < 0).mean()), 4), "ratio_to_Dstar": round(float(d.mean() / (np.log(1 / a) / I)), 3)}
            r["F2_bpi_wlc"] = {str(a): e2[name][lab]["F2_bpi_wlc"][str(a)] for a in ALPHAS}
            rec[lab] = r
        DL[name] = rec
        print("DL", name, {k: rec["cond_delay_nu240"][k]["0.05"]["mean_delay"] for k in rec["cond_delay_nu240"]}, flush=True)
    out["delay"] = DL
    return out


def e8():
    aq = D.aq_channels(); out = {}
    cases = [(k, None) for k in INSC] + [(f"AQ_{c}|boot168", v["pool"]) for c, v in aq.items()]
    for name, pool in cases:
        if pool is None:
            lrv, vb = C.psi_lrv_vb(**C.SYN_IN[name])
        else:
            lrv = float(np.mean([pool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(pool) - 167, 24)]))
            vb = float(np.mean([pool[i:i + 24].sum() ** 2 / 24 for i in range(0, len(pool) - 23, 6)]))
        target = max(lrv, vb); rec = {}
        for n0 in (336, 1000, 3000, 10000):
            under, over, eps, fa = [], [], [], []
            for s in SEEDS:
                rng = np.random.default_rng(s + 77)
                X = gen_null(name, rng, 400, n0 + W, pool, 168)
                f, st = Q.bpi_wlc(X, [0.05], n0=n0, logA_override=lambda a: np.log(1 / a), return_st=True)
                f = f[0.05]
                under.append((st < vb).any(1)); over.append((st > (1 + 3 * Q.KAPPA) * lrv).any(1))
                eps.append(np.abs(st / (1 + Q.KAPPA) - target).max(1) / target); fa.append((f >= 0) & (f <= W))
            u, o, e, h = map(np.concatenate, (under, over, eps, fa))
            rec[str(n0)] = {"P_lower_violation": wil(u), "P_upper_violation": wil(o) if pool is None else None,
                            "eps_max_p50": round(float(np.quantile(e, .5)), 3), "eps_max_p95": round(float(np.quantile(e, .95)), 3),
                            "F2_FA_A=1/alpha": wil(h)}
        out[name] = {"sigma2_LR": round(lrv, 3), "v_b": round(vb, 3), "by_n0": rec}
        print("E8", name, {n0: (r["P_lower_violation"][0], r["F2_FA_A=1/alpha"][0], r["eps_max_p95"]) for n0, r in rec.items()}, flush=True)
    return out


if __name__ == "__main__":
    t0 = time.time(); exp, outp = sys.argv[1], sys.argv[2]
    r = e7() if exp == "e7" else e8()
    json.dump({"experiment": exp, "result": r, "runtime_s": round(time.time() - t0, 1)}, open(outp, "w", encoding="utf-8"), indent=1)
