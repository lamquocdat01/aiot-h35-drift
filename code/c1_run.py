"""H35 S1-1c experiments E1-E5 per PREREG_H35 v1 + Amendments 1-2.
Usage: python c1_run.py <exp> <out.json>     exp in {e1_syn, e1_aq, e1_intel, e2_arl0, e2_delay, e3, e4, e5, real, gas}
"""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P
import q2_adaptive as Q
import c1_data as D

ALPHA = 0.05
SEEDS = [20261121 + k for k in range(5)]
N0, W = 336, 720
SYN_IN = {"AR1_0": dict(phi=()), "AR1_0.5": dict(phi=(0.5,)), "AR1_0.8": dict(phi=(0.8,)), "AR1_0.95": dict(phi=(0.95,)),
          "ARMA11_0.8_+0.4": dict(phi=(0.8,), theta=(0.4,)), "ARMA11_0.8_-0.5": dict(phi=(0.8,), theta=(-0.5,)),
          "ARMA21_1.2_-0.4_+0.3": dict(phi=(1.2, -0.4), theta=(0.3,))}
SYN_OUT = ["AR0.8_t5", "AR0.8_skew", "longmem"]
BOOTS = (24, 72, 168)


# ---------------------------------------------------------------- generators
def gen_out(name, rng, R, n):
    if name == "longmem":
        x = np.zeros((R, n)); u = rng.standard_normal(R) * 3
        for t in range(n):
            u = 0.999 * u + 0.13 * rng.standard_normal(R); x[:, t] = u + rng.standard_normal(R)
        return x
    e = rng.standard_t(5, (R, n + 500)) / np.sqrt(5 / 3) if name == "AR0.8_t5" else rng.exponential(size=(R, n + 500)) - 1.0
    x = np.zeros_like(e)
    for t in range(1, e.shape[1]): x[:, t] = 0.8 * x[:, t - 1] + e[:, t]
    return x[:, -n:]


def gen(name, rng, R, n, pool=None, blk=None):
    if pool is not None: return P.cbb(pool, R, n, blk, rng)
    if name in SYN_IN: return Q.gen_arma(rng, R, n, **SYN_IN[name])
    return gen_out(name, rng, R, n)


def psi_lrv_vb(phi=(), theta=(), b=24, H=5000):
    """true sigma2_LR and v_b for unit-innovation ARMA via impulse response."""
    psi = np.zeros(H); psi[0] = 1.0
    for j in range(1, H):
        v = theta[j - 1] if j - 1 < len(theta) else 0.0
        for i, ph in enumerate(phi):
            if j - 1 - i >= 0: v += ph * psi[j - 1 - i]
        psi[j] = v
    gam = np.array([np.dot(psi[: H - h], psi[h:]) for h in range(b)])
    lrv = psi.sum() ** 2
    vb = gam[0] + 2 * sum((1 - h / b) * gam[h] for h in range(1, b))
    return float(lrv), float(vb)


# ---------------------------------------------------------------- calibrated thresholds (iid N(0,1)), alpha -> h
def calib_cusum_h(alpha, n, k=0.5, sims=200000, seed=7):
    rng = np.random.default_rng(seed); mx = []
    for _ in range(sims // 20000):
        Z = rng.standard_normal((20000, n))
        sp = np.zeros(20000); sm = np.zeros(20000); m = np.zeros(20000)
        for t in range(n):
            sp = np.maximum(0, sp + Z[:, t] - k); sm = np.maximum(0, sm - Z[:, t] - k); m = np.maximum(m, np.maximum(sp, sm))
        mx.append(m)
    return float(np.quantile(np.concatenate(mx), 1 - alpha))


def calib_ph(alpha, n, sims=200000, seed=8):
    rng = np.random.default_rng(seed); mx = []
    for _ in range(sims // 20000):
        mx.append(P.ph_stat(rng.standard_normal((20000, n))).max(1))
    return float(np.quantile(np.concatenate(mx), 1 - alpha))


LB = (W - Q.B) // (Q.B + Q.G) + 1  # 24 blocks in the audit window
TH = {}


def thresholds(alphas=(0.05, 1e-3)):
    for a in alphas:
        if a in TH: continue
        TH[a] = {"cusum": calib_cusum_h(a, W), "ph": calib_ph(a, W) if a == 0.05 else None,
                 "cusum_blocks": calib_cusum_h(a, LB)}
    return TH


# ---------------------------------------------------------------- baselines B1c, B2c
def levinson_phi(gam, n_eff, pmax=30):
    R = gam.shape[0]; g0 = np.maximum(gam[:, 0], 1e-9)
    phi = np.zeros((R, pmax + 1)); e = g0.copy()
    best = n_eff * np.log(e); bphi = np.zeros((R, pmax + 1)); be = e.copy()
    for k in range(1, pmax + 1):
        acc = gam[:, k] - (phi[:, 1:k] * gam[:, k - 1:0:-1]).sum(1) if k > 1 else gam[:, 1]
        kap = np.clip(acc / np.maximum(e, 1e-12), -0.999, 0.999)
        nphi = phi.copy(); nphi[:, k] = kap
        if k > 1: nphi[:, 1:k] = phi[:, 1:k] - kap[:, None] * phi[:, k - 1:0:-1]
        phi = nphi; e = e * (1 - kap ** 2)
        aic = n_eff * np.log(np.maximum(e, 1e-12)) + 2 * k
        bt = aic < best
        best = np.where(bt, aic, best); bphi = np.where(bt[:, None], phi, bphi); be = np.where(bt, e, be)
    return bphi, be


def autocov0(Xb, pmax=30):
    n = Xb.shape[1]
    return np.stack([(Xb[:, : n - h] * Xb[:, h:]).sum(1) / (n - h) for h in range(pmax + 1)], 1)


def b1c_frozen(X, alphas, n0=N0, b=Q.B, g=Q.G):
    R, n = X.shape
    gam = autocov0(X[:, :n0]); st = Q.s2_tilde(gam, n0, b)
    step = b + g; L = (n - n0 - b) // step + 1
    sp = np.zeros(R); sm = np.zeros(R); out = {a: np.full(R, -1) for a in alphas}
    for l in range(L):
        s = n0 + l * step; z = X[:, s:s + b].sum(1) / np.sqrt(b * st)
        sp = np.maximum(0, sp + z - 0.5); sm = np.maximum(0, sm - z - 0.5); M = np.maximum(sp, sm)
        for a in alphas:
            hit = (M >= TH[a]["cusum_blocks"]) & (out[a] < 0); out[a][hit] = s + b - n0
    return out


def b2c_arres(X, alphas, n0=N0):
    R, n = X.shape
    gam = autocov0(X[:, :n0]); phi, e = levinson_phi(gam, n0)
    p = 30; pred = np.zeros((R, n - n0))
    for i in range(1, p + 1):
        pred += phi[:, i:i + 1] * X[:, n0 - i: n - i]
    Z = (X[:, n0:] - pred) / np.sqrt(e)[:, None]
    out = {}
    S = P.cusum_stat(Z)
    for a in alphas:
        hit = S >= TH[a]["cusum"]; out[a] = np.where(hit.any(1), hit.argmax(1) + 1, -1)
    return out


def m_monitors(X, rng, n0=N0, which=("M1", "M2", "M3", "M4", "M5", "M6")):
    sd = np.sqrt(np.maximum(X[:, :n0].var(1), 1e-9))[:, None]
    Z = X[:, n0:] / sd; o = {}
    if "M1" in which: o["M1_cusum"] = P.first_cross(P.cusum_stat(Z), TH[0.05]["cusum"])
    if "M2" in which: o["M2_pagehinkley"] = P.first_cross(P.ph_stat(Z), TH[0.05]["ph"])
    if "M3" in which: o["M3_adwin"] = P.mon_adwin(Z)
    if "M4" in which: o["M4_repeated_p"] = P.mon_repeated_p(Z)
    if "M5" in which: o["M5_eprocess"] = P.mon_eprocess(Z)
    if "M6" in which: o["M6_ctm"] = P.mon_ctm(Z, rng)
    return o


def all_fa(X, rng):
    o = m_monitors(X[:, : N0 + W], rng)
    Xw = X[:, : N0 + W]
    o["F1_sn_block"] = Q.sn_block(Xw, [ALPHA])[ALPHA]
    o["F2_bpi_wlc"] = Q.bpi_wlc(Xw, [ALPHA])[ALPHA]
    o["B1c_frozen_lrv_cusum"] = b1c_frozen(Xw, [ALPHA])[ALPHA]
    o["B2c_ar_residual_cusum"] = b2c_arres(Xw, [ALPHA])[ALPHA]
    return {k: (v >= 0) & (v <= W) for k, v in o.items()}


def fa_block(gen_fn, R):
    acc = {}; ps = {}
    for s in SEEDS:
        rng = np.random.default_rng(s); X = gen_fn(rng, R)
        for k, h in all_fa(X, rng).items():
            acc.setdefault(k, []).append(h); ps.setdefault(k, []).append(round(float(h.mean()), 4))
    return {k: {"rate_wilson95": P.wilson(int(np.concatenate(v).sum()), len(np.concatenate(v))), "per_seed": ps[k]} for k, v in acc.items()}


# ---------------------------------------------------------------- experiments
def e1_syn():
    thresholds(); out = {}
    for name in list(SYN_IN) + SYN_OUT:
        out[name] = fa_block(lambda rng, R: gen(name, rng, R, N0 + W), 400); print("E1", name, {k: v["rate_wilson95"][0] for k, v in out[name].items()}, flush=True)
    return out


def e1_real(chans, R=200):
    thresholds(); out = {}
    for ch, d in chans.items():
        for bb in BOOTS:
            out[f"{ch}|boot{bb}"] = fa_block(lambda rng, R_: P.cbb(d["pool"], R_, N0 + W, bb, rng), R)
            print("E1", ch, bb, {k: v["rate_wilson95"][0] for k, v in out[f"{ch}|boot{bb}"].items()}, flush=True)
    return out


def e2_arl0():
    thresholds(); out = {}; NL = 20000; R = 100
    aq = D.aq_channels()
    cases = [(k, None) for k in SYN_IN] + [(f"AQ_{k}|boot168", v["pool"]) for k, v in aq.items()]
    for name, pool in cases:
        acc = {}
        for s in SEEDS:
            rng = np.random.default_rng(s + 50)
            X = P.cbb(pool, R, N0 + NL, 168, rng) if pool is not None else Q.gen_arma(rng, R, N0 + NL, **SYN_IN[name])
            res = {"F1_sn_block": Q.sn_block(X, [ALPHA])[ALPHA], "F2_bpi_wlc": Q.bpi_wlc(X, [ALPHA])[ALPHA],
                   "B1c_frozen_lrv_cusum": b1c_frozen(X, [ALPHA])[ALPHA], "B2c_ar_residual_cusum": b2c_arres(X, [ALPHA])[ALPHA]}
            sd = np.sqrt(X[:, :N0].var(1))[:, None]; Z = X[:, N0:] / sd
            res["M1_cusum"] = P.first_cross(P.cusum_stat(Z), TH[ALPHA]["cusum"]); res["M5_eprocess"] = P.mon_eprocess(Z)
            for k, v in res.items(): acc.setdefault(k, []).append(v)
        out[name] = {}
        for k, v in acc.items():
            f = np.concatenate(v).astype(float); cens = f < 0; rl = np.where(cens, NL, f)
            out[name][k] = {"ARL0_lower_bound_mean": round(float(rl.mean()), 1), "censored_frac": round(float(cens.mean()), 3),
                            "median": round(float(np.median(rl)), 1)}
        print("ARL0", name, {k: (v["ARL0_lower_bound_mean"], v["censored_frac"]) for k, v in out[name].items()}, flush=True)
    return out


def delay_block(Xfn, S, I, nu, alphas, R=200):
    acc = {}
    for s in SEEDS:
        rng = np.random.default_rng(s + 1000); X = Xfn(rng, R); X[:, N0 + nu:] += S
        res = {"F1_sn_block": Q.sn_block(X, alphas), "F2_bpi_wlc": Q.bpi_wlc(X, alphas),
               "B1c_frozen_lrv_cusum": b1c_frozen(X, [a for a in alphas if a in TH]),
               "B2c_ar_residual_cusum": b2c_arres(X, [a for a in alphas if a in TH])}
        sd = np.sqrt(X[:, :N0].var(1))[:, None]; Z = X[:, N0:] / sd; Sx = P.cusum_stat(Z)
        res["M1_cusum"] = {a: P.first_cross(Sx, TH[a]["cusum"]) for a in alphas if a in TH}
        for m, d in res.items():
            for a, f in d.items(): acc.setdefault(m, {}).setdefault(a, []).append(f)
    out = {}
    for m, d in acc.items():
        out[m] = {}
        for a, v in d.items():
            f = np.concatenate(v); post = f >= nu; dd = (f[post] - nu).astype(float); Ds = np.log(1 / a) / I
            out[m][str(a)] = {"pre_change_alarm": round(float(((f >= 0) & (f < nu)).mean()), 4), "miss": round(float((f < 0).mean()), 4),
                              "mean_delay": round(float(dd.mean()), 1) if len(dd) else None,
                              "ci95": [round(float(dd.mean() - 1.96 * dd.std(ddof=1) / np.sqrt(len(dd))), 1), round(float(dd.mean() + 1.96 * dd.std(ddof=1) / np.sqrt(len(dd))), 1)] if len(dd) > 2 else None,
                              "ratio_to_Dstar": round(float(dd.mean() / Ds), 3) if len(dd) else None}
    return out


def pred_overhead(a):
    return round(float((1 + Q.KAPPA) * (1 + Q.G / Q.B) * (np.log((int(np.ceil(Q.HMON / (Q.B + Q.G))) + Q.WM) / a) + np.log(len(Q.TH))) / np.log(1 / a)), 3)


def e2_delay():
    thresholds(); out = {}; alphas = [0.05, 1e-3, 1e-6]; POST = 3000
    aq = D.aq_channels()
    cases = [(k, None) for k in SYN_IN] + [(f"AQ_{k}|boot168", v["pool"]) for k, v in aq.items()]
    for name, pool in cases:
        if pool is None: sig2 = psi_lrv_vb(**SYN_IN[name])[0]
        else: sig2 = float(np.mean([pool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(pool) - 167, 24)]))
        S = np.sqrt(2 * sig2 * np.log(1 / 0.05) / 60.0); I = S ** 2 / (2 * sig2)
        rec = {"S": round(float(S), 4), "sigma2_LR": round(sig2, 3), "pred_overhead": {str(a): pred_overhead(a) for a in alphas}}
        for nu, lab in ((0, "ARL1_nu0"), (240, "cond_delay_nu240")):
            fn = (lambda rng, R, nu=nu: P.cbb(pool, R, N0 + nu + POST, 168, rng)) if pool is not None else (lambda rng, R, nu=nu: Q.gen_arma(rng, R, N0 + nu + POST, **SYN_IN[name]))
            rec[lab] = delay_block(fn, S, I, nu, alphas)
        out[name] = rec
        print("DELAY", name, {m: rec["cond_delay_nu240"][m].get("1e-06", {}).get("ratio_to_Dstar") for m in rec["cond_delay_nu240"]}, flush=True)
    return out


def e3():
    thresholds(); out = {}; aq = D.aq_channels(); co = aq["CO(GT)"]["pool"]
    for b in (12, 24, 48):
        for g in (3, 6, 12):
            Q.B, Q.G = b, g
            rec = {}
            for name, pool in (("AR1_0.8", None), ("AQ_CO|boot168", co)):
                hs = []
                for s in SEEDS:
                    rng = np.random.default_rng(s)
                    X = P.cbb(pool, 400, N0 + W, 168, rng) if pool is not None else Q.gen_arma(rng, 400, N0 + W, phi=(0.8,))
                    f = Q.bpi_wlc(X, [ALPHA], b=b, g=g)[ALPHA]; hs.append((f >= 0) & (f <= W))
                h = np.concatenate(hs)
                sig2 = psi_lrv_vb(phi=(0.8,))[0] if pool is None else float(np.mean([pool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(pool) - 167, 24)]))
                S = np.sqrt(2 * sig2 * np.log(20) / 60.0); I = S ** 2 / (2 * sig2)
                ds = []
                for s in SEEDS:
                    rng = np.random.default_rng(s + 1000)
                    X = P.cbb(pool, 200, N0 + 240 + 3000, 168, rng) if pool is not None else Q.gen_arma(rng, 200, N0 + 240 + 3000, phi=(0.8,))
                    X[:, N0 + 240:] += S; f = Q.bpi_wlc(X, [1e-6], b=b, g=g)[1e-6]; ds.append(f[f >= 240] - 240)
                d = np.concatenate(ds).astype(float)
                rec[name] = {"FA": P.wilson(int(h.sum()), len(h)), "ratio_1e-6": round(float(d.mean() / (np.log(1e6) / I)), 3) if len(d) else None}
            out[f"F2_b{b}_g{g}"] = rec; print("E3", b, g, rec, flush=True)
    Q.B, Q.G = 24, 6
    for b in (24, 48):
        Q.BSN, Q.GSN = b, b; rec = {}
        for name, pool in (("AR1_0.8", None), ("AQ_CO|boot168", co)):
            hs = []
            for s in SEEDS:
                rng = np.random.default_rng(s)
                X = P.cbb(pool, 400, N0 + W, 168, rng) if pool is not None else Q.gen_arma(rng, 400, N0 + W, phi=(0.8,))
                f = Q.sn_block(X, [ALPHA], b=b, g=b)[ALPHA]; hs.append((f >= 0) & (f <= W))
            h = np.concatenate(hs); rec[name] = {"FA": P.wilson(int(h.sum()), len(h))}
        out[f"F1_b{b}_g{b}"] = rec; print("E3", "F1", b, rec, flush=True)
    Q.BSN, Q.GSN = 24, 24
    return out


def est_traj(X, n0=N0, b=Q.B, g=Q.G, wm=Q.WM):
    """pre-inflation estimator max(lrv, v_b) actually used at each block (lagged), as in bpi_wlc."""
    R, n = X.shape; step = b + g; L = (n - n0 - b) // step + 1
    hist = []; out = []
    st0 = Q.s2_tilde(autocov0(X[:, :n0]), n0, b) / (1 + Q.KAPPA)
    for l in range(L):
        s = n0 + l * step
        hist.append(Q.s2_tilde(autocov0(X[:, :s]), s, b) / (1 + Q.KAPPA))
        out.append(hist[l - wm] if l >= wm else st0)
    return np.stack(out, 1)


def e4():
    out = {}; NL = 3000; R = 200; aq = D.aq_channels(); intel = D.intel_residuals()
    cases = [(k, None, None) for k in SYN_IN] + [(f"AQ_{k}|boot168", v["pool"], None) for k, v in aq.items()] + \
            [(f"INTEL_{m}|boot168", v["pool"], "explore") for m, v in list(intel.items())[:5]]
    for name, pool, tag in cases:
        if pool is None:
            lrv, vb = psi_lrv_vb(**SYN_IN[name]); target = max(lrv, vb)
        else:
            lrv = float(np.mean([pool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(pool) - 167, 24)]))
            vb = float(np.mean([pool[i:i + 24].sum() ** 2 / 24 for i in range(0, len(pool) - 23, 6)])); target = max(lrv, vb)
        eps = []; under = []
        for s in SEEDS:
            rng = np.random.default_rng(s + 77)
            X = P.cbb(pool, R, N0 + NL, 168, rng) if pool is not None else Q.gen_arma(rng, R, N0 + NL, **SYN_IN[name])
            T = est_traj(X)
            eps.append(np.abs(T - target).max(1) / target); under.append(((1 + Q.KAPPA) * T < vb).any(1))
        e = np.concatenate(eps); u = np.concatenate(under)
        out[name] = {"target": round(target, 3), "sigma2_LR": round(lrv, 3), "v_b": round(vb, 3), "eps_max_p50": round(float(np.quantile(e, .5)), 3),
                     "eps_max_p95": round(float(np.quantile(e, .95)), 3), "P_any_under_vb": P.wilson(int(u.sum()), len(u)), "tag": tag or "decision"}
        print("E4", name, out[name], flush=True)
    return out


def e5():
    """R3 semi-synthetic 16 Gas channels; iid joint-row bootstrap (order not confirmed, Amendment 2)."""
    B, G, Rr = D.gas_residuals(); pool = Rr[B == 1]; pool = (pool - pool.mean(0)) / pool.std(0)
    N = 16; T = 4000; Rep = 200; b, g = 24, 24; LMIN = 50
    S = np.sqrt(2 * 1.0 * np.log(20) / 60.0)
    CG = np.array([0.1, 0.2, 0.35, 0.5, 0.75, 1.0]) / np.sqrt(b)
    res = {}
    for rule in ("R3a_eBH_reset", "R3b_discovery_investing"):
        acc = {"life": [], "sup": [], "cur_short": [], "nrec": []}
        for s in SEEDS:
            rng = np.random.default_rng(s + 500)
            start = np.zeros((Rep, N), int); onset = start + rng.geometric(1 / 300, (Rep, N))
            logE = np.zeros((Rep, N)); logM = np.zeros((Rep, N, 12)); bsum = np.zeros((Rep, N))
            Rc = np.zeros(Rep); Fc = np.zeros(Rep); sup = np.zeros(Rep); kidx = np.tile(np.arange(1, N + 1.0), (Rep, 1)); nst = np.full(Rep, float(N))
            tau = np.full(Rep, -1); cur = np.zeros(Rep)
            for t in range(1, T + 1):
                rows = pool[rng.integers(0, len(pool), Rep)]
                x = rows + S * (t >= onset)
                ph = (t - start - 1) % (b + g)  # phase within (gap, block)
                inblk = (ph >= g)
                bsum += np.where(inblk, x, 0.0)
                done = inblk & (ph == b + g - 1)
                if done.any():
                    Y = bsum[..., None]
                    inc = np.concatenate([Y * CG - Y ** 2 * CG ** 2 / 2, -Y * CG - Y ** 2 * CG ** 2 / 2], axis=2)
                    logM = np.where(done[..., None], logM + inc, logM)
                    logE = np.where(done, np.logaddexp.reduce(logM, axis=2) - np.log(12), logE)
                    bsum = np.where(done, 0.0, bsum)
                    if rule == "R3a_eBH_reset":
                        srt = -np.sort(-logE, 1); thr = np.log(N / (ALPHA * np.arange(1, N + 1)))
                        ok = srt >= thr; anyok = ok.any(1); kmax = N - np.argmax(ok[:, ::-1], 1)
                        cut = srt[np.arange(Rep), kmax - 1]; alarm = (logE >= cut[:, None]) & anyok[:, None] & done
                    else:
                        gam = np.where(kidx <= N, 1 / (2 * N), N / (2 * kidx * np.maximum(kidx - 1, 1)))
                        alarm = done & (logE >= -np.log(ALPHA * gam * (Rc[:, None] + 1)))
                    if alarm.any():
                        fa = alarm & (t < onset); nA = alarm.sum(1); nF = fa.sum(1)
                        short = (alarm & ((t - start) <= LMIN)).any(1); hit = short & (tau < 0)
                        tau[hit] = t; cur[hit] = nF[hit] / np.maximum(nA[hit], 1)
                        Rc += nA; Fc += nF
                        order = np.cumsum(alarm, 1); kidx = np.where(alarm, nst[:, None] + order, kidx); nst += nA
                        start = np.where(alarm, t, start); onset = np.where(alarm, t + rng.geometric(1 / 300, (Rep, N)), onset)
                        logE = np.where(alarm, 0.0, logE); logM = np.where(alarm[..., None], 0.0, logM); bsum = np.where(alarm, 0.0, bsum)
                sup = np.maximum(sup, np.where(Rc > 0, Fc / np.maximum(Rc, 1), 0))
            acc["life"].append(np.where(Rc > 0, Fc / np.maximum(Rc, 1), 0)); acc["sup"].append(sup)
            acc["cur_short"].append(cur); acc["nrec"].append(Rc)
        ci = lambda a: [round(float(a.mean()), 4), round(float(a.mean() - 1.96 * a.std(ddof=1) / np.sqrt(len(a))), 4), round(float(a.mean() + 1.96 * a.std(ddof=1) / np.sqrt(len(a))), 4)]
        res[rule] = {"lifetime_FDR_T": ci(np.concatenate(acc["life"])), "E_sup_FDP": ci(np.concatenate(acc["sup"])),
                     "current_FDR_at_tau_short": ci(np.concatenate(acc["cur_short"])), "mean_recal": round(float(np.concatenate(acc["nrec"]).mean()), 2)}
        print("E5", rule, res[rule], flush=True)
    res["note"] = "iid joint-row bootstrap of batch-1 oos residuals (16 sensors); sample order not confirmed (Amendment 2) -> illustrative only"
    return res


def real():
    """drift that: first alarm after calibration on real series (exploratory, no ground truth)."""
    thresholds(); out = {}
    def run(name, full, start):
        x = np.nan_to_num(full[start - N0:], nan=0.0)[None, :]
        o = {"F1_sn_block": Q.sn_block(x, [ALPHA])[ALPHA][0], "F2_bpi_wlc": Q.bpi_wlc(x, [ALPHA])[ALPHA][0],
             "B1c_frozen_lrv_cusum": b1c_frozen(x, [ALPHA])[ALPHA][0], "B2c_ar_residual_cusum": b2c_arres(x, [ALPHA])[ALPHA][0]}
        sd = np.sqrt(x[:, :N0].var()); Z = x[:, N0:] / sd
        o["M1_cusum"] = P.first_cross(P.cusum_stat(Z), TH[ALPHA]["cusum"])[0]; o["M5_eprocess"] = P.mon_eprocess(Z)[0]
        out[name] = {k: int(v) for k, v in o.items()}; out[name]["n_monitored"] = int(x.shape[1] - N0)
        print("REAL", name, out[name], flush=True)
    aq = D.aq_channels()
    for k, v in aq.items():
        st = int(np.argmax(v["index"] >= np.datetime64("2004-04-01")))
        run(f"AQ_{k}_from_2004-04-01_hours", v["full"], st)
    intel = D.intel_residuals()
    for m, v in intel.items():
        run(f"INTEL_{m}_from_day8_10min", v["full"], 7 * 144)
    return out


def gas():
    B, G, R = D.gas_residuals(); out = {"note": "batch-level only (Amendment 2): mean normalized residual per (sensor, batch); iid bootstrap CI within batch"}
    rng = np.random.default_rng(SEEDS[0]); tab = {}
    for s in range(16):
        row = {}
        for b in range(1, 11):
            v = R[B == b, s]; v = v[np.isfinite(v)]
            bs = np.array([rng.choice(v, len(v)).mean() for _ in range(500)])
            row[str(b)] = [round(float(v.mean()), 3), round(float(np.quantile(bs, .025)), 3), round(float(np.quantile(bs, .975)), 3)]
        tab[f"sensor{s + 1}"] = row
    out["mean_resid"] = tab
    out["n_per_batch"] = {str(b): int((B == b).sum()) for b in range(1, 11)}
    out["frac_sensors_CI_excl0_by_batch"] = {str(b): round(float(np.mean([(tab[k][str(b)][1] > 0) or (tab[k][str(b)][2] < 0) for k in tab])), 3) for b in range(1, 11)}
    print(out["frac_sensors_CI_excl0_by_batch"])
    return out


def main():
    exp, outp = sys.argv[1], sys.argv[2]; t0 = time.time()
    if exp == "e1_syn": r = e1_syn()
    elif exp == "e1_aq": r = e1_real({f"AQ_{k}": v for k, v in D.aq_channels().items()})
    elif exp == "e1_intel": r = e1_real({f"INTEL_{m}": v for m, v in D.intel_residuals().items()})
    elif exp == "e2_arl0": r = e2_arl0()
    elif exp == "e2_delay": r = e2_delay()
    elif exp == "e3": r = e3()
    elif exp == "e4": r = e4()
    elif exp == "e5": r = e5()
    elif exp == "real": r = real()
    elif exp == "gas": r = gas()
    r = {"experiment": exp, "result": r, "thresholds": {str(k): v for k, v in TH.items()}, "runtime_s": round(time.time() - t0, 1)}
    json.dump(r, open(outp, "w", encoding="utf-8"), indent=1)


if __name__ == "__main__":
    main()
