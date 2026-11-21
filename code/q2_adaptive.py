"""H35 Q1ROUTE Q-2: adaptive block plug-in window-limited e-CUSUM (BPI-WLC) vs blocked self-normalized (Thm 1)
vs original e-process.  Finite-sample checks of Theorem 3 (asymptotic) — false alarm and delay / Tartakovsky bound.

BPI-WLC(b, g, W_m, A, kappa):
  blocks of length b separated by gaps g after a burn-in of n0 in-control hours;
  long-run variance estimated online, LAGGED by W_m blocks (data before block l-W_m; burn-in only at start):
     AR-sieve: Levinson-Durbin on uncentred sample autocov up to PMAX, order by AIC,
               s2_LR = s2_p / (1 - sum phi)^2 ;   Bartlett block variance v_b = g0 + 2 sum_{h<b}(1-h/b) g(h)
     s2_tilde = (1+kappa) * max(s2_LR, v_b)
  z_l = Y_l / sqrt(b s2_tilde_l)  (Y_l = block sum);  for each start j in window W_m (blocks):
     M_{j,l} = mean_theta prod_{i=j..l} exp(theta z_i - theta^2/2), theta in +/- geometric grid
  alarm at first block end with max_j M_{j,l} >= A = (W + W_m)/alpha  (W = blocks in the 720-h audit window).
Benchmark (Tartakovsky 2017, Ex. 5.1): D* = |log alpha| / I,  I = S^2 / (2 sigma_LR^2).
Usage: python q2_adaptive.py <AirQualityUCI.csv> <out.json>
"""
import json, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P

SEEDS = [20261121 + k for k in range(5)]
N0 = 336          # burn-in (in-control, e.g. calibration fortnight)
HMON = 720        # audit window (h)
B, G = 24, 6      # BPI blocks
WM = 60           # window (blocks) of the CUSUM starts (round-1: 30 -> 60 to cover delays at small alpha)
PMAX = 30
KAPPA = 0.10
TH = np.geomspace(0.05, 4.0, 24); TH = np.concatenate([TH, -TH])
BSN, GSN = 24, 24 # Thm 1 self-normalized blocks
CG = np.array([0.1, 0.2, 0.35, 0.5, 0.75, 1.0])
LAMS2 = P.LAMS2


# ------------------------------------------------------------------ null generators
def gen_arma(rng, R, n, phi=(), theta=(), burn=500):
    p, q = len(phi), len(theta)
    e = rng.standard_normal((R, n + burn + q))
    x = np.zeros((R, n + burn + q))
    for t in range(max(p, q), n + burn + q):
        v = e[:, t].copy()
        for i, th in enumerate(theta): v += th * e[:, t - 1 - i]
        for i, ph in enumerate(phi): v += ph * x[:, t - 1 - i]
        x[:, t] = v
    return x[:, -n:]


def lrv_arma(phi=(), theta=()):
    return (1 + sum(theta)) ** 2 / (1 - sum(phi)) ** 2  # unit innovation variance


# ------------------------------------------------------------------ estimator
def levinson_aic(gam, n_eff):
    """gam (R, PMAX+1) autocov. returns s2_LR via AIC-chosen AR order (vectorised)."""
    R = gam.shape[0]
    g0 = np.maximum(gam[:, 0], 1e-9)
    phi = np.zeros((R, PMAX + 1)); e = g0.copy()
    best_aic = n_eff * np.log(e); best_lrv = g0.copy()  # order 0
    for k in range(1, PMAX + 1):
        acc = gam[:, k] - (phi[:, 1:k] * gam[:, k - 1:0:-1]).sum(1) if k > 1 else gam[:, 1]
        kap = np.clip(acc / np.maximum(e, 1e-12), -0.999, 0.999)
        newphi = phi.copy()
        newphi[:, k] = kap
        if k > 1:
            newphi[:, 1:k] = phi[:, 1:k] - kap[:, None] * phi[:, k - 1:0:-1]
        phi = newphi
        e = e * (1 - kap ** 2)
        aic = n_eff * np.log(np.maximum(e, 1e-12)) + 2 * k
        sphi = np.clip(phi[:, 1:k + 1].sum(1), -0.99, 0.995)
        lrv = e / (1 - sphi) ** 2
        better = aic < best_aic
        best_aic = np.where(better, aic, best_aic); best_lrv = np.where(better, lrv, best_lrv)
    return best_lrv


def s2_tilde(gam, n_eff, b):
    lrv = levinson_aic(gam, n_eff)
    h = np.arange(1, b)
    vb = gam[:, 0] + 2 * ((1 - h / b) * gam[:, 1:b]).sum(1)
    return (1 + KAPPA) * np.maximum(lrv, vb)


# ------------------------------------------------------------------ monitors
def bpi_wlc(X, alphas, n0=N0, b=B, g=G, wm=WM, wblocks=None, nmon=None, lag=True, logA_override=None):
    """X (R, n0+nmon). returns dict alpha -> alarm time (hours after monitoring start, -1 none)."""
    R, n = X.shape
    nmon = n - n0 if nmon is None else nmon
    step = b + g
    L = (nmon - b) // step + 1
    wblocks = int(np.ceil(HMON / step)) if wblocks is None else wblocks
    logA = {a: (np.log((wblocks + wm) / a) if logA_override is None else logA_override(a)) for a in alphas}
    first = {a: np.full(R, -1) for a in alphas}
    # running uncentred cross products C_h = sum_t x_t x_{t+h} over data seen so far
    C = np.zeros((R, PMAX + 1)); seen = 0

    def absorb(upto):
        nonlocal seen
        if upto <= seen: return
        lo = max(0, seen - PMAX)
        seg = X[:, lo:upto]
        for hh in range(PMAX + 1):
            # products x_t x_{t+h} with t+h in [seen, upto)
            a = seg[:, : seg.shape[1] - hh]; c = seg[:, hh:]
            tpos = np.arange(lo + hh, upto)
            mask = tpos >= seen
            C[:, hh] += (a[:, mask] * c[:, mask]).sum(1)
        seen = upto

    logP = np.full((R, wm, len(TH)), -np.inf)  # ring buffer of starts
    hist = []                                   # estimator after absorbing data before block l
    absorb(n0)
    st0 = s2_tilde(C / np.maximum(seen - np.arange(PMAX + 1), 1), seen, b)
    for l in range(L):
        s = n0 + l * step
        absorb(s)
        gam = C / np.maximum(seen - np.arange(PMAX + 1), 1)
        hist.append(s2_tilde(gam, seen, b))
        # LAGGED estimator: data older than the CUSUM window (round-1 fix: no post-change contamination
        # while a start inside the window is still collecting evidence)
        st = (hist[l - wm] if l >= wm else st0) if lag else hist[l]
        Y = X[:, s:s + b].sum(1)
        z = Y / np.sqrt(b * st)
        inc = z[:, None] * TH - TH ** 2 / 2
        logP = logP + inc[:, None, :]
        logP[:, l % wm, :] = inc                # new start (overwrites the oldest)
        M = np.logaddexp.reduce(logP, axis=2) - np.log(len(TH))
        Mmax = M.max(1)
        tend = (s + b) - n0                     # hour (after monitoring start) of the block end
        for a in alphas:
            hit = (Mmax >= logA[a]) & (first[a] < 0)
            first[a][hit] = tend
        absorb(s + step)
    return first


def sn_block(X, alphas, n0=N0, b=BSN, g=GSN):
    R, n = X.shape
    nmon = n - n0
    burn = X[:, :n0]
    sig0 = np.sqrt(b * np.maximum(burn.var(1), 1e-9))[:, None]
    step = b + g; K = (nmon - b) // step + 1
    S = np.stack([X[:, n0 + k * step: n0 + k * step + b].sum(1) for k in range(K)], 1)
    ends = n0 + np.arange(K) * step + b - n0
    logs = []
    for d in (1.0, -1.0):
        Yd = d * S
        lam = CG[None, None, :] / sig0[:, :, None]
        inc = Yd[..., None] * lam - (Yd[..., None] ** 2) * lam ** 2 / 2
        logs.append(np.cumsum(inc, axis=1))
    Lg = np.concatenate(logs, axis=2)
    logE = np.logaddexp.reduce(Lg, axis=2) - np.log(Lg.shape[2])
    out = {}
    for a in alphas:
        hit = logE >= np.log(1 / a)
        out[a] = np.where(hit.any(1), ends[hit.argmax(1)], -1)
    return out


def eproc_orig(X, alphas, n0=N0):
    R, n = X.shape
    sd = np.sqrt(np.maximum(X[:, :n0].var(1), 1e-9))[:, None]
    Z = X[:, n0:] / sd
    logw = np.cumsum(Z[..., None] * LAMS2 - LAMS2 ** 2 / 2, axis=1)
    logE = np.logaddexp.reduce(logw, axis=2) - np.log(len(LAMS2))
    out = {}
    for a in alphas:
        hit = logE >= np.log(1 / a)
        out[a] = np.where(hit.any(1), hit.argmax(1) + 1, -1)
    return out


# ------------------------------------------------------------------ experiments
def fa_rate(first, horizon=HMON):
    h = (first >= 0) & (first <= horizon)
    return h


def wil(h):
    return P.wilson(int(h.sum()), len(h))


def b1_pool(csv):
    df = P.load(csv)
    s1 = df["PT08.S1(CO)"].to_numpy(float); co = df["CO(GT)"].to_numpy(float); idx = df.index
    march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
    roos = np.full(len(s1), np.nan)
    for d in np.unique(days[march]):
        te = march & (days == d); tr = march & (days != d)
        c = P.fit(s1[tr], co[tr]); roos[te] = P.pred(c, s1[te]) - co[te]
    z = roos[march & np.isfinite(roos)]
    return (z - z.mean()) / z.std()


def main():
    t0 = time.time()
    csv, outp = sys.argv[1], sys.argv[2]
    zpool = b1_pool(csv)
    nulls = {"AR1_0": dict(phi=()), "AR1_0.5": dict(phi=(0.5,)), "AR1_0.8": dict(phi=(0.8,)),
             "AR1_0.95": dict(phi=(0.95,)), "ARMA11_0.8_+0.4": dict(phi=(0.8,), theta=(0.4,)),
             "ARMA11_0.8_-0.5": dict(phi=(0.8,), theta=(-0.5,)), "ARMA21_1.2_-0.4_+0.3": dict(phi=(1.2, -0.4), theta=(0.3,))}
    res = {"config": dict(seeds=SEEDS, N0=N0, HMON=HMON, b=B, g=G, WM=WM, PMAX=PMAX, kappa=KAPPA,
                          theta_grid=f"+/- geomspace(0.05,4,24)", A="(W+WM)/alpha, W=ceil(720/(b+g))=24",
                          sn=dict(b=BSN, g=GSN), alpha_FA=0.05)}
    # ---- false alarm in the 720-h window
    FA = {}
    RF = 400
    for name, kw in list(nulls.items()) + [("B1_boot24", 24), ("B1_boot72", 72), ("B1_boot168", 168)]:
        acc = {"bpi_wlc": [], "sn_block_thm1": [], "eprocess_orig": []}; perseed = {k: [] for k in acc}
        for s in SEEDS:
            rng = np.random.default_rng(s)
            if isinstance(kw, dict):
                X = gen_arma(rng, RF, N0 + HMON, **kw)
            else:
                X = P.cbb(zpool, RF, N0 + HMON, kw, rng)
            f = {"bpi_wlc": bpi_wlc(X, [0.05])[0.05], "sn_block_thm1": sn_block(X, [0.05])[0.05],
                 "eprocess_orig": eproc_orig(X, [0.05])[0.05]}
            for k, v in f.items():
                h = fa_rate(v); acc[k].append(h); perseed[k].append(round(float(h.mean()), 4))
        FA[name] = {k: {"rate_wilson95": wil(np.concatenate(v)), "per_seed": perseed[k]} for k, v in acc.items()}
        print("FA", name, {k: FA[name][k]["rate_wilson95"] for k in acc}, round(time.time() - t0), flush=True)
    res["false_alarm_720h"] = FA
    # ---- delay vs Tartakovsky bound
    alphas = [0.05, 1e-3, 1e-6, 1e-12]
    NU = 240; POST = 3000; RD = 200
    DL = {}
    lrv_b1 = float(np.mean([zpool[i:i + 168].sum() ** 2 / 168 for i in range(0, len(zpool) - 167, 24)]))
    dcases = [(k, nulls[k]) for k in ("AR1_0", "AR1_0.5", "AR1_0.8", "AR1_0.95", "ARMA11_0.8_+0.4")] + [("B1_boot72", 72)]
    for name, kw in dcases:
        sig2 = lrv_arma(**kw) if isinstance(kw, dict) else lrv_b1
        S = np.sqrt(2 * sig2 * np.log(1 / 0.05) / 60.0)  # D*(0.05) = 60 h
        I = S ** 2 / (2 * sig2)
        rec = {"S": round(float(S), 4), "sigma2_LR": round(float(sig2), 3), "I": round(float(I), 5),
               "Dstar": {str(a): round(float(np.log(1 / a) / I), 1) for a in alphas}}
        acc = {m: {a: [] for a in alphas} for m in ("bpi_wlc", "sn_block_thm1", "eprocess_orig")}
        for s in SEEDS:
            rng = np.random.default_rng(s + 1000)
            n = N0 + NU + POST
            X = gen_arma(rng, RD, n, **kw) if isinstance(kw, dict) else P.cbb(zpool, RD, n, kw, rng)
            X[:, N0 + NU:] += S
            outs = {"bpi_wlc": bpi_wlc(X, alphas), "sn_block_thm1": sn_block(X, alphas), "eprocess_orig": eproc_orig(X, alphas)}
            for m, o in outs.items():
                for a in alphas:
                    acc[m][a].append(o[a])
        for m in acc:
            rec[m] = {}
            for a in alphas:
                f = np.concatenate(acc[m][a])
                pre = (f >= 0) & (f < NU)
                post = f >= NU
                d = (f[post] - NU).astype(float)
                rec[m][str(a)] = {"pre_change_alarm": round(float(pre.mean()), 4),
                                  "miss": round(float((f < 0).mean()), 4),
                                  "mean_delay": round(float(d.mean()), 1) if len(d) else None,
                                  "ratio_to_Dstar": round(float(d.mean() / (np.log(1 / a) / I)), 3) if len(d) else None,
                                  "ratio_predicted_overhead": (round(float((1 + KAPPA) * (1 + G / B) * (np.log((int(np.ceil(HMON / (B + G))) + WM) / a) + np.log(len(TH))) / np.log(1 / a)), 3) if m == "bpi_wlc" else None),
                                  "ci95_delay": [round(float(d.mean() - 1.96 * d.std(ddof=1) / np.sqrt(len(d))), 1),
                                                 round(float(d.mean() + 1.96 * d.std(ddof=1) / np.sqrt(len(d))), 1)] if len(d) > 2 else None}
        DL[name] = rec
        print("DELAY", name, {m: {a: (rec[m][str(a)]["ratio_to_Dstar"], rec[m][str(a)]["pre_change_alarm"]) for a in alphas} for m in acc}, round(time.time() - t0), flush=True)
    res["delay"] = DL
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(outp, "w", encoding="utf-8") as fo:
        json.dump(res, fo, indent=1)


if __name__ == "__main__":
    main()
