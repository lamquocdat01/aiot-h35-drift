"""H35 S1-1a Pilot B: replay UCI Air Quality (CO(GT) ~ PT08.S1(CO)).
Part 1 (PA3): real alpha of 6 monitors on null residual (in-tolerance by construction) that keeps the
  real serial dependence: circular block bootstrap (block 24 h) of out-of-sample March-2004 residual,
  horizon 720 h; iid-bootstrap control. Plus sliding 168-h windows on the real March residual.
Part 2 (PA6): alarm-triggered recalibration over 12 months; refit on pre-alarm window (selection data)
  vs post-alarm window (correction rule) vs fixed 30-day schedule vs none; post-recal bias + MAE.
Usage: python pilot_B.py <AirQualityUCI.csv> <out.json>
"""
import json, sys, time
import numpy as np
import pandas as pd

ALPHA = 0.05
SEEDS = [20261121 + k for k in range(5)]
H = 720            # null-audit horizon (h)
BPS = 400          # bootstrap series per seed  (5 x 400 = 2000)
BLK = 24
W = 168            # recal window (h)
HEV = 168          # eval window after recal (h)
LAMS = np.array([0.1, 0.2, 0.35, 0.5, 0.75, 1.0])
LAMS2 = np.concatenate([LAMS, -LAMS])
EPS_GRID = np.linspace(0.025, 0.975, 20)


def load(path):
    df = pd.read_csv(path, sep=";", decimal=",")
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")].dropna(how="all")
    ts = pd.to_datetime(df["Date"] + " " + df["Time"], format="%d/%m/%Y %H.%M.%S")
    df.index = ts
    df = df.replace(-200, np.nan)
    df = df[~df.index.duplicated()].asfreq("h")
    return df


def fit(s1, co):
    m = np.isfinite(s1) & np.isfinite(co)
    if m.sum() < 24:
        return None
    A = np.c_[np.ones(m.sum()), s1[m]]
    coef, *_ = np.linalg.lstsq(A, co[m], rcond=None)
    return coef


def pred(coef, s1):
    return coef[0] + coef[1] * s1


# ---------------- monitors on standardized null series X (B, H), sigma = 1 assumed -------------
def mon_eprocess(X):
    B, n = X.shape
    logw = np.zeros((B, len(LAMS2)))
    first = np.full(B, -1)
    thr = np.log(1 / ALPHA)
    for t in range(n):
        logw += X[:, t:t + 1] * LAMS2 - LAMS2 ** 2 / 2
        le = np.logaddexp.reduce(logw, axis=1) - np.log(len(LAMS2))
        hit = (le >= thr) & (first < 0)
        first[hit] = t
    return first


def cusum_stat(X, k=0.5):
    B, n = X.shape
    sp = np.zeros(B); sm = np.zeros(B); out = np.zeros((B, n))
    for t in range(n):
        sp = np.maximum(0, sp + X[:, t] - k); sm = np.maximum(0, sm - X[:, t] - k)
        out[:, t] = np.maximum(sp, sm)
    return out


def ph_stat(X, delta=0.25):
    B, n = X.shape
    cs = np.cumsum(X, axis=1)
    mean = cs / np.arange(1, n + 1)
    mu = np.cumsum(X - mean - delta, axis=1); md = np.cumsum(X - mean + delta, axis=1)
    up = mu - np.minimum.accumulate(np.minimum(mu, 0), axis=1)
    dn = np.maximum.accumulate(np.maximum(md, 0), axis=1) - md
    return np.maximum(up, dn)


def first_cross(S, h):
    hit = S >= h
    return np.where(hit.any(1), hit.argmax(1), -1)


def mon_repeated_p(X, tmin=10):
    n = X.shape[1]
    z = np.cumsum(X, axis=1) / np.sqrt(np.arange(1, n + 1))
    hit = np.abs(z) >= 1.959964
    hit[:, :tmin] = False
    return np.where(hit.any(1), hit.argmax(1), -1)


def mon_adwin(X, delta=ALPHA, mmin=10, every=4):
    B, n = X.shape
    cs = np.cumsum(X, axis=1); cs2 = np.cumsum(X ** 2, axis=1)
    first = np.full(B, -1)
    for t in range(2 * mmin, n, every):
        nt = t + 1
        n0 = np.arange(mmin, nt - mmin + 1)
        n1 = nt - n0
        S0 = cs[:, n0 - 1]; S1 = cs[:, t:t + 1] - S0
        diff = np.abs(S0 / n0 - S1 / n1)
        var = np.maximum(cs2[:, t] / nt - (cs[:, t] / nt) ** 2, 1e-9)[:, None]
        m = 1.0 / (1.0 / n0 + 1.0 / n1)
        dp = delta / np.log(nt)
        eps = np.sqrt(2.0 / m * var * np.log(2.0 / dp)) + 2.0 / (3 * m) * np.log(2.0 / dp)
        hit = (diff > eps).any(1) & (first < 0)
        first[hit] = t
    return first


def ctm_one(A, rng):
    """simple-mixture power martingale on smoothed conformal p-values; score A (B,n), big = strange."""
    B, n = A.shape
    logM = np.zeros((B, len(EPS_GRID)))
    first = np.full(B, -1)
    thr = np.log(2 / ALPHA)
    for t in range(n):
        a = A[:, t:t + 1]
        past = A[:, :t + 1]
        gt = (past > a).sum(1); eq = (past == a).sum(1)
        p = (gt + rng.random(B) * eq) / (t + 1)
        p = np.clip(p, 1e-12, 1)
        logM += np.log(EPS_GRID) + (EPS_GRID - 1) * np.log(p)[:, None]
        le = np.logaddexp.reduce(logM, axis=1) - np.log(len(EPS_GRID))
        hit = (le >= thr) & (first < 0)
        first[hit] = t
    return first


def mon_ctm(X, rng):
    f1 = ctm_one(X, rng); f2 = ctm_one(-X, rng)
    f = np.where((f1 >= 0) & ((f2 < 0) | (f1 <= f2)), f1, f2)
    return f


def calib_thresholds(rng, n=H, B=20000):
    Z = rng.standard_normal((B, n))
    hc = float(np.quantile(cusum_stat(Z).max(1), 1 - ALPHA))
    hp = float(np.quantile(ph_stat(Z).max(1), 1 - ALPHA))
    return hc, hp


def run_monitors(X, rng, hc, hp):
    return {"cusum": first_cross(cusum_stat(X), hc), "page_hinkley": first_cross(ph_stat(X), hp),
            "adwin": mon_adwin(X), "repeated_p": mon_repeated_p(X), "eprocess": mon_eprocess(X),
            "ctm": mon_ctm(X, rng)}


def cbb(z, B, n, blk, rng):
    L = len(z)
    nb = int(np.ceil(n / blk))
    st = rng.integers(0, L, size=(B, nb))
    idx = (st[..., None] + np.arange(blk)) % L
    return z[idx.reshape(B, -1)[:, :n]]


def wilson(k, n, zc=1.96):
    p = k / n; d = 1 + zc ** 2 / n
    c = (p + zc ** 2 / (2 * n)) / d; h = zc * np.sqrt(p * (1 - p) / n + zc ** 2 / (4 * n ** 2)) / d
    return [round(p, 4), round(c - h, 4), round(c + h, 4)]


# ---------------- Part 2 replay ----------------
class EProc:
    def __init__(self): self.reset()
    def reset(self): self.logw = np.zeros(len(LAMS2))
    def step(self, x):
        self.logw += x * LAMS2 - LAMS2 ** 2 / 2
        return np.logaddexp.reduce(self.logw) - np.log(len(LAMS2)) >= np.log(1 / ALPHA)


class Cusum:
    def __init__(self, h): self.h = h; self.reset()
    def reset(self): self.sp = 0.0; self.sm = 0.0
    def step(self, x):
        self.sp = max(0, self.sp + x - 0.5); self.sm = max(0, self.sm - x - 0.5)
        return max(self.sp, self.sm) >= self.h


def replay(s1, co, coef0, sigma, start, policy, monitor=None):
    """returns predictions, list of events (alarm_t, refit_end_t, sign)."""
    n = len(s1)
    coef = coef0.copy()
    yhat = np.full(n, np.nan)
    events = []
    pending = None  # post-window: time at which refit happens
    t = start
    last_sched = start
    while t < n:
        yhat[t] = pred(coef, s1[t])
        if policy == "sched" and t - last_sched >= 720:
            c = fit(s1[t - W:t], co[t - W:t]); r_old = pred(coef, s1[t - W:t]) - co[t - W:t]
            sg = np.sign(np.nanmean(r_old)); coef = c if c is not None else coef; last_sched = t
            events.append((t, t, sg))
        elif policy.startswith("alarm"):
            if pending is not None:
                if t >= pending[0]:
                    tau, sg = pending[1], pending[2]
                    c = fit(s1[tau:t], co[tau:t]); coef = c if c is not None else coef; monitor.reset(); pending = None
                    events.append((tau, t, sg))
            else:
                r = yhat[t] - co[t]
                if np.isfinite(r) and monitor.step(r / sigma):
                    lo = max(start, t + 1 - W)
                    r_old = pred(coef, s1[lo:t + 1]) - co[lo:t + 1]
                    sg = np.sign(np.nanmean(r_old))
                    if policy == "alarm_pre":
                        c = fit(s1[lo:t + 1], co[lo:t + 1]); coef = c if c is not None else coef; monitor.reset()
                        events.append((t, t, sg))
                    else:
                        pending = (t + W, t, sg)
        t += 1
    return yhat, events


def event_bias(yhat, co, events, n):
    vals = []
    for k, (tau, tre, sg) in enumerate(events):
        a = tau + W; b = a + HEV
        if k + 1 < len(events): b = min(b, events[k + 1][1])
        if b - a < 24 or a >= n: continue
        r = yhat[a:b] - co[a:b]
        if np.isfinite(r).sum() < 12: continue
        vals.append((sg * np.nanmean(r), np.nanmean(r)))
    return np.array(vals) if vals else np.zeros((0, 2))


def block_boot_mean(v, rng, blk=3, B=2000):
    n = len(v)
    if n < 3: return [float(np.mean(v)) if n else None, None, None]
    nb = int(np.ceil(n / blk))
    st = rng.integers(0, n - blk + 1, size=(B, nb))
    idx = (st[..., None] + np.arange(blk)).reshape(B, -1)[:, :n]
    m = v[idx].mean(1)
    return [round(float(v.mean()), 4), round(float(np.quantile(m, .025)), 4), round(float(np.quantile(m, .975)), 4)]


def main():
    t0 = time.time()
    df = load(sys.argv[1])
    s1 = df["PT08.S1(CO)"].to_numpy(float); co = df["CO(GT)"].to_numpy(float)
    idx = df.index
    march = (idx >= "2004-03-10") & (idx < "2004-04-01")
    days = idx.normalize()
    # out-of-sample March residual (leave-one-day-out)
    roos = np.full(len(s1), np.nan)
    for d in np.unique(days[march]):
        te = march & (days == d); tr = march & (days != d)
        c = fit(s1[tr], co[tr]); roos[te] = pred(c, s1[te]) - co[te]
    coef0 = fit(s1[march], co[march])
    rin = pred(coef0, s1[march]) - co[march]
    sigma = float(np.nanstd(rin, ddof=2))
    z = roos[march]; z = z[np.isfinite(z)]
    z = (z - z.mean()) / z.std()
    acf = [round(float(np.corrcoef(z[:-k], z[k:])[0, 1]), 3) for k in (1, 2, 6, 12, 24)]
    out = {"config": dict(alpha=ALPHA, seeds=SEEDS, H=H, boot_series_per_seed=BPS, block=BLK, W=W, HEV=HEV,
                          calib="2004-03-10..2004-03-31, CO(GT)=a+b*PT08.S1", missing="-200 -> NaN, skipped"),
           "data": {"hours": int(len(df)), "co_missing_frac": round(float(np.isnan(co).mean()), 4),
                    "s1_missing_frac": round(float(np.isnan(s1).mean()), 4),
                    "march_valid_pairs": int(np.isfinite(rin).sum()), "coef0": coef0.round(6).tolist(),
                    "sigma_march": round(sigma, 4), "acf_oos_march_lags_1_2_6_12_24": acf}}
    rng0 = np.random.default_rng(SEEDS[0])
    hc, hp = calib_thresholds(rng0)
    out["thresholds_iid_calibrated"] = {"cusum_h(k=0.5)": round(hc, 3), "page_hinkley_lambda(delta=0.25)": round(hp, 3),
                                        "eprocess": "1/alpha", "ctm": "2/alpha per side", "adwin": "delta=alpha",
                                        "repeated_p": "|z|>=1.96 each hour"}
    # Part 1
    part1 = {}
    for mode in ("block24", "iid"):
        hits = {}
        per_seed = {}
        for s in SEEDS:
            rng = np.random.default_rng(s)
            X = cbb(z, BPS, H, BLK if mode == "block24" else 1, rng)
            f = run_monitors(X, rng, hc, hp)
            for k, v in f.items():
                hits.setdefault(k, []).append(v >= 0)
                per_seed.setdefault(k, []).append(round(float((v >= 0).mean()), 4))
        part1[mode] = {k: {"rate_wilson95": wilson(int(np.concatenate(v).sum()), len(np.concatenate(v))),
                           "per_seed": per_seed[k]} for k, v in hits.items()}
        print(mode, {k: part1[mode][k]["rate_wilson95"] for k in part1[mode]}, flush=True)
    # sliding real windows on actual March oos residual (168 h, step 12 h)
    zr = roos[march] / np.nanstd(roos[march])
    zr = np.where(np.isfinite(zr), zr, 0.0)  # NaN hour -> zero increment approx (documented)
    L = 168; starts = np.arange(0, len(zr) - L, 12)
    Xw = np.stack([zr[s:s + L] for s in starts])
    hcw, hpw = calib_thresholds(np.random.default_rng(SEEDS[1]), n=L)
    fw = {"cusum": first_cross(cusum_stat(Xw), hcw), "page_hinkley": first_cross(ph_stat(Xw), hpw),
          "adwin": mon_adwin(Xw), "repeated_p": mon_repeated_p(Xw), "eprocess": mon_eprocess(Xw),
          "ctm": mon_ctm(Xw, np.random.default_rng(SEEDS[2]))}
    part1["real_march_windows_168h"] = {"n_windows(overlapping)": int(len(starts)),
                                        **{k: round(float((v >= 0).mean()), 4) for k, v in fw.items()}}
    out["part1_null_audit"] = part1
    # Part 2 replay
    start = int(np.argmax(idx >= "2004-04-01"))
    n = len(s1)
    pol = {}
    variants = {"none": ("none", None), "sched30d": ("sched", None),
                "eproc_pre": ("alarm_pre", EProc()), "eproc_post": ("alarm_post", EProc()),
                "cusum_pre": ("alarm_pre", Cusum(hc)), "cusum_post": ("alarm_post", Cusum(hc))}
    for name, (p, m) in variants.items():
        yhat, ev = replay(s1, co, coef0, sigma, start, p, m)
        err = yhat[start:] - co[start:]
        ok = np.isfinite(err)
        eb = event_bias(yhat, co, ev, n)
        rec = {"n_recal": len(ev), "MAE": round(float(np.abs(err[ok]).mean()), 4),
               "mean_resid": round(float(err[ok].mean()), 4), "n_eval_events": int(len(eb))}
        if len(eb):
            rec["signed_bias_seeds"] = [block_boot_mean(eb[:, 0], np.random.default_rng(s)) for s in SEEDS]
            rec["raw_bias_seeds"] = [block_boot_mean(eb[:, 1], np.random.default_rng(s)) for s in SEEDS]
        # MAE difference vs sched via moving-block bootstrap (168 h)
        pol[name] = rec
        pol[name]["_err"] = err
        print(name, {k: v for k, v in rec.items() if k != "_err"}, flush=True)
    base = pol["sched30d"]["_err"]
    for name in pol:
        e = pol[name].pop("_err")
        if name == "sched30d": continue
        d = np.abs(e) - np.abs(base); okd = np.isfinite(d); d = np.where(okd, d, 0.0)
        cnt = okd.astype(float)
        nb = int(np.ceil(len(d) / 168)); res = []
        for s in SEEDS:
            rng = np.random.default_rng(s)
            st = rng.integers(0, len(d) - 168, size=(1000, nb))
            ii = (st[..., None] + np.arange(168)).reshape(1000, -1)[:, :len(d)]
            bm = d[ii].sum(1) / cnt[ii].sum(1)
            res.append([round(float(d.sum() / cnt.sum()), 4), round(float(np.quantile(bm, .025)), 4),
                        round(float(np.quantile(bm, .975)), 4)])
        pol[name]["dMAE_vs_sched_seeds"] = res
    out["part2_replay"] = pol
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("done", out["runtime_s"])


if __name__ == "__main__":
    main()
