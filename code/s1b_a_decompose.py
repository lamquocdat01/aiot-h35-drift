"""H35 S1-1b(a): decompose B1 false alarms of the e-process (block-24 bootstrap null, H=720 h, alpha=0.05).
Variants: raw (pilot), (1) diurnal cycle removed, (2) block-mean e-process (b=24 h, block sd from calibration
month), (3) online-whitened AR(1) with estimated phi (causal, burn-in 48 h), + combos.
Usage: python -I s1b_a_decompose.py <AirQualityUCI.csv> <out.json>
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P

ALPHA, SEEDS, H, BPS = P.ALPHA, P.SEEDS, P.H, P.BPS
LAMS2 = P.LAMS2


def eproc_first(X, burn=0):
    B, n = X.shape
    logw = np.zeros((B, len(LAMS2))); first = np.full(B, -1); thr = np.log(1 / ALPHA)
    for t in range(n):
        x = np.nan_to_num(X[:, t:t + 1])
        if t < burn:
            continue
        logw += x * LAMS2 - LAMS2 ** 2 / 2
        le = np.logaddexp.reduce(logw, axis=1) - np.log(len(LAMS2))
        hit = (le >= thr) & (first < 0); first[hit] = t
    return first


def whiten_online(X, burn=48):
    """x_t = (z_t - phi_hat_{t-1} z_{t-1}) / s_hat_{t-1}; phi_hat, s_hat from data up to t-1 (causal)."""
    B, n = X.shape
    out = np.full((B, n), np.nan)
    sxy = np.zeros(B); sxx = np.zeros(B)
    for t in range(1, n):
        if t >= burn:
            ph = np.clip(sxy / np.maximum(sxx, 1e-9), -0.99, 0.99)
            # innovation variance under stationarity: var(z)(1-phi^2), var(z) from past
            vz = X[:, :t].var(1)
            s = np.sqrt(np.maximum(vz * (1 - ph ** 2), 1e-6))
            out[:, t] = (X[:, t] - ph * X[:, t - 1]) / s
        sxy += X[:, t] * X[:, t - 1]; sxx += X[:, t - 1] ** 2
    return out


def block_means(X, b, sd_b):
    B, n = X.shape
    nb = n // b
    return X[:, :nb * b].reshape(B, nb, b).mean(2) / sd_b


def rate(first):
    return P.wilson(int((first >= 0).sum()), len(first))


def main():
    df = P.load(sys.argv[1])
    s1 = df["PT08.S1(CO)"].to_numpy(float); co = df["CO(GT)"].to_numpy(float); idx = df.index
    march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
    roos = np.full(len(s1), np.nan)
    for d in np.unique(days[march]):
        te = march & (days == d); tr = march & (days != d)
        c = P.fit(s1[tr], co[tr]); roos[te] = P.pred(c, s1[te]) - co[te]
    hr = idx.hour.to_numpy()
    m = march & np.isfinite(roos)
    z = roos[m]; h = hr[m]
    z = (z - z.mean()) / z.std()
    prof = np.array([z[h == k].mean() for k in range(24)])
    zds = z - prof[h]; zds = (zds - zds.mean()) / zds.std()
    acf = lambda v: [round(float(np.corrcoef(v[:-k], v[k:])[0, 1]), 3) for k in (1, 2, 6, 12, 24)]
    # AR(1) fit on calibration pool (reported only)
    phi_pool = float(np.sum(z[1:] * z[:-1]) / np.sum(z[:-1] ** 2))
    b = 24
    sdb_raw = float(np.array([z[i:i + b].mean() for i in range(0, len(z) - b + 1, b)]).std(ddof=1))
    sdb_ds = float(np.array([zds[i:i + b].mean() for i in range(0, len(zds) - b + 1, b)]).std(ddof=1))
    res = {"config": {"alpha": ALPHA, "seeds": SEEDS, "H": H, "series_per_seed": BPS, "boot_block": 24,
                      "block_mean_b": b, "whiten": "online causal AR(1), burn-in 48 h"},
           "pool": {"n": int(len(z)), "acf_raw": acf(z), "acf_deseason": acf(zds), "phi_pool_raw": round(phi_pool, 3),
                    "diurnal_profile_range_sd": round(float(prof.max() - prof.min()), 3),
                    "sd_blockmean24_raw": round(sdb_raw, 3), "sd_blockmean24_deseason": round(sdb_ds, 3),
                    "n_blocks_for_sd": int(len(z) // b)}}
    variants = {}
    for s in SEEDS:
        rng = np.random.default_rng(s)
        Xr = P.cbb(z, BPS, H, 24, rng)
        Xd = P.cbb(zds, BPS, H, 24, rng)
        Xiid = P.cbb(z, BPS, H, 1, rng)
        f = {"raw": eproc_first(Xr),
             "1_deseason": eproc_first(Xd),
             "2_block24_raw": eproc_first(block_means(Xr, b, sdb_raw)),
             "3_whiten_phihat_raw": eproc_first(whiten_online(Xr), burn=48),
             "ds+block24": eproc_first(block_means(Xd, b, sdb_ds)),
             "ds+whiten_phihat": eproc_first(whiten_online(Xd), burn=48),
             "iid_control_whiten_phihat": eproc_first(whiten_online(Xiid), burn=48)}
        for k, v in f.items():
            variants.setdefault(k, []).append(v)
    res["eprocess_false_alarm"] = {k: {"rate_wilson95": rate(np.concatenate(v)),
                                       "per_seed": [round(float((x >= 0).mean()), 4) for x in v]}
                                   for k, v in variants.items()}
    for k, v in res["eprocess_false_alarm"].items():
        print(k, v["rate_wilson95"], v["per_seed"])
    print(res["pool"])
    with open(sys.argv[2], "w", encoding="utf-8") as fo:
        json.dump(res, fo, indent=1)


if __name__ == "__main__":
    main()
