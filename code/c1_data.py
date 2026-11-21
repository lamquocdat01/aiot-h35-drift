"""H35 S1-1c data: UCI Air Quality (4 referenced channels), UCI Gas Drift #270 (batch level), Intel Lab (20 motes).
Rules exactly as PREREG_H35 v1 + Amendments 1-2."""
import os, gzip, json
import numpy as np
import pandas as pd

DS = os.environ.get("THS_DATASETS", r"D:\THS_Datasets")
AQ_CSV = os.path.join(DS, "UCI_AirQuality", "AirQualityUCI.csv")
GAS_DIR = os.path.join(DS, "UCI_GasDrift", "conc")
INTEL_DIR = os.path.join(DS, "IntelLab")
AQ_PAIRS = [("PT08.S1(CO)", "CO(GT)"), ("PT08.S2(NMHC)", "C6H6(GT)"), ("PT08.S3(NOx)", "NOx(GT)"), ("PT08.S4(NO2)", "NO2(GT)")]


def _fit(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 24: return None
    A = np.c_[np.ones(m.sum()), x[m]]
    return np.linalg.lstsq(A, y[m], rcond=None)[0]


def load_aq():
    df = pd.read_csv(AQ_CSV, sep=";", decimal=",")
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")].dropna(how="all")
    df.index = pd.to_datetime(df["Date"] + " " + df["Time"], format="%d/%m/%Y %H.%M.%S")
    df = df.replace(-200, np.nan)
    return df[~df.index.duplicated()].asfreq("h")


def aq_channels():
    """per channel: pool (standardized oos March residual, NaN dropped), full residual series (March fit,
    standardized by oos sd, NaN kept), index, and meta."""
    df = load_aq(); idx = df.index
    march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
    out = {}
    for sx, sy in AQ_PAIRS:
        x = df[sx].to_numpy(float); y = df[sy].to_numpy(float)
        roos = np.full(len(x), np.nan)
        for d in np.unique(days[march]):
            te = march & (days == d); tr = march & (days != d)
            c = _fit(x[tr], y[tr])
            if c is not None: roos[te] = c[0] + c[1] * x[te] - y[te]
        c0 = _fit(x[march], y[march])
        pool = roos[march & np.isfinite(roos)]
        mu, sd = pool.mean(), pool.std()
        full = (c0[0] + c0[1] * x - y) / sd
        out[sy] = {"pool": (pool - mu) / sd, "full": full, "index": idx, "march_valid": int(np.isfinite(roos[march]).sum()),
                   "coef": c0.tolist(), "sd": float(sd)}
    return out


# ---------------------------------------------------------------- Gas Drift #270 (batch level only)
def load_gas():
    rows = []
    for b in range(1, 11):
        with open(os.path.join(GAS_DIR, f"batch{b}.dat")) as f:
            for line in f:
                tok = line.split()
                if not tok: continue
                g, c = tok[0].split(";")
                feat = np.zeros(128)
                for t in tok[1:]:
                    k, v = t.split(":"); feat[int(k) - 1] = float(v)
                rows.append((b, int(g), float(c), feat))
    B = np.array([r[0] for r in rows]); G = np.array([r[1] for r in rows]); C = np.array([r[2] for r in rows])
    F = np.stack([r[3] for r in rows])
    return B, G, C, F


def gas_residuals(seed=20261121, ridge=1e-3):
    """residual r = c_hat - c, per (sensor, gas) ridge regression on the sensor's 8 features fit on batch 1;
    5-fold out-of-sample in batch 1; normalized by batch-1 oos sd of (sensor, gas). returns R (n,16)."""
    B, G, C, F = load_gas()
    n = len(B); R = np.full((n, 16), np.nan)
    rng = np.random.default_rng(seed)
    for s in range(16):
        X = F[:, 8 * s: 8 * s + 8]
        for g in np.unique(G):
            m1 = (B == 1) & (G == g); mall = (G == g)
            if m1.sum() < 10: continue
            mu, sd = X[m1].mean(0), X[m1].std(0) + 1e-12
            Z = (X - mu) / sd
            def fit(idx):
                A = np.c_[np.ones(idx.sum()), Z[idx]]
                return np.linalg.solve(A.T @ A + ridge * idx.sum() * np.eye(9), A.T @ C[idx])
            pred = lambda w, idx: np.c_[np.ones(idx.sum()), Z[idx]] @ w
            ids = np.where(m1)[0]; folds = rng.permutation(len(ids)) % 5
            roos = np.full(n, np.nan)
            for k in range(5):
                te = np.zeros(n, bool); te[ids[folds == k]] = True
                tr = m1 & ~te
                w = fit(tr); roos[te] = pred(w, te) - C[te]
            sdg = np.nanstd(roos[m1])
            w = fit(m1)
            other = mall & (B != 1)
            R[m1, s] = roos[m1] / sdg
            R[other, s] = (pred(w, other) - C[other]) / sdg
    return B, G, R


# ---------------------------------------------------------------- Intel Lab
def intel_residuals():
    df = pd.read_csv(os.path.join(INTEL_DIR, "data.txt.gz"), sep=r"\s+", header=None,
                     names=["date", "time", "epoch", "mote", "temp", "hum", "light", "volt"], on_bad_lines="skip",
                     engine="python")
    df = df.dropna(subset=["date", "time", "mote", "temp"])
    df["ts"] = pd.to_datetime(df["date"] + " " + df["time"], errors="coerce", format="mixed")
    df = df.dropna(subset=["ts"]); df["mote"] = df["mote"].astype(int)
    df.loc[(df.temp < -10) | (df.temp > 60), "temp"] = np.nan
    t0 = pd.Timestamp("2004-02-28 00:00:00")
    df = df[df.ts >= t0]
    df["bin"] = ((df.ts - t0).dt.total_seconds() // 600).astype(int)
    piv = df.pivot_table(index="bin", columns="mote", values="temp", aggfunc="mean")
    piv = piv.reindex(range(int(piv.index.max()) + 1))
    locs = pd.read_csv(os.path.join(INTEL_DIR, "mote_locs.txt"), sep=r"\s+", header=None, names=["mote", "x", "y"])
    locs = locs[locs.mote.isin(piv.columns)]
    first7 = piv.iloc[: 7 * 144]
    cnt = first7.notna().sum().reindex(locs.mote).fillna(0)
    order = sorted(cnt.index, key=lambda m: (-cnt[m], m))  # Amendment 1(4)
    chosen = order[:20]
    P = locs.set_index("mote")[["x", "y"]]
    res = {}
    for m in chosen:
        d = np.hypot(P.x - P.loc[m, "x"], P.y - P.loc[m, "y"]); d = d.drop(m)
        nb = list(d.sort_values().index[:5])
        nbv = piv[nb]
        med = nbv.median(axis=1, skipna=True).where(nbv.notna().sum(axis=1) >= 3)
        r = (piv[m] - med).to_numpy(float)
        pool = r[: 7 * 144]; pool = pool[np.isfinite(pool)]
        mu, sd = pool.mean(), pool.std()
        res[int(m)] = {"pool": (pool - mu) / sd, "full": (r - mu) / sd, "n_valid_first7": int(cnt[m]), "neighbors": [int(v) for v in nb]}
    return res
