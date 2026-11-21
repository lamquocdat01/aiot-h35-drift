"""H35 S1-1b: Monte Carlo checks for Lemma 1 / Theorem 1 (self-normalized blocked e-process).
E-process: block sums S_k over blocks of length b separated by gaps g; for each direction s in {+1,-1}
and lambda in grid c/sigma0 (sigma0 = a-priori scale, only affects power):
   M_n = prod_k exp(lam*Y_k - lam^2*Y_k^2/2),  Y_k = s*S_k - b*delta0.
E_n = average of the 12 products; alarm if E_n >= 1/alpha.  Valid (Thm 1) if block sums are independent
across blocks (gap >= m for m-dependent data) and symmetric about a centre within the tolerance.
Usage: python -I s1b_thm1_check.py <AirQualityUCI.csv> <out.json>
"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pilot_B as P

ALPHA = 0.05
SEEDS = [20261121 + k for k in range(5)]
REPS = 2000
T = 720
CGRID = np.array([0.1, 0.2, 0.35, 0.5, 0.75, 1.0])


def block_sums(X, b, g):
    B, n = X.shape
    step = b + g
    K = (n - b) // step + 1
    starts = np.arange(K) * step
    idx = starts[:, None] + np.arange(b)
    return X[:, idx].sum(2)  # (B, K)


def sn_eprocess_alarm(S, b, sigma0, delta0=0.0):
    lam = CGRID / sigma0
    logs = []
    for s in (1.0, -1.0):
        Y = s * S - b * delta0  # (B,K)
        inc = Y[..., None] * lam - (Y[..., None] ** 2) * lam ** 2 / 2  # (B,K,L)
        logs.append(np.cumsum(inc, axis=1))
    L = np.concatenate(logs, axis=2)
    logE = np.logaddexp.reduce(L, axis=2) - np.log(L.shape[2])
    return (logE >= np.log(1 / ALPHA)).any(1)


def lemma_grid():
    """Lemma 1 numeric: sup over u<=0, v of u-(u^2+v^2)/2+log cosh(v(1-u)) (two-point symmetric case)."""
    u = np.linspace(-20, 0, 2001)[:, None]; v = np.linspace(-40, 40, 8001)[None, :]
    a = v * (1 - u)
    lc = np.abs(a) + np.log1p(np.exp(-2 * np.abs(a))) - np.log(2)
    f = u - (u ** 2 + v ** 2) / 2 + lc
    i = np.unravel_index(np.argmax(f), f.shape)
    return {"max": float(f.max()), "argmax_u": float(u[i[0], 0]), "argmax_v": float(v[0, i[1]])}


def lemma_mc(rng):
    """E exp(lam Y - lam^2 Y^2/2) for several symmetric Z and shifts s<=0, and asymmetric Z (counterexample)."""
    out = {}
    n = 2_000_000
    Zs = {"gauss": rng.standard_normal(n), "rademacher": rng.choice([-1.0, 1.0], n),
          "t3": rng.standard_t(3, n), "laplace": rng.laplace(size=n)}
    for name, Z in Zs.items():
        best = -1
        for lam in (0.3, 1.0, 2.0, 4.0):
            for s in (0.0, -0.2, -0.5, -1.0, -3.0):
                Y = s + Z
                best = max(best, float(np.mean(np.exp(lam * Y - lam ** 2 * Y ** 2 / 2))))
        out[f"sym_{name}_max_over_lam_s"] = round(best, 4)
    # asymmetric mean-zero: centred exponential, lam chosen adversarially
    Z = rng.exponential(size=n) - 1.0
    vals = {lam: float(np.mean(np.exp(lam * (-Z) - lam ** 2 * Z ** 2 / 2))) for lam in (0.3, 0.6, 1.0, 1.5, 2.0)}
    out["asym_neg_centred_exp_by_lam"] = {str(k): round(v, 4) for k, v in vals.items()}
    # mean-zero two-point asymmetric: Z = 9 w.p. 0.1, -1 w.p. 0.9 ; direction -Z
    Z = np.where(rng.random(n) < 0.1, 9.0, -1.0)
    vals = {lam: float(np.mean(np.exp(lam * (-Z) - lam ** 2 * Z ** 2 / 2))) for lam in (0.3, 0.6, 1.0, 1.5)}
    out["asym_two_point_by_lam"] = {str(k): round(v, 4) for k, v in vals.items()}
    return out


def gen(kind, rng, B, n):
    if kind.startswith("ma"):
        m = int(kind[2:].split("_")[0])
        e = rng.standard_normal((B, n + m))
        w = np.ones(m + 1) / np.sqrt(m + 1)
        X = np.stack([np.convolve(r, w, mode="valid") for r in e])
        return X
    if kind.startswith("ar"):
        phi = float(kind[2:])
        X = np.zeros((B, n)); x = rng.standard_normal(B) / np.sqrt(1 - phi ** 2)
        for t in range(n):
            x = phi * x + rng.standard_normal(B); X[:, t] = x
        return X
    if kind == "skew_iid":
        return rng.exponential(size=(B, n)) - 1.0
    if kind == "skew_ma23":
        e = rng.exponential(size=(B, n + 23)) - 1.0
        return np.stack([np.convolve(r, np.ones(24), mode="valid") for r in e]) / np.sqrt(24)
    if kind == "longmem":  # sum of AR(0.999) + noise -> very slow mixing
        X = np.zeros((B, n)); x = rng.standard_normal(B) * 3
        for t in range(n):
            x = 0.999 * x + 0.13 * rng.standard_normal(B); X[:, t] = x + rng.standard_normal(B)
        return X
    raise ValueError(kind)


def rate(h):
    return P.wilson(int(h.sum()), len(h))


def main():
    res = {"config": dict(alpha=ALPHA, seeds=SEEDS, reps_per_seed=REPS, T=T, cgrid=CGRID.tolist())}
    res["lemma1_grid_two_point"] = lemma_grid()
    res["lemma1_mc"] = lemma_mc(np.random.default_rng(SEEDS[0]))
    print(res["lemma1_grid_two_point"], res["lemma1_mc"], flush=True)
    cases = [  # (kind, b, g, delta0, centre_shift_per_obs)
        ("ma23", 24, 24, 0, 0), ("ma23", 24, 0, 0, 0), ("ma23", 24, 12, 0, 0),
        ("ar0.8", 24, 24, 0, 0), ("ar0.8", 24, 0, 0, 0), ("ar0.95", 24, 24, 0, 0), ("ar0.95", 24, 0, 0, 0),
        ("skew_iid", 1, 0, 0, 0), ("skew_iid", 24, 0, 0, 0), ("skew_ma23", 24, 24, 0, 0),
        ("longmem", 24, 24, 0, 0),
        ("ma23", 24, 24, 0.1, 0.1),   # centre at tolerance boundary
        ("ma23", 24, 24, 0.1, 0.0),   # centre inside tolerance
        ("ma23", 24, 24, 0.0, 0.1),   # centre outside tolerance (H1-like, not a null)
    ]
    out = {}
    for kind, b, g, d0, shift in cases:
        hits, ps = [], []
        for s in SEEDS:
            rng = np.random.default_rng(s)
            X = gen(kind, rng, REPS, T) + shift
            sig0 = np.sqrt(b)  # a-priori scale (iid unit variance); power only
            h = sn_eprocess_alarm(block_sums(X, b, g), b, sig0, d0)
            hits.append(h); ps.append(round(float(h.mean()), 4))
        key = f"{kind}|b{b}|g{g}|d0={d0}|shift={shift}"
        out[key] = {"rate_wilson95": rate(np.concatenate(hits)), "per_seed": ps}
        print(key, out[key], flush=True)
    res["synthetic"] = out
    # UCI B1 null pool
    df = P.load(sys.argv[1])
    s1 = df["PT08.S1(CO)"].to_numpy(float); co = df["CO(GT)"].to_numpy(float); idx = df.index
    march = (idx >= "2004-03-10") & (idx < "2004-04-01"); days = idx.normalize()
    roos = np.full(len(s1), np.nan)
    for d in np.unique(days[march]):
        te = march & (days == d); tr = march & (days != d)
        c = P.fit(s1[tr], co[tr]); roos[te] = P.pred(c, s1[te]) - co[te]
    z = roos[march & np.isfinite(roos)]; z = (z - z.mean()) / z.std()
    sdb = float(np.array([z[i:i + 24].sum() for i in range(0, len(z) - 23, 24)]).std(ddof=1))
    uci = {}
    for bb in (24, 72, 168):
        for b, g in ((24, 0), (24, 24), (48, 0), (48, 48)):
            hits = []
            for s in SEEDS:
                rng = np.random.default_rng(s)
                X = P.cbb(z, 400, T, bb, rng)
                hits.append(sn_eprocess_alarm(block_sums(X, b, g), b, sdb * np.sqrt(b / 24)))
            key = f"bootblk{bb}|b{b}|g{g}"
            uci[key] = rate(np.concatenate(hits)); print(key, uci[key], flush=True)
    res["uci_b1_selfnormalized"] = uci
    with open(sys.argv[2], "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
