"""H35 Q-2 self-check: counterexamples for Theorem 3 (break one assumption at a time).
Usage: python q2_selfcheck.py <AirQualityUCI.csv> <out.json>"""
import json, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import q2_adaptive as Q, pilot_B as P

def ar1_only(gam, n_eff, b):          # misspecified AR(1) phi-hat, no Bartlett max, no inflation
    phi = np.clip(gam[:, 1] / np.maximum(gam[:, 0], 1e-9), -0.99, 0.99)
    return gam[:, 0] * (1 - phi ** 2) / (1 - phi) ** 2

def gen_skew(rng, R, n, phi=0.8):
    e = rng.exponential(size=(R, n + 500)) - 1.0
    x = np.zeros_like(e)
    for t in range(1, e.shape[1]): x[:, t] = phi * x[:, t - 1] + e[:, t]
    return x[:, -n:]

def gen_longmem(rng, R, n):
    x = np.zeros((R, n)); u = rng.standard_normal(R) * 3
    for t in range(n):
        u = 0.999 * u + 0.13 * rng.standard_normal(R); x[:, t] = u + rng.standard_normal(R)
    return x

def main():
    zpool = Q.b1_pool(sys.argv[1]); out = {}
    R = 400
    def fa(name, mk, **kw):
        hs = []
        for s in Q.SEEDS:
            rng = np.random.default_rng(s); X = mk(rng)
            f = Q.bpi_wlc(X, [0.05], **kw)[0.05]; hs.append((f >= 0) & (f <= Q.HMON))
        out[name] = P.wilson(int(np.concatenate(hs).sum()), R * len(Q.SEEDS)); print(name, out[name], flush=True)
    n = Q.N0 + Q.HMON
    ar95 = lambda rng: Q.gen_arma(rng, R, n, phi=(0.95,))
    b168 = lambda rng: P.cbb(zpool, R, n, 168, rng)
    # (a)-1 no union bound: A = 1/alpha
    fa("a1_A=1/alpha|AR0.95", ar95, logA_override=lambda a: np.log(1 / a))
    fa("a1_A=1/alpha|B1boot168", b168, logA_override=lambda a: np.log(1 / a))
    # (a)-2 misspecified AR(1) estimator, kappa=0, no Bartlett max
    orig = Q.s2_tilde; Q.s2_tilde = ar1_only
    fa("a2_AR1only_k0|B1boot168", b168); fa("a2_AR1only_k0|B1boot168_A=1/alpha", b168, logA_override=lambda a: np.log(1 / a))
    Q.s2_tilde = orig
    # (a)-3 non-Gaussian skewed innovations
    fa("a3_skew_exp_AR0.8", lambda rng: gen_skew(rng, R, n)); fa("a3_skew_exp_AR0.8_A=1/alpha", lambda rng: gen_skew(rng, R, n), logA_override=lambda a: np.log(1 / a))
    # (a)-4 long memory / slow mixing
    fa("a4_longmem", lambda rng: gen_longmem(rng, R, n)); fa("a4_longmem_A=1/alpha", lambda rng: gen_longmem(rng, R, n), logA_override=lambda a: np.log(1 / a))
    # (b) delay counterexamples on AR(0.8), D*(0.05)=60 h
    kw = dict(phi=(0.8,)); s2 = Q.lrv_arma(**kw); S = np.sqrt(2 * s2 * np.log(20) / 60); I = S ** 2 / (2 * s2)
    NU, POST, RD = 240, 3000, 200
    def dl(name, alpha=1e-6, **kwb):
        ds, miss = [], []
        for s in Q.SEEDS:
            rng = np.random.default_rng(s + 1000); X = Q.gen_arma(rng, RD, Q.N0 + NU + POST, **kw); X[:, Q.N0 + NU:] += S
            f = Q.bpi_wlc(X, [alpha], **kwb)[alpha]; ds.append(f[f >= NU] - NU); miss.append(f < 0)
        d = np.concatenate(ds).astype(float)
        out[name] = {"ratio": round(float(d.mean() / (np.log(1 / alpha) / I)), 3), "miss": round(float(np.concatenate(miss).mean()), 4)}
        print(name, out[name], flush=True)
    dl("b0_reference_alpha1e-6")
    dl("b1_nonlagged_estimator", lag=False)
    dl("b2_short_window_WM10", wm=10)
    Q.KAPPA = 1.0; dl("b3_kappa=1"); Q.KAPPA = 0.10
    json.dump(out, open(sys.argv[2], "w", encoding="utf-8"), indent=1)

if __name__ == "__main__":
    main()
