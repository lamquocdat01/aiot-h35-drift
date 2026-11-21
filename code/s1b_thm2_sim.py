"""H35 S1-1b: Monte Carlo for Theorem 2 (discovery-investing alpha allocation over sensor-epochs).
Same model as pilot_A (N sensors, AR(1) whitened with known phi, common factor, null>=Lmin, ramp drift).
Hypothesis k = k-th started sensor-epoch (global order, ties by sensor index). Level fixed at start sigma_k:
    alpha_k = alpha * gamma_k * (R(sigma_k) + 1),
    gamma_k = 1/(2N) for k<=N,  N/(2k(k-1)) for k>N   (sum = 1).
Alarm when e-process of the epoch >= 1/alpha_k.
Rules: iii (pilot epoch-Bonferroni), iv (Thm 2), and assumption-breaking variants:
  iv_nonpred : level recomputed every step with current R(t-1) (not fixed at epoch start)
  iv_const   : gamma_k = 1/N for all k (not summable)
  iv_raw     : Thm-2 levels but e-process on raw AR(1) residual (invalid e-process)
Metrics: lifetime FDR at fixed t and stopping times, E[sup_t FDP_t], delay by epoch j.
Usage: python -I s1b_thm2_sim.py <out.json> [N ...]
"""
import json, sys, time
import numpy as np

ALPHA = 0.05
SEEDS = [20261121 + k for k in range(5)]
REPS = 1000
T = 4000
PHI, RHO, LMIN, LEXTRA, SLOPE = 0.5, 0.3, 50, 150, 0.02
LAMS = np.array([0.1, 0.2, 0.35, 0.5, 0.75, 1.0])
JMAX = 12
FIXED_T = [1000, 4000]


def gamma(k, n, const=False):
    k = np.asarray(k, float)
    if const:
        return np.full_like(k, 1.0 / n)
    return np.where(k <= n, 1.0 / (2 * n), n / (2 * k * np.maximum(k - 1, 1)))


def run(rule, seed, n, reps=REPS):
    rng = np.random.default_rng(seed)
    K = len(LAMS)
    logw = np.zeros((reps, n, K))
    j = np.ones((reps, n), dtype=int)
    last = np.zeros((reps, n), dtype=int)
    onset = LMIN + rng.geometric(1.0 / LEXTRA, size=(reps, n))
    prev_r = np.zeros((reps, n))
    noise_prev = rng.standard_normal((reps, n)) / np.sqrt(1 - PHI ** 2)
    R = np.zeros(reps, int); F = np.zeros(reps, int)
    kidx = np.tile(np.arange(1, n + 1), (reps, 1)).astype(float)  # hypothesis index of current epoch
    nstarted = np.full(reps, n, dtype=float)
    const = rule == "iv_const"
    lvl = ALPHA * gamma(kidx, n, const) * 1.0  # R(0)+1 = 1
    supFDP = np.zeros(reps)
    fixed = {}
    tau_short = np.full(reps, -1); tauFDP = np.zeros(reps)
    dsum = np.zeros(JMAX + 1); dn = np.zeros(JMAX + 1)
    lam2 = LAMS ** 2 / 2
    for t in range(1, T + 1):
        c = rng.standard_normal((reps, 1)); u = rng.standard_normal((reps, n))
        noise = PHI * noise_prev + np.sqrt(RHO) * c + np.sqrt(1 - RHO) * u
        d = SLOPE * np.maximum(0, t - onset)
        r = d + noise
        x = r if rule in ("iv_raw", "iii_raw") else r - PHI * prev_r
        logw += x[..., None] * LAMS - lam2
        logE = np.logaddexp.reduce(logw, axis=2) - np.log(K)
        if rule in ("iii", "iii_raw"):
            thr = np.log(n / (ALPHA * (1.0 / (j * (j + 1)))))
        elif rule == "iii_harm":  # g_j = 1/j, not summable
            thr = np.log(n * j / ALPHA)
        elif rule == "iv_nonpred":
            thr = -np.log(ALPHA * gamma(kidx, n) * (R[:, None] + 1))
        else:
            thr = -np.log(lvl)
        alarm = logE >= thr
        isfalse = alarm & (t < onset); istrue = alarm & (t >= onset)
        nA = alarm.sum(1); nF = isfalse.sum(1)
        R += nA; F += nF
        if alarm.any():
            jj = np.minimum(j, JMAX)
            np.add.at(dsum, jj[istrue], (t - onset)[istrue]); np.add.at(dn, jj[istrue], 1)
            short = (alarm & ((t - last) <= LMIN)).any(1)
            hit = short & (tau_short < 0)
            tau_short[hit] = t; tauFDP[hit] = F[hit] / np.maximum(R[hit], 1)
            # new hypotheses: indices assigned in sensor order among alarmed sensors
            order = np.cumsum(alarm, axis=1)
            newk = nstarted[:, None] + order
            kidx = np.where(alarm, newk, kidx)
            nstarted += nA
            newlvl = ALPHA * gamma(kidx, n, const) * (R[:, None] + 1)  # R(sigma_k) includes alarms at t
            lvl = np.where(alarm, newlvl, lvl)
            logw[alarm] = 0.0; j[alarm] += 1; last[alarm] = t
            onset[alarm] = t + 1 + LMIN + rng.geometric(1.0 / LEXTRA, size=alarm.sum())
        prev_r = np.where(alarm, noise, r)
        noise_prev = noise
        fdp = np.where(R > 0, F / np.maximum(R, 1), 0.0)
        supFDP = np.maximum(supFDP, fdp)
        if t in FIXED_T:
            fixed[t] = fdp.copy()
    miss = tau_short < 0; tauFDP[miss] = F[miss] / np.maximum(R[miss], 1)
    return dict(fixed=fixed, sup=supFDP, tau=tauFDP, dsum=dsum, dn=dn, R=R, F=F)


def ci(a):
    a = np.asarray(a, float); m = a.mean(); se = a.std(ddof=1) / np.sqrt(len(a))
    return [round(float(m), 4), round(float(m - 1.96 * se), 4), round(float(m + 1.96 * se), 4)]


def main():
    t0 = time.time()
    args = sys.argv[2:]
    rules = ("iii", "iv", "iv_nonpred", "iv_const", "iv_raw")
    if args and args[0] == "--rules":
        rules = tuple(args[1].split(",")); args = args[2:]
    Ns = [int(x) for x in args] or [16]
    res = {"config": dict(alpha=ALPHA, seeds=SEEDS, reps_per_seed=REPS, T=T, phi=PHI, rho=RHO, Lmin=LMIN,
                          slope=SLOPE, gamma="1/(2N) k<=N; N/(2k(k-1)) k>N")}
    for n in Ns:
        for rule in rules:
            outs = [run(rule, s, n) for s in SEEDS]
            cat = lambda key: np.concatenate([o[key] for o in outs])
            ds = sum(o["dsum"] for o in outs); dn = sum(o["dn"] for o in outs)
            rr = {"fixed_FDR": {str(t): ci(np.concatenate([o["fixed"][t] for o in outs])) for t in FIXED_T},
                  "E_sup_FDP": ci(cat("sup")), "tau_short_FDR": ci(cat("tau")),
                  "per_seed_E_sup_FDP": [round(float(o["sup"].mean()), 4) for o in outs],
                  "mean_alarms": round(float(cat("R").mean()), 2), "mean_false": round(float(cat("F").mean()), 3),
                  "delay_by_epoch": {str(k): (round(float(ds[k] / dn[k]), 1) if dn[k] else None) for k in range(1, JMAX + 1)}}
            res[f"N{n}_{rule}"] = rr
            print(n, rule, rr["E_sup_FDP"], rr["tau_short_FDR"], rr["fixed_FDR"], rr["mean_alarms"],
                  {k: rr["delay_by_epoch"][k] for k in ("1", "2", "5", "10")}, flush=True)
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(sys.argv[1], "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
