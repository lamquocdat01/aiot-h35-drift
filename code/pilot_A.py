"""H35 S1-1a Pilot A: lifetime FDR with alarm-triggered recalibration (numpy simulation).

Model (per sensor i, epoch j): after each (re)calibration the sensor is in tolerance (null, drift=0)
for L = Lmin + Geom(1/Lextra) steps, then a ramp drift d_t = s*(t-onset) starts (alternative).
Residual r_t = d_t + AR(1) noise (phi known); whitened x_t = r_t - phi*r_{t-1}.
Cross-sensor dependence: innovations share a common factor (corr rho).
E-process: mixture over lambda grid of prod exp(lambda*x - lambda^2/2)  (test martingale under null).
Alarm = recalibration: wealth reset to 1, drift reset to 0, new epoch j+1.
Alarm is FALSE if t < onset (sensor still in tolerance).

Rules: (i) const 1/alpha + reset; (ii) e-BH over N sensors each step + reset;
       (iii) epoch threshold N/(alpha*g_j), g_j = 1/(j(j+1));  (iii') per-sensor 1/(alpha*g_j) (no /N).
Stopping times: fixed t; tau_first (first alarm); tau_short (first alarm within Lmin of that sensor's
last recal); tau_short5 (tau_short if R<=5 at that time, else T).
"""
import json, sys, time
import numpy as np

ALPHA = 0.05
SEEDS = [20261121 + k for k in range(5)]
REPS = 1000
N = 16
T = 4000
PHI = 0.5
RHO = 0.3
LMIN = 50
LEXTRA = 150
SLOPE = 0.02
LAMS = np.array([0.1, 0.2, 0.35, 0.5, 0.75, 1.0])
FIXED_T = [500, 1000, 2000, 4000]
JMAX = 12


def run(rule, seed, n=N, reps=REPS):
    rng = np.random.default_rng(seed)
    K = len(LAMS)
    logw = np.zeros((reps, n, K))
    j = np.ones((reps, n), dtype=int)
    last = np.zeros((reps, n), dtype=int)
    onset = LMIN + rng.geometric(1.0 / LEXTRA, size=(reps, n))
    prev_r = np.zeros((reps, n))
    noise_prev = rng.standard_normal((reps, n)) / np.sqrt(1 - PHI ** 2)
    prev_d = np.zeros((reps, n))
    R = np.zeros(reps, int); F = np.zeros(reps, int)
    fixed = {t: None for t in FIXED_T}
    taus = {k: np.full(reps, -1) for k in ("first", "short", "short5")}
    tauFDP = {k: np.zeros(reps) for k in taus}
    tauCur = {k: np.zeros(reps) for k in taus}  # FDP of alarms issued at tau (current decision)
    delay_sum = np.zeros(JMAX + 1); delay_n = np.zeros(JMAX + 1)
    false_by_j = np.zeros(JMAX + 1); epochs_by_j = np.zeros(JMAX + 1)
    lam2 = LAMS ** 2 / 2
    for t in range(1, T + 1):
        c = rng.standard_normal((reps, 1))
        u = rng.standard_normal((reps, n))
        eps = np.sqrt(RHO) * c + np.sqrt(1 - RHO) * u
        noise = PHI * noise_prev + eps
        d = SLOPE * np.maximum(0, t - onset)
        r = d + noise
        x = r - PHI * prev_r
        # after a recal at t-1 the previous drift was reset; r_{t-1} as observed stays (noise continuity)
        logw += x[..., None] * LAMS - lam2
        logE = np.logaddexp.reduce(logw, axis=2) - np.log(K)
        if rule == "i":
            alarm = logE >= np.log(1 / ALPHA)
        elif rule == "iii":
            g = 1.0 / (j * (j + 1))
            alarm = logE >= np.log(n / (ALPHA * g))
        elif rule == "iii_p":
            g = 1.0 / (j * (j + 1))
            alarm = logE >= np.log(1 / (ALPHA * g))
        elif rule == "ii":
            s = -np.sort(-logE, axis=1)
            thr = np.log(n / (ALPHA * np.arange(1, n + 1)))
            ok = s >= thr
            any_ok = ok.any(axis=1)
            kmax = n - np.argmax(ok[:, ::-1], axis=1)  # largest k with ok
            cut = s[np.arange(reps), kmax - 1]
            alarm = (logE >= cut[:, None]) & any_ok[:, None]
        else:
            raise ValueError(rule)
        isfalse = alarm & (t < onset)
        istrue = alarm & (t >= onset)
        nA = alarm.sum(1); nF = isfalse.sum(1)
        R += nA; F += nF
        # stats by epoch
        if alarm.any():
            jj = np.minimum(j, JMAX)
            np.add.at(epochs_by_j, jj[alarm], 1)
            np.add.at(false_by_j, jj[isfalse], 1)
            np.add.at(delay_sum, jj[istrue], (t - onset)[istrue])
            np.add.at(delay_n, jj[istrue], 1)
            short = (alarm & ((t - last) <= LMIN)).any(1)
            first_any = nA > 0
            for k, cond in (("first", first_any), ("short", short), ("short5", short & (R <= 5))):
                hit = cond & (taus[k] < 0)
                taus[k][hit] = t
                tauFDP[k][hit] = F[hit] / np.maximum(R[hit], 1)
                tauCur[k][hit] = nF[hit] / np.maximum(nA[hit], 1)
            # reset
            logw[alarm] = 0.0
            j[alarm] += 1
            last[alarm] = t
            onset[alarm] = t + 1 + LMIN + rng.geometric(1.0 / LEXTRA, size=alarm.sum())
        prev_r = np.where(alarm, noise, r)
        noise_prev = noise
        if t in fixed:
            fixed[t] = (F / np.maximum(R, 1)).copy()
    for k in taus:  # not stopped -> stop at T
        miss = taus[k] < 0
        tauFDP[k][miss] = F[miss] / np.maximum(R[miss], 1)
        tauCur[k][miss] = 0.0
    out = {"fixed": {t: fixed[t] for t in FIXED_T},
           "tau": {k: tauFDP[k] for k in taus}, "tau_cur": {k: tauCur[k] for k in taus},
           "tau_hit": {k: (taus[k] > 0) for k in taus},
           "delay_sum": delay_sum, "delay_n": delay_n, "false_by_j": false_by_j, "epochs_by_j": epochs_by_j,
           "R": R, "F": F}
    return out


def ci(a):
    a = np.asarray(a, float)
    m = a.mean(); se = a.std(ddof=1) / np.sqrt(len(a))
    return [round(float(m), 4), round(float(m - 1.96 * se), 4), round(float(m + 1.96 * se), 4)]


def counterexample_A(seed, reps=200000):
    """N=1. Step1 of each epoch: e=1/alpha w.p. alpha/2 (false alarm), else c; step2 drift -> true alarm.
    tau = first step-1 alarm. Current-decision FDP at tau = 1; lifetime FDP at tau = 1/K."""
    rng = np.random.default_rng(seed)
    p = ALPHA / 2
    c = 0.5 / (1 - p)
    Kep = rng.geometric(p, size=reps)  # epoch index of the first step-1 alarm
    life = 1.0 / Kep  # F=1, R = (Kep-1) true + 1 false = Kep
    return {"E_e_step1": p * (1 / ALPHA) + (1 - p) * c, "tau_finite_frac": 1.0,
            "FDR_current_at_tau": 1.0, "lifetime_FDP_at_tau_mc": ci(life),
            "lifetime_FDP_at_tau_exact": float(p * (-np.log(p)) / (1 - p))}


def counterexample_B(seed, reps=200000):
    rng = np.random.default_rng(seed)
    z1 = np.where(rng.random(reps) < 0.5, 2.0, 0.0)
    ztau = np.where(z1 == 2.0, 2.0, 1.0)
    return {"E_Z_tau_mc": ci(ztau), "E_Z_tau_exact": 1.5}


def main():
    t0 = time.time()
    res = {"config": dict(alpha=ALPHA, seeds=SEEDS, reps_per_seed=REPS, N=N, T=T, phi=PHI, rho=RHO, Lmin=LMIN,
                          Lextra_mean=LEXTRA, slope=SLOPE, lambdas=LAMS.tolist(), g_j="1/(j(j+1))"),
           "bound_alpha_1plog": ALPHA * (1 + np.log(1 / ALPHA)), "bound_alpha_over_1m": ALPHA / (1 - ALPHA)}
    rules = {}
    for rule in ("i", "ii", "iii", "iii_p"):
        outs = [run(rule, s) for s in SEEDS]
        cat = lambda f: np.concatenate([f(o) for o in outs])
        rr = {"fixed_t_lifetimeFDR": {str(t): ci(cat(lambda o: o["fixed"][t])) for t in FIXED_T},
              "tau_lifetimeFDR": {k: ci(cat(lambda o: o["tau"][k])) for k in ("first", "short", "short5")},
              "tau_currentFDR": {k: ci(cat(lambda o: o["tau_cur"][k])) for k in ("first", "short", "short5")},
              "tau_hit_frac": {k: round(float(cat(lambda o: o["tau_hit"][k]).mean()), 4) for k in ("first", "short", "short5")},
              "per_seed_tau_short_lifetime": [round(float(o["tau"]["short"].mean()), 4) for o in outs],
              "per_seed_tau_first_lifetime": [round(float(o["tau"]["first"].mean()), 4) for o in outs],
              "mean_alarms_per_rep": round(float(cat(lambda o: o["R"]).mean()), 2)}
        ds = sum(o["delay_sum"] for o in outs); dn = sum(o["delay_n"] for o in outs)
        fb = sum(o["false_by_j"] for o in outs); eb = sum(o["epochs_by_j"] for o in outs)
        rr["delay_by_epoch"] = {str(k): (round(float(ds[k] / dn[k]), 1) if dn[k] > 0 else None) for k in range(1, JMAX + 1)}
        rr["n_true_by_epoch"] = {str(k): int(dn[k]) for k in range(1, JMAX + 1)}
        rr["false_frac_by_epoch"] = {str(k): (round(float(fb[k] / eb[k]), 4) if eb[k] > 0 else None) for k in range(1, JMAX + 1)}
        rules[rule] = rr
        print(rule, rr["fixed_t_lifetimeFDR"], rr["tau_lifetimeFDR"], rr["tau_currentFDR"], rr["delay_by_epoch"], flush=True)
    res["rules"] = rules
    # single-stream check (N=1) for rule i to separate multiplicity from reset
    outs = [run("i", s, n=1) for s in SEEDS]
    cat = lambda f: np.concatenate([f(o) for o in outs])
    res["rule_i_N1"] = {"fixed_t_lifetimeFDR": {str(t): ci(cat(lambda o: o["fixed"][t])) for t in FIXED_T},
                        "tau_lifetimeFDR": {k: ci(cat(lambda o: o["tau"][k])) for k in ("first", "short", "short5")},
                        "tau_currentFDR": {k: ci(cat(lambda o: o["tau_cur"][k])) for k in ("first", "short", "short5")}}
    print("i N=1", res["rule_i_N1"], flush=True)
    res["counterexample_A"] = [counterexample_A(s) for s in SEEDS]
    res["counterexample_B"] = [counterexample_B(s) for s in SEEDS]
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(sys.argv[1], "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, ensure_ascii=False)
    print("done", res["runtime_s"])


if __name__ == "__main__":
    main()
