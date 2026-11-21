# PREREG_H35 — English translation (not the registered document)

> **The original pre-registration is `prereg/PREREG_H35.md`, written in Vietnamese.** Its SHA-256 is recorded in `prereg/PREREG_H35.sha256` and in the tag `prereg-H35-v1` (commit dated 21-11-2026 11:09 +0700, before any confirmatory run). This English file is a translation made on 21-11-2026 during release preparation; it is **not** covered by the hash. If the two versions differ, the Vietnamese original prevails. Paths are given as in the original (Windows separators). Venue names in the original reflect the submission plan at the time of registration.

# PREREG_H35 — S1-1c (pre-registration, 21-11-2026, version v1)
Topic H35. Tentative venues: A = Technometrics (Q1, conditional), B = Measurement Science and Technology (Q2).
Written **before running** any 1c experiment. SHA-256 hash in `prereg\PREREG_H35.sha256`; commit + tag `prereg-H35-v1` in `%USERPROFILE%\Repos\aiot-h35-drift` (local, not pushed).
Every deviation from this file is recorded in `prereg\AMENDMENT_LOG.md` **before** re-running.
Earlier S1 results (pilot, 1b, Q-0, Q-2) are **exploratory**. 1c is the confirmatory test, uses the stated seeds, and parameters are not adjusted after seeing results.

## 0. Common constants (fixed in advance)
| Item | Value |
|---|---|
| α (window-wise false alarm) | 0.05. Delay additionally measured at α ∈ {10⁻³, 10⁻⁶} |
| Seeds | 20261121, 20261122, 20261123, 20261124, 20261125 (5 seeds) |
| Replications | Synthetic nulls: 5 × 400 = 2 000 series. Delay: 5 × 200 = 1 000. Real pools: 5 × 200 = 1 000 per channel. ARL₀: 5 × 100 = 500 long series |
| Burn-in (training, in control) | n₀ = 336 steps |
| Audit window | W = 720 steps after burn-in |
| Real-data bootstrap null | circular block, block length ∈ {24, 72, 168} steps; all three reported, conclusions based on **168** (most conservative for autocorrelation) |
| CI | Wilson 95 % for proportions; mean ± 1.96·SE for delay/ARL; block bootstrap (block 168) for differences on real series |
| Time unit | UCI AQ: hours. Gas Drift: sample index within batch (assumes file order = measurement order — **an assumption, to be stated in the report**). Intel Lab: 10-minute means |

## 1. Monitors (parameters fixed in advance)
| Symbol | Monitor | Parameters | Threshold at α |
|---|---|---|---|
| M1 | Two-sided CUSUM on residual standardised by burn-in σ | k = 0.5 | h calibrated with 20 000 iid N(0,1) series so that P(alarm in W) = α |
| M2 | Page-Hinkley | δ = 0.25 | as M1 |
| M3 | ADWIN | δ = α, sub-window ≥ 10, check every 4 steps | — |
| M4 | Repeated p < 0.05 (z of cumulative mean) | from step 10 | \|z\| ≥ 1.96 |
| M5 | Betting e-process (mixture of 12 Gaussian λ) | σ = burn-in sd | 1/α |
| M6 | Conformal test martingale (simple mixture, two-sided) | — | 2/α per side |
| F1 | Self-normalised blocks (Theorem 1) | b = 24, g = 24 | 1/α |
| **F2** | **BPI-WLC (Theorem 3)** | b = 24, g = 6, W_m = 60, κ = 0.1, p_max = 30 (AIC), θ grid ±geomspace(0.05; 4; 24) | A = (⌈W/(b+g)⌉ + W_m)/α |
| B1c | CUSUM with long-run variance **frozen** from training (CSW/Fremdt line) | blocks b = 24, g = 6; σ̂²_LR = (1+κ)·max(AR-sieve, Bartlett) from burn-in, not updated; Page CUSUM on z_l, k = 0.5 | h calibrated on iid N(0,1) blocks so that P(alarm in W) = α |
| B2c | Control chart on AR residuals (Alwan–Roberts) | AR(p̂) (AIC, p ≤ 30) fitted on burn-in, frozen; two-sided CUSUM on residuals standardised by σ̂_e, k = 0.5 | h calibrated on iid as M1 |
| R3a | Instantaneous e-BH + reset (16 channels) | F1 e-process per channel | e-BH at level α |
| R3b | Discovery-based allocation (Theorem 2, e-LOND type) | γ_k of THEORY §3 | — |

## 2. Data and channels (fixed in advance)
| Set | Channels | Residual | Null (in tolerance) | Notes |
|---|---|---|---|---|
| UCI Air Quality (available) | 4 referenced pairs: PT08.S1→CO(GT), PT08.S2→C6H6(GT), PT08.S3→NOx(GT), PT08.S4→NO2(GT). PT08.S5 (O3) **excluded** because it has no reference | linear regression on 10–31/03/2004; leave-one-day-out out-of-sample residual | out-of-sample residuals of March | −200 → NaN; NaN dropped from pool |
| UCI Gas Sensor Array Drift (downloaded in 1c) | 16 sensors | per sensor: regression of concentration on 8 features, **separately per gas**, fitted on batch 1; 5-fold out-of-sample residuals within batch 1 → null pool; batches 2–10 → real drift | 5-fold out-of-sample within batch 1 | no real time; sample index only. Quote licence verbatim |
| Intel Berkeley Lab (downloaded in 1c) | 20 motes with most valid temperature values in the first 7 days | r_i = x_i − median of neighbours (5 nearest motes by `mote_locs.txt`), 10-minute means | first 7 days | no licence → report aggregate numbers only, no redistribution; temperature outside [−10; 60] °C → NaN |
| Synthetic within the scope of Theorem 3 | AR(1) φ ∈ {0; 0.5; 0.8; 0.95}; ARMA(1,1) (0.8; +0.4), (0.8; −0.5); ARMA(2,1) (1.2, −0.4; +0.3) | — | — | Gaussian |
| Synthetic **outside** scope | AR(1) 0.8 with standardised t(5) noise; AR(1) 0.8 with exp − 1 noise; long memory (AR 0.999 + white noise) | — | — | |

## 3. Experiments
- **E1 (R1, α audit):** M1–M6, F1, F2, B1c, B2c on all real channels (4 + 16 + 20) and all synthetic nulls. Metric: false alarm within W.
- **E2 (corrections vs baselines):**
  - false alarm within W;
  - **ARL₀** (null series of 20 000 steps, censored at 20 000, report mean + censored fraction);
  - **ARL₁** (step change right after burn-in);
  - **conditional delay** E[τ − ν | τ ≥ ν] with ν = 240;
  - miss rate within 3 000 steps.
  - S chosen so that D*(0.05) = 60 steps using the true σ²_LR (synthetic) or σ²_LR estimated from the pool by batch means with block 168 (real).
  - Real drift: number of alarms and time of first alarm after calibration, on the real series after the null period (UCI 04-2004 → 04-2005; Gas batches 2–10; Intel day 8 → end).
- **E3 (sensitivity, grid fixed in advance):** F2 with b ∈ {12, 24, 48} × g ∈ {3, 6, 12}; F1 with b ∈ {24, 48}, g = b. Nulls: AR 0.8, real B1 (CO). Metrics: false alarm, delay ratio at α = 10⁻⁶.
- **E4 (assumption (E)):** on each null series, the largest relative error over time
  ε_max = max_l \|σ̂²_l − σ²_ref\| / σ²_ref  (σ̂² before κ inflation).
  - σ²_ref: synthetic uses the true σ²_LR **and** v_b; real uses batch means with block 168 on the pool.
  - Report the 50th/95th percentiles of ε_max, and P(∃l: σ̃²_l < v̂_b,ref).
- **E5 (R3):** semi-synthetic, 16 Gas Drift channels.
  - Each channel: bootstrap series from the null pool; step drift of size S injected at a time ~ Geom(1/300) after each calibration; calibration = reset of the series to null.
  - Compare R3a with R3b: instantaneous FDR at τ_short, lifetime FDR, E[sup FDP], number of calibrations.
  - T = 4 000, 5 × 200 replications.
- **E6 (ESP32-S3 + INA219):** block of 10⁴ F2 updates, measure RAM, µJ, latency. Board not connected → **manual task**, not blocking.

## 4. Hypotheses and PASS/FAIL criteria (fixed in advance)
| # | Hypothesis | PASS when | FAIL when |
|---|---|---|---|
| H1 (R1) | ≥ 1 "guaranteed" monitor (M3, M5, M6) has actual α > α on real data | Lower CI bound > 0.05 (bootstrap 168) on ≥ 50 % of real channels of **each** dataset | Not met on any dataset |
| H2 (Theorem 3 (a)) | F2 keeps false alarm ≤ α within the scope of the theorem and on real pools | Lower CI bound ≤ 0.05 on **every** in-scope null and every real channel (bootstrap 168) | Some null/channel with lower CI bound > 0.05 |
| H3 (Theorem 3 (b)) | F2 delay follows the predicted overhead and is faster than F1 | Delay/D* ratio at α = 10⁻⁶ ≤ 1.2 × predicted overhead on every in-scope synthetic null; F2 delay < F1 delay in ≥ 80 % of configurations (null × α) | Either is violated |
| H4 (baselines) | B2c (AR residual) breaks false alarm on real data; B1c keeps false alarm but is slower than F2 at small α | Descriptive report; **no** PASS/FAIL (exploratory) | — |
| H5 (out of scope) | F2 breaks under long memory; t(5) and skew may still hold | Report verbatim; exploratory | — |
| H6 (assumption (E)) | Estimator uniformly consistent over time in finite samples | 95th percentile of ε_max ≤ κ = 0.10 on AR φ ≤ 0.8 and ARMA | 95th percentile > 0.10 |
| H7 (R3) | e-BH + reset breaks instantaneous FDR but not lifetime FDR; R3b keeps E[sup FDP] ≤ α | Instantaneous FDR (R3a) has lower CI bound > α **and** E[sup FDP] (R3b) has upper CI bound ≤ α | Otherwise |
Venue conclusion: H2 and H3 PASS are necessary to keep A = Technometrics. If either FAILs → propose B = MST (decided at S2).

## 5. Reporting
- Negative results reported verbatim.
- Summary numbers go to `results\numbers.json` via `code\build_numbers.py`; figures to `results\fig\*.pdf`.
- Final report: `S1\S1_REPORT_H35.md` (≤ 3 pages). *(Translator's note: internal working report, not part of this repository.)*
