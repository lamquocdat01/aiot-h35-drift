# PREREG_H35 — S1-1c (đăng ký trước, 21-11-2026, phiên bản v1)
Đề tài H35. Venue tạm: A = Technometrics (Q1 có điều kiện), B = Measurement Science and Technology (Q2).
Viết **trước khi chạy** mọi thí nghiệm 1c. Hash SHA-256 ở `prereg\PREREG_H35.sha256`; commit + tag `prereg-H35-v1` ở `%USERPROFILE%\Repos\aiot-h35-drift` (cục bộ, không push).
Mọi lệch khỏi file này ghi vào `prereg\AMENDMENT_LOG.md` **trước** khi chạy lại.
Kết quả S1 trước đó (pilot, 1b, Q-0, Q-2) là **thăm dò**. 1c là kiểm định xác nhận, dùng seed đã nêu, không chỉnh tham số sau khi xem kết quả.

## 0. Hằng số chung (chốt trước)
| Mục | Giá trị |
|---|---|
| α (báo giả theo cửa sổ) | 0,05. Trễ đo thêm ở α ∈ {10⁻³, 10⁻⁶} |
| Seed | 20261121, 20261122, 20261123, 20261124, 20261125 (5 seed) |
| Lần lặp | Null tổng hợp: 5 × 400 = 2 000 chuỗi. Trễ: 5 × 200 = 1 000. Pool thật: 5 × 200 = 1 000 / kênh. ARL₀: 5 × 100 = 500 chuỗi dài |
| Burn-in (huấn luyện, trong kiểm soát) | n₀ = 336 bước |
| Cửa sổ kiểm toán | W = 720 bước sau burn-in |
| Bootstrap null thật | khối vòng, độ dài khối ∈ {24, 72, 168} bước; báo cả ba, kết luận theo **168** (bảo thủ nhất cho tự tương quan) |
| CI | Wilson 95 % cho tỷ lệ; mean ± 1,96·SE cho trễ/ARL; bootstrap khối (khối 168) cho chênh lệch trên chuỗi thật |
| Đơn vị thời gian | UCI AQ: giờ. Gas Drift: chỉ số mẫu trong batch (giả định thứ tự file = thứ tự đo — **giả định, ghi vào báo cáo**). Intel Lab: trung bình 10 phút |

## 1. Monitor (tham số chốt trước)
| Ký hiệu | Monitor | Tham số | Ngưỡng ở α |
|---|---|---|---|
| M1 | CUSUM hai phía trên residual chuẩn hoá theo σ burn-in | k = 0,5 | h hiệu chuẩn bằng 20 000 chuỗi iid N(0,1) để P(báo trong W) = α |
| M2 | Page-Hinkley | δ = 0,25 | như M1 |
| M3 | ADWIN | δ = α, cửa sổ con ≥ 10, kiểm mỗi 4 bước | — |
| M4 | p < 0,05 lặp (z của trung bình tích luỹ) | từ bước 10 | \|z\| ≥ 1,96 |
| M5 | e-process cược (hỗn hợp 12 λ Gauss) | σ = sd burn-in | 1/α |
| M6 | Conformal test martingale (simple mixture, hai phía) | — | 2/α mỗi phía |
| F1 | Khối tự chuẩn hoá (Định lý 1) | b = 24, g = 24 | 1/α |
| **F2** | **BPI-WLC (Định lý 3)** | b = 24, g = 6, W_m = 60, κ = 0,1, p_max = 30 (AIC), lưới θ ±geomspace(0,05; 4; 24) | A = (⌈W/(b+g)⌉ + W_m)/α |
| B1c | CUSUM với phương sai dài hạn **đóng băng** từ huấn luyện (dòng CSW/Fremdt) | khối b = 24, g = 6; σ̂²_LR = (1+κ)·max(AR-sieve, Bartlett) từ burn-in, không cập nhật; Page CUSUM trên z_l, k = 0,5 | h hiệu chuẩn bằng iid N(0,1) khối cho P(báo trong W) = α |
| B2c | Biểu đồ kiểm soát trên residual AR (Alwan–Roberts) | AR(p̂) (AIC, p ≤ 30) fit trên burn-in, đóng băng; CUSUM hai phía trên residual chuẩn hoá theo σ̂_e, k = 0,5 | h hiệu chuẩn iid như M1 |
| R3a | e-BH tức thời + reset (16 kênh) | e-process F1 từng kênh | e-BH mức α |
| R3b | Phân bổ theo phát hiện (Định lý 2, kiểu e-LOND) | γ_k của THEORY §3 | — |

## 2. Dữ liệu và kênh (chốt trước)
| Bộ | Kênh | Residual | Null (trong dung sai) | Ghi chú |
|---|---|---|---|---|
| UCI Air Quality (đã có) | 4 cặp có tham chiếu: PT08.S1→CO(GT), PT08.S2→C6H6(GT), PT08.S3→NOx(GT), PT08.S4→NO2(GT). PT08.S5 (O3) **loại** vì không có tham chiếu | hồi quy tuyến tính trên 10–31/03/2004; residual ngoài mẫu bỏ-một-ngày | residual ngoài mẫu tháng 3 | −200 → NaN; NaN bỏ khỏi pool |
| UCI Gas Sensor Array Drift (tải 1c) | 16 cảm biến | mỗi cảm biến: hồi quy nồng độ từ 8 đặc trưng, **riêng từng khí**, fit trên batch 1; residual ngoài mẫu 5-fold trong batch 1 → pool null; batch 2–10 → drift thật | 5-fold ngoài mẫu trong batch 1 | không có thời gian thật; chỉ số mẫu. Trích license nguyên văn |
| Intel Berkeley Lab (tải 1c) | 20 mote có nhiều giá trị nhiệt độ hợp lệ nhất trong 7 ngày đầu | r_i = x_i − median láng giềng (5 mote gần nhất theo `mote_locs.txt`), trung bình 10 phút | 7 ngày đầu | không license → chỉ báo số tổng hợp, không phát hành lại; nhiệt độ ngoài [−10; 60] °C → NaN |
| Tổng hợp trong phạm vi Định lý 3 | AR(1) φ ∈ {0; 0,5; 0,8; 0,95}; ARMA(1,1) (0,8; +0,4), (0,8; −0,5); ARMA(2,1) (1,2, −0,4; +0,3) | — | — | Gauss |
| Tổng hợp **ngoài** phạm vi | AR(1) 0,8 với nhiễu t(5) chuẩn hoá; AR(1) 0,8 với nhiễu exp − 1; nhớ dài (AR 0,999 + nhiễu trắng) | — | — | |

## 3. Thí nghiệm
- **E1 (R1, kiểm toán α):** M1–M6, F1, F2, B1c, B2c trên mọi kênh thật (4 + 16 + 20) và mọi null tổng hợp. Chỉ số: báo giả trong W.
- **E2 (sửa vs baseline):**
  - báo giả trong W;
  - **ARL₀** (chuỗi null dài 20 000 bước, kiểm duyệt tại 20 000, báo mean + tỷ lệ bị kiểm duyệt);
  - **ARL₁** (thay đổi bậc ngay sau burn-in);
  - **trễ có điều kiện** E[τ − ν | τ ≥ ν] với ν = 240;
  - tỷ lệ bỏ sót trong 3 000 bước.
  - S chọn để D*(0,05) = 60 bước theo σ²_LR thật (tổng hợp) hoặc σ²_LR ước từ pool bằng batch means khối 168 (thật).
  - Drift thật: số lần báo và thời điểm báo đầu sau hiệu chuẩn, trên chuỗi thật sau giai đoạn null (UCI 04-2004 → 04-2005; Gas batch 2–10; Intel ngày 8 → cuối).
- **E3 (độ nhạy, lưới chốt trước):** F2 với b ∈ {12, 24, 48} × g ∈ {3, 6, 12}; F1 với b ∈ {24, 48}, g = b. Null: AR 0,8, B1 thật (CO). Chỉ số: báo giả, tỷ số trễ ở α = 10⁻⁶.
- **E4 (giả thiết (E)):** trên mỗi chuỗi null, sai số tương đối lớn nhất theo thời gian
  ε_max = max_l \|σ̂²_l − σ²_ref\| / σ²_ref  (σ̂² trước thổi κ).
  - σ²_ref: tổng hợp dùng σ²_LR thật **và** v_b; thật dùng batch means khối 168 trên pool.
  - Báo phân vị 50/95 của ε_max, và P(∃l: σ̃²_l < v̂_b,ref).
- **E5 (R3):** bán tổng hợp 16 kênh Gas Drift.
  - Mỗi kênh: chuỗi bootstrap từ pool null; drift tiêm dạng bậc S tại thời điểm ~ Geom(1/300) sau mỗi lần hiệu chuẩn; hiệu chuẩn = reset chuỗi về null.
  - So R3a với R3b: FDR tức thời tại τ_short, FDR vòng đời, E[sup FDP], số lần hiệu chuẩn.
  - T = 4 000, 5 × 200 lần lặp.
- **E6 (ESP32-S3 + INA219):** khối 10⁴ cập nhật của F2, đo RAM, µJ, latency. Bo chưa cắm → **việc tay**, không chặn.

## 4. Giả thuyết và tiêu chí PASS/FAIL (chốt trước)
| # | Giả thuyết | PASS khi | FAIL khi |
|---|---|---|---|
| H1 (R1) | ≥ 1 monitor "có bảo đảm" (M3, M5, M6) có α thực > α trên dữ liệu thật | Cận dưới CI > 0,05 (bootstrap 168) trên ≥ 50 % số kênh thật của **mỗi** bộ dữ liệu | Không đạt ở bộ nào |
| H2 (Định lý 3 (a)) | F2 giữ báo giả ≤ α trong phạm vi định lý và trên pool thật | Cận dưới CI ≤ 0,05 ở **mọi** null trong phạm vi và mọi kênh thật (bootstrap 168) | Có null/kênh với cận dưới CI > 0,05 |
| H3 (Định lý 3 (b)) | Trễ F2 theo hao phí dự đoán, nhanh hơn F1 | Tỷ số trễ/D* ở α = 10⁻⁶ ≤ 1,2 × hao phí dự đoán ở mọi null tổng hợp trong phạm vi; trễ F2 < trễ F1 ở ≥ 80 % cấu hình (null × α) | Vi phạm một trong hai |
| H4 (baseline) | B2c (residual AR) vỡ báo giả trên dữ liệu thật; B1c giữ báo giả nhưng chậm hơn F2 ở α nhỏ | Báo cáo mô tả; **không** có PASS/FAIL (thăm dò) | — |
| H5 (ngoài phạm vi) | F2 vỡ ở nhớ dài; t(5) và lệch có thể vẫn giữ | Báo nguyên văn; thăm dò | — |
| H6 (giả thiết (E)) | Ước lượng nhất quán đều theo thời gian ở mẫu hữu hạn | Phân vị 95 của ε_max ≤ κ = 0,10 trên AR φ ≤ 0,8 và ARMA | Phân vị 95 > 0,10 |
| H7 (R3) | e-BH + reset vỡ FDR tức thời nhưng không vỡ FDR vòng đời; R3b giữ cả E[sup FDP] ≤ α | FDR tức thời (R3a) có cận dưới CI > α **và** E[sup FDP] (R3b) có cận trên CI ≤ α | Ngược lại |
Kết luận venue: H2 và H3 PASS là điều kiện cần để giữ A = Technometrics. Một trong hai FAIL → đề xuất B = MST (S2 quyết).

## 5. Báo cáo
- Kết quả âm báo nguyên văn.
- Số tổng hợp vào `results\numbers.json` qua `code\build_numbers.py`; hình vào `results\fig\*.pdf`.
- Báo cáo cuối: `S1\S1_REPORT_H35.md` (≤ 3 trang).
