# Giai đoạn 2 & 3 — BLAST/Mash baseline + Robustness-to-noise

Bộ script này đã được **chạy thử và xác nhận hoạt động đúng** trên dữ liệu giả lập
(synthetic FASTA). Anh chỉ cần thay bằng dữ liệu thật của T0–T3 và closed-world corpus.

## 0. Cài đặt môi trường

```bash
apt-get update && apt-get install -y ncbi-blast+ mash
pip install biopython numpy pandas scikit-learn matplotlib
```

## 1. Chuẩn bị dữ liệu đầu vào

### Trường hợp A — dữ liệu của anh có cấu trúc thư mục kiểu NCBI Virus download

Nếu dữ liệu của anh có dạng:

```
sars-cov2-time-forward-dataset/
├── T0_train/
│   ├── sequences_..._Alpha_B.1.1.7_2000_....fasta
│   ├── sequences_..._Beta_B.1.351_2000_....fasta
│   └── ...
├── T1_test_pre_omicron/...
├── T2_test_transition/...
└── T3_test_post_omicron/...
```

(nhiều file FASTA/period, một file/lineage, chưa có metadata.csv) — đặt nguyên thư mục
`sars-cov2-time-forward-dataset/` vào `data/`, rồi chạy script gộp tự động:

```bash
# bước 1: chạy dry-run trước để kiểm tra parsing đúng, KHÔNG ghi file
python3 scripts/00a_prepare_dataset.py \
  --root data/sars-cov2-time-forward-dataset --dry-run

# đọc kỹ phần "declared vs actual" — nếu có dòng MISMATCH, nghĩa là số sequence
# thực tế trong file khác số ghi trong tên file (có thể do download thiếu/OK,
# không phải lỗi script — chỉ cần soát lại nếu chênh lệch quá lớn hoặc bất thường)

# bước 2: chạy thật, sinh ra data/T0.fasta...T3.fasta + data/metadata.csv
python3 scripts/00a_prepare_dataset.py \
  --root data/sars-cov2-time-forward-dataset
```

Script tự động:
- Lấy `period` từ tên thư mục con (T0/T1/T2/T3).
- Lấy `lineage` từ tên file (token ngay sau `__`, ví dụ `Alpha`, `Omicron`).
- Lấy `accession` từ FASTA header (token đầu tiên sau `>`, bỏ version `.1`/`.2`).
- Gộp tất cả file cùng period thành 1 file `data/T<n>.fasta`, sinh `data/metadata.csv`.
- Cảnh báo nếu có accession trùng lặp hoặc file không parse được tên đúng quy ước.

Sau bước này, đi thẳng tới **Mục 2 — Chạy Giai đoạn 2** bên dưới.

### Trường hợp B — anh tự chuẩn bị sẵn 1 file/period

Nếu không dùng script trên, đặt thủ công vào `data/`:

- `data/T0.fasta`, `data/T1.fasta`, `data/T2.fasta`, `data/T3.fasta` — mỗi file là các
  genome tương ứng period đó (đúng như corpus đã dùng để train/eval encoder, Bảng 1
  trong bài báo).
- `data/metadata.csv` — bắt buộc có 4 cột:

  ```csv
  accession,lineage,period,collection_date
  MZ123456,Alpha,T0,2021-03-01
  MZ123457,Omicron,T2,2021-11-15
  ...
  ```

  `period` phải là một trong `T0`, `T1`, `T2`, `T3`. `accession` phải khớp với ID
  trong header FASTA (script tự cắt bỏ version suffix `.1` và mọi thứ sau `|`, nên
  `>MZ123456.1 Severe acute...` vẫn khớp được với `MZ123456`).

Nếu dữ liệu closed-world corpus muốn benchmark thêm (Bước 2.4, tùy chọn), chuẩn bị
tương tự với period tự đặt tên khác (ví dụ `CW_train`, `CW_test`) và điều chỉnh script
`03_compute_metrics.py` cho phù hợp (phần AUROC hiện đang hard-code IND=T1/OOD=Omicron
T2+T3 theo đúng Section 3.6 của bài báo).

## 2. Chạy Giai đoạn 2 — BLAST/Mash baseline

```bash
cd pipeline

# Bước 2.2: chọn reference set đại diện từ T0 (mặc định 30 seq/lineage, đổi bằng --n-per-class)
python3 scripts/01_build_reference_set.py \
  --fasta data/T0.fasta --meta data/metadata.csv --n-per-class 30 --outdir results

# Bước 2.3 (phần 1): tách T0 thành calibration set (T0 trừ reference) để fit threshold
python3 scripts/00_split_t0_calib.py \
  --fasta data/T0.fasta --ref-meta results/reference_metadata.csv --out results/T0_calib.fasta

# Bước 2.3 (phần 2): chạy BLAST + Mash cho calibration set và T1/T2/T3
bash scripts/02_run_baselines.sh T0_calib results/T0_calib.fasta results/reference.fasta
bash scripts/02_run_baselines.sh T1 data/T1.fasta results/reference.fasta
bash scripts/02_run_baselines.sh T2 data/T2.fasta results/reference.fasta
bash scripts/02_run_baselines.sh T3 data/T3.fasta results/reference.fasta

# Bước 2.4/2.5: fit threshold trên T0_calib (95th percentile, giống theta của paper),
# tính rejection rate theo lineage/period (giống Bảng 6) và AUROC/FPR95 (giống Bảng 4)
python3 scripts/03_compute_metrics.py --meta data/metadata.csv --method blast
python3 scripts/03_compute_metrics.py --meta data/metadata.csv --method mash
```

**Kết quả sau bước này:**
- `results/table6_blast.csv`, `results/table6_mash.csv` — rejection rate theo
  lineage/period, format thẳng vào Bảng 6 mở rộng của bài báo.
- `results/table7_comparison.csv` — AUROC, FPR95, runtime, threshold của cả 2
  phương pháp — đây là **Bảng 7 (BLAST/Mash baseline comparison)** đề xuất trong kế
  hoạch, đặt cạnh Bảng 4 (4 baseline OOD hiện có: MSP/ODIN/energy/k-NN) để so sánh
  đầy đủ 6 phương pháp.

So sánh 2 dòng cuối `table7_comparison.csv` với AUROC/FPR95 của encoder (0.990/4.36%
theo Section 4.11) là con số chính cần đưa vào Discussion.

## 3. Chạy Giai đoạn 3 — Robustness-to-noise

Chạy trên period chứa Omicron để đo trực tiếp khả năng phát hiện biến thể mới dưới
nhiễu (khuyến nghị T2, có thể lặp lại cho T3):

```bash
# Bước 3.1 + 3.2: sinh nhiễu (substitution 0.5/1/2/5%, truncation 5/10/20%) và
# chạy BLAST+Mash ở từng mức
bash scripts/05_run_robustness.sh T2 data/T2.fasta

# Bước 3.3: tổng hợp accuracy theo từng mức nhiễu + vẽ Figure 7 (đường BLAST, Mash)
python3 scripts/06_aggregate_robustness.py --meta data/metadata.csv --period T2 --method blast
python3 scripts/06_aggregate_robustness.py --meta data/metadata.csv --period T2 --method mash
```

**Quan trọng:** để Figure 7 có đủ 2-3 đường so sánh (encoder vs BLAST vs Mash) như kế
hoạch đề ra, cần chạy **chính pipeline encoder + Mahalanobis hiện có** của anh trên các
file nhiễu đã sinh sẵn ở `results/noise_T2/*.fasta` (cùng seed, cùng accession, nên so
sánh hoàn toàn công bằng), rồi tự thêm dòng vào `results/robustness_summary.csv` theo
đúng format:

```csv
method,period,level,accuracy,n
encoder,T2,level_0,<số của anh>,21
encoder,T2,sub_0.005,<số của anh>,21
...
```

Sau đó chạy lại để vẽ:

```bash
python3 scripts/06_aggregate_robustness.py --plot-only
```

**Kết quả sau bước này:**
- `results/robustness_summary.csv` — accuracy theo mức nhiễu cho từng method.
- `figures/figure7_robustness.png` — hình so sánh độ suy giảm hiệu năng theo nhiễu,
  đúng như Figure 7 đề xuất trong kế hoạch.

## 4. Đã kiểm thử

Toàn bộ pipeline (bước 0 → 3) đã được chạy thử thành công trên dữ liệu synthetic
(160 seq T0 / 33 T1 / 21 T2 / 23 T3, genome dài 2.000 nt) để xác nhận không có lỗi cú
pháp/logic. Với dữ liệu thật (genome ~29.903 nt, hàng nghìn sequence mỗi period),
BLAST sẽ chạy chậm hơn — cân nhắc thêm `-num_threads` (script đã tự dùng `nproc`) và
chạy nền (`nohup ... &`) nếu makeblastdb/blastn mất nhiều phút.

## 5. Các điểm cần anh quyết định trước khi chạy dữ liệu thật

1. **`--n-per-class` của reference set** (mặc định 30): cân bằng giữa "đại diện đủ tốt"
   và "vẫn đúng tinh thần reference nhỏ, hiệu quả" mà reviewer yêu cầu — có thể thử vài
   giá trị (10/30/50) và báo cáo độ nhạy trong bài nếu muốn chắc chắn hơn.
2. **BLAST hay Mash là baseline chính trong bài?** Nên báo cáo cả hai (đã có sẵn) —
   BLAST trả lời trực tiếp tên công cụ reviewer nêu, Mash củng cố luận điểm về tốc độ.
3. **Chạy robustness trên T2, T3, hay cả hai?** T2 (giai đoạn Omicron mới xuất hiện,
   n nhỏ) là quan trọng nhất về mặt câu chuyện early-warning; T3 (Omicron đã áp đảo)
   có thể chạy thêm nếu còn thời gian để củng cố.
