# Hướng dẫn chi tiết sử dụng mã nguồn — Giai đoạn 2 & 3

Hướng dẫn này đi **từng bước một**, kèm lệnh chạy cụ thể, output mong đợi ở mỗi bước,
và cách xử lý các lỗi thường gặp. Mã nguồn nằm trong thư mục `pipeline/` với cấu trúc:

```
pipeline/
├── README.md                       (tóm tắt nhanh)
├── HUONG_DAN_SU_DUNG.md            (file này — hướng dẫn chi tiết)
├── data/
│   └── sars-cov2-time-forward-dataset/   ← đặt dữ liệu thật vào đây
│       ├── T0_train/
│       │   ├── sequences_..._Alpha_B.1.1.7_2000_....fasta
│       │   ├── sequences_..._Beta_B.1.351_2000_....fasta
│       │   ├── sequences_..._Gamma_P.1_2000_....fasta
│       │   └── sequences_..._Delta_B.1.617.2_2000_....fasta
│       ├── T1_test_pre_omicron/
│       ├── T2_test_transition/
│       └── T3_test_post_omicron/
├── results/                        (sẽ được sinh ra khi chạy)
├── figures/                        (sẽ được sinh ra khi chạy)
└── scripts/
    ├── 00a_prepare_dataset.py      ← BƯỚC MỚI: tiền xử lý dữ liệu thật
    ├── 00_split_t0_calib.py
    ├── 01_build_reference_set.py
    ├── 02_run_baselines.sh
    ├── 03_compute_metrics.py
    ├── 04_generate_noise.py
    ├── 05_run_robustness.sh
    └── 06_aggregate_robustness.py
```

> **Lý do cần bước tiền xử lý mới:** dữ liệu thật có cấu trúc *nhiều file FASTA
> mỗi period* (một file/lineage) và *chưa có sẵn `metadata.csv`* — khác với giả định
> cũ (1 file gộp/period + metadata.csv có sẵn). Script `00a_prepare_dataset.py` giải
> quyết toàn bộ phần này tự động, sau đó các bước Giai đoạn 2 & 3 (01→06) chạy y hệt
> như trước — không cần sửa thêm gì.

---

## Bước 0 — Cài đặt môi trường

```bash
apt-get update
apt-get install -y ncbi-blast+ mash
pip install biopython numpy pandas scikit-learn matplotlib
```

Kiểm tra cài đặt thành công:
```bash
blastn -version
mash --version
python3 -c "import Bio, numpy, pandas, sklearn, matplotlib; print('OK')"
```
Nếu cả 3 lệnh chạy không lỗi, môi trường đã sẵn sàng.

---

## Bước 1 — Tiền xử lý: gộp dữ liệu thật thành đúng định dạng

> **Đây là bước bắt buộc mới** — không có trong phiên bản hướng dẫn cũ vì trước đó
> giả định `T0.fasta`…`T3.fasta` và `metadata.csv` đã có sẵn. Với dữ liệu thật, cần
> chạy bước này trước tất cả các bước còn lại.

### 1.1. Kiểm tra cấu trúc dữ liệu đầu vào

Đảm bảo thư mục `data/sars-cov2-time-forward-dataset/` đặt đúng vị trí và có 4 thư mục
con `T0_train/`, `T1_test_pre_omicron/`, `T2_test_transition/`, `T3_test_post_omicron/`.
Mỗi thư mục chứa một số file `.fasta`, mỗi file một lineage, tên file theo định dạng:

```
sequences_<date>_<id>__<Lineage>_<PangoLineage>_<count>_<start>_<end>.fasta
```

Ví dụ thực tế trong `T0_train/`:
```
sequences_20260705_8911960__Alpha_B.1.1.7_2000_31Dec2019_30Sep2021.fasta
sequences_20260705_5268294__Beta_B.1.351_2000_31Dec2019_30Sep2021.fasta
sequences_20260705_4057718__Gamma_P.1_2000_31Dec2019_30Sep2021.fasta
sequences_20260705_7446879__Delta_B.1.617.2_2000_31Dec2019_30Sep2021.fasta
```

Script sẽ tự động đọc:
- `period` từ tên thư mục (`T0_train` → `T0`, `T1_test_pre_omicron` → `T1`, v.v.)
- `lineage` từ tên file (token ngay sau `__`, ví dụ `Alpha`, `Beta`, `Gamma`, `Delta`)
- `accession` từ FASTA header (token đầu tiên sau `>`, bỏ phần version `.1`/`.2`)

### 1.2. Chạy dry-run (bắt buộc trước khi ghi file)

Luôn chạy `--dry-run` trước để xác nhận script parse đúng, **không ghi bất kỳ file nào**:

```bash
python3 scripts/00a_prepare_dataset.py \
  --root data/sars-cov2-time-forward-dataset --dry-run
```

**Output mong đợi (dạng minh họa):**
```
[T0/Alpha] sequences_..._Alpha_B.1.1.7_2000_31Dec2019_30Sep2021.fasta: declared=2000 actual=2000 [OK]
[T0/Beta]  sequences_..._Beta_B.1.351_2000_31Dec2019_30Sep2021.fasta:  declared=2000 actual=2000 [OK]
[T0/Gamma] sequences_..._Gamma_P.1_2000_31Dec2019_30Sep2021.fasta:     declared=2000 actual=2000 [OK]
[T0/Delta] sequences_..._Delta_B.1.617.2_2000_31Dec2019_30Sep2021.fasta: declared=2000 actual=2000 [OK]
...

=== Summary ===
  T0: 8000 sequences written
  T1: ...
  T2: ...
  T3: ...

[dry-run] no files written. Re-run without --dry-run once the summary above looks correct.
```

**Những gì cần kiểm tra trong output:**

| Dấu hiệu | Ý nghĩa | Hành động |
|---|---|---|
| `[OK]` trên tất cả dòng | Số sequence thực tế = số ghi trong tên file | Tiếp tục Bước 1.3 |
| `[MISMATCH]` trên một số dòng | Số sequence thực tế ≠ số ghi trong tên file | Kiểm tra lại nếu chênh nhiều — có thể file download bị thiếu |
| `[WARN] could not parse lineage from filename` | Tên file không khớp format chuẩn | Kiểm tra tên file thủ công |
| `[WARNING] duplicate accessions` | Có accession trùng giữa các file trong cùng period | Bình thường nếu số lượng nhỏ; kiểm tra nếu nhiều |

### 1.3. Chạy thật — sinh file gộp và metadata

Sau khi dry-run trông ổn, chạy thật:

```bash
python3 scripts/00a_prepare_dataset.py \
  --root data/sars-cov2-time-forward-dataset
```

Script sẽ tự sinh vào `data/`:

| File được tạo | Nội dung |
|---|---|
| `data/T0.fasta` | Toàn bộ genome T0 (Alpha + Beta + Gamma + Delta gộp lại) |
| `data/T1.fasta` | Toàn bộ genome T1 |
| `data/T2.fasta` | Toàn bộ genome T2 (chứa Omicron mới xuất hiện) |
| `data/T3.fasta` | Toàn bộ genome T3 |
| `data/metadata.csv` | Bảng ánh xạ `accession → lineage / period / date range` |

**Output mong đợi:**
```
...
=== Summary ===
  T0: 8000 sequences written -> data/T0.fasta
  T1: ...
  T2: ...
  T3: ...

[done] wrote 10000 rows to data/metadata.csv
[done] period FASTA files written to: data/T0.fasta, data/T1.fasta, data/T2.fasta, data/T3.fasta
```

**Định dạng `metadata.csv` được sinh ra** (4 cột bắt buộc):
```csv
accession,lineage,period,collection_date
MZ123456,Alpha,T0,31Dec2019-30Sep2021
MZ123457,Delta,T0,31Dec2019-30Sep2021
MZ123458,Omicron,T2,01Nov2021-31Jan2022
```

**Lưu ý về accession:** script tự động cắt bỏ phần version (`.1`, `.2`...) và mọi thứ
sau ký tự `|` trong FASTA header, nên header dạng `>OQ123456.1 Severe acute...` sẽ cho
accession `OQ123456`. Để kiểm tra nhanh header thật:
```bash
grep "^>" data/T0.fasta | head -3
head -3 data/metadata.csv
```

---

## Bước 2 — Giai đoạn 2: BLAST/Mash baseline

Chạy toàn bộ từ thư mục `pipeline/`. **Đảm bảo Bước 1 đã hoàn tất** và 4 file
`data/T*.fasta` + `data/metadata.csv` đã tồn tại trước khi bắt đầu.

### 2.1. Xây reference set đại diện từ T0

```bash
python3 scripts/01_build_reference_set.py \
  --fasta data/T0.fasta \
  --meta data/metadata.csv \
  --n-per-class 30 \
  --outdir results
```

**Output mong đợi:**
```
[reference] Alpha: picked 30/2000
[reference] Beta: picked 30/2000
[reference] Gamma: picked 30/2000
[reference] Delta: picked 30/2000
[reference] wrote 120 sequences to results/reference.fasta
```

Kiểm tra: `results/reference.fasta` (120 sequence) và `results/reference_metadata.csv`
đã được tạo trước khi chạy tiếp.

### 2.2. Tách tập calibration (T0 trừ reference set)

```bash
python3 scripts/00_split_t0_calib.py \
  --fasta data/T0.fasta \
  --ref-meta results/reference_metadata.csv \
  --out results/T0_calib.fasta
```

**Output mong đợi:**
```
[calib] wrote 7880 T0 sequences (excluding 120 reference accessions) to results/T0_calib.fasta
```

Tập này dùng để fit threshold một cách công bằng (không cho reference set "tự khớp
với chính nó").

### 2.3. Chạy BLAST + Mash cho từng giai đoạn

```bash
bash scripts/02_run_baselines.sh T0_calib results/T0_calib.fasta results/reference.fasta
bash scripts/02_run_baselines.sh T1 data/T1.fasta results/reference.fasta
bash scripts/02_run_baselines.sh T2 data/T2.fasta results/reference.fasta
bash scripts/02_run_baselines.sh T3 data/T3.fasta results/reference.fasta
```

Mỗi lệnh in ra thời gian chạy, ví dụ:
```
[baseline] done -> results/blast_T2.tsv , results/mash_T2.tsv
blast_seconds=12.34
mash_seconds=0.45
```

Với dữ liệu thật (~29.903 nt/genome, hàng nghìn sequence), BLAST cho T0_calib
(~7.880 genome) có thể mất vài phút — nếu lâu, chạy nền:
```bash
nohup bash scripts/02_run_baselines.sh T0_calib results/T0_calib.fasta results/reference.fasta > log_t0calib.txt 2>&1 &
```

**Kiểm tra output:** mỗi period phải có đủ 2 file `blast_<period>.tsv` và
`mash_<period>.tsv` trong `results/` trước khi sang bước tiếp theo — nếu thiếu, xem
lại log lỗi (`makeblastdb`/`blastn`/`mash` in lỗi ra `stderr`).

### 2.4. Tính threshold + rejection rate + AUROC/FPR95

```bash
python3 scripts/03_compute_metrics.py --meta data/metadata.csv --method blast
python3 scripts/03_compute_metrics.py --meta data/metadata.csv --method mash
```

**Output mong đợi (ví dụ minh họa, số thật sẽ khác):**
```
[blast] threshold (novelty score, 95.0th pct of T0_calib) = 2.15
  T1 Alpha: n=200 rejection_rate=0.020
  T1 Delta: n=1000 rejection_rate=0.015
  T2 Omicron: n=27 rejection_rate=0.740
  T3 Omicron: n=1000 rejection_rate=0.910
[blast] AUROC (T1 IND vs pooled Omicron OOD) = 0.9350
[blast] FPR95 = 0.0620
[blast] total wall-clock over T1+T2+T3 = 45.20s
```

**File kết quả cần lấy ra để đưa vào bài:**
- `results/table6_blast.csv`, `results/table6_mash.csv` → nội dung mở rộng cho Bảng 6.
- `results/table7_comparison.csv` → 2 dòng mới (blast, mash) để đặt cạnh 4 baseline
  OOD hiện có trong Bảng 4, cộng với dòng encoder (AUROC=0.990, FPR95=4.36% theo
  Section 4.11) tự thêm vào để có bảng so sánh đầy đủ 6 phương pháp.

Mở file để xem:
```bash
cat results/table7_comparison.csv
column -s, -t results/table6_blast.csv | less
```

---

## Bước 3 — Giai đoạn 3: Robustness-to-noise

### 3.1–3.2. Sinh nhiễu + chạy BLAST/Mash trên từng mức nhiễu

```bash
bash scripts/05_run_robustness.sh T2 data/T2.fasta
```

Lệnh này tự động:
- Sinh 8 phiên bản của `T2.fasta`: `level_0` (sạch), `sub_0.005/0.01/0.02/0.05`
  (substitution 0.5–5%), `trunc_0.05/0.1/0.2` (xóa liên tục 5–20% chiều dài).
- Chạy BLAST + Mash cho từng phiên bản, lưu vào `results/robustness_T2/`.

**Output mong đợi:**
```
[noise] wrote substitution level 0.005 -> results/noise_T2/sub_0.005.fasta
...
=== level: sub_0.05 ===
[robustness] done -> results/robustness_T2
```

Với genome dài thật (29.903 nt) và nhiều sequence, bước này có thể mất 5–15 phút tùy
số lượng — nên chạy nền như gợi ý ở Bước 2.3 nếu cần.

Có thể lặp lại cho T3 nếu muốn thêm dữ liệu củng cố:
```bash
bash scripts/05_run_robustness.sh T3 data/T3.fasta
```

### 3.3. Tổng hợp accuracy theo mức nhiễu (cho BLAST và Mash)

```bash
python3 scripts/06_aggregate_robustness.py --meta data/metadata.csv --period T2 --method blast
python3 scripts/06_aggregate_robustness.py --meta data/metadata.csv --period T2 --method mash
```

**Output mong đợi:**
```
[blast] level_0: accuracy=0.9950 (n=1027)
[blast] sub_0.005: accuracy=0.9820 (n=1027)
[blast] sub_0.01: accuracy=0.9410 (n=1027)
...
[plot] wrote figures/figure7_robustness.png
```

Ở bước này, `figures/figure7_robustness.png` đã có 2 đường (blast, mash) — đủ để xem
sơ bộ, nhưng **chưa đầy đủ** vì thiếu đường của encoder.

### 3.4. Thêm kết quả encoder vào cùng biểu đồ (bắt buộc để hoàn thiện Figure 7)

Đây là bước duy nhất **KHÔNG** chạy bằng script có sẵn — cần dùng chính pipeline
encoder + Mahalanobis hiện tại:

1. Chạy encoder trên đúng 8 file nhiễu đã sinh sẵn:
   ```
   results/noise_T2/level_0.fasta
   results/noise_T2/sub_0.005.fasta
   results/noise_T2/sub_0.01.fasta
   results/noise_T2/sub_0.02.fasta
   results/noise_T2/sub_0.05.fasta
   results/noise_T2/trunc_0.05.fasta
   results/noise_T2/trunc_0.1.fasta
   results/noise_T2/trunc_0.2.fasta
   ```
   (accession ID giữ nguyên như T2.fasta gốc nên khớp thẳng với `metadata.csv`,
   không cần xử lý gì thêm).

2. Với mỗi file, dùng đúng threshold θ=13.91 đã fix từ T0 (không fit lại), tính tỷ lệ
   phân loại đúng (known → không bị reject; Omicron → bị reject) — chính là "accuracy"
   cùng định nghĩa với 2 baseline kia.

3. Thêm 8 dòng vào `results/robustness_summary.csv` theo đúng format:
   ```csv
   method,period,level,accuracy,n
   encoder,T2,level_0,0.9960,1027
   encoder,T2,sub_0.005,0.9930,1027
   encoder,T2,sub_0.01,0.9880,1027
   encoder,T2,sub_0.02,0.9790,1027
   encoder,T2,sub_0.05,0.9510,1027
   encoder,T2,trunc_0.05,0.9940,1027
   encoder,T2,trunc_0.1,0.9820,1027
   encoder,T2,trunc_0.2,0.9340,1027
   ```
   (số liệu trên chỉ là ví dụ minh họa định dạng — thay bằng số thật).

4. Vẽ lại biểu đồ với đủ 3 đường:
   ```bash
   python3 scripts/06_aggregate_robustness.py --plot-only
   ```

**Output cuối cùng:** `figures/figure7_robustness.png` với 3 đường (encoder, blast,
mash) — đây chính là Figure 7 cần chèn vào Section 4.13 của bài báo.

---

## Tổng quan toàn bộ quy trình

```
Bước 1: Tiền xử lý dữ liệu thật
  └─ 00a_prepare_dataset.py --dry-run   (kiểm tra)
  └─ 00a_prepare_dataset.py             (ghi file)
      ↓
      data/T0.fasta  data/T1.fasta  data/T2.fasta  data/T3.fasta
      data/metadata.csv

Bước 2: Giai đoạn 2 — BLAST/Mash baseline
  └─ 01_build_reference_set.py          → results/reference.fasta
  └─ 00_split_t0_calib.py              → results/T0_calib.fasta
  └─ 02_run_baselines.sh (×4 period)   → results/blast_T*.tsv, mash_T*.tsv
  └─ 03_compute_metrics.py (×2 method) → results/table6_*.csv, table7_comparison.csv

Bước 3: Giai đoạn 3 — Robustness-to-noise
  └─ 05_run_robustness.sh T2           → results/noise_T2/*.fasta, robustness_T2/
  └─ 06_aggregate_robustness.py (×2)   → results/robustness_summary.csv
  └─ [manual] thêm dòng encoder        → results/robustness_summary.csv
  └─ 06_aggregate_robustness.py --plot-only → figures/figure7_robustness.png
```

## Bảng tổng hợp: file nào dùng cho phần nào của bài báo

| File output | Dùng cho |
|---|---|
| `results/table6_blast.csv`, `table6_mash.csv` | Mở rộng Bảng 6 (rejection rate theo lineage/period) |
| `results/table7_comparison.csv` | Bảng 7 mới — so sánh AUROC/FPR95/runtime giữa encoder, BLAST, Mash |
| `results/robustness_summary.csv` | Số liệu nguồn cho Figure 7 |
| `figures/figure7_robustness.png` | Figure 7 — đường cong suy giảm theo mức nhiễu |

---

## Xử lý lỗi thường gặp

| Triệu chứng | Nguyên nhân khả dĩ | Cách sửa |
|---|---|---|
| `[WARN] could not parse lineage from filename` (Bước 1.2) | Tên file không đúng format chuẩn `sequences_...__Lineage_Pango_count_start_end.fasta` | Kiểm tra tên file thủ công; nếu cần đổi tên hàng loạt dùng `rename` hoặc script nhỏ |
| `[MISMATCH]` nhiều dòng (Bước 1.2) | File download thiếu sequence so với số ghi trong tên | Download lại file đó; nếu chênh ít (1–2) có thể bỏ qua |
| `[WARNING] ... accessions were not found in the FASTA file` (Bước 2.1) | Định dạng accession trong `metadata.csv` không khớp header FASTA | So sánh header thật (`grep "^>" data/T0.fasta \| head`) với cột `accession` trong metadata; thường do header GISAID dùng định dạng EPI_ISL — cần re-run `00a_prepare_dataset.py` với logic parse phù hợp |
| `FileNotFoundError: results/blast_T*.tsv` (Bước 2.4) | Bước 2.3 chưa chạy xong hoặc bị lỗi giữa chừng | Chạy lại từng lệnh `bash scripts/02_run_baselines.sh ...` và đọc kỹ log |
| BLAST chạy rất lâu | Genome dài (~30kb) × hàng nghìn sequence × T0_calib lớn | Chạy nền (`nohup ... &`), hoặc giảm tạm `--n-per-class` để test nhanh trước khi chạy full |
| Mash cho kết quả kỳ lạ / distance toàn giống nhau | Quên cờ `-i` khi sketch | Đã được fix sẵn trong `02_run_baselines.sh`/`05_run_robustness.sh`, không cần sửa nếu dùng đúng script gửi kèm |
| AUROC = 1.0 hoặc 0.5 bất thường | Dữ liệu test quá nhỏ (vài chục sequence) | Bình thường khi test với dữ liệu mẫu nhỏ; sẽ mượt hơn với dữ liệu thật đầy đủ (T2 n≈1.027) |

---

## Lưu ý quan trọng trước khi chạy dữ liệu thật

1. **Xác nhận định dạng FASTA header:** nếu file gốc dùng GISAID EPI_ISL thay vì
   accession NCBI, hãy upload thử 1 file FASTA thật để kiểm tra header trước khi chạy
   full — tránh phải chạy lại từ đầu nếu parse sai.

2. **`--n-per-class` của reference set** (mặc định 30): cân bằng giữa đại diện đủ tốt
   và tinh thần reference nhỏ — có thể thử vài giá trị (10/30/50) và báo cáo độ nhạy
   nếu muốn chắc chắn hơn.

3. **BLAST hay Mash là baseline chính?** Nên báo cáo cả hai — BLAST trả lời trực tiếp
   tên công cụ reviewer nêu, Mash củng cố luận điểm về tốc độ.

4. **Robustness: T2, T3, hay cả hai?** T2 (Omicron mới xuất hiện, n nhỏ) là quan trọng
   nhất về mặt câu chuyện early-warning; T3 (Omicron áp đảo) có thể chạy thêm nếu còn
   thời gian để củng cố.

---

Nếu trong lúc chạy với dữ liệu thật gặp lỗi cụ thể, gửi đoạn log lỗi (và nếu tiện,
upload trực tiếp 1 file FASTA thật vào đây), tôi sẽ debug trực tiếp cùng.
