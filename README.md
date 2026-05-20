# ViralGen AI - Hệ Thống Giám Sát Virus Thế Giới Mở

Hệ thống giám sát virus sử dụng trí tuệ nhân tạo để phân tích và định danh biến thể SARS-CoV-2, với khả năng phát hiện biến thể lạ (Out-of-Distribution Detection).

## Cấu trúc thư mục

```
viralgen-ai-app/
├── app/                          # Ứng dụng Web Flask
│   ├── app.py                    # Flask application
│   ├── templates/                # HTML templates
│   ├── static/css/style.css      # CSS
│   └── uploads/                  # Thư mục upload
├── dataset_test/                 # Dữ liệu test mẫu
│   ├── IND/                      # Biến thể đã biết (Alpha, Beta, Delta, Gamma, Omicron)
│   └── OOD/                      # Biến thể lạ (Lambda, Zeta, Mu, Eta)
├── weight/                       # Model weights & artifacts (tạo bởi train_pipeline.py)
│   ├── open_world_encoder.pth    # PyTorch model weights
│   ├── tfidf.pkl                 # TF-IDF transformer
│   ├── scaler.pkl                # StandardScaler
│   ├── label_encoder.pkl         # Label Encoder
│   ├── class_stats.pkl           # Mahalanobis statistics
│   └── metadata.pkl              # Config metadata
├── train_pipeline.py             # Script train model (chạy 1 lần)
├── predict.py                    # Module dự đoán (load từ weight/)
├── test_predict.py               # Script test
├── run.py                        # Script chạy web app
└── requirements.txt              # Dependencies
```

## Cách sử dụng

### Bước 1: Cài đặt dependencies

```bash
pip install -r requirements.txt
```

### Bước 2: Train model (chạy 1 lần)

```bash
python train_pipeline.py
```

Script sẽ:
- Đọc dữ liệu từ `dataset_test/IND/`
- Train TF-IDF, StandardScaler, Label Encoder
- Train model OpenWorldVirusEncoder (20 epochs)
- Tính Mahalanobis statistics
- Lưu tất cả vào folder `weight/`

### Bước 3: Chạy ứng dụng Web

```bash
python run.py
```

Truy cập: http://localhost:5000

### Bước 4: Test nhanh (command line)

```bash
# Test tất cả file
python test_predict.py

# Test 1 file cụ thể
python predict.py --input dataset_test/IND/sequences_20260517_1817244_B.1.1.7_Alpha.fasta
```

## API Endpoints

| Endpoint | Method | Mô tả |
|----------|--------|-------|
| `/` | GET | Trang chủ (upload file) |
| `/upload` | POST | Upload file để phân tích |
| `/analyze/<filename>` | GET | Kết quả phân tích |
| `/api/analyze` | POST | API phân tích (JSON) |
| `/api/test-files` | GET | Danh sách file test |
| `/test/<path>` | GET | Test file mẫu |

## Kết quả phân tích

Hệ thống trả về:

- **Biến thể dự đoán**: Alpha, Beta, Delta, Gamma, Omicron hoặc UNKNOWN_VARIANT
- **Độ tự tin**: Phần trăm tin cậy (0-100%)
- **Khoảng cách Mahalanobis**: Đo lường sự khác biệt với các biến thể đã biết
- **MC Variance**: Độ bất định của mô hình
- **Similarity Scores**: Chỉ số tương đồng với từng biến thể
- **Mutation Sites**: Các vị trí đột biến phát hiện

## Công nghệ

- **Backend**: Python, Flask
- **ML/DL**: PyTorch, scikit-learn, SciPy
- **Bioinformatics**: Biopython
- **Frontend**: Bootstrap 5, Chart.js