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

## Hướng dẫn cài đặt

### Yêu cầu hệ thống

- **Python**: 3.9 trở lên (khuyến nghị 3.10+)
- **RAM**: Tối thiểu 4GB (khuyến nghị 8GB)
- **Disk**: Tối thiểu 2GB trống
- **Git**: Để clone repository

### Windows

```bash
# 1. Cài Python từ https://www.python.org/downloads/
#    Đánh dấu chọn "Add Python to PATH" khi cài đặt

# 2. Clone project
git clone <repository-url>
cd viralgen-ai-app

# 3. Tạo virtual environment (khuyến nghị)
python -m venv venv
venv\Scripts\activate

# 4. Cài đặt dependencies
pip install -r requirements.txt

# 5. Chạy ứng dụng (script tự động)
start.bat
```

### macOS

```bash
# 1. Cài Homebrew (nếu chưa có)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# 2. Cài Python
brew install python@3.10

# 3. Clone project
git clone <repository-url>
cd viralgen-ai-app

# 4. Tạo virtual environment (khuyến nghị)
python3 -m venv venv
source venv/bin/activate

# 5. Cài đặt dependencies
pip install -r requirements.txt

# 6. Chạy ứng dụng
chmod +x start.sh
./start.sh
```

### Ubuntu / Debian

```bash
# 1. Cập nhật package và cài Python
sudo apt update
sudo apt install python3 python3-pip python3-venv git -y

# 2. Clone project
git clone <repository-url>
cd viralgen-ai-app

# 3. Tạo virtual environment (khuyến nghị)
python3 -m venv venv
source venv/bin/activate

# 4. Cài đặt dependencies
pip install -r requirements.txt

# 5. Chạy ứng dụng
chmod +x start.sh
./start.sh
```

---

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