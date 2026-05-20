@echo off
REM ViralGen AI - Start Script cho Windows
REM Script khởi động ứng dụng

echo ============================================================
echo ViralGen AI - Hệ Thống Giám Sát Virus Thế Giới Mở
echo ============================================================

REM Kiểm tra Python
python --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Python không được cài đặt
    pause
    exit /b 1
)

REM Cài đặt dependencies nếu cần
echo 📦 Kiểm tra dependencies...
pip install -r requirements.txt -q

REM Tạo thư mục uploads nếu chưa tồn tại
if not exist "app\uploads" mkdir app\uploads

REM Khởi động ứng dụng
echo.
echo 🚀 Khởi động ứng dụng...
echo 📋 Hướng dẫn:
echo   1. Mở trình duyệt và truy cập http://localhost:5000
echo   2. Upload file FASTA/CSV hoặc chọn file test mẫu
echo   3. Xem kết quả phân tích biến thể virus
echo.
echo ============================================================

REM Chạy ứng dụng
python run.py

pause