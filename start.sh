#!/bin/bash
# ViralGen AI - Start Script
# Script khởi động ứng dụng

echo "============================================================"
echo "🧬 ViralGen AI - Hệ Thống Giám Sát Virus Thế Giới Mở"
echo "============================================================"

# Kiểm tra Python
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 không được cài đặt"
    exit 1
fi

# Kiểm tra pip
if ! command -v pip3 &> /dev/null; then
    echo "❌ pip3 không được cài đặt"
    exit 1
fi

# Cài đặt dependencies nếu cần
echo "📦 Kiểm tra dependencies..."
pip3 install -r requirements.txt -q

# Tạo thư mục uploads nếu chưa tồn tại
mkdir -p app/uploads

# Khởi động ứng dụng
echo ""
echo "🚀 Khởi động ứng dụng..."
echo "📋 Hướng dẫn:"
echo "  1. Mở trình duyệt và truy cập http://localhost:5000"
echo "  2. Upload file FASTA/CSV hoặc chọn file test mẫu"
echo "  3. Xem kết quả phân tích biến thể virus"
echo ""
echo "============================================================"

# Chạy ứng dụng
python3 run.py