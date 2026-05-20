#!/usr/bin/env python3
"""
ViralGen AI - Run Script
Khởi động ứng dụng Web Flask (load model từ weight/)
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from app.app import app, init_analyzer


def main():
    print("=" * 60)
    print("🧬 ViralGen AI - Hệ Thống Giám Sát Virus Thế Giới Mở")
    print("=" * 60)

    print("\n📦 Đang load model từ weight/...")
    if init_analyzer():
        print("✅ Analyzer đã sẵn sàng")
    else:
        print("❌ Không thể load analyzer.")
        print("Hãy chạy: python train_pipeline.py")
        sys.exit(1)

    host = os.environ.get('HOST', '0.0.0.0')
    port = int(os.environ.get('PORT', 5000))

    print(f"\n🌐 Server: http://{host}:{port}")
    print("=" * 60)
    print("📋 Hướng dẫn:")
    print("  1. Mở trình duyệt → http://localhost:5000")
    print("  2. Upload file .fasta/.csv hoặc chọn file test mẫu")
    print("  3. Xem kết quả phân tích biến thể virus")
    print("=" * 60)

    try:
        app.run(host=host, port=port, debug=True)
    except KeyboardInterrupt:
        print("\n🛑 Đã dừng server.")
    except Exception as e:
        print(f"\n❌ Lỗi: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()