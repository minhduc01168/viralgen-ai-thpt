#!/usr/bin/env python3
"""
ViralGen AI - Ứng dụng Web Flask
Load model từ weight/ folder (không cần huấn luyện lại)
"""

import os
import sys
import json
from datetime import datetime
from flask import Flask, render_template, request, jsonify, redirect, url_for, flash
from werkzeug.utils import secure_filename

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from predict import VirusAnalyzer

app = Flask(__name__)
app.secret_key = 'viralgen-ai-secret-key-2026'

UPLOAD_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uploads')
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'csv'}
MAX_CONTENT_LENGTH = 50 * 1024 * 1024

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = MAX_CONTENT_LENGTH

os.makedirs(UPLOAD_FOLDER, exist_ok=True)

analyzer = None


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def init_analyzer():
    """Load analyzer từ weight/ folder"""
    global analyzer
    weight_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'weight')
    try:
        analyzer = VirusAnalyzer(weight_dir=weight_dir)
        print(f"✅ Đã load analyzer từ: {weight_dir}")
        return True
    except FileNotFoundError as e:
        print(f"❌ {e}")
        print("Hãy chạy: python train_pipeline.py")
        return False
    except Exception as e:
        print(f"❌ Lỗi load analyzer: {e}")
        import traceback
        traceback.print_exc()
        return False


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        flash('Không tìm thấy file', 'error')
        return redirect(request.url)

    file = request.files['file']
    if file.filename == '':
        flash('Chưa chọn file', 'error')
        return redirect(request.url)

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        filename = f"{timestamp}_{filename}"
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)
        flash('Upload file thành công!', 'success')
        return redirect(url_for('analyze', filename=filename))

    flash('Định dạng file không hỗ trợ. Chỉ chấp nhận .fasta, .fa, .csv', 'error')
    return redirect(url_for('index'))


@app.route('/analyze/<filename>')
def analyze(filename):
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    if not os.path.exists(filepath):
        flash('File không tồn tại', 'error')
        return redirect(url_for('index'))

    try:
        result = analyzer.predict_variant(filepath)
        summary = analyzer.get_analysis_summary(result)
        return render_template('result.html',
                             filename=filename,
                             result=result,
                             summary=summary,
                             has_model=True)
    except Exception as e:
        flash(f'Lỗi khi phân tích: {str(e)}', 'error')
        return redirect(url_for('index'))


@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    if 'file' not in request.files:
        return jsonify({'error': 'Không tìm thấy file'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'Chưa chọn file'}), 400

    if not allowed_file(file.filename):
        return jsonify({'error': 'Định dạng file không hỗ trợ'}), 400

    try:
        filename = secure_filename(file.filename)
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        temp_filename = f"{timestamp}_{filename}"
        temp_path = os.path.join(app.config['UPLOAD_FOLDER'], temp_filename)
        file.save(temp_path)

        result = analyzer.predict_variant(temp_path)
        os.remove(temp_path)

        return jsonify({'success': True, 'result': result})
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/test-files')
def get_test_files():
    test_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'dataset_test')
    files = []
    for root, dirs, filenames in os.walk(test_dir):
        for filename in filenames:
            if filename.endswith(('.fasta', '.fa', '.csv')):
                filepath = os.path.join(root, filename)
                rel_path = os.path.relpath(filepath, test_dir)
                files.append({
                    'name': filename,
                    'path': rel_path,
                    'size': os.path.getsize(filepath),
                    'type': 'IND' if 'IND' in root else 'OOD'
                })
    return jsonify({'files': files})


@app.route('/test/<path:filepath>')
def test_file(filepath):
    test_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'dataset_test')
    full_path = os.path.join(test_dir, filepath)

    if not os.path.exists(full_path):
        flash('File test không tồn tại', 'error')
        return redirect(url_for('index'))

    try:
        result = analyzer.predict_variant(full_path)
        summary = analyzer.get_analysis_summary(result)
        return render_template('result.html',
                             filename=os.path.basename(full_path),
                             result=result,
                             summary=summary,
                             has_model=True)
    except Exception as e:
        flash(f'Lỗi khi phân tích: {str(e)}', 'error')
        return redirect(url_for('index'))


@app.errorhandler(413)
def too_large(e):
    flash('File quá lớn. Kích thước tối đa là 50MB', 'error')
    return redirect(url_for('index'))


@app.errorhandler(404)
def not_found(e):
    return render_template('404.html'), 404


if __name__ == '__main__':
    print("🚀 Khởi động ViralGen AI Web Application...")

    if init_analyzer():
        print("✅ Analyzer đã sẵn sàng")
    else:
        print("❌ Không thể load analyzer. Hãy chạy: python train_pipeline.py")
        sys.exit(1)

    app.run(host='0.0.0.0', port=5000, debug=True)