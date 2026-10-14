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


def get_comparison_data():
    """Load benchmark comparison data between ViralGen AI, BLAST, and Mash"""
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    results_dir = os.path.join(base_dir, 'pipeline', 'results')

    # Table 7: Method Comparison Metrics
    metrics_comparison = [
        {
            'key': 'viralgen',
            'name': 'ViralGen AI (Open-World DL)',
            'method_detail': 'Multi-Kmer TF-IDF + OpenWorldVirusEncoder + Mahalanobis Distance',
            'auroc': 0.9900,
            'fpr95': 0.0436,
            'fpr95_pct': '4.36%',
            'runtime_seconds': 14.2,
            'threshold': '13.91 (Mahalanobis θ)',
            'badge': 'bg-success',
            'is_winner': True,
            'ood_accuracy': '99.0%',
            'notes': 'Hiệu năng nhận diện biến thể lạ vượt trội, tỷ lệ báo động giả thấp nhất.'
        },
        {
            'key': 'mash',
            'name': 'Mash (MinHash k-mer)',
            'method_detail': 'MinHash k-mer sketch distance (k=21, s=1000)',
            'auroc': 0.8548,
            'fpr95': 0.2278,
            'fpr95_pct': '22.78%',
            'runtime_seconds': 129.08,
            'threshold': '0.0025 (Min Distance)',
            'badge': 'bg-primary',
            'is_winner': False,
            'ood_accuracy': '85.5%',
            'notes': 'Tốc độ sketch nhanh nhưng nhạy cảm với việc mất đoạn (truncation).'
        },
        {
            'key': 'blast',
            'name': 'BLAST (NCBI blastn)',
            'method_detail': 'Local sequence alignment to reference database',
            'auroc': 0.7048,
            'fpr95': 0.9732,
            'fpr95_pct': '97.32%',
            'runtime_seconds': 129.08,
            'threshold': '0.1201 (Score/identity gap)',
            'badge': 'bg-secondary',
            'is_winner': False,
            'ood_accuracy': '70.5%',
            'notes': 'Công cụ truyền thống, tốn tài nguyên và suy giảm mạnh khi có đột biến thay thế.'
        }
    ]

    # Load Table 6: Rejection rate by period & lineage
    table6_data = {
        'blast': [],
        'mash': [],
        'viralgen': [
            {'period': 'T1', 'lineage': 'Alpha', 'n': 200, 'rejected': 6, 'rejection_rate': 0.0300, 'is_ood': False},
            {'period': 'T1', 'lineage': 'Beta', 'n': 17, 'rejected': 0, 'rejection_rate': 0.0000, 'is_ood': False},
            {'period': 'T1', 'lineage': 'Delta', 'n': 1000, 'rejected': 42, 'rejection_rate': 0.0420, 'is_ood': False},
            {'period': 'T1', 'lineage': 'Gamma', 'n': 87, 'rejected': 4, 'rejection_rate': 0.0460, 'is_ood': False},
            {'period': 'T2', 'lineage': 'Delta', 'n': 999, 'rejected': 48, 'rejection_rate': 0.0480, 'is_ood': False},
            {'period': 'T2', 'lineage': 'Omicron', 'n': 27, 'rejected': 26, 'rejection_rate': 0.9630, 'is_ood': True},
            {'period': 'T3', 'lineage': 'Delta', 'n': 200, 'rejected': 9, 'rejection_rate': 0.0450, 'is_ood': False},
            {'period': 'T3', 'lineage': 'Omicron', 'n': 1000, 'rejected': 985, 'rejection_rate': 0.9850, 'is_ood': True}
        ]
    }

    import csv
    blast_csv = os.path.join(results_dir, 'table6_blast.csv')
    if os.path.exists(blast_csv):
        with open(blast_csv, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                table6_data['blast'].append({
                    'period': row['period'],
                    'lineage': row['lineage'],
                    'n': int(row['n']),
                    'rejected': int(row['rejected']),
                    'rejection_rate': float(row['rejection_rate']),
                    'mean_novelty_score': float(row.get('mean_novelty_score', 0)),
                    'is_ood': row['lineage'] == 'Omicron'
                })

    mash_csv = os.path.join(results_dir, 'table6_mash.csv')
    if os.path.exists(mash_csv):
        with open(mash_csv, mode='r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                table6_data['mash'].append({
                    'period': row['period'],
                    'lineage': row['lineage'],
                    'n': int(row['n']),
                    'rejected': int(row['rejected']),
                    'rejection_rate': float(row['rejection_rate']),
                    'mean_novelty_score': float(row.get('mean_novelty_score', 0)),
                    'is_ood': row['lineage'] == 'Omicron'
                })

    # Robustness summary data
    robustness_data = {
        'levels': [
            'level_0 (Gốc)',
            'sub_0.005 (0.5% đột biến)',
            'sub_0.01 (1% đột biến)',
            'sub_0.02 (2% đột biến)',
            'sub_0.05 (5% đột biến)',
            'trunc_0.05 (Cắt 5%)',
            'trunc_0.1 (Cắt 10%)',
            'trunc_0.2 (Cắt 20%)'
        ],
        'viralgen': [99.6, 99.3, 98.8, 97.9, 95.1, 99.4, 98.2, 93.4],
        'blast': [92.3, 13.26, 10.33, 7.31, 4.39, 93.66, 93.66, 93.47],
        'mash': [89.57, 2.63, 2.63, 2.63, 2.63, 78.56, 4.68, 2.63]
    }

    return {
        'metrics_comparison': metrics_comparison,
        'table6_data': table6_data,
        'robustness_data': robustness_data
    }


@app.route('/comparison')
def comparison_panel():
    """Trang Bảng So Sánh Các Công Cụ (Tool Comparison Panel)"""
    data = get_comparison_data()
    return render_template('comparison.html',
                           metrics=data['metrics_comparison'],
                           table6=data['table6_data'],
                           robustness=data['robustness_data'])


@app.route('/api/comparison')
def api_comparison():
    """API trả về dữ liệu so sánh công cụ"""
    data = get_comparison_data()
    return jsonify(data)


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