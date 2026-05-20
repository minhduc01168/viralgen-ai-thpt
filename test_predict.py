#!/usr/bin/env python3
"""
ViralGen AI - Test Script
Test module predict với các file trong dataset_test (load từ weight/)
"""

import os
import sys
import json
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from predict import VirusAnalyzer


def test_predict():
    """Test predict với tất cả file trong dataset_test"""

    print("=" * 60)
    print("🧪 ViralGen AI - Test Module Predict")
    print("=" * 60)

    # Load analyzer từ weight/
    weight_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'weight')
    try:
        analyzer = VirusAnalyzer(weight_dir=weight_dir)
    except FileNotFoundError as e:
        print(f"❌ {e}")
        print("Hãy chạy: python train_pipeline.py")
        return

    # Thu thập file test
    test_files = []
    base = os.path.dirname(os.path.abspath(__file__))

    ind_dir = os.path.join(base, 'dataset_test', 'IND')
    ood_dir = os.path.join(base, 'dataset_test', 'OOD')

    for folder, ftype in [(ind_dir, 'IND'), (ood_dir, 'OOD')]:
        if os.path.exists(folder):
            for filename in sorted(os.listdir(folder)):
                if filename.endswith(('.fasta', '.fa')):
                    test_files.append({
                        'path': os.path.join(folder, filename),
                        'name': filename,
                        'type': ftype
                    })

    print(f"\n📁 Tìm thấy {len(test_files)} file test:")
    for i, fi in enumerate(test_files, 1):
        print(f"  {i}. {fi['name']} ({fi['type']})")

    # Test từng file
    results = []
    print("\n" + "=" * 60)
    print("🔬 Bắt đầu test phân tích...")
    print("=" * 60)

    for i, fi in enumerate(test_files, 1):
        print(f"\n📄 Test {i}/{len(test_files)}: {fi['name']}")
        print("-" * 40)
        try:
            result = analyzer.predict_variant(fi['path'])
            print(f"  Biến thể: {result['variant']}")
            print(f"  Độ tự tin: {result['confidence']}%")
            print(f"  Mahalanobis: {result['mahalanobis_distance']}")
            print(f"  UNKNOWN: {result['is_unknown']}")

            result['test_type'] = fi['type']
            result['expected_variant'] = fi['name'].split('_')[-1].replace('.fasta', '').replace('.fa', '')
            results.append(result)
        except Exception as e:
            print(f"❌ Lỗi: {e}")
            import traceback
            traceback.print_exc()

    # Tóm tắt
    print("\n" + "=" * 60)
    print("📊 TÓM TẮT KẾT QUẢ TEST")
    print("=" * 60)

    if results:
        ind_results = [r for r in results if r['test_type'] == 'IND']
        ood_results = [r for r in results if r['test_type'] == 'OOD']

        print(f"\n📈 Thống kê:")
        print(f"  - Tổng: {len(results)}")
        print(f"  - IND: {len(ind_results)}")
        print(f"  - OOD: {len(ood_results)}")

        if ind_results:
            correct_ind = sum(1 for r in ind_results if not r['is_unknown'])
            print(f"  - IND đúng: {correct_ind}/{len(ind_results)} ({correct_ind/len(ind_results)*100:.1f}%)")

        if ood_results:
            correct_ood = sum(1 for r in ood_results if r['is_unknown'])
            print(f"  - OOD đúng: {correct_ood}/{len(ood_results)} ({correct_ood/len(ood_results)*100:.1f}%)")

        # Chi tiết từng file
        print(f"\n📋 Chi tiết:")
        for r in results:
            status = "✅" if (r['test_type'] == 'IND' and not r['is_unknown']) or (r['test_type'] == 'OOD' and r['is_unknown']) else "❌"
            print(f"  {status} [{r['test_type']}] {r['file_name']} → {r['variant']} (dist={r['mahalanobis_distance']})")

        # Lưu JSON
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        result_file = f"test_results_{timestamp}.json"
        with open(result_file, 'w', encoding='utf-8') as f:
            json.dump({'timestamp': timestamp, 'total': len(results), 'results': results}, f, indent=2, ensure_ascii=False)
        print(f"\n💾 Đã lưu: {result_file}")

    print("\n" + "=" * 60)
    print("✅ Hoàn thành test!")
    print("=" * 60)


if __name__ == '__main__':
    test_predict()