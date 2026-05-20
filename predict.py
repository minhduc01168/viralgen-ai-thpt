#!/usr/bin/env python3
"""
ViralGen AI - Module dự đoán biến thể virus SARS-CoV-2
Load toàn bộ artifacts từ folder weight/ (không cần huấn luyện lại)
"""

import os
import itertools
import pickle
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO
import scipy.linalg
from scipy.spatial.distance import mahalanobis

import torch
import torch.nn as nn

import warnings
warnings.filterwarnings('ignore')


# ============================================================
# 1. MÔ HÌNH
# ============================================================
class OpenWorldVirusEncoder(nn.Module):
    def __init__(self, input_dim, num_classes):
        super(OpenWorldVirusEncoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(0.4),
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 128),
            nn.ReLU()
        )
        self.classifier = nn.Linear(128, num_classes)

    def forward(self, x):
        embedding = self.encoder(x)
        logits = self.classifier(embedding)
        return logits, embedding


# ============================================================
# 2. HÀM TIỆN ÍCH
# ============================================================
NUCLEOTIDES = ['A', 'C', 'G', 'T']
KMERS_3 = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=3)]
KMERS_4 = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=4)]
KMERS_5 = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=5)]
ALL_MULTI_KMERS = KMERS_3 + KMERS_4 + KMERS_5


def count_multi_kmers(sequence):
    seq_str = str(sequence).upper()
    kmer_counts = Counter()
    for k in [3, 4, 5]:
        kmers_in_seq = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
        kmer_counts.update(kmers_in_seq)
    return [kmer_counts.get(kmer, 0) for kmer in ALL_MULTI_KMERS]


def parse_fasta_description(description):
    """Parse mô tả FASTA thành thông tin dễ hiểu cho người dùng."""
    info = {
        'virus_name': 'SARS-CoV-2',
        'virus_name_vi': 'Virus Corona gây bệnh COVID-19',
        'sample_id': '',
        'sample_type': 'Mẫu vật',
        'origin': '',
        'year': '',
        'simple_desc': '',
    }

    desc = description.strip()

    # Tách ID (phần đầu trước khoảng trắng)
    parts = desc.split(None, 1)
    if parts:
        info['sample_id'] = parts[0]

    # Tìm "isolate" hoặc "strain"
    desc_lower = desc.lower()
    if 'isolate' in desc_lower:
        idx = desc_lower.index('isolate')
        after = desc[idx+7:].strip().split()
        if after:
            info['sample_type'] = f'Mẫu phân lập {after[0]}'
    elif 'strain' in desc_lower:
        idx = desc_lower.index('strain')
        after = desc[idx+6:].strip().split()
        if after:
            info['sample_type'] = f'Chủng {after[0]}'

    # Tìm năm (20xx)
    import re
    year_match = re.search(r'20\d{2}', desc)
    if year_match:
        info['year'] = year_match.group()

    # Tìm quốc gia/từ khóa địa lý
    countries = {
        'USA': 'Mỹ', 'UK': 'Anh', 'India': 'Ấn Độ', 'Brazil': 'Brazil',
        'China': 'Trung Quốc', 'Japan': 'Nhật Bản', 'Korea': 'Hàn Quốc',
        'Germany': 'Đức', 'France': 'Pháp', 'Italy': 'Ý', 'Spain': 'Tây Ban Nha',
        'Australia': 'Úc', 'Canada': 'Canada', 'Mexico': 'Mexico',
        'South Africa': 'Nam Phi', 'Nigeria': 'Nigeria', 'Thailand': 'Thái Lan',
        'Vietnam': 'Việt Nam', 'Singapore': 'Singapore', 'Malaysia': 'Malaysia',
        'Indonesia': 'Indonesia', 'Philippines': 'Philippines',
        'Switzerland': 'Thụy Sĩ', 'Sweden': 'Thụy Điển', 'Netherlands': 'Hà Lan',
        'Argentina': 'Argentina', 'Chile': 'Chile', 'Colombia': 'Colombia',
        'Peru': 'Peru', 'Turkey': 'Thổ Nhĩ Kỳ', 'Russia': 'Nga',
        'Israel': 'Israel', 'Egypt': 'Ai Cập', 'Japan': 'Nhật Bản',
    }
    for eng, viet in countries.items():
        if eng.lower() in desc_lower:
            info['origin'] = viet
            break

    # Tạo mô tả đơn giản
    simple_parts = [info['virus_name_vi']]
    if info['sample_id']:
        simple_parts.append(f"Mã mẫu: {info['sample_id']}")
    if info['origin']:
        simple_parts.append(f"Nguồn gốc: {info['origin']}")
    if info['year']:
        simple_parts.append(f"Năm: {info['year']}")
    info['simple_desc'] = ' | '.join(simple_parts)

    return info


# Tham chiếu Spike protein SARS-CoV-2 (đoạn đầu, ~200 aa)
# Dùng để so sánh và phát hiện đột biến amino acid
SPIKE_REFERENCE_START = 21563  # Vị trí bắt đầu Spike trên genome
SPIKE_SITES_OF_INTEREST = {
    # Wuhan-Hu-1 reference positions (trên genome)
    23403: {'ref': 'A', 'alt': 'G', 'gene': 'Spike', 'codon': 'D614G', 'note': 'Spike D614G - gia tăng lây nhiễm'},
    23604: {'ref': 'C', 'alt': 'G', 'gene': 'Spike', 'codon': 'P681R', 'note': 'Spike P681R (Alpha/Beta/Gamma)'},
    24469: {'ref': 'T', 'alt': 'C', 'gene': 'Spike', 'codon': 'N501Y', 'note': 'Spike N501Y - tăng ái lực ACE2'},
    21563: {'ref': 'A', 'alt': 'T', 'gene': 'Spike', 'codon': 'L5F', 'note': 'Spike signal peptide'},
    22813: {'ref': 'G', 'alt': 'T', 'gene': 'Spike', 'codon': 'K417N', 'note': 'Spike K417N (Beta/Omicron)'},
    22812: {'ref': 'A', 'alt': 'C', 'gene': 'Spike', 'codon': 'K417T', 'note': 'Spike K417T (Gamma)'},
    23012: {'ref': 'G', 'alt': 'A', 'gene': 'Spike', 'codon': 'E484K', 'note': 'Spike E484K - escape kháng thể'},
    23063: {'ref': 'A', 'alt': 'G', 'gene': 'Spike', 'codon': 'D484A', 'note': 'Spike E484A (Omicron)'},
    23525: {'ref': 'C', 'alt': 'T', 'gene': 'Spike', 'codon': 'T478K', 'note': 'Spike T478K (Delta/Omicron)'},
    23599: {'ref': 'G', 'alt': 'T', 'gene': 'Spike', 'codon': 'D614G', 'note': 'Spike D614G'},
}


# ============================================================
# 3. VIRUS ANALYZER
# ============================================================
class VirusAnalyzer:
    """Phân tích & dự đoán biến thể virus. Load tất cả từ weight/ folder."""

    def __init__(self, weight_dir=None):
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        if weight_dir is None:
            weight_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'weight')

        self.weight_dir = weight_dir
        self.model = None
        self.tfidf = None
        self.scaler = None
        self.label_encoder = None
        self.class_stats = None
        self.metadata = None

        self._load_all()

    # ----------------------------------------------------------
    # Load
    # ----------------------------------------------------------
    def _load_all(self):
        if not os.path.isdir(self.weight_dir):
            raise FileNotFoundError(
                f"Không tìm thấy folder weight/: {self.weight_dir}\n"
                "Hãy chạy train_pipeline.py trước để tạo weight/"
            )

        with open(os.path.join(self.weight_dir, 'metadata.pkl'), 'rb') as f:
            self.metadata = pickle.load(f)
        with open(os.path.join(self.weight_dir, 'tfidf.pkl'), 'rb') as f:
            self.tfidf = pickle.load(f)
        with open(os.path.join(self.weight_dir, 'scaler.pkl'), 'rb') as f:
            self.scaler = pickle.load(f)
        with open(os.path.join(self.weight_dir, 'label_encoder.pkl'), 'rb') as f:
            self.label_encoder = pickle.load(f)
        with open(os.path.join(self.weight_dir, 'class_stats.pkl'), 'rb') as f:
            self.class_stats = pickle.load(f)

        model_path = os.path.join(self.weight_dir, 'open_world_encoder.pth')
        self.model = OpenWorldVirusEncoder(
            input_dim=self.metadata['input_dim'],
            num_classes=self.metadata['num_classes']
        ).to(self.device)
        self.model.load_state_dict(torch.load(model_path, map_location=self.device))
        self.model.eval()

        print(f"✅ VirusAnalyzer loaded từ: {self.weight_dir}")

    # ----------------------------------------------------------
    # Tiền xử lý FASTA - trả về TẤT CẢ trình tự
    # ----------------------------------------------------------
    def preprocess_fasta(self, file_path):
        """Đọc tất cả trình tự từ file FASTA, trả về features + thông tin từng trình tự."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Không tìm thấy file: {file_path}")

        records = list(SeqIO.parse(file_path, "fasta"))
        if not records:
            raise ValueError("File FASTA không chứa trình tự nào")

        # Thu thập thông tin từng trình tự
        sequence_infos = []
        all_features = []

        for i, record in enumerate(records):
            seq_str = str(record.seq).upper()
            seq_len = len(seq_str)

            # Preview: 80 nucleotide đầu + ... + 80 nucleotide cuối
            if seq_len > 160:
                preview = seq_str[:80] + "..." + seq_str[-80:]
            else:
                preview = seq_str

            # Parse mô tả thành thông tin dễ hiểu
            parsed = parse_fasta_description(record.description if record.description else record.id)

            sequence_infos.append({
                'index': i + 1,
                'id': record.id,
                'description': record.description[:150] if record.description else record.id,
                'simple_desc': parsed['simple_desc'],
                'sample_type': parsed['sample_type'],
                'origin': parsed['origin'],
                'year': parsed['year'],
                'virus_name': parsed['virus_name_vi'],
                'length': seq_len,
                'preview': preview,
                'gc_content': round((seq_str.count('G') + seq_str.count('C')) / max(seq_len, 1) * 100, 2),
                'n_counts': {
                    'A': seq_str.count('A'),
                    'C': seq_str.count('C'),
                    'G': seq_str.count('G'),
                    'T': seq_str.count('T'),
                    'N': seq_str.count('N'),
                }
            })

            if seq_len >= 25000:
                all_features.append(count_multi_kmers(record.seq))

        if not all_features:
            # Nếu không trình tự nào đủ dài, vẫn lấy trình tự dài nhất
            longest = max(records, key=lambda r: len(r.seq))
            all_features.append(count_multi_kmers(longest.seq))

        # Trung bình features của tất cả trình tự đủ dài
        features_array = np.array(all_features)
        mean_features = np.mean(features_array, axis=0, keepdims=True)

        features_tfidf = self.tfidf.transform(mean_features).toarray()
        features_scaled = self.scaler.transform(features_tfidf)

        return features_scaled, sequence_infos

    def preprocess_csv(self, file_path):
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Không tìm thấy file: {file_path}")

        df = pd.read_csv(file_path)
        seq_column = None
        for col in df.columns:
            if any(kw in col.lower() for kw in ['sequence', 'seq', 'genome']):
                seq_column = col
                break
        if seq_column is None:
            seq_column = df.columns[0]

        results = []
        sequence_infos = []
        for idx, row in df.iterrows():
            seq = str(row[seq_column]).upper()
            seq_len = len(seq)
            preview = seq[:80] + "..." + seq[-80:] if seq_len > 160 else seq
            sequence_infos.append({
                'index': idx + 1,
                'id': f"row_{idx}",
                'description': f"CSV row {idx}",
                'length': seq_len,
                'preview': preview,
                'gc_content': round((seq.count('G') + seq.count('C')) / max(seq_len, 1) * 100, 2),
            })
            if seq_len >= 25000:
                results.append(count_multi_kmers(seq))

        if not results:
            raise ValueError("Không tìm thấy trình tự nào đủ dài trong file CSV")

        features_array = np.array(results)
        mean_features = np.mean(features_array, axis=0, keepdims=True)
        features_tfidf = self.tfidf.transform(mean_features).toarray()
        features_scaled = self.scaler.transform(features_tfidf)
        return features_scaled, sequence_infos

    # ----------------------------------------------------------
    # MC Dropout helper
    # ----------------------------------------------------------
    @staticmethod
    def _enable_dropout(m):
        if type(m) == nn.Dropout:
            m.train()

    # ----------------------------------------------------------
    # Phát hiện đột biến vị trí cụ thể trên genome
    # ----------------------------------------------------------
    def detect_mutations_detailed(self, file_path):
        """So sánh trình tự với reference để tìm đột biến tại các vị trí quan trọng."""
        try:
            records = list(SeqIO.parse(file_path, "fasta"))
            if not records:
                return [], {}

            # Lấy trình tự dài nhất để so sánh
            record = max(records, key=lambda r: len(r.seq))
            seq = str(record.seq).upper()

            # Thông tin trình tự được dùng để detect mutations
            source_seq = {
                'id': record.id,
                'description': record.description[:150] if record.description else record.id,
                'length': len(seq),
                'index': next((i+1 for i, r in enumerate(records) if r.id == record.id), 0)
            }

            found_mutations = []
            for pos, info in SPIKE_SITES_OF_INTEREST.items():
                if pos < len(seq):
                    actual = seq[pos]
                    if actual != info['ref']:
                        found_mutations.append({
                            'position': pos,
                            'ref_base': info['ref'],
                            'actual_base': actual,
                            'gene': info['gene'],
                            'codon_change': info['codon'],
                            'note': info['note'],
                            'type': 'SNP',
                            'confidence': 95.0
                        })

            # Tìm thêm các vùng bất thường tần suất nucleotide
            nucleotide_freq = Counter(seq)
            total = len(seq)
            window_size = 1000
            for i in range(0, len(seq) - window_size, window_size):
                window = seq[i:i+window_size]
                window_freq = Counter(window)
                for nuc in ['A', 'C', 'G', 'T']:
                    expected_freq = nucleotide_freq[nuc] / total
                    actual_freq = window_freq.get(nuc, 0) / window_size
                    if abs(actual_freq - expected_freq) > 0.1:
                        found_mutations.append({
                            'position': i + window_size // 2,
                            'ref_base': f'{nuc}:{expected_freq:.2%}',
                            'actual_base': f'{nuc}:{actual_freq:.2%}',
                            'gene': 'Genome',
                            'codon_change': f'Frequency shift',
                            'note': f'Tần suất {nuc} bất thường tại vùng {i}-{i+window_size}',
                            'type': 'Frequency',
                            'confidence': min(95, abs(actual_freq - expected_freq) * 500)
                        })

            return found_mutations[:20], source_seq
        except Exception as e:
            print(f"Lỗi detect_mutations_detailed: {e}")
            return [], {}

    # ----------------------------------------------------------
    # Dự đoán
    # ----------------------------------------------------------
    def predict_variant(self, file_path, T=3.0, mc_passes=30, mahalanobis_threshold=150.0):
        file_ext = os.path.splitext(file_path)[1].lower()
        if file_ext in ['.fasta', '.fa']:
            features_scaled, sequence_infos = self.preprocess_fasta(file_path)
        elif file_ext == '.csv':
            features_scaled, sequence_infos = self.preprocess_csv(file_path)
        else:
            raise ValueError(f"Định dạng file không hỗ trợ: {file_ext}")

        sample_tensor = torch.FloatTensor(features_scaled).to(self.device)

        # MC Dropout
        self.model.eval()
        self.model.apply(self._enable_dropout)

        mc_logits, mc_embeddings = [], []
        with torch.no_grad():
            for _ in range(mc_passes):
                logits, emb = self.model(sample_tensor)
                mc_logits.append(logits.cpu().numpy()[0])
                mc_embeddings.append(emb.cpu().numpy()[0])

        # Temperature Scaling
        mean_logits = np.mean(mc_logits, axis=0)
        scaled_logits = mean_logits / T
        exp_logits = np.exp(scaled_logits - np.max(scaled_logits))
        probs = exp_logits / exp_logits.sum()

        pred_class = int(np.argmax(probs))
        confidence = float(probs[pred_class])
        prediction_variance = float(np.var(mc_logits, axis=0).mean())

        # Mahalanobis
        mean_embedding = np.mean(mc_embeddings, axis=0)
        variant_names = list(self.label_encoder.classes_)

        min_dist = float('inf')
        best_class = pred_class
        for label, stats in self.class_stats.items():
            dist = float(mahalanobis(mean_embedding, stats['centroid'], stats['cov_inv']))
            if dist < min_dist:
                min_dist = dist
                best_class = int(stats['index'])

        is_unknown = bool(min_dist > mahalanobis_threshold or prediction_variance > 1.5)

        # Similarity scores
        similarity_scores = {}
        for label, stats in self.class_stats.items():
            dist = float(mahalanobis(mean_embedding, stats['centroid'], stats['cov_inv']))
            similarity = max(0.0, 100.0 - (dist / mahalanobis_threshold * 100.0))
            similarity_scores[label] = round(similarity, 2)

        # Mutation sites (chi tiết hơn)
        if file_ext in ['.fasta', '.fa']:
            mutation_sites, mutation_source = self.detect_mutations_detailed(file_path)
        else:
            mutation_sites, mutation_source = [], {}

        result = {
            'file_name': os.path.basename(file_path),
            'variant': variant_names[best_class] if not is_unknown else 'UNKNOWN_VARIANT',
            'confidence': round(confidence * 100, 2),
            'mahalanobis_distance': round(min_dist, 2),
            'mahalanobis_threshold': mahalanobis_threshold,
            'mc_variance': round(prediction_variance, 4),
            'is_unknown': is_unknown,
            'mutation_sites': mutation_sites,
            'mutation_source': mutation_source,
            'similarity_scores': similarity_scores,
            'probabilities': {variant_names[i]: round(float(probs[i]) * 100, 2) for i in range(len(probs))},
            'sequences': sequence_infos,
            'total_sequences': len(sequence_infos),
            'sequences_analyzed': sum(1 for s in sequence_infos if s['length'] >= 25000),
        }
        return result

    # ----------------------------------------------------------
    # Tóm tắt
    # ----------------------------------------------------------
    @staticmethod
    def get_analysis_summary(result):
        lines = [
            f"📁 File: {result['file_name']}",
            f"🧬 Biến thể: {result['variant']}",
            f"📊 Độ tự tin: {result['confidence']}%",
            f"📏 Mahalanobis: {result['mahalanobis_distance']} (Ngưỡng: {result['mahalanobis_threshold']})",
            f"📈 MC Variance: {result['mc_variance']}",
            f"🔢 Số trình tự: {result['total_sequences']} (phân tích: {result['sequences_analyzed']})",
        ]
        if result['is_unknown']:
            lines.append("🚨 UNKNOWN_VARIANT - Cần giám sát thêm")
        else:
            lines.append(f"✅ Biến thể đã biết: {result['variant']}")

        lines.append("\n📊 Tương đồng:")
        for variant, score in result['similarity_scores'].items():
            lines.append(f"  - {variant}: {score}%")

        if result['mutation_sites']:
            lines.append(f"\n🔬 Đột biến ({len(result['mutation_sites'])}):")
            for mut in result['mutation_sites'][:5]:
                lines.append(f"  - {mut['position']}: {mut.get('codon_change', mut['type'])}")

        return "\n".join(lines)


# ============================================================
# 4. CLI
# ============================================================
def main():
    import argparse
    parser = argparse.ArgumentParser(description='ViralGen AI - Dự đoán biến thể virus')
    parser.add_argument('--input', '-i', required=True, help='File FASTA/CSV cần phân tích')
    parser.add_argument('--weight', '-w', default=None, help='Thư mục weight')
    args = parser.parse_args()

    analyzer = VirusAnalyzer(weight_dir=args.weight)
    result = analyzer.predict_variant(args.input)
    summary = analyzer.get_analysis_summary(result)

    print("\n" + "=" * 60)
    print("🌍 KẾT QUẢ PHÂN TÍCH VIRALGEN AI")
    print("=" * 60)
    print(summary)
    print("=" * 60)

    import json
    print("\n📋 JSON:")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()