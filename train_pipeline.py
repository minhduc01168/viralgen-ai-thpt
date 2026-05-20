#!/usr/bin/env python3
"""
ViralGen AI - Training Pipeline
Chạy 1 lần để train model + tfidf + scaler + class_stats, lưu vào folder weight/
"""

import os
import sys
import itertools
import pickle
from collections import Counter

import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.preprocessing import StandardScaler
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
import scipy.linalg

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 1. HẰNG SỐ
# ============================================================
NUCLEOTIDES = ['A', 'C', 'G', 'T']
kmers_3 = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=3)]
kmers_4 = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=4)]
kmers_5 = [''.join(p) for p in itertools.product(NUCLEOTIDES, repeat=5)]
ALL_MULTI_KMERS = kmers_3 + kmers_4 + kmers_5  # 1344 chiều

VARIANT_NAMES = ['Alpha', 'Beta', 'Delta', 'Gamma', 'Omicron']
NUM_CLASSES = len(VARIANT_NAMES)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
WEIGHT_DIR = os.path.join(BASE_DIR, 'weight')
TRAIN_DIR = os.path.join(BASE_DIR, 'dataset_test', 'IND')

INPUT_DIM = len(ALL_MULTI_KMERS)  # 1344
EMBEDDING_DIM = 128
EPOCHS = 20
BATCH_SIZE = 32
LEARNING_RATE = 0.001
LABEL_SMOOTHING = 0.1


# ============================================================
# 2. MÔ HÌNH
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
# 3. HÀM TIỆN ÍCH
# ============================================================
def count_multi_kmers(sequence):
    seq_str = str(sequence).upper()
    kmer_counts = Counter()
    for k in [3, 4, 5]:
        kmers_in_seq = [seq_str[i:i+k] for i in range(len(seq_str) - k + 1)]
        kmer_counts.update(kmers_in_seq)
    return [kmer_counts.get(kmer, 0) for kmer in ALL_MULTI_KMERS]


# ============================================================
# 4. PIPELINE
# ============================================================
def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"🖥️  Device: {device}")
    print(f"📂 Train dir: {TRAIN_DIR}")
    print(f"💾 Weight dir: {WEIGHT_DIR}")

    os.makedirs(WEIGHT_DIR, exist_ok=True)

    # ----- Bước 1: Đọc dữ liệu & trích xuất Multi-Kmer -----
    print("\n🔬 Bước 1: Đọc dữ liệu & trích xuất Multi-Kmer (1344 chiều)...")
    data_rows, labels = [], []
    fasta_files = [
        os.path.join(dp, f)
        for dp, dn, filenames in os.walk(TRAIN_DIR)
        for f in filenames if f.lower().endswith(('.fasta', '.fa'))
    ]
    print(f"   Tìm thấy {len(fasta_files)} file FASTA")

    for file_path in fasta_files:
        label = os.path.splitext(os.path.basename(file_path))[0].split('_')[-1].capitalize()
        for record in SeqIO.parse(file_path, "fasta"):
            if len(record.seq) >= 25000:
                data_rows.append(count_multi_kmers(record.seq))
                labels.append(label)

    df_raw = pd.DataFrame(data_rows, columns=ALL_MULTI_KMERS)
    y_raw = np.array(labels)
    print(f"   ✅ {df_raw.shape[0]} mẫu, {df_raw.shape[1]} đặc trưng")

    # ----- Bước 2: TF-IDF + StandardScaler -----
    print("\n📊 Bước 2: Huấn luyện TF-IDF & StandardScaler...")
    tfidf = TfidfTransformer()
    X_tfidf = tfidf.fit_transform(df_raw.values).toarray()

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_tfidf)
    print(f"   ✅ X_scaled shape: {X_scaled.shape}")

    # ----- Bước 3: Label Encoding & Train/Test Split -----
    print("\n🏷️  Bước 3: Label Encoding & Train/Test Split...")
    le = LabelEncoder()
    y_encoded = le.fit_transform(y_raw)

    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y_encoded, test_size=0.2, random_state=42, stratify=y_encoded
    )
    print(f"   Train: {X_train.shape[0]} mẫu, Test: {X_test.shape[0]} mẫu")
    print(f"   Classes: {le.classes_}")

    train_dataset = TensorDataset(torch.FloatTensor(X_train), torch.LongTensor(y_train))
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    test_dataset = TensorDataset(torch.FloatTensor(X_test), torch.LongTensor(y_test))
    test_loader = DataLoader(test_dataset, batch_size=BATCH_SIZE, shuffle=False)

    # ----- Bước 4: Huấn luyện mô hình -----
    print(f"\n🚀 Bước 4: Huấn luyện OpenWorldVirusEncoder ({EPOCHS} epochs)...")
    model = OpenWorldVirusEncoder(input_dim=INPUT_DIM, num_classes=NUM_CLASSES).to(device)
    criterion = nn.CrossEntropyLoss(label_smoothing=LABEL_SMOOTHING)
    optimizer = optim.Adam(model.parameters(), lr=LEARNING_RATE)

    best_acc = 0.0
    for epoch in range(EPOCHS):
        model.train()
        running_loss = 0.0
        for inputs, labels_batch in train_loader:
            inputs, labels_batch = inputs.to(device), labels_batch.to(device)
            optimizer.zero_grad()
            logits, _ = model(inputs)
            loss = criterion(logits, labels_batch)
            loss.backward()
            optimizer.step()
            running_loss += loss.item()

        model.eval()
        correct, total = 0, 0
        with torch.no_grad():
            for inputs, labels_batch in test_loader:
                inputs, labels_batch = inputs.to(device), labels_batch.to(device)
                logits, _ = model(inputs)
                _, predicted = torch.max(logits.data, 1)
                total += labels_batch.size(0)
                correct += (predicted == labels_batch).sum().item()

        epoch_acc = 100 * correct / total
        avg_loss = running_loss / len(train_loader)
        print(f"   Epoch [{epoch+1}/{EPOCHS}] - Loss: {avg_loss:.4f} - Test Acc: {epoch_acc:.2f}%")

        if epoch_acc > best_acc:
            best_acc = epoch_acc

    print(f"   ✅ Best Test Acc: {best_acc:.2f}%")

    # ----- Bước 5: Tính class_stats (Mahalanobis) -----
    print("\n📐 Bước 5: Tính Centroid & Covariance cho Mahalanobis...")
    model.eval()
    train_embeddings, train_labels_list = [], []

    with torch.no_grad():
        for inputs, labels_batch in train_loader:
            inputs = inputs.to(device)
            _, emb = model(inputs)
            train_embeddings.append(emb.cpu().numpy())
            train_labels_list.append(labels_batch.numpy())

    X_emb = np.vstack(train_embeddings)
    y_emb = np.concatenate(train_labels_list)

    class_stats = {}
    for c in range(NUM_CLASSES):
        X_c = X_emb[y_emb == c]
        centroid = np.mean(X_c, axis=0)
        cov_matrix = np.cov(X_c, rowvar=False) + np.eye(EMBEDDING_DIM) * 1e-6
        cov_inv = scipy.linalg.pinvh(cov_matrix)
        class_stats[le.classes_[c]] = {
            'centroid': centroid,
            'cov_inv': cov_inv,
            'index': c
        }

    print(f"   ✅ Đã tính class_stats cho {NUM_CLASSES} biến thể")

    # ----- Bước 6: Lưu tất cả vào weight/ -----
    print(f"\n💾 Bước 6: Lưu artifacts vào {WEIGHT_DIR}/...")

    # 6a. Model weights
    model_path = os.path.join(WEIGHT_DIR, 'open_world_encoder.pth')
    torch.save(model.state_dict(), model_path)
    print(f"   ✅ Model → {model_path}")

    # 6b. TF-IDF
    tfidf_path = os.path.join(WEIGHT_DIR, 'tfidf.pkl')
    with open(tfidf_path, 'wb') as f:
        pickle.dump(tfidf, f)
    print(f"   ✅ TF-IDF → {tfidf_path}")

    # 6c. Scaler
    scaler_path = os.path.join(WEIGHT_DIR, 'scaler.pkl')
    with open(scaler_path, 'wb') as f:
        pickle.dump(scaler, f)
    print(f"   ✅ Scaler → {scaler_path}")

    # 6d. Label Encoder
    le_path = os.path.join(WEIGHT_DIR, 'label_encoder.pkl')
    with open(le_path, 'wb') as f:
        pickle.dump(le, f)
    print(f"   ✅ Label Encoder → {le_path}")

    # 6e. Class stats (Mahalanobis)
    stats_path = os.path.join(WEIGHT_DIR, 'class_stats.pkl')
    with open(stats_path, 'wb') as f:
        pickle.dump(class_stats, f)
    print(f"   ✅ Class Stats → {stats_path}")

    # 6f. Config metadata
    metadata = {
        'input_dim': INPUT_DIM,
        'num_classes': NUM_CLASSES,
        'embedding_dim': EMBEDDING_DIM,
        'variant_names': list(le.classes_),
        'kmer_sizes': [3, 4, 5],
        'label_smoothing': LABEL_SMOOTHING,
        'best_test_acc': best_acc,
    }
    meta_path = os.path.join(WEIGHT_DIR, 'metadata.pkl')
    with open(meta_path, 'wb') as f:
        pickle.dump(metadata, f)
    print(f"   ✅ Metadata → {meta_path}")

    print("\n" + "=" * 60)
    print("🎉 HOÀN TẤT PIPELINE!")
    print(f"   Tất cả artifacts đã được lưu trong: {WEIGHT_DIR}/")
    print(f"   Best Test Accuracy: {best_acc:.2f}%")
    print("=" * 60)


if __name__ == "__main__":
    main()