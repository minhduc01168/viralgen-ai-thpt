#!/usr/bin/env python3
"""
ViralGen AI - Configuration
"""

import os

class Config:
    """Cấu hình cơ sở"""
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'viralgen-ai-secret-key-2026'
    UPLOAD_FOLDER = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uploads')
    MAX_CONTENT_LENGTH = 50 * 1024 * 1024  # 50MB
    ALLOWED_EXTENSIONS = {'fasta', 'fa', 'csv'}
    
    # Model configuration
    MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models', 'open_world_encoder.pth')
    TRAIN_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'dataset_test', 'IND')
    
    # Virus analysis parameters
    KMER_SIZES = [3, 4, 5]
    MIN_SEQUENCE_LENGTH = 25000
    MAHALANOBIS_THRESHOLD = 25.0
    TEMPERATURE_SCALING = 3.0
    MC_DROPOUT_PASSES = 30
    
    # Variant names
    VARIANT_NAMES = ['Alpha', 'Beta', 'Delta', 'Gamma', 'Omicron']


class DevelopmentConfig(Config):
    """Cấu hình development"""
    DEBUG = True


class ProductionConfig(Config):
    """Cấu hình production"""
    DEBUG = False


# Dictionary ánh xạ tên cấu hình sang class
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'default': DevelopmentConfig
}