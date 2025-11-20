# Radial Shell Encoding Project

## Overview
Implementation of a multi-resolution radial shell encoding system for spatial transcriptomics similarity search. Enables region-to-region comparison across ~600 Xenium samples with 5000 genes at multiple spatial resolutions (10-160 μm).

## Key Components

### 1. Data Structure
- **Input**: Binned Xenium data at 5 resolutions (10, 20, 40, 80, 160 μm)
- **Storage**: Zarr-compressed spatial expression matrices
- **Genes**: 300-5000 per sample, ~70% spatially variable
- **Query Size**: Up to 500 bins

### 2. Encoding Pipeline
```
Raw Bins → Patch Extraction → Radial Shell Encoding → PCA Compression → FAISS Index
```

### 3. Core Algorithms

#### Radial Shell Encoding
- Divides circular patches into 5 concentric shells
- Computes average gene expression per shell
- Captures center-periphery spatial organization
- Rotation/translation invariant

#### PCA Compression
- Reduces 17,500 dimensions (3500 genes × 5 shells) to 256
- Retains 95-99% variance
- Enables 64× storage reduction and faster search

#### Multi-Resolution Architecture
- Resolution-specific patch radii (e.g., 10μm: 10-100 bins, 160μm: 3-8 bins)
- Physical patch sizes comparable across resolutions (200-1000 μm)
- Auto-detection of appropriate radius based on query size

## Database Structure
```
database/
├── 10um/
│   ├── patches_metadata.parquet      # Sample ID, center coords, radius, n_bins
│   ├── embeddings/
│   │   └── embeddings_r*.h5          # PCA-compressed (256D)
│   ├── faiss_indices/
│   │   └── faiss_r*.index            # Inner product search
│   └── pca_model_10um.pkl            # Fitted transformer
├── 20um/
├── 40um/
├── 80um/
└── 160um/
```

## Performance Targets
- **Storage**: ~17 GB total (from 1.1 TB naive)
- **Query Time**: ~2 seconds including encoding
- **Memory**: ~3 GB per resolution
- **Scalability**: Works for 10,000+ samples

## Integration Points
- Uses existing `XeniumProcessor` for data loading
- Compatible with current zarr/bin structure
- Leverages haystack results for spatially variable genes
- Outputs compatible with existing visualization pipeline

## Implementation Status
- [x] Core encoding module (radial_shell_encoder.py)
- [x] Database builder (database_builder.py)
- [x] Search engine (similarity_search.py)
- [x] CLI interface (radial_shell_system.py)
- [x] Integration with XeniumProcessor (radial_shell_integration_example.py)
- [x] Validation scripts (test_radial_shell_system.py)
- [x] Documentation (RADIAL_SHELL_README.md, RADIAL_SHELL_QUICKSTART.md)

## Files Created

### Core Implementation
1. `radial_shell_encoder.py` - Radial shell encoding logic and patch generation
2. `database_builder.py` - Multi-resolution database building pipeline
3. `similarity_search.py` - Search engine with FAISS indexing
4. `radial_shell_system.py` - Main CLI interface

### Integration & Testing
5. `test_radial_shell_system.py` - Comprehensive validation tests
6. `radial_shell_integration_example.py` - Flask API integration examples

### Documentation
7. `RADIAL_SHELL_README.md` - Complete system documentation
8. `RADIAL_SHELL_QUICKSTART.md` - Quick start guide
9. `requirements_radial_shell.txt` - Python dependencies
10. `.claude/commands/radial-search.md` - Claude command integration

## Quick Start

```bash
# Install dependencies
pip install -r requirements_radial_shell.txt

# Run tests
python test_radial_shell_system.py

# Build database
python radial_shell_system.py build --data-dir /path/to/data --output-dir spatial_db

# Search
python radial_shell_system.py search --database spatial_db --sample /path/to/sample \
  --x-center 5000 --y-center 5000 --radius 500 --resolution 40
```
