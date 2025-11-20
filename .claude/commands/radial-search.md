---
description: Build radial shell encoding database and perform similarity searches
---

# Radial Shell Encoding - Spatial Similarity Search

This command helps build and query the radial shell encoding database for spatial transcriptomics similarity search.

## Usage Examples

### Build Database
```bash
python radial_shell_system.py build --data-dir /path/to/xenium/data --output-dir spatial_database
```

### Search for Similar Regions
```bash
python radial_shell_system.py search \
  --database spatial_database \
  --sample MySample01 \
  --x-center 5000 --y-center 5000 --radius 50 \
  --resolution 40 \
  --top-k 100
```

### Build Database for Specific Resolutions
```bash
python radial_shell_system.py build \
  --data-dir /path/to/xenium/data \
  --output-dir spatial_database \
  --resolutions 20 40 80
```

## Features
- Multi-resolution support (10-160 μm)
- PCA-compressed storage (~17 GB vs 1.1 TB naive)
- Fast FAISS-based similarity search (~2 seconds)
- Spatially variable gene filtering
- Rotation/translation invariant encoding

## Key Parameters
- `--n-shells`: Number of concentric shells (default: 5)
- `--pca-dims`: PCA compression dimensions (default: 256)
- `--use-variable-genes`: Use only spatially variable genes (default: True)
