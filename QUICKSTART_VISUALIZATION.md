# Quick Start: Visualize Your Database

## Step 1: Make sure your database is built

```bash
python radial_shell_system.py build \
  --data-dir /home/vlad/xenium_mundus/data_full/data/GPL33762/GSE243168 \
  --output-dir spatial_database \
  --resolutions 160 \
  --use-variable-genes False
```

Expected output:
```
[Phase 0/3] Discovering global gene coordinate system...
Total unique genes across all samples: 321
Using all 321 genes

[Phase 1/3] Generating patches and raw embeddings...
Processing samples: 100%|████████████████| 3/3

[Phase 2/3] Fitting PCA...
Original dimensions: 1605
Compressed dimensions: 256

[Phase 3/3] Building FAISS indices...
✓ DATABASE BUILD COMPLETE
```

## Step 2: Visualize the database

```bash
python visualize_database.py \
  --database spatial_database \
  --resolution 160 \
  --gene-name TP53 \
  --output-dir database_plots
```

This will create 4 visualizations in `database_plots/`:

### 1. `radial_shell_concept.png`
Shows how the radial shell encoding works:
- Left: Concentric shells dividing a circular patch
- Right: Feature matrix (genes × shells)

### 2. `database_summary.png`
Database statistics:
- Number of patches per radius
- Number of patches per sample
- Patch size distribution
- Configuration summary

### 3. `feature_statistics.png`
Feature analysis across different radii:
- Feature value distributions (violin plots)
- Average gene activity per gene

### 4. `gene_encoding_TP53.png`
Detailed view of how TP53 is encoded:
- Top row: Heatmap of TP53 expression across shells and patches
- Bottom row: Average expression per shell with error bars

## Step 3: Try different genes

```bash
# Visualize a different gene
python visualize_database.py \
  --database spatial_database \
  --resolution 160 \
  --gene-name ACTB \
  --output-dir database_plots
```

## Step 4: Check what's stored

The database contains:

```
spatial_database/160um/
├── config.json                  # Configuration (genes, parameters)
├── patches_metadata.parquet     # Patch locations and info
├── pca_models.pkl              # PCA transformation models
├── faiss_index_r3.bin          # FAISS index for radius 3
├── faiss_index_r5.bin          # FAISS index for radius 5
├── faiss_index_r8.bin          # FAISS index for radius 8
├── raw_embeddings_r3.h5        # Raw features (before PCA)
├── raw_embeddings_r5.h5
└── raw_embeddings_r8.h5
```

### Check the configuration:

```bash
cat spatial_database/160um/config.json
```

```json
{
  "global_genes": ["ACTB", "TP53", ...],
  "use_variable_genes": false,
  "n_shells": 5,
  "pca_dims": 256,
  "resolution_um": 160
}
```

### Check the metadata:

```python
import pandas as pd

metadata = pd.read_parquet('spatial_database/160um/patches_metadata.parquet')
print(metadata.head())
```

```
   patch_id     sample_id  resolution_um  radius  center_x  center_y  n_bins
0  GSM778...  GSM7780155            160       3     160.0     160.0      29
1  GSM778...  GSM7780155            160       3     160.0     320.0      29
2  GSM778...  GSM7780155            160       3     160.0     480.0      29
...
```

### Check raw embeddings:

```python
import h5py

with h5py.File('spatial_database/160um/raw_embeddings_r5.h5', 'r') as f:
    embeddings = f['embeddings'][:]
    patch_ids = f['patch_ids'][:]

print(f"Shape: {embeddings.shape}")
print(f"Patches: {len(patch_ids)}")
print(f"Features per patch: {embeddings.shape[1]}")
print(f"Expected: {321 genes} × {5 shells} = {321*5} features")
```

```
Shape: (710, 1605)
Patches: 710
Features per patch: 1605
Expected: 321 genes × 5 shells = 1605 features
```

## Step 5: Understand the heatmap

When you look at `gene_encoding_TP53.png`:

**Top row (heatmap):**
- **X-axis:** Patch index (0-100, showing first 100 patches)
- **Y-axis:** Shell index (1-5, from center to edge)
- **Color:** TP53 expression level in that shell for that patch

**Patterns to look for:**
- Horizontal bands → Gene shows consistent radial pattern across patches
- Vertical bands → Some patches have high expression, others low
- Gradients → Gene expression changes from center to edge

**Bottom row (bar chart):**
- **X-axis:** Shell index (1-5)
- **Y-axis:** Average TP53 expression across all patches
- **Error bars:** Standard deviation

**Interpretations:**
- Increasing from center → Expression higher at edges
- Decreasing from center → Expression higher in center
- Flat → No radial pattern

## Step 6: Run a search

```bash
python radial_shell_system.py search \
  --database spatial_database \
  --sample /home/vlad/xenium_mundus/data_full/data/GPL33762/GSE243168/GSM7780155 \
  --resolution 160 \
  --x-center 100 --y-center 100 \
  --radius 800 \
  --top-k 20 \
  --output search_results.csv
```

View results:
```bash
head -5 search_results.csv
```

```
patch_id,sample_id,similarity,n_bins,center_x,center_y,radius,resolution_um
GSM7780155_r5_p0,GSM7780155,1.0000,81,100.0,100.0,5,160
GSM7780154_r5_p234,GSM7780154,0.9523,81,112.0,98.0,5,160
GSM7780155_r5_p156,GSM7780155,0.9401,79,80.0,120.0,5,160
```

## Troubleshooting

### Database not found
```
Error: Resolution directory not found: spatial_database/160um
```
**Solution:** Build the database first (Step 1)

### Gene not found
```
Error: 'TP53' is not in list
```
**Solution:** Check available genes in config.json:
```bash
python -c "import json; print(json.load(open('spatial_database/160um/config.json'))['global_genes'][:10])"
```

### No raw embeddings
```
Warning: Raw embeddings not found
```
**Solution:** Raw embeddings are saved automatically during build. If missing, rebuild the database.

## What the plots show you

### Radial Shell Concept
- **Purpose:** Understand the encoding method
- **Use:** Explain to others how it works

### Database Summary
- **Purpose:** Check database quality and coverage
- **Use:** Verify patches are distributed across samples and radii

### Feature Statistics
- **Purpose:** Understand feature distributions
- **Use:** Diagnose issues (e.g., all zeros, outliers)

### Gene Encoding
- **Purpose:** See what the system "sees" for a gene
- **Use:** Validate that spatial patterns are captured

## Advanced: Custom visualization

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load embeddings for radius 5
with h5py.File('spatial_database/160um/raw_embeddings_r5.h5', 'r') as f:
    embeddings = f['embeddings'][:]

# Extract features for gene 0 (ACTB)
gene_idx = 0
n_shells = 5
patch_idx = 0

gene_features = embeddings[patch_idx, gene_idx*n_shells:(gene_idx+1)*n_shells]

# Plot
plt.figure(figsize=(8, 5))
plt.bar(range(n_shells), gene_features)
plt.xlabel('Shell Index')
plt.ylabel('Expression')
plt.title(f'ACTB Expression in Patch {patch_idx}')
plt.xticks(range(n_shells), [f'Shell {i+1}' for i in range(n_shells)])
plt.savefig('actb_patch0.png')
```

## Next: Integrate with Frontend

See `SIMILARITY_SEARCH_INTEGRATION.md` for adding a search button to your web interface.
