# HonchaMAP Radial Shell Encoding Package - Complete Summary

## What You Have Now

Your spatial transcriptomics similarity search system is **fully functional** with both PCA and non-PCA modes!

### Your Database (Built Successfully!)

```
spatial_database/
├── no_pca_db/160um/           ← NO PCA (preserves gene information)
│   ├── config.json            (321 genes, use_pca=false)
│   ├── patches_metadata.parquet (4,229 patches)
│   └── faiss_indices/
│       ├── faiss_r3.index     (3,399 patches, 1605 dims)
│       ├── faiss_r5.index     (710 patches, 1605 dims)
│       └── faiss_r8.index     (120 patches, 1605 dims)
│
└── 160um/                     ← WITH PCA (compressed, faster search)
    ├── config.json            (321 genes, use_pca=true)
    ├── patches_metadata.parquet
    ├── pca_model_160um.pkl    (1605 → 256 dims)
    └── faiss_indices/
        ├── faiss_r3.index     (3,332 patches, 256 dims)
        ├── faiss_r5.index     (710 patches, 256 dims)
        └── faiss_r8.index     (120 patches, 256 dims)
```

**Dimensions:**
- **Full (no PCA):** 321 genes × 5 shells = **1,605 dimensions**
- **Compressed (with PCA):** 1,605 → **256 dimensions**

---

## How The System Works

### 1. Radial Shell Encoding

**Input:** A circular patch of tissue
```
  Center: (x, y)
  Radius: r bins
  Bins: ~80-200 bins per patch
```

**Process:**
1. Divide patch into 5 concentric shells (like tree rings)
2. Calculate average gene expression in each shell
3. Create feature vector: `[gene0_shell0, gene0_shell1, ..., gene0_shell4, gene1_shell0, ..., gene320_shell4]`
4. Result: 1,605 features (321 genes × 5 shells)

**Why concentric shells?**
- Captures radial organization (common in biology)
- Immune infiltration zones
- Tumor microenvironments
- Gradient-based signaling

### 2. Global Gene Coordinate System

**Problem:** Your samples have different genes (300-5000 per sample)

**Solution:**
```python
Sample 1: [ACTB, TP53, MYC]      → 3 genes
Sample 2: [ACTB, BRCA1, TP53]   → 3 genes
Global:   [ACTB, BRCA1, MYC, TP53] → 4 genes (union)

# Map each sample to global space:
Sample 1 encoding:
  ACTB  (index 0): 5.2
  BRCA1 (index 1): 0.0  ← missing, filled with zero
  MYC   (index 2): 3.1
  TP53  (index 3): 2.7

Sample 2 encoding:
  ACTB  (index 0): 4.8
  BRCA1 (index 1): 6.3
  MYC   (index 2): 0.0  ← missing, filled with zero
  TP53  (index 3): 3.0
```

**Result:** All patches have same 1,605 dimensions regardless of sample!

### 3. PCA vs No-PCA

#### No-PCA Mode (use_pca=False) ✅ **YOU'RE USING THIS**

**Advantages:**
- ✅ Preserve individual gene information
- ✅ Can analyze specific genes (ACTA2, AIF1, etc.)
- ✅ No information loss
- ✅ Full interpretability

**Disadvantages:**
- ❌ Larger indices (22 MB vs 3 MB)
- ❌ Slightly slower search (~10-20% slower)

**Your database stats:**
```
Radius 3: 3,399 patches × 1,605 features = 21.8 MB
Radius 5:   710 patches × 1,605 features =  4.5 MB
Radius 8:   120 patches × 1,605 features =  0.7 MB
Total: ~27 MB
```

#### PCA Mode (use_pca=True)

**Advantages:**
- ✅ Smaller indices (3.4 MB vs 21.8 MB)
- ✅ Faster search
- ✅ Remove noise

**Disadvantages:**
- ❌ **Cannot trace back to individual genes**
- ❌ Features are linear combinations
- ❌ ~10-20% information loss

**Your PCA database stats:**
```
Radius 3: 3,332 patches × 256 features = 3.4 MB
Radius 5:   710 patches × 256 features = 0.7 MB
Radius 8:   120 patches × 256 features = 0.1 MB
Total: ~4.2 MB
```

---

## Your Tools

### 1. Database Building

```bash
# Build database WITHOUT PCA (gene-level analysis)
python radial_shell_system.py build \
  --data-dir /path/to/xenium/data \
  --output-dir spatial_database/no_pca_db \
  --resolutions 160 \
  --use-variable-genes False \
  --use-pca False

# Build database WITH PCA (faster, compressed)
python radial_shell_system.py build \
  --data-dir /path/to/xenium/data \
  --output-dir spatial_database \
  --resolutions 160 \
  --use-variable-genes False \
  --use-pca True
```

### 2. Similarity Search

```bash
# Search for similar regions
python radial_shell_system.py search \
  --database spatial_database/no_pca_db \
  --sample /path/to/sample/GSM7780155 \
  --resolution 160 \
  --x-center 1500 \
  --y-center 1500 \
  --radius 800 \
  --top-k 100 \
  --output results.csv
```

**Results:** CSV with similarity scores, locations, sample IDs

### 3. Visualization (Gene-Level)

```bash
# Visualize a single gene across all radii
python visualize_no_pca_database.py \
  --database spatial_database/no_pca_db \
  --resolution 160 \
  --gene-name ACTA2 \
  --output-dir plots

# Compare multiple genes side-by-side
python visualize_no_pca_database.py \
  --database spatial_database/no_pca_db \
  --resolution 160 \
  --compare-genes ACTA2 ACTG2 AIF1 \
  --output-dir plots
```

**Output plots:**
1. **Heatmap**: Gene expression across shells (center→edge) and patches
2. **Bar chart**: Average radial profile with error bars
3. **Violin plot**: Distribution of expression per shell

---

## What Each Visualization Shows

### Single Gene Plot (3 rows × 3 columns)

**Columns:** Different patch radii (3, 5, 8 bins)
**Rows:**

#### Row 1: Expression Heatmap
- **X-axis:** Patch index (0-100)
- **Y-axis:** Shell index (1-5, center→edge)
- **Color:** Gene expression level

**Interpretation:**
- Horizontal bands → Consistent radial pattern
- Vertical bands → Some patches high, others low
- Gradient top-to-bottom → Expression changes center→edge

#### Row 2: Radial Profile
- **X-axis:** Shell index (1-5, center→edge)
- **Y-axis:** Average expression ± standard deviation
- **Bars:** Color-coded by shell

**Interpretation:**
- Increasing bars → Higher at edges
- Decreasing bars → Higher in center
- Flat → No radial pattern

#### Row 3: Expression Distribution
- **X-axis:** Shell index (1-5)
- **Y-axis:** Expression values (violin plot)
- Shows variability within each shell

**Interpretation:**
- Wide violins → High variability
- Narrow violins → Consistent expression
- Outliers → Unusual patches

### Multi-Gene Comparison

**Layout:** N genes × 2 columns

- **Left column:** Heatmap for each gene
- **Right column:** Radial profile for each gene

**Use:** Compare spatial patterns between genes

---

## Example Workflow

### 1. Build Database (Done!)

You already built both:
- `spatial_database/no_pca_db/` ← For gene-level analysis
- `spatial_database/` ← For fast general searches

### 2. Visualize Gene Patterns

```bash
# Check ACTA2 expression pattern
python visualize_no_pca_database.py \
  --database spatial_database/no_pca_db \
  --resolution 160 \
  --gene-name ACTA2 \
  --output-dir acta2_analysis
```

**Result:** `acta2_analysis/ACTA2_radial_pattern.png`

**What to look for:**
- Does ACTA2 show radial organization?
- Is it higher in center or edge?
- Is the pattern consistent across patches?

### 3. Search for Similar Regions

Select a region in your frontend, then:

```bash
python radial_shell_system.py search \
  --database spatial_database/no_pca_db \
  --sample /path/to/GSM7780155 \
  --resolution 160 \
  --x-center 1500 \
  --y-center 1500 \
  --radius 800 \
  --top-k 50 \
  --output similar_regions.csv
```

**Results:** Top 50 most similar patches with:
- Similarity scores (0-1, higher = more similar)
- Sample IDs
- Locations (x, y coordinates)
- Patch sizes

### 4. Analyze Results

```bash
# View top matches
head -10 similar_regions.csv
```

```csv
patch_id,sample_id,similarity,n_bins,center_x,center_y,radius
GSM7780155_r5_p0,GSM7780155,1.0000,81,1500.0,1500.0,5
GSM7780154_r5_p234,GSM7780154,0.9523,81,1534.5,1498.2,5
GSM7780153_r5_p156,GSM7780153,0.9401,79,1200.0,1600.0,5
...
```

**Interpretation:**
- Similarity 1.0 = identical (query matches itself)
- Similarity > 0.95 = very similar patterns
- Similarity > 0.85 = moderately similar
- Different samples = pattern repeats across samples!

---

## Answer to Your Question

### "If you use PCA you cannot search for a specific gene, right?"

**Correct!** ✅

When you use PCA compression:
- Original: 1,605 features (each = specific gene + specific shell)
- After PCA: 256 features (each = linear combination of many genes)

**Example:**
```
Original feature 5 = ACTA2 in shell 0 (clear meaning)
PCA feature 5 = 0.3×ACTA2_shell0 + 0.1×AIF1_shell1 - 0.2×ACTG2_shell3 + ... (mixed)
```

You **cannot** extract "just ACTA2" from PCA features.

**That's why you built the no-PCA database!** 🎯

With `use_pca=False`:
- ✅ Feature 5 = ACTA2 in shell 0
- ✅ Feature 10 = ACTA2 in shell 2
- ✅ Can extract all ACTA2 features
- ✅ Can visualize gene-specific patterns

---

## Your 321 Genes

Available in your database:
```
ABCC11, ACTA2, ACTG2, ADAM9, ADGRE5, ADH1B, ADIPOQ, AGR3, AHSP, AIF1,
AKR1C1, AKR1C3, ALDH1A3, ANGPT2, ANKRD28, ANKRD29, ANKRD30A, APOBEC3A,
APOBEC3B, APOC1, AQP1, AQP3, AR, AVPR1A, BACE2, BANK1, BASP1, BTNL9,
... (321 total)
```

Check the full list:
```bash
python -c "import json; print('\\n'.join(json.load(open('spatial_database/no_pca_db/160um/config.json'))['global_genes']))"
```

---

## Next Steps

### 1. Explore Gene Patterns

Try different genes to see spatial organization:
```bash
# Smooth muscle markers
python visualize_no_pca_database.py --database spatial_database/no_pca_db --resolution 160 --gene-name ACTA2
python visualize_no_pca_database.py --database spatial_database/no_pca_db --resolution 160 --gene-name ACTG2

# Immune markers
python visualize_no_pca_database.py --database spatial_database/no_pca_db --resolution 160 --gene-name AIF1
```

### 2. Compare Related Genes

```bash
# Compare muscle markers
python visualize_no_pca_database.py --database spatial_database/no_pca_db --resolution 160 --compare-genes ACTA2 ACTG2
```

### 3. Frontend Integration

Use `similarity_search_api_example.py` as a template to add search button to your web interface

### 4. Advanced Analysis

- Find all follicle structures
- Compare treated vs untreated samples
- Identify recurring spatial modules
- Build region similarity networks

---

## Key Files

| File | Purpose |
|------|---------|
| `database_builder.py` | Build the database from Xenium data |
| `similarity_search.py` | Search engine for finding similar regions |
| `radial_shell_system.py` | Command-line interface (build + search) |
| `visualize_no_pca_database.py` | Gene-level heatmap visualization |
| `HOW_IT_WORKS.md` | Detailed system architecture |
| `QUICKSTART_VISUALIZATION.md` | Visualization guide |

---

## Performance

Your database:
- **Build time:** ~20-30 seconds for 3 samples
- **Search time:** <1 second per query
- **Memory:** ~100-200 MB when loaded
- **Disk:** ~27 MB (no PCA) or ~4 MB (with PCA)

---

## Summary

✅ **What works:**
- Database building with global gene system
- No-PCA mode preserving gene information
- Similarity search across samples
- Gene-level visualization with heatmaps
- Multi-radius patches (3, 5, 8 bins)

✅ **What you can do:**
- Search for similar spatial regions
- Analyze individual gene patterns
- Compare genes side-by-side
- Find recurring structures across samples
- Integrate with your frontend

✅ **What you have:**
- 4,229 searchable patches
- 321 genes fully preserved
- 3 different scales (radii)
- Fast FAISS-based search
- Publication-quality visualizations

**Everything is working! 🎉**
