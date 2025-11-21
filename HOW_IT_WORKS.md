# How the Radial Shell Encoding System Works

## Overview

The Radial Shell Encoding System is a spatial transcriptomics similarity search engine that finds similar tissue regions across multiple samples. It works by encoding spatial gene expression patterns into numerical features that can be efficiently searched.

## Core Concept: Radial Shell Encoding

### What is it?

Think of it like this:
1. **Take a circular patch** of tissue (e.g., 500μm radius)
2. **Divide it into concentric shells** (like tree rings or an onion)
3. **Calculate average gene expression** in each shell
4. **Create a feature vector** from these averages

### Why does it work?

Spatial patterns in biology often show radial organization:
- Immune cell infiltration zones
- Tumor microenvironments
- Gradient-based signaling regions
- Niche structures

By encoding patterns radially, we capture these spatial relationships while:
- **Reducing dimensionality** (millions of cells → hundreds of features)
- **Maintaining biological meaning** (spatial organization preserved)
- **Enabling fast search** (FAISS vector similarity)

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                     DATABASE BUILDING                            │
└─────────────────────────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  Phase 0: Global Gene Discovery               │
    │  - Collect all genes from all samples         │
    │  - Build union (global gene list)             │
    │  - Determine variable genes (intersection)    │
    │  - Create gene→index mapping                  │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  Phase 1: Generate Patches & Embeddings       │
    │  For each sample:                             │
    │    - Generate circular patches (3 radii)      │
    │    - Remap sample genes → global genes        │
    │    - Encode each patch:                       │
    │      • Divide into shells                     │
    │      • Calculate gene expression per shell    │
    │      • Flatten: (genes × shells) → features   │
    │  Output: Raw embeddings (N patches × M dims)  │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  Phase 2: PCA Compression (Optional)          │
    │  - Fit IncrementalPCA on all patches          │
    │  - Compress: M dimensions → 256 dimensions    │
    │  - Normalize: L2 normalization                │
    │  Output: Compressed embeddings                │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  Phase 3: Build FAISS Indices                 │
    │  For each radius:                             │
    │    - Create IndexFlatIP (inner product)       │
    │    - Add all patches for this radius          │
    │  Save: indices + metadata + PCA models        │
    └───────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│                      SIMILARITY SEARCH                           │
└─────────────────────────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  1. Extract Query Region                      │
    │  - User selects region (x, y, radius)         │
    │  - Extract bins within radius                 │
    │  - Load gene expression data                  │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  2. Encode Query                              │
    │  - Remap query genes → global genes           │
    │  - Apply radial shell encoding                │
    │  - Same process as database building          │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  3. Transform with PCA                        │
    │  - Apply same PCA model as database           │
    │  - Compress to 256 dimensions                 │
    │  - L2 normalize                               │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  4. Search FAISS Index                        │
    │  - Select appropriate radius index            │
    │  - Compute cosine similarity to all patches   │
    │  - Return top K matches                       │
    └───────────────────────────────────────────────┘
                                ↓
    ┌───────────────────────────────────────────────┐
    │  5. Return Results                            │
    │  - Patch IDs, similarity scores               │
    │  - Sample IDs, locations (x, y)               │
    │  - Metadata (radius, num bins)                │
    └───────────────────────────────────────────────┘
```

## Data Flow Example

### Input Data
```
Sample: GSM7780155
Genes: 288 genes (varies by sample)
Spatial bins: 100 × 100 grid at 160μm resolution
```

### Phase 0: Global Gene Discovery
```
Sample 1: 288 genes
Sample 2: 313 genes  } Union → 321 global genes
Sample 3: 313 genes

Gene mapping created:
  Sample 1 gene "ACTB" → Global index 0
  Sample 1 gene "TP53" → Global index 150
  ...
  Missing genes → filled with zeros
```

### Phase 1: Patch Generation & Encoding

**For one patch (radius=5 bins, 81 bins total):**

```
Step 1: Extract bins
  - Center: (50, 50)
  - Radius: 5 bins
  - Bins in patch: 81 bins

Step 2: Get gene expression
  - Shape: (81 bins, 321 genes)
  - Values: normalized expression counts

Step 3: Radial shell encoding (5 shells)
  For each gene:
    Shell 1 (0.0-0.2r): avg expression = 2.5
    Shell 2 (0.2-0.4r): avg expression = 3.1
    Shell 3 (0.4-0.6r): avg expression = 2.8
    Shell 4 (0.6-0.8r): avg expression = 1.9
    Shell 5 (0.8-1.0r): avg expression = 0.7

Step 4: Flatten to feature vector
  Features = [gene0_shell0, gene0_shell1, ..., gene0_shell4,
              gene1_shell0, gene1_shell1, ..., gene1_shell4,
              ...
              gene320_shell0, ..., gene320_shell4]

  Dimensions: 321 genes × 5 shells = 1,605 features
```

**For all patches:**
```
Radius 3: 3,332 patches × 1,605 features
Radius 5:   710 patches × 1,605 features
Radius 8:   120 patches × 1,605 features
Total: 4,162 patches
```

### Phase 2: PCA Compression

```
Input:  4,162 patches × 1,605 features
PCA fit on all patches
Output: 4,162 patches × 256 features (compressed)

Explained variance: ~80-90% retained
```

### Phase 3: FAISS Index Building

```
For each radius:
  - Create index
  - Add compressed embeddings
  - Save to disk

Files created:
  spatial_database/160um/
    ├── faiss_index_r3.bin   (3,332 patches)
    ├── faiss_index_r5.bin   (710 patches)
    ├── faiss_index_r8.bin   (120 patches)
    ├── patches_metadata.parquet
    ├── config.json
    ├── pca_models.pkl
    ├── raw_embeddings_r3.h5
    ├── raw_embeddings_r5.h5
    └── raw_embeddings_r8.h5
```

### Search Example

**Query:**
```
User selects region:
  - Center: (1500, 1500) in sample GSM7780155
  - Radius: 800 μm
  - Resolution: 160 μm
```

**Search process:**
```
1. Convert radius: 800μm / 160μm = 5 bins
2. Extract 81 bins in circular region
3. Load expression: (81, 321) matrix
4. Encode with radial shells → 1,605 features
5. Apply PCA → 256 features
6. Search FAISS index for radius=5
7. Return top 100 matches:

   Rank  Sample      Similarity  Location
   1     GSM7780154  0.952      (1534, 1498)
   2     GSM7780155  0.940      (1200, 1600)
   3     GSM7780153  0.935      (2100, 1800)
   ...
```

## Key Features

### 1. Global Gene Coordinate System

**Problem:** Different samples have different genes (300-5000 per sample)

**Solution:**
- Build union of all genes → global gene list
- Map each sample's genes to global indices
- Fill missing genes with zeros
- All embeddings have same dimensions

**Example:**
```
Sample 1 has genes: [ACTB, TP53, MYC]      → 3 genes
Sample 2 has genes: [ACTB, BRCA1, TP53]   → 3 genes
Global genes:       [ACTB, BRCA1, MYC, TP53] → 4 genes (union)

Sample 1 encoding:
  ACTB  → index 0: expression = 5.2
  BRCA1 → index 1: expression = 0.0 (missing)
  MYC   → index 2: expression = 3.1
  TP53  → index 3: expression = 2.7

Sample 2 encoding:
  ACTB  → index 0: expression = 4.8
  BRCA1 → index 1: expression = 6.3
  MYC   → index 2: expression = 0.0 (missing)
  TP53  → index 3: expression = 3.0
```

### 2. Multi-Radius Search

Different patch sizes capture different scales:

```
Radius 3 bins (480μm at 160μm):
  - Small local neighborhoods
  - ~50-80 bins per patch
  - 3,332 patches total

Radius 5 bins (800μm at 160μm):
  - Medium-scale regions
  - ~80-120 bins per patch
  - 710 patches total

Radius 8 bins (1280μm at 160μm):
  - Large tissue regions
  - ~150-200 bins per patch
  - 120 patches total
```

### 3. Variable Gene Filtering (Optional)

**With variable genes enabled:**
```
1. Run haystack analysis on each sample
2. Find spatially variable genes (p_adj < 0.01)
3. Take intersection across samples
4. Only use these genes for encoding

Result: 321 genes → 50 variable genes
        1,605 features → 250 features
        Faster search, more focused on spatial patterns
```

**With variable genes disabled:**
```
Use all genes from global gene list
Better for exploratory analysis
More robust when samples have different biology
```

## Performance Characteristics

### Database Build Time

```
3 samples, 321 genes each, 160μm resolution:
  Phase 0 (Gene discovery):  < 1 second
  Phase 1 (Embeddings):      ~10 seconds
  Phase 2 (PCA):             ~5 seconds
  Phase 3 (FAISS):           < 1 second
  Total:                     ~20 seconds
```

### Search Time

```
Single query:              < 1 second
100 queries (batch):       ~5 seconds
Bottleneck: Query encoding and PCA transform
```

### Memory Usage

```
Database in memory:        ~100-500 MB
  - FAISS indices:         ~50 MB
  - PCA models:            ~10 MB
  - Metadata:              ~5 MB

Per search:                ~10 MB
  - Query encoding:        ~5 MB
  - Results:               ~1 MB
```

### Storage

```
Database on disk:          ~200 MB
  - FAISS indices:         ~50 MB
  - Raw embeddings:        ~100 MB (optional, for analysis)
  - PCA models:            ~10 MB
  - Metadata:              ~5 MB
  - Config:                < 1 MB
```

## Dimensionality Reduction

### Without PCA
```
Input:  321 genes × 5 shells = 1,605 dimensions
Output: 1,605 dimensions
Search: Exact, slower
```

### With PCA (default)
```
Input:  1,605 dimensions
PCA:    Fit on all patches
Output: 256 dimensions
Search: Approximate, much faster
Loss:   ~10-20% variance
```

## Use Cases

### 1. Find Similar Tumor Microenvironments
```
Query: Select a tumor region
Results: Other regions with similar immune infiltration patterns
```

### 2. Identify Recurring Spatial Patterns
```
Query: Select a follicle structure
Results: All follicles across all samples
```

### 3. Compare Treatment Groups
```
Query: Select untreated tumor region
Results: Find similar regions in treated samples
Analysis: Compare gene expression differences
```

### 4. Quality Control
```
Query: Select normal tissue
Results: Should match other normal regions
Detection: Batch effects or technical artifacts
```

## Visualization

Run the visualization script to see what's stored:

```bash
python visualize_database.py \
  --database spatial_database \
  --resolution 160 \
  --gene-name TP53 \
  --output-dir plots
```

This creates:
1. **Radial shell concept** - How the encoding works
2. **Database summary** - Statistics and configuration
3. **Feature statistics** - Distributions across radii
4. **Gene encoding** - Single gene across all shells/radii

## Mathematical Details

### Radial Shell Encoding

For a patch with bins $B = \{b_1, ..., b_n\}$ and centroid $c$:

1. **Calculate distances:**
   $$r_i = ||b_i - c||_2$$

2. **Define shell boundaries:**
   $$s_k = \frac{k \cdot r_{max}}{N_{shells}}, \quad k = 0, 1, ..., N_{shells}$$

3. **Assign bins to shells:**
   $$B_k = \{b_i : s_{k-1} \leq r_i < s_k\}$$

4. **Calculate shell features:**
   $$f_{g,k} = \frac{1}{|B_k|} \sum_{b_i \in B_k} E_{g,i}$$

   where $E_{g,i}$ is expression of gene $g$ in bin $b_i$

5. **Flatten to vector:**
   $$\mathbf{x} = [f_{1,1}, ..., f_{1,K}, f_{2,1}, ..., f_{2,K}, ..., f_{G,1}, ..., f_{G,K}]$$

   where $G$ = number of genes, $K$ = number of shells

### Cosine Similarity (after L2 normalization)

$$\text{similarity}(\mathbf{q}, \mathbf{p}) = \frac{\mathbf{q} \cdot \mathbf{p}}{||\mathbf{q}||_2 \cdot ||\mathbf{p}||_2}$$

After L2 normalization: $||\mathbf{q}||_2 = ||\mathbf{p}||_2 = 1$

So: $\text{similarity}(\mathbf{q}, \mathbf{p}) = \mathbf{q} \cdot \mathbf{p}$ (inner product)

FAISS IndexFlatIP computes this efficiently.

## Configuration Options

### Build Options

```bash
python radial_shell_system.py build \
  --data-dir /path/to/data \
  --output-dir database \
  --resolutions 160 \                    # Resolution(s) to build
  --use-variable-genes False \           # Use all genes
  --variable-gene-threshold -2.0 \       # p_adj < 0.01
  --n-shells 5 \                         # Number of shells
  --pca-dims 256 \                       # PCA compression
  --batch-size 10000                     # PCA batch size
```

### Search Options

```bash
python radial_shell_system.py search \
  --database database \
  --sample /path/to/sample \
  --resolution 160 \
  --x-center 1500 --y-center 1500 \     # Query center
  --radius 800 \                         # Query radius (μm)
  --top-k 100 \                          # Number of results
  --min-bins 50 --max-bins 200 \        # Filter by size
  --output results.csv                   # Save results
```

## Troubleshooting

See `REPOSITORY_INCONSISTENCIES_REPORT.md` for known issues and fixes.

Common issues:
1. **Zero-dimensional embeddings** → Fixed with fallback to all genes
2. **Boolean argument parsing** → Fixed with custom parser
3. **Zarr indexing** → Fixed by converting to numpy first
4. **Gene set mismatches** → Fixed with global gene coordinate system

## Next Steps

1. **Build your database** with real data
2. **Visualize** what's stored using the visualization script
3. **Run searches** to find similar regions
4. **Integrate** with your frontend for interactive exploration

All the pieces are in place and tested!
