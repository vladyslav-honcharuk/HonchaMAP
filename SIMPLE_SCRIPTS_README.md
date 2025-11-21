# Simple Scripts - Easy to Understand Version

These scripts do the same thing as the full package, but in simple, readable code.

## What They Do

**The Concept in Plain English:**
1. Take a circular patch of tissue (like a cookie cutter)
2. Divide it into 5 rings (like a bullseye target)
3. Calculate average gene expression in each ring
4. Compare patterns between patches using similarity

**Why rings?**
- Many biological patterns are radial (center → edge)
- Tumor infiltration, immune zones, gradient signals
- Captures spatial organization

## Scripts

### 1. `simple_search.py` - Complete search in one file (263 lines)

Does everything:
- ✅ Loads samples
- ✅ Extracts circular patches
- ✅ Encodes into radial features
- ✅ Builds searchable database
- ✅ Finds similar patches

**Usage:**
```bash
python simple_search.py --data-dir /path/to/samples
```

**What it does:**
```
1. Scans all samples
2. Generates patches every 10 bins
3. For each patch:
   - Divides into 5 concentric rings
   - Calculates average gene expression per ring
   - Stores as feature vector
4. Creates a query patch
5. Finds 10 most similar patches
```

**Output:**
```
Top 10 matches:
 1. GSM7780155 ( 50,  50)  similarity: 1.000
 2. GSM7780154 (120, 140)  similarity: 0.952
 3. GSM7780153 ( 80, 100)  similarity: 0.928
...

Detailed comparison table showing query vs match gene expression in each ring

Spatial visualization: spatial_matches_GENENAME.png
  - Shows query + top 3 matches
  - Original tissue coordinates (not circular)
  - Gene expression as color
  - Blue circle = patch radius
```

### 2. `show_gene_pattern.py` - Show what's stored for one gene (126 lines)

Shows exactly what the database contains for a specific gene.

**Usage:**
```bash
python show_gene_pattern.py \
  --database spatial_database/no_pca_db/160um \
  --gene ACTA2 \
  --radius 5
```

**Output:**
- Left plot: Heatmap of gene across shells and patches
- Right plot: Average pattern (center → edge)
- Text: Is gene higher in center or edge?

## How The Encoding Works

### Input: One circular patch
```
Patch at (100, 100), radius 5 bins
Contains 81 bins (cells/spots)
Each bin has expression for 321 genes
```

### Processing: Divide into rings
```
Ring 1 (center):     distance 0.0 - 1.0 from center
Ring 2:              distance 1.0 - 2.0
Ring 3:              distance 2.0 - 3.0
Ring 4:              distance 3.0 - 4.0
Ring 5 (edge):       distance 4.0 - 5.0
```

### Output: Feature vector
```
For gene ACTA2:
  Ring 1: avg expression = 2.5
  Ring 2: avg expression = 3.1
  Ring 3: avg expression = 2.8
  Ring 4: avg expression = 1.9
  Ring 5: avg expression = 0.7

For gene AIF1:
  Ring 1: avg expression = 0.3
  Ring 2: avg expression = 0.5
  ...

Final vector: [2.5, 3.1, 2.8, 1.9, 0.7, 0.3, 0.5, ...]
             gene0_ring0-4, gene1_ring0-4, ...
             = 321 genes × 5 rings = 1,605 numbers
```

### Search: Find similar patterns
```
Query patch encoding:    [2.5, 3.1, 2.8, 1.9, 0.7, ...]
Database patch 1:        [2.4, 3.0, 2.9, 2.0, 0.8, ...]  ← similarity 0.98
Database patch 2:        [1.2, 1.5, 1.8, 2.1, 2.4, ...]  ← similarity 0.65
...

Use cosine similarity to compare
Higher score = more similar pattern
```

## Example: Understanding a Gene Pattern

Let's say you run:
```bash
python show_gene_pattern.py --database spatial_database/no_pca_db/160um --gene ACTA2
```

**Possible outputs:**

### Pattern 1: Higher in center
```
Ring values: [5.2, 4.1, 3.0, 1.8, 0.5]
→ ACTA2 is HIGHER in the CENTER

Interpretation: Gene enriched in center of patches
```

### Pattern 2: Higher at edges
```
Ring values: [0.8, 1.5, 2.3, 3.4, 4.9]
→ ACTA2 is HIGHER at the EDGES

Interpretation: Gene enriched at periphery
```

### Pattern 3: No pattern
```
Ring values: [2.1, 2.0, 2.2, 1.9, 2.1]
→ ACTA2 has NO CLEAR radial pattern

Interpretation: Uniform distribution
```

## Code Structure (simple_search.py)

```python
# STEP 1: Load sample data (20 lines)
def load_sample(sample_path):
    # Read genes.csv
    # Load zarr file
    # Return genes and expression

# STEP 2: Extract and encode patches (50 lines)
def extract_patch(expression, center_x, center_y, radius):
    # Find bins in circle
    # Extract their expression

def encode_radial(x_coords, y_coords, patch_data, n_shells=5):
    # Calculate distances from center
    # Divide into rings
    # Average expression per ring per gene
    # Return feature vector

# STEP 3: Build database (60 lines)
def build_database(data_dir, radius):
    # Load all samples
    # Generate patches (every 10 bins)
    # Encode each patch
    # Return database

# STEP 4: Search (20 lines)
def search(query_encoding, database, top_k=10):
    # Calculate similarity to all patches
    # Return top K matches

# MAIN: Put it all together (60 lines)
def main():
    # Build database
    # Create query
    # Search
    # Print results
```

## Differences from Full Package

### Full Package (database_builder.py, etc.)
- **Pros:** Production-ready, optimized, FAISS indices, PCA compression
- **Cons:** Complex, many files, harder to understand

### Simple Scripts
- **Pros:** Easy to read, all in one file, shows exactly what happens
- **Cons:** Slower (no FAISS), loads everything in memory, no optimization

**When to use simple scripts:**
- Learning how it works
- Small datasets (< 10 samples)
- Quick prototyping
- Understanding the algorithm

**When to use full package:**
- Production use
- Large datasets (100+ samples)
- Need fast search
- Integration with web app

## Example Session

```bash
# Run the simple search
python simple_search.py --data-dir /path/to/samples

# Output shows:
# - How many samples found
# - How many patches generated
# - Query location
# - Top 10 similar patches with coordinates

# Then visualize a gene
python show_gene_pattern.py \
  --database spatial_database/no_pca_db/160um \
  --gene ACTA2

# Creates ACTA2_simple.png showing:
# - Left: Expression in 50 patches across 5 rings
# - Right: Average pattern
# - Text: Interpretation
```

## Key Takeaways

1. **Radial encoding** = divide patch into rings, average expression per ring
2. **Feature vector** = one number per gene per ring (321 genes × 5 rings = 1,605 features)
3. **Similarity** = cosine similarity between feature vectors
4. **Search** = find patches with most similar feature vectors

That's it! The whole system boils down to this simple concept.
