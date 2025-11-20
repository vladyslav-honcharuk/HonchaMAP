# Radial Shell Encoding - Quick Start Guide

## Installation

### 1. Install Dependencies

```bash
pip install -r requirements_radial_shell.txt
```

### 2. Verify Installation

```bash
python test_radial_shell_system.py
```

This will run validation tests on synthetic data. All tests should pass.

## Building Your First Database

### Step 1: Prepare Your Data

Ensure your Xenium samples have:
- Binned zarr files: `zarr/bins_size_{resolution}.zarr.zip`
- Haystack results: `haystack_results.csv` (for spatially variable genes)
- Gene list: `genes.csv`

Directory structure:
```
/path/to/xenium/data/
├── Sample_01/
│   ├── zarr/
│   │   ├── bins_size_10.zarr.zip
│   │   ├── bins_size_20.zarr.zip
│   │   ├── bins_size_40.zarr.zip
│   │   ├── bins_size_80.zarr.zip
│   │   └── bins_size_160.zarr.zip
│   ├── haystack_results.csv
│   └── genes.csv
├── Sample_02/
│   └── ...
└── Sample_03/
    └── ...
```

### Step 2: Build Database

**For all resolutions:**
```bash
python radial_shell_system.py build \
  --data-dir /path/to/xenium/data \
  --output-dir spatial_database \
  --verbose
```

**For specific resolutions (recommended for testing):**
```bash
python radial_shell_system.py build \
  --data-dir /path/to/xenium/data \
  --output-dir spatial_database \
  --resolutions 40 80 \
  --verbose
```

**Build time estimates:**
- 100 samples, 1 resolution: ~30-60 minutes
- 600 samples, 5 resolutions: ~5-8 hours (16 cores)

### Step 3: Verify Database

```bash
python radial_shell_system.py info --database spatial_database
```

Expected output:
```
==================================================================
DATABASE INFO: spatial_database
==================================================================

40μm Resolution:
  Total patches: 300,000
  Samples: 100
    Radius 5: 45,000 patches
    Radius 10: 85,000 patches
    Radius 15: 110,000 patches
    Radius 20: 60,000 patches
  PCA dimensions: 256
  Variance explained: 96.8%
  Storage: 2.15 GB

...
```

## Using the Database

### Example 1: Search from Coordinates

```bash
python radial_shell_system.py search \
  --database spatial_database \
  --sample /path/to/xenium/data/Sample_01 \
  --x-center 5000 \
  --y-center 5000 \
  --radius 500 \
  --resolution 40 \
  --top-k 50 \
  --output search_results.csv
```

### Example 2: Search from Selection (Box)

Create `selection.json`:
```json
{
  "type": "box",
  "xRange": [4500, 5500],
  "yRange": [4500, 5500]
}
```

Then search:
```bash
python radial_shell_system.py search \
  --database spatial_database \
  --sample /path/to/xenium/data/Sample_01 \
  --selection-file selection.json \
  --resolution 40 \
  --output results.csv
```

### Example 3: Search from Selection (Lasso)

Create `lasso_selection.json`:
```json
{
  "type": "lasso",
  "lassoPoints": {
    "x": [4500, 5000, 5500, 5500, 4500],
    "y": [4500, 4500, 5000, 5500, 5500]
  }
}
```

```bash
python radial_shell_system.py search \
  --database spatial_database \
  --sample /path/to/xenium/data/Sample_01 \
  --selection-file lasso_selection.json \
  --resolution 40
```

### Example 4: Filter Results by Size

```bash
python radial_shell_system.py search \
  --database spatial_database \
  --sample /path/to/xenium/data/Sample_01 \
  --x-center 5000 --y-center 5000 --radius 500 \
  --resolution 40 \
  --min-bins 100 \
  --max-bins 800 \
  --top-k 100
```

## Python API Examples

### Search from User Selection

```python
from similarity_search import MultiResolutionSpatialSearch
from pathlib import Path

# Initialize search engine
searcher = MultiResolutionSpatialSearch(Path("spatial_database"))

# Define selection (from your web app)
selection = {
    'type': 'box',
    'xRange': [4000, 6000],
    'yRange': [4000, 6000]
}

# Search
results = searcher.search_from_selection(
    sample_path=Path("/path/to/sample"),
    selection_bounds=selection,
    resolution_um=40,
    k=100,
    min_bins=50,
    max_bins=500
)

# Display
print(f"Found {len(results)} similar patches\n")
for idx, row in results.head(10).iterrows():
    print(f"{idx+1}. {row['sample_id']}")
    print(f"   Similarity: {row['similarity']:.4f}")
    print(f"   Location: ({row['center_x']:.0f}, {row['center_y']:.0f})")
    print(f"   Radius: {row['radius']} bins ({row['n_bins']} bins total)")
    print()
```

### Integration with XeniumProcessor

```python
from xenium_processor import XeniumProcessor
from similarity_search import MultiResolutionSpatialSearch
from pathlib import Path

# Your existing code to get selection
with XeniumProcessor("Sample_01", base_folder="/path/to/data") as processor:
    # User makes a selection in your web app
    selection_bounds = {
        'type': 'box',
        'xRange': [4000, 6000],
        'yRange': [4000, 6000]
    }

    # Search for similar regions
    searcher = MultiResolutionSpatialSearch(Path("spatial_database"))

    results = searcher.search_from_selection(
        sample_path=processor.zarr_path.parent,
        selection_bounds=selection_bounds,
        resolution_um=40,
        k=50
    )

    # Now visualize the results
    for _, patch in results.head(5).iterrows():
        print(f"Similar patch in {patch['sample_id']} "
              f"at ({patch['center_x']:.0f}, {patch['center_y']:.0f})")

        # You can load and visualize this patch
        # ... your visualization code ...
```

## Common Issues and Solutions

### Issue: "Database not found"

**Solution:** Make sure you've built the database first:
```bash
python radial_shell_system.py build --data-dir /path/to/data --output-dir spatial_database
```

### Issue: "Sample not found"

**Solution:** Provide the full path to the sample directory (not just the name):
```bash
--sample /full/path/to/xenium/data/Sample_01
```

### Issue: "No spatially variable genes found"

**Solution:** Either:
1. Generate `haystack_results.csv` for your samples first
2. Or disable variable gene filtering: `--use-variable-genes False`

### Issue: "Out of memory during PCA fitting"

**Solution:** Reduce batch size:
```bash
python radial_shell_system.py build \
  --data-dir /path/to/data \
  --output-dir spatial_database \
  --batch-size 5000
```

### Issue: "Search is slow"

**Solutions:**
1. Use GPU-accelerated FAISS: `pip install faiss-gpu`
2. Reduce number of patches by increasing stride in patch generation
3. Search at lower resolution (e.g., 80μm instead of 10μm)

## Performance Tuning

### For Faster Building

```bash
# Use fewer PCA dimensions (trades quality for speed)
python radial_shell_system.py build \
  --data-dir /path/to/data \
  --output-dir spatial_database \
  --pca-dims 128

# Or build specific resolutions only
python radial_shell_system.py build \
  --data-dir /path/to/data \
  --output-dir spatial_database \
  --resolutions 40 80
```

### For More Detailed Encoding

```bash
# Use more shells for finer spatial detail
python radial_shell_system.py build \
  --data-dir /path/to/data \
  --output-dir spatial_database \
  --n-shells 7 \
  --pca-dims 512
```

### For Better Search Quality

```python
# Use higher k and filter afterward
results = searcher.search_from_coordinates(
    sample_path=sample_path,
    x_center=5000,
    y_center=5000,
    radius_physical=500,
    resolution_um=40,
    k=500,  # Get more candidates
    min_bins=100,
    max_bins=800
)

# Then filter to top 50 by similarity
top_results = results.head(50)
```

## Next Steps

1. **Build database** for your samples
2. **Test searches** on known patterns
3. **Integrate with your app** using the Flask endpoints
4. **Validate results** against manual annotations
5. **Optimize parameters** based on your data characteristics

For more details, see:
- `RADIAL_SHELL_README.md` - Complete documentation
- `radial_shell_integration_example.py` - Flask API examples
- `test_radial_shell_system.py` - Validation tests

## Support

If you encounter issues:
1. Check the log file: `radial_shell_system.log`
2. Run validation tests: `python test_radial_shell_system.py`
3. Enable verbose mode: `--verbose` flag
4. Check database info: `python radial_shell_system.py info --database spatial_database`
