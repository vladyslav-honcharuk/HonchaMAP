# Similarity Search Integration for XeniumProcessor

This guide explains how to add similarity search functionality to your web frontend, allowing users to select regions and find similar spatial patterns.

## Overview

The similarity search integration allows users to:
1. Select a region on the tissue visualization (box or lasso selection)
2. Click a button to search for similar regions
3. View top matching regions from the database
4. Download full results as CSV

## Files Created

1. `xenium_similarity_search_addon.py` - Methods to add to XeniumProcessor class
2. `similarity_search_api_example.py` - Flask API endpoint example
3. This documentation file

## Installation Steps

### Step 1: Add Methods to XeniumProcessor

Add the two methods from `xenium_similarity_search_addon.py` to your `XeniumProcessor` class:

```python
# In your xenium_processor.py file

class XeniumProcessor:
    # ... your existing methods ...

    def run_similarity_search_on_selection(
        self,
        selection_bounds: dict,
        database_path: str,
        output_folder: str = "search_results",
        bin_size: int = 160,
        resolution_um: int = 160,
        top_k: int = 100,
        min_bins: Optional[int] = None,
        max_bins: Optional[int] = None
    ) -> dict:
        """Run similarity search for a selected region."""
        # Copy implementation from xenium_similarity_search_addon.py
        ...

    def _calculate_query_from_selection(self, selection_bounds: dict) -> dict:
        """Calculate query center and radius from selection bounds."""
        # Copy implementation from xenium_similarity_search_addon.py
        ...
```

Or simply import and attach them:

```python
from xenium_similarity_search_addon import (
    run_similarity_search_on_selection,
    _calculate_query_from_selection
)

class XeniumProcessor:
    # ... your existing code ...

# Add the methods to the class
XeniumProcessor.run_similarity_search_on_selection = run_similarity_search_on_selection
XeniumProcessor._calculate_query_from_selection = _calculate_query_from_selection
```

### Step 2: Add Flask API Endpoint

Add the API endpoint to your Flask backend (see `similarity_search_api_example.py`):

```python
@app.route('/api/similarity-search', methods=['POST'])
def run_similarity_search():
    """API endpoint to run similarity search on selected region."""
    # See similarity_search_api_example.py for full implementation
    ...
```

### Step 3: Add Frontend Button and JavaScript

Add a button to your HTML interface:

```html
<!-- Add this near your other analysis buttons -->
<button onclick="runSimilaritySearch()"
        id="similarity-search-btn"
        class="btn btn-primary">
    🔍 Find Similar Regions
</button>

<!-- Results display area -->
<div id="search-results" class="search-results-panel"></div>
```

Add the JavaScript function (see `similarity_search_api_example.py` for full code):

```javascript
async function runSimilaritySearch() {
    const selectionBounds = getCurrentSelection(); // Your existing selection

    const response = await fetch('/api/similarity-search', {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({
            sample_id: currentSampleId,
            selection_bounds: selectionBounds,
            database_path: '/path/to/spatial_database',
            top_k: 100
        })
    });

    const results = await response.json();
    if (results.success) {
        displaySearchResults(results.top_matches);
    }
}
```

## Usage Example

### 1. Build the Database (One-time)

First, build the similarity search database using your Xenium data:

```bash
python radial_shell_system.py build \
  --data-dir /path/to/xenium/data \
  --output-dir spatial_database \
  --resolutions 160 \
  --use-variable-genes False
```

### 2. User Workflow

1. **User selects a region** on the tissue visualization (box or lasso)
2. **User clicks "Find Similar Regions"** button
3. **Backend processes the search:**
   - Calculates center/radius from selection
   - Encodes the selected region using radial shell encoding
   - Searches the database using FAISS for similar patterns
   - Returns top matches sorted by similarity

4. **Frontend displays results:**
   - Table showing top matches with similarity scores
   - Sample IDs, locations, and sizes
   - Option to download full results as CSV

### 3. Example API Call

```json
POST /api/similarity-search

{
  "sample_id": "GSM7780155",
  "selection_bounds": {
    "type": "box",
    "xRange": [1000, 2000],
    "yRange": [1000, 2000],
    "binSize": 160
  },
  "database_path": "/home/user/spatial_database",
  "top_k": 100,
  "resolution_um": 160
}
```

### 4. Example Response

```json
{
  "success": true,
  "results_file": "search_results/GSM7780155_search_20251120_143052.csv",
  "num_results": 87,
  "top_matches": [
    {
      "rank": 1,
      "sample_id": "GSM7780154",
      "similarity": 0.9523,
      "center_x": 1534.5,
      "center_y": 1498.2,
      "radius": 5,
      "n_bins": 81,
      "patch_id": "GSM7780154_r5_p0234"
    },
    {
      "rank": 2,
      "sample_id": "GSM7780155",
      "similarity": 0.9401,
      "center_x": 1200.0,
      "center_y": 1600.0,
      "radius": 5,
      "n_bins": 81,
      "patch_id": "GSM7780155_r5_p0156"
    }
    // ... more results
  ],
  "query_info": {
    "center_x": 1500,
    "center_y": 1500,
    "radius_um": 707.1,
    "selection_type": "box",
    "width": 1000,
    "height": 1000
  },
  "message": "Found 87 similar regions"
}
```

## Configuration

### Database Path

Set the database path in your configuration:

```python
# In your Flask app config
SPATIAL_DATABASE_PATH = "/home/vlad/spatial_database"
```

### Output Folder

Results are saved to `search_results/` by default. Configure as needed:

```python
results = processor.run_similarity_search_on_selection(
    selection_bounds=selection_bounds,
    database_path=database_path,
    output_folder="custom_output_folder",  # Change this
    top_k=100
)
```

### Search Parameters

- `top_k`: Number of results to return (default: 100)
- `resolution_um`: Resolution to search at (default: 160)
- `bin_size`: Bin size for query encoding (default: 160)
- `min_bins`, `max_bins`: Filter by patch size

## Results Format

### CSV File

Results are saved as CSV with these columns:
- `patch_id` - Unique identifier for the patch
- `sample_id` - Sample name (e.g., "GSM7780155")
- `similarity` - Cosine similarity score (0-1, higher = more similar)
- `n_bins` - Number of bins in the patch
- `center_x`, `center_y` - Patch center coordinates (bin units)
- `radius` - Patch radius (bins)
- `resolution_um` - Resolution

### In-Memory Results

The method returns a dictionary with:
- `success`: Boolean indicating success
- `results_file`: Path to saved CSV
- `num_results`: Total number of matches
- `top_matches`: List of top 10 matches with all details
- `query_info`: Information about the query region
- `message`: Human-readable status message

## Error Handling

The method handles various error cases:

```python
# Example error responses
{
    "success": false,
    "error": "No bins found in selection area"
}

{
    "success": false,
    "error": "Database not found: /path/to/database"
}

{
    "success": true,
    "num_results": 0,
    "message": "No results found"
}
```

## Advanced Usage

### Filtering by Patch Size

```python
results = processor.run_similarity_search_on_selection(
    selection_bounds=selection_bounds,
    database_path=database_path,
    min_bins=50,   # Only patches with >= 50 bins
    max_bins=200,  # Only patches with <= 200 bins
    top_k=50
)
```

### Multi-Resolution Search

Search at different resolutions:

```python
# Search at 80μm resolution
results_80 = processor.run_similarity_search_on_selection(
    selection_bounds=selection_bounds,
    database_path=database_path,
    resolution_um=80,
    bin_size=80
)

# Search at 160μm resolution
results_160 = processor.run_similarity_search_on_selection(
    selection_bounds=selection_bounds,
    database_path=database_path,
    resolution_um=160,
    bin_size=160
)
```

## Performance Notes

- **Database Build**: One-time operation, can take hours for large datasets
- **Search Speed**: Typically <1 second for 100 results from thousands of patches
- **Memory**: Database is loaded once and cached in the search engine
- **Disk Space**: Results CSVs are small (~10-50 KB per search)

## Troubleshooting

### "Database not found" Error

Ensure you've built the database first:
```bash
python radial_shell_system.py build --data-dir /path/to/data --output-dir spatial_database
```

### "No bins found in selection area" Error

The selected region is too small or outside the tissue bounds. Try:
- Selecting a larger region
- Checking coordinate ranges match your data

### "Failed to open zarr file" Error

Ensure the bin_size parameter matches your data:
```python
# Check available bin sizes in your data
ls /path/to/sample/zarr/
# bins_size_10.zarr.zip, bins_size_40.zarr.zip, bins_size_160.zarr.zip, etc.
```

## Next Steps

After integration, you can:

1. **Visualize similar patches** on the tissue
2. **Compare gene expression** between query and matches
3. **Export match coordinates** for downstream analysis
4. **Build region similarity networks**
5. **Identify recurring spatial patterns** across samples

## Support

For issues or questions, check:
- `REPOSITORY_INCONSISTENCIES_REPORT.md` - Known issues
- `RADIAL_SHELL_README.md` - Radial shell encoding details
- Test the system with synthetic data first
