"""
Example Flask API endpoint for running similarity search from frontend selections.

Add this to your Flask backend (wherever you handle the XeniumProcessor API calls).
"""

from flask import Flask, request, jsonify
from pathlib import Path
import json

# Your existing imports
# from xenium_processor import XeniumProcessor

# Example API endpoint
@app.route('/api/similarity-search', methods=['POST'])
def run_similarity_search():
    """
    API endpoint to run similarity search on a selected region.

    Expected POST JSON:
    {
        "sample_id": "GSM7780155",
        "selection_bounds": {
            "type": "box",  // or "lasso"
            "xRange": [1000, 2000],  // for box
            "yRange": [1000, 2000],   // for box
            // OR for lasso:
            // "lassoPoints": {"x": [...], "y": [...]}
            "binSize": 160
        },
        "database_path": "/path/to/spatial_database",
        "output_folder": "search_results",  // optional
        "top_k": 100,  // optional
        "resolution_um": 160  // optional
    }

    Returns JSON:
    {
        "success": true,
        "results_file": "/path/to/results.csv",
        "num_results": 50,
        "top_matches": [
            {
                "rank": 1,
                "sample_id": "GSM7780154",
                "similarity": 0.95,
                "center_x": 1500.5,
                "center_y": 1500.5,
                "radius": 5,
                "n_bins": 81,
                "patch_id": "GSM7780154_r5_p123"
            },
            ...
        ],
        "query_info": {
            "center_x": 1500,
            "center_y": 1500,
            "radius_um": 500,
            "selection_type": "box"
        },
        "message": "Found 50 similar regions"
    }
    """
    try:
        data = request.get_json()

        # Extract parameters
        sample_id = data.get('sample_id')
        selection_bounds = data.get('selection_bounds')
        database_path = data.get('database_path')
        output_folder = data.get('output_folder', 'search_results')
        top_k = data.get('top_k', 100)
        resolution_um = data.get('resolution_um', 160)
        bin_size = data.get('bin_size', 160)

        # Validate required parameters
        if not sample_id or not selection_bounds or not database_path:
            return jsonify({
                'success': False,
                'error': 'Missing required parameters: sample_id, selection_bounds, database_path'
            }), 400

        # Initialize XeniumProcessor
        # Replace XENIUM_BASE_FOLDER with your actual config
        XENIUM_BASE_FOLDER = "/home/vlad/xenium_mundus/data_full/data"

        with XeniumProcessor(sample_id, XENIUM_BASE_FOLDER) as processor:
            # Run similarity search
            results = processor.run_similarity_search_on_selection(
                selection_bounds=selection_bounds,
                database_path=database_path,
                output_folder=output_folder,
                bin_size=bin_size,
                resolution_um=resolution_um,
                top_k=top_k
            )

        if results['success']:
            return jsonify(results), 200
        else:
            return jsonify(results), 500

    except Exception as e:
        import traceback
        return jsonify({
            'success': False,
            'error': str(e),
            'traceback': traceback.format_exc()
        }), 500


# Frontend JavaScript example
"""
// Example JavaScript code for your frontend button

async function runSimilaritySearch() {
    // Get current selection bounds from your Plotly chart
    const selectionBounds = getCurrentSelection(); // Your existing selection logic

    // Show loading indicator
    showLoadingSpinner();

    try {
        const response = await fetch('/api/similarity-search', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                sample_id: currentSampleId,  // Your current sample
                selection_bounds: selectionBounds,
                database_path: '/path/to/spatial_database',  // Configure this
                output_folder: 'search_results',
                top_k: 100,
                resolution_um: 160
            })
        });

        const results = await response.json();

        if (results.success) {
            console.log(`Found ${results.num_results} similar regions`);
            console.log(`Results saved to: ${results.results_file}`);

            // Display top matches
            displaySearchResults(results.top_matches);

            // Optionally download the CSV
            if (confirm('Download full results as CSV?')) {
                downloadResults(results.results_file);
            }
        } else {
            alert('Search failed: ' + results.error);
        }
    } catch (error) {
        console.error('Search error:', error);
        alert('Search failed: ' + error.message);
    } finally {
        hideLoadingSpinner();
    }
}

function displaySearchResults(topMatches) {
    // Create a table or list showing top matches
    const resultsDiv = document.getElementById('search-results');

    let html = '<h3>Top Similar Regions</h3>';
    html += '<table class="results-table">';
    html += '<tr><th>Rank</th><th>Sample</th><th>Similarity</th><th>Location</th><th>Size</th></tr>';

    topMatches.forEach(match => {
        html += `<tr>
            <td>${match.rank}</td>
            <td>${match.sample_id}</td>
            <td>${(match.similarity * 100).toFixed(1)}%</td>
            <td>(${match.center_x.toFixed(0)}, ${match.center_y.toFixed(0)})</td>
            <td>${match.n_bins} bins</td>
        </tr>`;
    });

    html += '</table>';
    resultsDiv.innerHTML = html;
}

// Add button to your HTML
<button onclick="runSimilaritySearch()" id="similarity-search-btn">
    🔍 Find Similar Regions
</button>
"""
