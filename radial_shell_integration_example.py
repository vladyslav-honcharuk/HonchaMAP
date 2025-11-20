#!/usr/bin/env python3
"""
Integration Example - How to integrate Radial Shell Search with your existing app.

This shows how to add similarity search endpoints to your Flask application.

Author: XeniumMundus Project
"""

from flask import Flask, request, jsonify
from pathlib import Path
from typing import Optional
import logging

from xenium_processor import XeniumProcessor
from similarity_search import MultiResolutionSpatialSearch


# Global search engine (initialized once)
SEARCH_ENGINE: Optional[MultiResolutionSpatialSearch] = None
DATABASE_PATH = Path("spatial_database")  # Configure this


def initialize_search_engine():
    """Initialize the search engine on startup."""
    global SEARCH_ENGINE
    if DATABASE_PATH.exists():
        SEARCH_ENGINE = MultiResolutionSpatialSearch(DATABASE_PATH)
        logging.info(f"Initialized search engine with resolutions: {SEARCH_ENGINE.resolutions}")
    else:
        logging.warning(f"Database not found: {DATABASE_PATH}")


def create_app():
    """Create Flask app with radial shell search endpoints."""
    app = Flask(__name__)

    # Initialize search engine
    initialize_search_engine()

    @app.route('/api/radial-search/from-selection', methods=['POST'])
    def radial_search_from_selection():
        """
        Search for similar regions using selection bounds.

        Request JSON:
        {
            "dsid": "sample_id",
            "base_folder": "/path/to/data",
            "selection": {
                "type": "box" | "lasso",
                "xRange": [min, max],  // for box
                "yRange": [min, max],  // for box
                "lassoPoints": {"x": [...], "y": [...]}  // for lasso
            },
            "resolution": 40,
            "binSize": 40,
            "topK": 100,
            "minBins": 50,
            "maxBins": 500
        }

        Returns:
        {
            "success": true,
            "results": [
                {
                    "patch_id": "...",
                    "sample_id": "...",
                    "similarity": 0.95,
                    "center_x": 5000,
                    "center_y": 5000,
                    "radius": 15,
                    "n_bins": 314
                },
                ...
            ],
            "query_info": {
                "n_query_bins": 500,
                "resolution": 40,
                "auto_radius": 15
            }
        }
        """
        if SEARCH_ENGINE is None:
            return jsonify({
                'success': False,
                'error': 'Search engine not initialized. Database not found.'
            }), 503

        try:
            data = request.json

            # Extract parameters
            dsid = data['dsid']
            base_folder = data.get('base_folder', '/path/to/default')
            selection = data['selection']
            resolution_um = data.get('resolution', 40)
            bin_size = data.get('binSize', resolution_um)
            top_k = data.get('topK', 100)
            min_bins = data.get('minBins')
            max_bins = data.get('maxBins')

            # Find sample path
            # You may need to adjust this based on your directory structure
            from glob import glob
            search_patterns = [
                f"{base_folder}/*/*/{dsid}",
                f"{base_folder}/*/{dsid}",
            ]

            sample_path = None
            for pattern in search_patterns:
                matches = glob(pattern)
                if matches:
                    sample_path = Path(matches[0])
                    break

            if sample_path is None:
                return jsonify({
                    'success': False,
                    'error': f'Sample not found: {dsid}'
                }), 404

            # Perform search
            results_df = SEARCH_ENGINE.search_from_selection(
                sample_path=sample_path,
                selection_bounds=selection,
                resolution_um=resolution_um,
                k=top_k,
                min_bins=min_bins,
                max_bins=max_bins,
                bin_size=bin_size
            )

            # Convert to JSON-friendly format
            results = results_df.to_dict('records')

            # Add query info
            n_query_bins = len(results_df) if len(results_df) > 0 else 0

            return jsonify({
                'success': True,
                'results': results,
                'query_info': {
                    'n_results': len(results),
                    'resolution': resolution_um,
                    'top_k': top_k
                }
            })

        except Exception as e:
            logging.error(f"Search error: {e}")
            import traceback
            traceback.print_exc()
            return jsonify({
                'success': False,
                'error': str(e)
            }), 500

    @app.route('/api/radial-search/from-coordinates', methods=['POST'])
    def radial_search_from_coordinates():
        """
        Search for similar regions using center coordinates and radius.

        Request JSON:
        {
            "dsid": "sample_id",
            "base_folder": "/path/to/data",
            "xCenter": 5000,
            "yCenter": 5000,
            "radius": 500,  // in micrometers
            "resolution": 40,
            "binSize": 40,
            "topK": 100,
            "minBins": 50,
            "maxBins": 500
        }
        """
        if SEARCH_ENGINE is None:
            return jsonify({
                'success': False,
                'error': 'Search engine not initialized.'
            }), 503

        try:
            data = request.json

            # Extract parameters
            dsid = data['dsid']
            base_folder = data.get('base_folder', '/path/to/default')
            x_center = data['xCenter']
            y_center = data['yCenter']
            radius = data['radius']
            resolution_um = data.get('resolution', 40)
            bin_size = data.get('binSize', resolution_um)
            top_k = data.get('topK', 100)
            min_bins = data.get('minBins')
            max_bins = data.get('maxBins')

            # Find sample path (same as above)
            from glob import glob
            search_patterns = [
                f"{base_folder}/*/*/{dsid}",
                f"{base_folder}/*/{dsid}",
            ]

            sample_path = None
            for pattern in search_patterns:
                matches = glob(pattern)
                if matches:
                    sample_path = Path(matches[0])
                    break

            if sample_path is None:
                return jsonify({
                    'success': False,
                    'error': f'Sample not found: {dsid}'
                }), 404

            # Perform search
            results_df = SEARCH_ENGINE.search_from_coordinates(
                sample_path=sample_path,
                x_center=x_center,
                y_center=y_center,
                radius_physical=radius,
                resolution_um=resolution_um,
                k=top_k,
                min_bins=min_bins,
                max_bins=max_bins,
                bin_size=bin_size
            )

            # Convert to JSON
            results = results_df.to_dict('records')

            return jsonify({
                'success': True,
                'results': results,
                'query_info': {
                    'n_results': len(results),
                    'x_center': x_center,
                    'y_center': y_center,
                    'radius_um': radius,
                    'resolution': resolution_um
                }
            })

        except Exception as e:
            logging.error(f"Search error: {e}")
            return jsonify({
                'success': False,
                'error': str(e)
            }), 500

    @app.route('/api/radial-search/database-info', methods=['GET'])
    def radial_search_database_info():
        """
        Get information about the radial shell database.

        Returns:
        {
            "success": true,
            "resolutions": [10, 20, 40, 80, 160],
            "total_patches": 7300000,
            "by_resolution": {
                "10": {"patches": 1200000, "samples": 600},
                "20": {"patches": 480000, "samples": 600},
                ...
            }
        }
        """
        if SEARCH_ENGINE is None:
            return jsonify({
                'success': False,
                'error': 'Search engine not initialized.'
            }), 503

        try:
            info = {
                'success': True,
                'resolutions': SEARCH_ENGINE.resolutions,
                'database_path': str(DATABASE_PATH),
                'by_resolution': {}
            }

            # Get stats for each resolution
            for res in SEARCH_ENGINE.resolutions:
                if res in SEARCH_ENGINE.metadata:
                    metadata = SEARCH_ENGINE.metadata[res]
                    info['by_resolution'][str(res)] = {
                        'patches': len(metadata),
                        'samples': metadata['sample_id'].nunique(),
                        'radii': sorted(metadata['radius'].unique().tolist())
                    }

            return jsonify(info)

        except Exception as e:
            logging.error(f"Info error: {e}")
            return jsonify({
                'success': False,
                'error': str(e)
            }), 500

    return app


# Example JavaScript client code
EXAMPLE_JS_CLIENT = """
// Example: Call radial search from JavaScript

async function searchSimilarRegions(selectionData) {
    const response = await fetch('/api/radial-search/from-selection', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({
            dsid: currentSampleId,
            base_folder: XENIUM_FOLDER,
            selection: selectionData,  // From your selection tool
            resolution: currentBinSize,
            binSize: currentBinSize,
            topK: 100,
            minBins: 50,
            maxBins: 500
        })
    });

    const data = await response.json();

    if (data.success) {
        console.log(`Found ${data.results.length} similar patches`);

        // Display results
        data.results.forEach((result, idx) => {
            console.log(`${idx+1}. ${result.sample_id} - Similarity: ${result.similarity.toFixed(4)}`);

            // You can visualize these patches by loading them
            visualizeSimilarPatch(result);
        });
    } else {
        console.error('Search failed:', data.error);
    }
}

// Example: Visualize a similar patch
async function visualizeSimilarPatch(patchInfo) {
    // Load the sample
    const processor = new XeniumProcessor(patchInfo.sample_id, XENIUM_FOLDER);

    // Get expression for the patch region
    const response = await fetch('/api/get-patch-expression', {
        method: 'POST',
        body: JSON.stringify({
            dsid: patchInfo.sample_id,
            center_x: patchInfo.center_x,
            center_y: patchInfo.center_y,
            radius: patchInfo.radius,
            resolution: currentBinSize
        })
    });

    // Visualize on your plot
    // ... your visualization code ...
}
"""


if __name__ == "__main__":
    # Example: Run the Flask app
    app = create_app()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    print("=" * 70)
    print("RADIAL SHELL SEARCH - FLASK INTEGRATION")
    print("=" * 70)
    print("\nEndpoints:")
    print("  POST /api/radial-search/from-selection")
    print("  POST /api/radial-search/from-coordinates")
    print("  GET  /api/radial-search/database-info")
    print("\nStarting server...")

    app.run(debug=True, host='0.0.0.0', port=5001)
