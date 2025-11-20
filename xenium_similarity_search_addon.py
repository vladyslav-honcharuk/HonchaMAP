"""
Add this method to your XeniumProcessor class to enable similarity search from frontend selections.
"""

import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional
import pandas as pd

# Import the similarity search engine
from similarity_search import MultiResolutionSpatialSearch


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
    """
    Run similarity search for a selected region and save results to a folder.

    Parameters:
    -----------
    selection_bounds : dict
        Selection bounds from frontend with keys:
        - 'type': 'box' or 'lasso'
        - For box: 'xRange' [min, max], 'yRange' [min, max]
        - For lasso: 'lassoPoints' {x: [...], y: [...]}
        - 'binSize': bin size used for visualization
    database_path : str
        Path to the radial shell database directory
    output_folder : str
        Folder to save search results (default: "search_results")
    bin_size : int
        Bin size to use for the search (default: 160)
    resolution_um : int
        Resolution to search in micrometers (default: 160)
    top_k : int
        Number of top results to return (default: 100)
    min_bins, max_bins : int, optional
        Filter results by number of bins in patches

    Returns:
    --------
    dict: Results containing:
        - 'success': bool
        - 'results_file': path to saved CSV
        - 'num_results': number of matches found
        - 'top_matches': list of top matches with details
        - 'query_info': information about the query region
        - 'error': error message if failed
    """
    try:
        # Calculate query center and radius from selection bounds
        query_info = self._calculate_query_from_selection(selection_bounds)

        if 'error' in query_info:
            return {'success': False, 'error': query_info['error']}

        x_center = query_info['center_x']
        y_center = query_info['center_y']
        radius_physical = query_info['radius_um']

        print(f"Running similarity search:")
        print(f"  Center: ({x_center:.1f}, {y_center:.1f})")
        print(f"  Radius: {radius_physical:.1f} μm")
        print(f"  Database: {database_path}")

        # Initialize the search engine
        searcher = MultiResolutionSpatialSearch(database_dir=Path(database_path))

        # Get the sample path (parent of zarr folder)
        sample_path = self.zarr_path.parent

        # Run the search
        results = searcher.search_from_coordinates(
            sample_path=sample_path,
            x_center=x_center,
            y_center=y_center,
            radius_physical=radius_physical,
            resolution_um=resolution_um,
            k=top_k,
            min_bins=min_bins,
            max_bins=max_bins,
            bin_size=bin_size
        )

        if len(results) == 0:
            return {
                'success': True,
                'num_results': 0,
                'top_matches': [],
                'query_info': query_info,
                'message': 'No results found'
            }

        # Create output folder if it doesn't exist
        output_dir = Path(output_folder)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate filename with timestamp
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        sample_name = sample_path.name
        results_filename = f"{sample_name}_search_{timestamp}.csv"
        results_path = output_dir / results_filename

        # Save results to CSV
        results.to_csv(results_path, index=False)
        print(f"Results saved to: {results_path}")

        # Prepare top matches for response
        top_matches = []
        for idx, row in results.head(min(10, top_k)).iterrows():
            top_matches.append({
                'rank': int(idx + 1),
                'sample_id': str(row['sample_id']),
                'similarity': float(row['similarity']),
                'center_x': float(row['center_x']),
                'center_y': float(row['center_y']),
                'radius': int(row['radius']),
                'n_bins': int(row['n_bins']),
                'patch_id': str(row['patch_id'])
            })

        return {
            'success': True,
            'results_file': str(results_path),
            'num_results': len(results),
            'top_matches': top_matches,
            'query_info': query_info,
            'message': f'Found {len(results)} similar regions'
        }

    except Exception as e:
        import traceback
        error_msg = f"Similarity search failed: {str(e)}"
        print(error_msg)
        print(traceback.format_exc())
        return {
            'success': False,
            'error': error_msg,
            'traceback': traceback.format_exc()
        }


def _calculate_query_from_selection(self, selection_bounds: dict) -> dict:
    """
    Calculate query center and radius from frontend selection bounds.

    Parameters:
    -----------
    selection_bounds : dict
        Selection data from frontend

    Returns:
    --------
    dict: Query parameters with center_x, center_y, radius_um
    """
    selection_type = selection_bounds.get('type')

    if selection_type == 'box':
        # For box selection, use center of box
        x_range = selection_bounds['xRange']
        y_range = selection_bounds['yRange']

        center_x = (x_range[0] + x_range[1]) / 2.0
        center_y = (y_range[0] + y_range[1]) / 2.0

        # Calculate radius as half the diagonal
        width = x_range[1] - x_range[0]
        height = y_range[1] - y_range[0]
        radius_um = np.sqrt(width**2 + height**2) / 2.0

        return {
            'center_x': center_x,
            'center_y': center_y,
            'radius_um': radius_um,
            'selection_type': 'box',
            'width': width,
            'height': height
        }

    elif selection_type == 'lasso':
        # For lasso selection, calculate centroid and bounding radius
        lasso_points = selection_bounds['lassoPoints']
        x_coords = np.array(lasso_points['x'])
        y_coords = np.array(lasso_points['y'])

        # Calculate centroid
        center_x = np.mean(x_coords)
        center_y = np.mean(y_coords)

        # Calculate maximum distance from centroid to any lasso point
        distances = np.sqrt((x_coords - center_x)**2 + (y_coords - center_y)**2)
        radius_um = np.max(distances)

        return {
            'center_x': center_x,
            'center_y': center_y,
            'radius_um': radius_um,
            'selection_type': 'lasso',
            'num_points': len(x_coords)
        }

    else:
        return {'error': f'Invalid selection type: {selection_type}'}


# Add these methods to your XeniumProcessor class:
#
# class XeniumProcessor:
#     ...
#
#     # Add the two methods above to your class
#     run_similarity_search_on_selection = run_similarity_search_on_selection
#     _calculate_query_from_selection = _calculate_query_from_selection
