#!/usr/bin/env python3
"""
Radial Shell Encoder - Core encoding logic for spatial transcriptomics patches.

This module provides functions to compute radial shell encodings that capture
the spatial organization of gene expression from center to periphery.

Author: XeniumMundus Project
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, List
from pathlib import Path
from scipy.spatial import cKDTree
import logging

logger = logging.getLogger(__name__)


class RadialShellEncoder:
    """
    Encodes spatial gene expression patterns using concentric shell averaging.

    This creates a rotation/translation invariant representation that captures
    center-periphery organization of gene expression.
    """

    def __init__(self, n_shells: int = 5):
        """
        Initialize the radial shell encoder.

        Parameters:
        -----------
        n_shells : int
            Number of concentric shells to divide patches into
        """
        self.n_shells = n_shells
        self.logger = logging.getLogger(__name__)

    def encode_patch(
        self,
        bin_coords: np.ndarray,
        gene_expression: np.ndarray,
        variable_gene_mask: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Compute radial shell encoding for a spatial patch.

        Parameters:
        -----------
        bin_coords : np.ndarray
            Coordinates of bins, shape (n_bins, 2) with columns [x, y]
        gene_expression : np.ndarray
            Gene expression matrix, shape (n_bins, n_genes)
        variable_gene_mask : np.ndarray, optional
            Boolean mask for spatially variable genes, shape (n_genes,)

        Returns:
        --------
        np.ndarray: Encoded features, shape (n_genes_filtered * n_shells,)
        """
        # Filter to spatially variable genes if mask provided
        if variable_gene_mask is not None:
            gene_expression = gene_expression[:, variable_gene_mask]

        n_bins, n_genes = gene_expression.shape

        if n_bins == 0:
            return np.zeros(n_genes * self.n_shells, dtype=np.float32)

        # Compute centroid
        centroid = bin_coords.mean(axis=0)

        # Calculate radial distances from centroid
        radii = np.linalg.norm(bin_coords - centroid, axis=1)
        max_radius = radii.max()

        if max_radius == 0:
            # All bins at same location - return mean expression repeated
            mean_expr = gene_expression.mean(axis=0)
            return np.tile(mean_expr, self.n_shells).astype(np.float32)

        # Define shell boundaries
        shell_boundaries = np.linspace(0, max_radius, self.n_shells + 1)

        # Compute average expression per shell
        features = []
        for gene_idx in range(n_genes):
            gene_expr = gene_expression[:, gene_idx]

            for shell_idx in range(self.n_shells):
                # Find bins in this shell
                in_shell = (
                    (radii >= shell_boundaries[shell_idx]) &
                    (radii < shell_boundaries[shell_idx + 1])
                )

                if in_shell.sum() > 0:
                    # Average expression in this shell
                    features.append(gene_expr[in_shell].mean())
                else:
                    # No bins in this shell
                    features.append(0.0)

        return np.array(features, dtype=np.float32)

    def encode_from_zarr(
        self,
        zarr_array: np.ndarray,
        bin_indices: np.ndarray,
        bin_size: int,
        variable_gene_mask: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Encode a patch directly from zarr array representation.

        Parameters:
        -----------
        zarr_array : np.ndarray
            Zarr array, shape (n_genes, width, height)
        bin_indices : np.ndarray
            Array of (x_bin, y_bin) indices, shape (n_bins, 2)
        bin_size : int
            Size of each bin in micrometers
        variable_gene_mask : np.ndarray, optional
            Boolean mask for spatially variable genes

        Returns:
        --------
        np.ndarray: Encoded features
        """
        # Extract expression values for selected bins
        x_bins = bin_indices[:, 0]
        y_bins = bin_indices[:, 1]

        # Get expression matrix: (n_genes, n_bins) -> transpose to (n_bins, n_genes)
        gene_expression = zarr_array[:, x_bins, y_bins].T

        # Create coordinates from bin indices
        bin_coords = bin_indices.astype(float) * bin_size

        return self.encode_patch(bin_coords, gene_expression, variable_gene_mask)


class PatchGenerator:
    """
    Generates circular patches from spatial transcriptomics data at multiple radii.
    """

    # Resolution-specific radii (in bin units) targeting 200-1000 μm physical size
    RADII_BY_RESOLUTION = {
        10:  [10, 20, 30, 50, 75, 100],
        20:  [8, 15, 25, 40],
        40:  [5, 10, 15, 20],
        80:  [4, 7, 12],
        160: [3, 5, 8]
    }

    def __init__(self, resolution_um: int, stride_factor: float = 0.5):
        """
        Initialize patch generator.

        Parameters:
        -----------
        resolution_um : int
            Bin size in micrometers (10, 20, 40, 80, or 160)
        stride_factor : float
            Fraction of radius to use as stride between patch centers
        """
        if resolution_um not in self.RADII_BY_RESOLUTION:
            raise ValueError(
                f"Invalid resolution {resolution_um}. "
                f"Must be one of {list(self.RADII_BY_RESOLUTION.keys())}"
            )

        self.resolution_um = resolution_um
        self.radii = self.RADII_BY_RESOLUTION[resolution_um]
        self.stride_factor = stride_factor
        self.logger = logging.getLogger(__name__)

    def generate_patches(
        self,
        bin_coords: np.ndarray,
        sample_id: str,
        min_bins: int = 10
    ) -> List[Dict]:
        """
        Generate patches at multiple radii for a sample.

        Parameters:
        -----------
        bin_coords : np.ndarray
            Coordinates of all bins, shape (n_bins, 2) with columns [x_bin, y_bin]
        sample_id : str
            Identifier for the sample
        min_bins : int
            Minimum number of bins required in a patch

        Returns:
        --------
        List[Dict]: List of patch metadata dictionaries
        """
        # Build KD-tree for efficient neighbor search
        tree = cKDTree(bin_coords)

        # Get coordinate ranges
        x_min, x_max = bin_coords[:, 0].min(), bin_coords[:, 0].max()
        y_min, y_max = bin_coords[:, 1].min(), bin_coords[:, 1].max()

        all_patches = []

        for radius in self.radii:
            stride = max(1, int(radius * self.stride_factor))

            # Grid of potential patch centers
            x_centers = np.arange(x_min + radius, x_max - radius + 1, stride)
            y_centers = np.arange(y_min + radius, y_max - radius + 1, stride)

            patch_id = 0
            for x_center in x_centers:
                for y_center in y_centers:
                    center = np.array([x_center, y_center])

                    # Find bins within radius
                    indices = tree.query_ball_point(center, r=radius)

                    if len(indices) < min_bins:
                        continue

                    all_patches.append({
                        'patch_id': f"{sample_id}_r{radius}_p{patch_id}",
                        'sample_id': sample_id,
                        'radius': radius,
                        'center': center,
                        'bin_indices': np.array(indices),
                        'n_bins': len(indices)
                    })
                    patch_id += 1

        self.logger.info(
            f"Generated {len(all_patches)} patches for sample {sample_id} "
            f"at resolution {self.resolution_um}μm"
        )

        return all_patches


def load_variable_genes(haystack_file: Path, threshold_log_padj: float = -2.0) -> np.ndarray:
    """
    Load spatially variable genes from haystack results.

    Parameters:
    -----------
    haystack_file : Path
        Path to haystack_results.csv file
    threshold_log_padj : float
        Log10 adjusted p-value threshold (default: -2.0 = p_adj < 0.01)

    Returns:
    --------
    np.ndarray: Boolean mask for spatially variable genes
    """
    if not haystack_file.exists():
        logger.warning(f"Haystack file not found: {haystack_file}. Using all genes.")
        return None

    try:
        df = pd.read_csv(haystack_file, index_col=0)

        # Check for required column
        if 'logpval_adj' not in df.columns:
            logger.warning("'logpval_adj' column not found in haystack results. Using all genes.")
            return None

        # Create mask for spatially variable genes
        mask = df['logpval_adj'] <= threshold_log_padj

        n_variable = mask.sum()
        n_total = len(mask)
        logger.info(
            f"Loaded {n_variable}/{n_total} ({100*n_variable/n_total:.1f}%) "
            f"spatially variable genes"
        )

        return mask.values

    except Exception as e:
        logger.error(f"Failed to load haystack results: {e}")
        return None


def auto_select_radius(query_size: int, resolution_um: int) -> int:
    """
    Auto-select appropriate patch radius based on query size.

    Parameters:
    -----------
    query_size : int
        Number of bins in query region
    resolution_um : int
        Resolution in micrometers

    Returns:
    --------
    int: Recommended radius in bin units
    """
    # Estimate radius from area (assuming circular)
    estimated_radius = int(np.sqrt(query_size / np.pi))

    # Get available radii for this resolution
    available_radii = PatchGenerator.RADII_BY_RESOLUTION.get(resolution_um, [10])

    # Find closest radius
    closest = min(available_radii, key=lambda r: abs(r - estimated_radius))

    return closest


if __name__ == "__main__":
    # Example usage
    print("Radial Shell Encoder Module")
    print("=" * 50)

    # Example: encode a simple patch
    encoder = RadialShellEncoder(n_shells=5)

    # Create synthetic data
    n_bins = 100
    n_genes = 50

    # Random circular patch
    angles = np.random.uniform(0, 2*np.pi, n_bins)
    radii = np.random.uniform(0, 50, n_bins)
    bin_coords = np.column_stack([
        radii * np.cos(angles),
        radii * np.sin(angles)
    ])

    # Random expression
    gene_expression = np.random.rand(n_bins, n_genes)

    # Encode
    encoding = encoder.encode_patch(bin_coords, gene_expression)

    print(f"Input: {n_bins} bins, {n_genes} genes")
    print(f"Output encoding shape: {encoding.shape}")
    print(f"Output encoding dimension: {len(encoding)}")
    print(f"Expected: {n_genes} genes × {encoder.n_shells} shells = {n_genes * encoder.n_shells}")
