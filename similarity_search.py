#!/usr/bin/env python3
"""
Similarity Search - Search engine for radial shell encoded spatial patterns.

This module provides a search interface to find similar spatial regions
across the multi-resolution database.

Author: XeniumMundus Project
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union
import pickle
import faiss
import h5py
import logging
import zarr

from radial_shell_encoder import (
    RadialShellEncoder,
    PatchGenerator,
    load_variable_genes,
    auto_select_radius
)


class MultiResolutionSpatialSearch:
    """
    Search engine supporting all resolutions with lazy loading.
    """

    def __init__(self, database_path: Path):
        """
        Initialize search engine.

        Parameters:
        -----------
        database_path : Path
            Path to database directory
        """
        self.db_path = Path(database_path)
        self.n_shells = 5  # Default, will be inferred from data

        if not self.db_path.exists():
            raise FileNotFoundError(f"Database not found: {self.db_path}")

        # Lazy-loaded per resolution
        self.metadata = {}
        self.indices = {}
        self.pca_models = {}
        self.embeddings = {}

        # Available resolutions
        self.resolutions = self._detect_resolutions()
        self.logger = logging.getLogger(__name__)

        self.logger.info(f"Initialized search engine with resolutions: {self.resolutions}")

    def _detect_resolutions(self) -> List[int]:
        """Detect available resolutions in database."""
        resolutions = []
        for item in self.db_path.iterdir():
            if item.is_dir() and item.name.endswith('um'):
                try:
                    res = int(item.name[:-2])
                    resolutions.append(res)
                except ValueError:
                    continue
        return sorted(resolutions)

    def _load_resolution(self, resolution_um: int, radius: int):
        """
        Load data for specific resolution and radius.

        Parameters:
        -----------
        resolution_um : int
            Resolution in micrometers
        radius : int
            Patch radius in bin units
        """
        key = (resolution_um, radius)
        if key in self.indices:
            return  # Already loaded

        resolution_dir = self.db_path / f"{resolution_um}um"

        if not resolution_dir.exists():
            raise FileNotFoundError(
                f"Resolution {resolution_um}μm not found in database"
            )

        # Load metadata (once per resolution)
        if resolution_um not in self.metadata:
            metadata_path = resolution_dir / "patches_metadata.parquet"
            if not metadata_path.exists():
                raise FileNotFoundError(f"Metadata not found: {metadata_path}")
            self.metadata[resolution_um] = pd.read_parquet(metadata_path)
            self.logger.info(
                f"Loaded metadata for {resolution_um}μm: "
                f"{len(self.metadata[resolution_um]):,} patches"
            )

        # Load PCA model (once per resolution)
        if resolution_um not in self.pca_models:
            pca_path = resolution_dir / f"pca_model_{resolution_um}um.pkl"
            if not pca_path.exists():
                raise FileNotFoundError(f"PCA model not found: {pca_path}")
            with open(pca_path, 'rb') as f:
                self.pca_models[resolution_um] = pickle.load(f)
            self.logger.info(f"Loaded PCA model for {resolution_um}μm")

        # Load FAISS index
        index_path = resolution_dir / f"faiss_indices/faiss_r{radius}.index"
        if not index_path.exists():
            raise FileNotFoundError(f"FAISS index not found: {index_path}")
        self.indices[key] = faiss.read_index(str(index_path))
        self.logger.info(
            f"Loaded FAISS index for {resolution_um}μm, radius {radius}: "
            f"{self.indices[key].ntotal:,} vectors"
        )

        # Load embeddings for reference
        emb_path = resolution_dir / f"embeddings/embeddings_r{radius}_compressed.h5"
        if emb_path.exists():
            with h5py.File(emb_path, 'r') as f:
                self.embeddings[key] = {
                    'embeddings': f['embeddings'][:],
                    'patch_ids': [pid.decode() for pid in f['patch_ids'][:]]
                }

    def search_from_coordinates(
        self,
        sample_path: Path,
        x_center: float,
        y_center: float,
        radius_physical: float,
        resolution_um: int,
        k: int = 100,
        min_bins: Optional[int] = None,
        max_bins: Optional[int] = None,
        bin_size: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Search for similar regions using spatial coordinates.

        Parameters:
        -----------
        sample_path : Path
            Path to sample directory
        x_center, y_center : float
            Center coordinates of query region (in bin units)
        radius_physical : float
            Physical radius in micrometers
        resolution_um : int
            Resolution to search in
        k : int
            Number of results to return
        min_bins, max_bins : int, optional
            Filter results by number of bins
        bin_size : int, optional
            Bin size for zarr file (default: same as resolution)

        Returns:
        --------
        pd.DataFrame: Search results with columns:
            - patch_id, sample_id, similarity, n_bins
            - center_x, center_y, radius
        """
        if bin_size is None:
            bin_size = resolution_um

        # Convert physical radius to bin units
        radius_bins = int(radius_physical / resolution_um)

        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"SPATIAL PATTERN SEARCH ({resolution_um}μm)")
        self.logger.info(f"{'='*60}")
        self.logger.info(f"Query center: ({x_center}, {y_center})")
        self.logger.info(f"Query radius: {radius_physical}μm ({radius_bins} bins)")

        # Load sample data
        zarr_path = sample_path / "zarr" / f"bins_size_{bin_size}.zarr.zip"
        if not zarr_path.exists():
            raise FileNotFoundError(f"Zarr file not found: {zarr_path}")

        store = zarr.storage.ZipStore(str(zarr_path), mode='r')
        z = zarr.open_array(store, mode='r')

        # Get bins within radius
        n_genes, width, height = z.shape
        y_indices, x_indices = np.meshgrid(
            np.arange(height),
            np.arange(width),
            indexing='ij'
        )

        # Calculate distances from center
        distances = np.sqrt(
            (x_indices - x_center)**2 + (y_indices - y_center)**2
        )
        in_radius = distances <= radius_bins

        # Get selected bin coordinates
        selected_bins = np.column_stack([
            x_indices[in_radius],
            y_indices[in_radius]
        ])

        n_query_bins = len(selected_bins)
        self.logger.info(f"Query bins: {n_query_bins}")

        if n_query_bins == 0:
            self.logger.warning("No bins in query region")
            store.close()
            return pd.DataFrame()

        # Auto-select search radius
        search_radius = auto_select_radius(n_query_bins, resolution_um)
        self.logger.info(f"Auto-selected search radius: {search_radius} bins")

        # Load spatially variable genes
        variable_gene_mask = None
        haystack_file = sample_path / "haystack_results.csv"
        if haystack_file.exists():
            variable_gene_mask = load_variable_genes(haystack_file)

        # Compute query embedding
        self.logger.info("\n[1/4] Computing query embedding...")
        encoder = RadialShellEncoder(n_shells=self.n_shells)
        query_embedding = encoder.encode_from_zarr(
            z,
            selected_bins,
            bin_size,
            variable_gene_mask
        )

        store.close()

        # Search
        return self._search_with_embedding(
            query_embedding,
            resolution_um,
            search_radius,
            k,
            min_bins,
            max_bins
        )

    def search_from_selection(
        self,
        sample_path: Path,
        selection_bounds: Dict,
        resolution_um: int,
        k: int = 100,
        min_bins: Optional[int] = None,
        max_bins: Optional[int] = None,
        bin_size: Optional[int] = None
    ) -> pd.DataFrame:
        """
        Search for similar regions using selection bounds (box or lasso).

        Parameters:
        -----------
        sample_path : Path
            Path to sample directory
        selection_bounds : Dict
            Selection bounds with keys:
            - 'type': 'box' or 'lasso'
            - For box: 'xRange', 'yRange'
            - For lasso: 'lassoPoints' with 'x' and 'y' lists
        resolution_um : int
            Resolution to search in
        k : int
            Number of results
        min_bins, max_bins : int, optional
            Filter results by number of bins
        bin_size : int, optional
            Bin size for zarr file

        Returns:
        --------
        pd.DataFrame: Search results
        """
        if bin_size is None:
            bin_size = resolution_um

        # Load sample data
        zarr_path = sample_path / "zarr" / f"bins_size_{bin_size}.zarr.zip"
        if not zarr_path.exists():
            raise FileNotFoundError(f"Zarr file not found: {zarr_path}")

        store = zarr.storage.ZipStore(str(zarr_path), mode='r')
        z = zarr.open_array(store, mode='r')

        # Get all bin coordinates
        n_genes, width, height = z.shape
        y_indices, x_indices = np.meshgrid(
            np.arange(height),
            np.arange(width),
            indexing='ij'
        )
        x_coords = x_indices.flatten() * bin_size
        y_coords = y_indices.flatten() * bin_size

        # Find bins in selection
        selection_type = selection_bounds.get('type')

        if selection_type == 'box':
            x_range = selection_bounds['xRange']
            y_range = selection_bounds['yRange']
            mask = (
                (x_coords >= x_range[0]) & (x_coords <= x_range[1]) &
                (y_coords >= y_range[0]) & (y_coords <= y_range[1])
            )
        elif selection_type == 'lasso':
            # Use matplotlib path for lasso
            from matplotlib.path import Path as MplPath
            lasso_points = list(zip(
                selection_bounds['lassoPoints']['x'],
                selection_bounds['lassoPoints']['y']
            ))
            path = MplPath(lasso_points)
            points = np.column_stack((x_coords, y_coords))
            mask = path.contains_points(points)
        else:
            raise ValueError(f"Invalid selection type: {selection_type}")

        # Get selected bins (convert back to bin coordinates)
        selected_flat = np.where(mask)[0]
        selected_bins = np.column_stack([
            x_indices.flatten()[selected_flat],
            y_indices.flatten()[selected_flat]
        ])

        n_query_bins = len(selected_bins)
        self.logger.info(f"Query bins: {n_query_bins}")

        if n_query_bins == 0:
            self.logger.warning("No bins in selection")
            store.close()
            return pd.DataFrame()

        # Auto-select search radius
        search_radius = auto_select_radius(n_query_bins, resolution_um)
        self.logger.info(f"Auto-selected search radius: {search_radius} bins")

        # Load spatially variable genes
        variable_gene_mask = None
        haystack_file = sample_path / "haystack_results.csv"
        if haystack_file.exists():
            variable_gene_mask = load_variable_genes(haystack_file)

        # Compute query embedding
        self.logger.info("\n[1/4] Computing query embedding...")
        encoder = RadialShellEncoder(n_shells=self.n_shells)
        query_embedding = encoder.encode_from_zarr(
            z,
            selected_bins,
            bin_size,
            variable_gene_mask
        )

        store.close()

        # Search
        return self._search_with_embedding(
            query_embedding,
            resolution_um,
            search_radius,
            k,
            min_bins,
            max_bins
        )

    def _search_with_embedding(
        self,
        query_embedding: np.ndarray,
        resolution_um: int,
        radius: int,
        k: int,
        min_bins: Optional[int],
        max_bins: Optional[int]
    ) -> pd.DataFrame:
        """
        Perform search with pre-computed embedding.

        Parameters:
        -----------
        query_embedding : np.ndarray
            Query embedding vector
        resolution_um : int
            Resolution to search
        radius : int
            Patch radius to search
        k : int
            Number of results
        min_bins, max_bins : int, optional
            Filter by number of bins

        Returns:
        --------
        pd.DataFrame: Search results
        """
        # Load data for this resolution/radius
        self._load_resolution(resolution_um, radius)

        # Transform with PCA
        self.logger.info("[2/4] Transforming with PCA...")
        pca = self.pca_models[resolution_um]
        query_compressed = pca.transform(
            query_embedding.reshape(1, -1)
        ).astype('float32')
        faiss.normalize_L2(query_compressed)

        # Search
        self.logger.info(f"[3/4] Searching...")
        key = (resolution_um, radius)
        similarities, indices = self.indices[key].search(query_compressed, k * 2)

        # Retrieve metadata
        self.logger.info("[4/4] Retrieving results...")
        metadata = self.metadata[resolution_um]
        results = metadata[metadata['radius'] == radius].iloc[indices[0]].copy()
        results['similarity'] = similarities[0]

        # Filter by bin count
        if min_bins is not None:
            results = results[results['n_bins'] >= min_bins]
        if max_bins is not None:
            results = results[results['n_bins'] <= max_bins]

        # Sort by similarity (highest first)
        results = results.sort_values('similarity', ascending=False).head(k)

        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"RESULTS: {len(results)} patches")
        self.logger.info(f"{'='*60}\n")

        return results.reset_index(drop=True)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    print("Similarity Search Module")
    print("=" * 50)
    print("Use radial_shell_system.py for searching databases")
