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
import json
import zarr

from .radial_shell_encoder import (
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
        self.configs = {}  # Store config.json per resolution

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

        # Load configuration (once per resolution) - MUST load config before PCA
        if resolution_um not in self.configs:
            config_path = resolution_dir / "config.json"
            if config_path.exists():
                with open(config_path, 'r') as f:
                    self.configs[resolution_um] = json.load(f)
                n_genes = len(self.configs[resolution_um].get('global_genes', []))
                use_pca = self.configs[resolution_um].get('use_pca', True)
                self.logger.info(f"Loaded config for {resolution_um}μm: {n_genes} global genes, use_pca={use_pca}")
            else:
                self.logger.warning(f"No config.json found for {resolution_um}μm - using legacy mode")
                self.configs[resolution_um] = None

        # Load PCA model (once per resolution) - only if config says to use PCA
        if resolution_um not in self.pca_models:
            config = self.configs.get(resolution_um, {})
            use_pca = config.get('use_pca', True) if config else True  # Default to True for legacy

            if use_pca:
                pca_path = resolution_dir / f"pca_model_{resolution_um}um.pkl"
                if pca_path.exists():
                    with open(pca_path, 'rb') as f:
                        self.pca_models[resolution_um] = pickle.load(f)
                    self.logger.info(f"Loaded PCA model for {resolution_um}μm")
                else:
                    raise FileNotFoundError(f"PCA model not found: {pca_path}")
            else:
                self.logger.info(f"Skipping PCA model load for {resolution_um}μm (use_pca=False)")

        # Load FAISS index
        index_path = resolution_dir / f"faiss_indices/faiss_r{radius}.index"
        if not index_path.exists():
            raise FileNotFoundError(f"FAISS index not found: {index_path}")
        self.indices[key] = faiss.read_index(str(index_path))
        self.logger.info(
            f"Loaded FAISS index for {resolution_um}μm, radius {radius}: "
            f"{self.indices[key].ntotal:,} vectors"
        )

        # Load embeddings for reference (check both compressed and full)
        emb_path_compressed = resolution_dir / f"embeddings/embeddings_r{radius}_compressed.h5"
        emb_path_full = resolution_dir / f"embeddings/embeddings_r{radius}_full.h5"

        if emb_path_compressed.exists():
            with h5py.File(emb_path_compressed, 'r') as f:
                self.embeddings[key] = {
                    'embeddings': f['embeddings'][:],
                    'patch_ids': [pid.decode() for pid in f['patch_ids'][:]]
                }
        elif emb_path_full.exists():
            with h5py.File(emb_path_full, 'r') as f:
                self.embeddings[key] = {
                    'embeddings': f['embeddings'][:],
                    'patch_ids': [pid.decode() for pid in f['patch_ids'][:]]
                }
        else:
            self.logger.warning(f"No embeddings file found for r{radius}")

    def _encode_query_global(
        self,
        sample_path: Path,
        z: zarr.Array,
        selected_bins: np.ndarray,
        bin_size: int,
        resolution_um: int
    ) -> np.ndarray:
        """
        Encode query patch using global gene coordinate system.

        Parameters:
        -----------
        sample_path : Path
            Path to sample directory
        z : zarr.Array
            Zarr expression array for sample
        selected_bins : np.ndarray
            Selected bin indices, shape (n_bins, 2)
        bin_size : int
            Bin size in micrometers
        resolution_um : int
            Resolution to use (for loading config)

        Returns:
        --------
        np.ndarray: Encoded query embedding
        """
        # Ensure config is loaded for this resolution
        if resolution_um not in self.configs:
            # Trigger loading by calling _load_resolution with any radius
            # (config is loaded per resolution, not per radius)
            available_radii = PatchGenerator(resolution_um).radii
            if available_radii:
                self._load_resolution(resolution_um, available_radii[0])

        config = self.configs.get(resolution_um)

        if config is None:
            # Legacy mode - no global genes, use old method
            self.logger.warning("Using legacy encoding (no global gene system)")
            encoder = RadialShellEncoder(n_shells=self.n_shells)
            return encoder.encode_from_zarr(z, selected_bins, bin_size, None)

        # Get global genes from config
        global_genes = config['global_genes']
        n_global_genes = len(global_genes)

        # Get global gene mask if using variable genes
        global_gene_mask = None
        if 'variable_gene_indices' in config:
            global_gene_mask = np.zeros(n_global_genes, dtype=bool)
            global_gene_mask[config['variable_gene_indices']] = True

        # Load sample genes
        genes_file = sample_path / "genes.csv"
        if not genes_file.exists():
            raise FileNotFoundError(f"genes.csv not found in {sample_path}")

        with open(genes_file, 'r') as f:
            sample_genes = [line.strip() for line in f if line.strip()]

        # Create mapping from sample genes to global genes
        global_gene_to_idx = {gene: idx for idx, gene in enumerate(global_genes)}
        gene_mapping = np.array([global_gene_to_idx.get(g, -1) for g in sample_genes])

        # Convert zarr to numpy and remap to global coordinate system
        zarr_np = np.asarray(z)
        n_genes, width, height = zarr_np.shape
        global_z = np.zeros((n_global_genes, width, height), dtype=zarr_np.dtype)

        # Map sample genes to global positions
        for sample_idx, global_idx in enumerate(gene_mapping):
            if global_idx >= 0:  # Gene exists in global set
                global_z[global_idx] = zarr_np[sample_idx]

        # Extract expression for selected bins
        x_bins = selected_bins[:, 0]
        y_bins = selected_bins[:, 1]
        gene_expression = global_z[:, x_bins, y_bins].T  # (n_bins, n_global_genes)

        # Create coordinates from bin indices
        bin_coords = selected_bins.astype(float) * bin_size

        # Encode using radial shell encoder
        encoder = RadialShellEncoder(n_shells=self.n_shells)
        return encoder.encode_patch(bin_coords, gene_expression, global_gene_mask)

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
        bin_size: Optional[int] = None,
        gene_id: Optional[Union[int, str]] = None
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
        gene_id : int or str, optional
            Gene index or name to search for. If provided, searches only
            for patterns of this specific gene. If None, searches across all genes.

        Returns:
        --------
        pd.DataFrame: Search results with columns:
            - patch_id, sample_id, similarity, n_bins
            - center_x, center_y, radius
        """
        if bin_size is None:
            bin_size = resolution_um

        # Convert physical coordinates to bin units
        x_center_bins = x_center / bin_size
        y_center_bins = y_center / bin_size
        radius_bins = int(radius_physical / resolution_um)

        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"SPATIAL PATTERN SEARCH ({resolution_um}μm)")
        self.logger.info(f"{'='*60}")
        self.logger.info(f"Query center (physical): ({x_center:.1f}, {y_center:.1f}) μm")
        self.logger.info(f"Query center (bins): ({x_center_bins:.1f}, {y_center_bins:.1f})")
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

        # Calculate distances from center (using bin coordinates)
        distances = np.sqrt(
            (x_indices - x_center_bins)**2 + (y_indices - y_center_bins)**2
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

        # Compute query embedding using global gene coordinate system
        self.logger.info("\n[1/4] Computing query embedding...")
        query_embedding = self._encode_query_global(
            sample_path,
            z,
            selected_bins,
            bin_size,
            resolution_um
        )

        store.close()

        # Search
        return self._search_with_embedding(
            query_embedding,
            resolution_um,
            search_radius,
            k,
            min_bins,
            max_bins,
            gene_id
        )

    def search_from_selection(
        self,
        sample_path: Path,
        selection_bounds: Dict,
        resolution_um: int,
        k: int = 100,
        min_bins: Optional[int] = None,
        max_bins: Optional[int] = None,
        bin_size: Optional[int] = None,
        gene_id: Optional[Union[int, str]] = None
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

        # Compute query embedding using global gene coordinate system
        self.logger.info("\n[1/4] Computing query embedding...")
        query_embedding = self._encode_query_global(
            sample_path,
            z,
            selected_bins,
            bin_size,
            resolution_um
        )

        store.close()

        # Search
        return self._search_with_embedding(
            query_embedding,
            resolution_um,
            search_radius,
            k,
            min_bins,
            max_bins,
            gene_id
        )

    def _search_with_embedding(
        self,
        query_embedding: np.ndarray,
        resolution_um: int,
        radius: int,
        k: int,
        min_bins: Optional[int],
        max_bins: Optional[int],
        gene_id: Optional[Union[int, str]] = None
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

        # Get configuration for this resolution
        config = self.configs.get(resolution_um, {})
        use_pca = config.get('use_pca', True)  # Default to True for backward compatibility

        # Handle single-gene search
        gene_idx = None
        if gene_id is not None:
            self.logger.info(f"[2/4] Extracting single gene embedding (gene_id={gene_id})...")

            # Get gene index
            if isinstance(gene_id, str):
                # Gene name provided - find its index
                global_genes = config.get('global_genes', [])
                if gene_id not in global_genes:
                    raise ValueError(f"Gene '{gene_id}' not found in global gene list")
                gene_idx = global_genes.index(gene_id)
                self.logger.info(f"Gene '{gene_id}' -> index {gene_idx}")
            else:
                # Gene index provided
                gene_idx = gene_id

            # Extract only the shells for this gene
            # Embedding format: [gene0_shell0, gene0_shell1, ..., gene0_shell4, gene1_shell0, ...]
            # For n_shells=5, gene i occupies indices [i*5:(i+1)*5]
            n_shells = config.get('n_shells', 5)
            start_idx = gene_idx * n_shells
            end_idx = start_idx + n_shells
            query_gene_embedding = query_embedding[start_idx:end_idx]

            self.logger.info(f"Gene index: {gene_idx}, shells extracted: {start_idx}:{end_idx}")
        else:
            # Use full embedding
            query_gene_embedding = query_embedding

        # Transform with PCA if database uses it
        if use_pca and resolution_um in self.pca_models:
            self.logger.info("[2.5/4] Transforming with PCA...")
            pca = self.pca_models[resolution_um]
            query_compressed = pca.transform(
                query_gene_embedding.reshape(1, -1)
            ).astype('float32')
            faiss.normalize_L2(query_compressed)
        else:
            # No PCA - use raw embedding
            if gene_id is not None:
                self.logger.info("[2.5/4] Using single-gene embedding (no PCA)...")
            else:
                self.logger.info("[2.5/4] Using full embedding (no PCA)...")
            query_compressed = query_gene_embedding.reshape(1, -1).astype('float32')
            faiss.normalize_L2(query_compressed)

        # Search
        self.logger.info(f"[3/4] Searching...")
        key = (resolution_um, radius)

        # For single-gene search, need to extract gene embeddings from stored data
        if gene_id is not None and key in self.embeddings:
            self.logger.info(f"Building temporary single-gene index for gene {gene_idx}...")

            # Get stored embeddings
            stored_embeddings = self.embeddings[key]['embeddings']  # shape: (n_patches, n_genes * n_shells)

            # Extract gene-specific embeddings
            n_shells = config.get('n_shells', 5)
            start_idx = gene_idx * n_shells
            end_idx = start_idx + n_shells
            gene_embeddings = stored_embeddings[:, start_idx:end_idx].astype('float32')

            # Normalize
            faiss.normalize_L2(gene_embeddings)

            # Create temporary index
            temp_index = faiss.IndexFlatIP(gene_embeddings.shape[1])
            temp_index.add(gene_embeddings)

            # Search using temporary index
            similarities, indices = temp_index.search(query_compressed, k * 2)
            self.logger.info(f"Searched {temp_index.ntotal} patches for gene {gene_id}")
        else:
            # Normal search using pre-built index
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
