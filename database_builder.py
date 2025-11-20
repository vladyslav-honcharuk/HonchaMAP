#!/usr/bin/env python3
"""
Database Builder - Build multi-resolution radial shell encoding database.

This module handles the complete pipeline:
1. Generate patches from samples
2. Compute radial shell encodings
3. Fit PCA compression
4. Build FAISS indices
5. Save metadata

Author: XeniumMundus Project
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Set
import h5py
import pickle
import logging
import json
from sklearn.decomposition import IncrementalPCA
import faiss
from tqdm import tqdm
import zarr

from radial_shell_encoder import (
    RadialShellEncoder,
    PatchGenerator,
    load_variable_genes
)


class MultiResolutionDatabaseBuilder:
    """
    Builds searchable database across multiple resolutions with PCA compression.
    """

    def __init__(
        self,
        output_dir: Path,
        n_shells: int = 5,
        pca_dims: int = 256,
        use_variable_genes: bool = True,
        batch_size: int = 10000
    ):
        """
        Initialize database builder.

        Parameters:
        -----------
        output_dir : Path
            Output directory for database
        n_shells : int
            Number of concentric shells
        pca_dims : int
            PCA compression dimensions
        use_variable_genes : bool
            Whether to use only spatially variable genes
        batch_size : int
            Batch size for PCA fitting
        """
        self.output_dir = Path(output_dir)
        self.n_shells = n_shells
        self.pca_dims = pca_dims
        self.use_variable_genes = use_variable_genes
        self.batch_size = batch_size

        self.encoder = RadialShellEncoder(n_shells=n_shells)
        self.logger = logging.getLogger(__name__)

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def _discover_global_genes(
        self,
        sample_paths: List[Path]
    ) -> Tuple[List[str], Optional[np.ndarray]]:
        """
        Discover all unique genes across samples and determine global gene set.

        Returns:
        --------
        Tuple of:
        - List of gene names in sorted order (global gene coordinate system)
        - Optional boolean mask indicating which genes are spatially variable
        """
        self.logger.info("\n[Gene Discovery] Building global gene coordinate system...")

        # Step 1: Collect all genes from all samples
        all_genes_per_sample = []
        variable_genes_per_sample = []

        for sample_path in tqdm(sample_paths, desc="Reading gene lists"):
            # Read genes.csv
            genes_file = sample_path / "genes.csv"
            if genes_file.exists():
                with open(genes_file, 'r') as f:
                    genes = [line.strip() for line in f if line.strip()]
                all_genes_per_sample.append(set(genes))
            else:
                self.logger.warning(f"No genes.csv found in {sample_path.name}, skipping")
                continue

            # If using variable genes, read haystack results
            if self.use_variable_genes:
                haystack_file = sample_path / "haystack_results.csv"
                if haystack_file.exists():
                    try:
                        df = pd.read_csv(haystack_file, index_col=0)
                        if 'logpval_adj' in df.columns:
                            # Get variable genes (default threshold -2.0 = p_adj < 0.01)
                            var_genes = set(df[df['logpval_adj'] <= -2.0].index.tolist())
                            variable_genes_per_sample.append(var_genes)
                        else:
                            self.logger.warning(f"No logpval_adj in haystack for {sample_path.name}")
                    except Exception as e:
                        self.logger.warning(f"Failed to read haystack for {sample_path.name}: {e}")

        if not all_genes_per_sample:
            raise ValueError("No genes.csv files found in any sample")

        # Step 2: Build union of all genes
        all_genes_union = set.union(*all_genes_per_sample)
        self.logger.info(f"Total unique genes across all samples: {len(all_genes_union)}")

        # Step 3: Determine global gene set
        global_genes = sorted(list(all_genes_union))  # Sort for consistent ordering
        global_gene_mask = None

        if self.use_variable_genes and variable_genes_per_sample:
            # Use intersection of variable genes (genes that are variable in ALL samples)
            # This is conservative but ensures consistency
            variable_genes_intersection = set.intersection(*variable_genes_per_sample)

            # If intersection is too small, use union (genes variable in ANY sample)
            if len(variable_genes_intersection) < 50:
                self.logger.warning(
                    f"Variable gene intersection only {len(variable_genes_intersection)} genes. "
                    f"Using union instead."
                )
                variable_genes_union = set.union(*variable_genes_per_sample)
                variable_genes_global = variable_genes_union
            else:
                variable_genes_global = variable_genes_intersection

            # Create boolean mask
            global_gene_mask = np.array([g in variable_genes_global for g in global_genes])
            n_variable = global_gene_mask.sum()
            self.logger.info(
                f"Using {n_variable}/{len(global_genes)} "
                f"({100*n_variable/len(global_genes):.1f}%) spatially variable genes"
            )
        else:
            self.logger.info(f"Using all {len(global_genes)} genes")

        return global_genes, global_gene_mask

    def _create_gene_mapping(
        self,
        sample_genes: List[str],
        global_genes: List[str]
    ) -> np.ndarray:
        """
        Create mapping from sample gene indices to global gene indices.

        Returns:
        --------
        Array of shape (n_sample_genes,) where each element is the index
        in global_genes, or -1 if gene not in global set.
        """
        global_gene_to_idx = {gene: idx for idx, gene in enumerate(global_genes)}
        mapping = np.array([global_gene_to_idx.get(g, -1) for g in sample_genes])
        return mapping

    def _encode_patch_global(
        self,
        global_z: np.ndarray,
        selected_bins: np.ndarray,
        bin_size: int,
        global_gene_mask: Optional[np.ndarray]
    ) -> np.ndarray:
        """
        Encode a patch using global gene coordinate system.

        Parameters:
        -----------
        global_z : np.ndarray
            Expression array in global coordinate system, shape (n_global_genes, width, height)
        selected_bins : np.ndarray
            Bin indices for this patch, shape (n_bins, 2)
        bin_size : int
            Bin size in micrometers
        global_gene_mask : Optional[np.ndarray]
            Boolean mask for spatially variable genes in global space

        Returns:
        --------
        np.ndarray: Encoded features
        """
        # Extract expression values for selected bins
        x_bins = selected_bins[:, 0]
        y_bins = selected_bins[:, 1]

        # Get expression matrix: (n_global_genes, n_bins) -> transpose to (n_bins, n_global_genes)
        gene_expression = global_z[:, x_bins, y_bins].T

        # Create coordinates from bin indices
        bin_coords = selected_bins.astype(float) * bin_size

        # Encode using radial shell encoder
        return self.encoder.encode_patch(bin_coords, gene_expression, global_gene_mask)

    def build_database(
        self,
        samples_by_resolution: Dict[int, List[Path]],
        bin_size_map: Optional[Dict[int, int]] = None
    ):
        """
        Build complete multi-resolution database.

        Parameters:
        -----------
        samples_by_resolution : Dict[int, List[Path]]
            Dictionary mapping resolution (μm) to list of sample directories
            Example: {10: [sample1_dir, sample2_dir], 20: [...], ...}
        bin_size_map : Dict[int, int], optional
            Custom mapping from resolution to bin size for zarr files
            Default uses resolution as bin size
        """
        if bin_size_map is None:
            bin_size_map = {res: res for res in samples_by_resolution.keys()}

        for resolution_um, sample_paths in samples_by_resolution.items():
            self.logger.info("=" * 70)
            self.logger.info(f"Processing {resolution_um}μm resolution")
            self.logger.info("=" * 70)

            bin_size = bin_size_map.get(resolution_um, resolution_um)

            self._build_resolution_database(
                resolution_um,
                sample_paths,
                bin_size
            )

    def _build_resolution_database(
        self,
        resolution_um: int,
        sample_paths: List[Path],
        bin_size: int
    ):
        """Build database for one resolution."""
        resolution_dir = self.output_dir / f"{resolution_um}um"
        resolution_dir.mkdir(parents=True, exist_ok=True)

        # Phase 0: Discover global gene coordinate system
        self.logger.info("\n[Phase 0/3] Discovering global gene coordinate system...")
        global_genes, global_gene_mask = self._discover_global_genes(sample_paths)

        # Phase 1: Generate patches and raw embeddings
        self.logger.info("\n[Phase 1/3] Generating patches and raw embeddings...")
        raw_embeddings_by_radius, metadata = self._generate_raw_embeddings(
            sample_paths,
            resolution_um,
            bin_size,
            global_genes,
            global_gene_mask
        )

        if not raw_embeddings_by_radius:
            self.logger.warning(f"No patches generated for {resolution_um}μm. Skipping.")
            return

        # Save metadata
        metadata_df = pd.DataFrame(metadata)
        metadata_path = resolution_dir / "patches_metadata.parquet"
        metadata_df.to_parquet(metadata_path, index=False)
        self.logger.info(f"Saved metadata: {len(metadata_df):,} patches")

        # Save global gene list and configuration
        config = {
            'global_genes': global_genes,
            'use_variable_genes': self.use_variable_genes,
            'n_shells': self.n_shells,
            'pca_dims': self.pca_dims,
            'resolution_um': resolution_um
        }
        if global_gene_mask is not None:
            config['variable_gene_indices'] = np.where(global_gene_mask)[0].tolist()

        config_path = resolution_dir / "config.json"
        with open(config_path, 'w') as f:
            json.dump(config, f, indent=2)
        self.logger.info(f"Saved configuration with {len(global_genes)} global genes")

        # Phase 2: Fit PCA
        self.logger.info("\n[Phase 2/3] Fitting PCA...")
        pca = self._fit_pca(raw_embeddings_by_radius)

        # Save PCA model
        pca_path = resolution_dir / f"pca_model_{resolution_um}um.pkl"
        with open(pca_path, 'wb') as f:
            pickle.dump(pca, f)
        self.logger.info(f"Saved PCA model: {pca_path}")

        # Phase 3: Compress and index
        self.logger.info("\n[Phase 3/3] Compressing embeddings and building FAISS indices...")
        for radius in raw_embeddings_by_radius.keys():
            self._compress_and_index(
                raw_embeddings_by_radius[radius],
                pca,
                resolution_dir,
                radius
            )

        self.logger.info(f"\n✓ Completed {resolution_um}μm resolution")

    def _generate_raw_embeddings(
        self,
        sample_paths: List[Path],
        resolution_um: int,
        bin_size: int,
        global_genes: List[str],
        global_gene_mask: Optional[np.ndarray]
    ) -> Tuple[Dict[int, Dict], List[Dict]]:
        """
        Generate patches and compute raw embeddings for all samples.

        Parameters:
        -----------
        global_genes : List[str]
            Ordered list of all genes in global coordinate system
        global_gene_mask : Optional[np.ndarray]
            Boolean mask for spatially variable genes in global space

        Returns:
        --------
        Tuple of:
        - Dict mapping radius to {'embeddings': array, 'patch_ids': array}
        - List of metadata dictionaries
        """
        patch_generator = PatchGenerator(resolution_um)
        embeddings_by_radius = {r: [] for r in patch_generator.radii}
        patch_ids_by_radius = {r: [] for r in patch_generator.radii}
        all_metadata = []

        for sample_idx, sample_path in enumerate(tqdm(sample_paths, desc="Processing samples")):
            try:
                # Load sample data
                sample_id = sample_path.name
                zarr_path = sample_path / "zarr" / f"bins_size_{bin_size}.zarr.zip"

                if not zarr_path.exists():
                    self.logger.warning(f"Zarr file not found: {zarr_path}")
                    continue

                # Open zarr array
                store = zarr.storage.ZipStore(str(zarr_path), mode='r')
                z = zarr.open_array(store, mode='r')

                # Get bin coordinates
                n_genes, width, height = z.shape
                y_indices, x_indices = np.meshgrid(
                    np.arange(height),
                    np.arange(width),
                    indexing='ij'
                )
                all_bin_coords = np.column_stack([
                    x_indices.flatten(),
                    y_indices.flatten()
                ])

                # Filter to non-empty bins
                # Sample middle gene to find non-empty bins
                middle_gene_idx = n_genes // 2
                sample_expr = z[middle_gene_idx, :, :].T.flatten()
                non_empty_mask = sample_expr > 0
                bin_coords = all_bin_coords[non_empty_mask]

                if len(bin_coords) == 0:
                    self.logger.warning(f"No non-empty bins in {sample_id}")
                    store.close()
                    continue

                # Load sample genes and create mapping to global genes
                genes_file = sample_path / "genes.csv"
                if not genes_file.exists():
                    self.logger.warning(f"No genes.csv for {sample_id}, skipping")
                    store.close()
                    continue

                with open(genes_file, 'r') as f:
                    sample_genes = [line.strip() for line in f if line.strip()]

                # Create mapping from sample genes to global genes
                gene_mapping = self._create_gene_mapping(sample_genes, global_genes)

                # Convert zarr to numpy and remap to global coordinate system
                # z shape: (n_sample_genes, width, height)
                # We need to create global_z: (n_global_genes, width, height)
                zarr_np = np.asarray(z)
                n_global_genes = len(global_genes)
                global_z = np.zeros((n_global_genes, width, height), dtype=zarr_np.dtype)

                # Map sample genes to global positions
                for sample_idx, global_idx in enumerate(gene_mapping):
                    if global_idx >= 0:  # Gene exists in global set
                        global_z[global_idx] = zarr_np[sample_idx]

                # Generate patches
                patches = patch_generator.generate_patches(
                    bin_coords,
                    sample_id
                )

                # Compute embeddings for each patch
                for patch in patches:
                    # Get bin indices for this patch
                    selected_bins = bin_coords[patch['bin_indices']]

                    # Encode patch using global gene space
                    embedding = self._encode_patch_global(
                        global_z,
                        selected_bins,
                        bin_size,
                        global_gene_mask
                    )

                    radius = patch['radius']
                    embeddings_by_radius[radius].append(embedding)
                    patch_ids_by_radius[radius].append(patch['patch_id'])

                    all_metadata.append({
                        'patch_id': patch['patch_id'],
                        'sample_id': patch['sample_id'],
                        'resolution_um': resolution_um,
                        'radius': patch['radius'],
                        'center_x': float(patch['center'][0]),
                        'center_y': float(patch['center'][1]),
                        'n_bins': patch['n_bins']
                    })

                store.close()

            except Exception as e:
                self.logger.error(f"Failed to process {sample_path}: {e}")
                continue

        # Convert to arrays
        for radius in patch_generator.radii:
            if len(embeddings_by_radius[radius]) == 0:
                self.logger.warning(f"No embeddings for radius {radius}")
                del embeddings_by_radius[radius]
            else:
                embeddings_by_radius[radius] = {
                    'embeddings': np.vstack(embeddings_by_radius[radius]),
                    'patch_ids': np.array(patch_ids_by_radius[radius])
                }
                self.logger.info(
                    f"  Radius {radius}: {len(embeddings_by_radius[radius]['embeddings']):,} patches"
                )

        return embeddings_by_radius, all_metadata

    def _fit_pca(self, embeddings_by_radius: Dict[int, Dict]) -> IncrementalPCA:
        """
        Fit PCA on all embeddings for this resolution.

        Parameters:
        -----------
        embeddings_by_radius : Dict
            Dictionary mapping radius to embeddings

        Returns:
        --------
        IncrementalPCA: Fitted PCA transformer
        """
        # Collect all embeddings
        all_embeddings = np.vstack([
            data['embeddings']
            for data in embeddings_by_radius.values()
        ])

        self.logger.info(f"  Fitting PCA on {len(all_embeddings):,} patches")
        self.logger.info(f"  Original dimensions: {all_embeddings.shape[1]}")

        # Fit PCA in batches
        pca = IncrementalPCA(n_components=self.pca_dims, batch_size=self.batch_size)

        for i in tqdm(range(0, len(all_embeddings), self.batch_size), desc="Fitting PCA"):
            batch = all_embeddings[i:i + self.batch_size]
            pca.partial_fit(batch)

        variance = pca.explained_variance_ratio_.sum()
        self.logger.info(f"  Compressed dimensions: {self.pca_dims}")
        self.logger.info(f"  Variance explained: {variance:.1%}")

        return pca

    def _compress_and_index(
        self,
        embedding_data: Dict,
        pca: IncrementalPCA,
        resolution_dir: Path,
        radius: int
    ):
        """
        Compress embeddings with PCA and build FAISS index.

        Parameters:
        -----------
        embedding_data : Dict
            Dictionary with 'embeddings' and 'patch_ids' keys
        pca : IncrementalPCA
            Fitted PCA transformer
        resolution_dir : Path
            Resolution output directory
        radius : int
            Patch radius
        """
        embeddings = embedding_data['embeddings']
        patch_ids = embedding_data['patch_ids']

        self.logger.info(f"\n  Radius {radius}:")
        self.logger.info(f"    Original: {embeddings.shape}")

        # Transform with PCA
        embeddings_compressed = pca.transform(embeddings).astype('float32')
        self.logger.info(f"    Compressed: {embeddings_compressed.shape}")

        # Save compressed embeddings
        emb_dir = resolution_dir / "embeddings"
        emb_dir.mkdir(exist_ok=True)
        emb_path = emb_dir / f"embeddings_r{radius}_compressed.h5"

        with h5py.File(emb_path, 'w') as f:
            f.create_dataset('embeddings', data=embeddings_compressed, compression='gzip')
            f.create_dataset('patch_ids', data=patch_ids.astype('S'))

        self.logger.info(f"    Saved embeddings: {emb_path}")

        # Build FAISS index
        # Normalize for inner product (cosine similarity)
        faiss.normalize_L2(embeddings_compressed)

        # Create index
        index = faiss.IndexFlatIP(embeddings_compressed.shape[1])
        index.add(embeddings_compressed)

        # Save index
        index_dir = resolution_dir / "faiss_indices"
        index_dir.mkdir(exist_ok=True)
        index_path = index_dir / f"faiss_r{radius}.index"
        faiss.write_index(index, str(index_path))

        self.logger.info(f"    FAISS index: {index.ntotal:,} vectors")


def discover_samples(
    base_dir: Path,
    min_resolutions: int = 1
) -> Dict[int, List[Path]]:
    """
    Discover Xenium samples organized by resolution.

    Parameters:
    -----------
    base_dir : Path
        Base directory containing sample subdirectories
    min_resolutions : int
        Minimum number of resolutions required per sample

    Returns:
    --------
    Dict[int, List[Path]]: Mapping from resolution to list of sample paths
    """
    logger = logging.getLogger(__name__)

    # Find all sample directories (those with zarr/ subdirectory)
    sample_dirs = []
    for item in base_dir.iterdir():
        if item.is_dir() and (item / "zarr").exists():
            sample_dirs.append(item)

    logger.info(f"Found {len(sample_dirs)} potential samples in {base_dir}")

    # Organize by resolution
    samples_by_resolution = {res: [] for res in [10, 20, 40, 80, 160]}

    for sample_dir in sample_dirs:
        zarr_dir = sample_dir / "zarr"
        available_resolutions = []

        for res in [10, 20, 40, 80, 160]:
            zarr_file = zarr_dir / f"bins_size_{res}.zarr.zip"
            if zarr_file.exists():
                samples_by_resolution[res].append(sample_dir)
                available_resolutions.append(res)

        if len(available_resolutions) >= min_resolutions:
            logger.debug(
                f"  {sample_dir.name}: {len(available_resolutions)} resolutions "
                f"{available_resolutions}"
            )

    # Remove empty resolutions
    samples_by_resolution = {
        res: samples
        for res, samples in samples_by_resolution.items()
        if len(samples) > 0
    }

    # Summary
    for res, samples in samples_by_resolution.items():
        logger.info(f"  {res}μm: {len(samples)} samples")

    return samples_by_resolution


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    print("Database Builder Module")
    print("=" * 50)
    print("Use radial_shell_system.py for building databases")
