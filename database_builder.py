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
from typing import List, Dict, Optional, Tuple
import h5py
import pickle
import logging
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

        # Phase 1: Generate patches and raw embeddings
        self.logger.info("\n[Phase 1/3] Generating patches and raw embeddings...")
        raw_embeddings_by_radius, metadata = self._generate_raw_embeddings(
            sample_paths,
            resolution_um,
            bin_size
        )

        if not raw_embeddings_by_radius:
            self.logger.warning(f"No patches generated for {resolution_um}μm. Skipping.")
            return

        # Save metadata
        metadata_df = pd.DataFrame(metadata)
        metadata_path = resolution_dir / "patches_metadata.parquet"
        metadata_df.to_parquet(metadata_path, index=False)
        self.logger.info(f"Saved metadata: {len(metadata_df):,} patches")

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
        bin_size: int
    ) -> Tuple[Dict[int, Dict], List[Dict]]:
        """
        Generate patches and compute raw embeddings for all samples.

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

                # Load spatially variable genes if requested
                variable_gene_mask = None
                if self.use_variable_genes:
                    haystack_file = sample_path / "haystack_results.csv"
                    variable_gene_mask = load_variable_genes(haystack_file)

                # Generate patches
                patches = patch_generator.generate_patches(
                    bin_coords,
                    sample_id
                )

                # Compute embeddings for each patch
                for patch in patches:
                    # Get bin indices for this patch
                    selected_bins = bin_coords[patch['bin_indices']]

                    # Encode patch
                    embedding = self.encoder.encode_from_zarr(
                        z,
                        selected_bins,
                        bin_size,
                        variable_gene_mask
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
