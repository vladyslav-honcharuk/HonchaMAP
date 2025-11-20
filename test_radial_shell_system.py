#!/usr/bin/env python3
"""
Test and Validation Script for Radial Shell Encoding System.

This script creates synthetic data and validates the complete pipeline:
- Radial shell encoding
- PCA compression
- FAISS indexing
- Similarity search

Author: XeniumMundus Project
"""

import numpy as np
import tempfile
import shutil
from pathlib import Path
import zarr
import pandas as pd
import logging

from radial_shell_encoder import (
    RadialShellEncoder,
    PatchGenerator,
    auto_select_radius
)
from database_builder import MultiResolutionDatabaseBuilder
from similarity_search import MultiResolutionSpatialSearch


def create_synthetic_sample(
    output_dir: Path,
    sample_id: str,
    resolution_um: int = 40,
    n_genes: int = 100,
    width: int = 100,
    height: int = 100,
    pattern: str = 'gradient'
) -> Path:
    """
    Create a synthetic Xenium sample with spatial patterns.

    Parameters:
    -----------
    output_dir : Path
        Output directory for sample
    sample_id : str
        Sample identifier
    resolution_um : int
        Bin resolution
    n_genes : int
        Number of genes
    width, height : int
        Spatial dimensions in bins
    pattern : str
        Type of spatial pattern ('gradient', 'spots', 'random')

    Returns:
    --------
    Path: Path to created sample directory
    """
    sample_dir = output_dir / sample_id
    zarr_dir = sample_dir / "zarr"
    zarr_dir.mkdir(parents=True, exist_ok=True)

    # Create expression data with spatial pattern
    expr_data = np.zeros((n_genes, width, height), dtype=np.float32)

    # Create coordinate grids
    x_grid, y_grid = np.meshgrid(np.arange(width), np.arange(height), indexing='ij')

    if pattern == 'gradient':
        # Gradient from center
        center_x, center_y = width // 2, height // 2
        distances = np.sqrt((x_grid - center_x)**2 + (y_grid - center_y)**2)
        max_dist = distances.max()

        for g in range(n_genes):
            # Each gene has gradient with some noise
            gradient = 1 - (distances / max_dist)
            noise = np.random.normal(0, 0.1, (width, height))
            expr_data[g] = np.maximum(0, gradient + noise)

    elif pattern == 'spots':
        # Random spots of expression
        for g in range(n_genes):
            # Create 5-10 spots
            n_spots = np.random.randint(5, 11)
            for _ in range(n_spots):
                spot_x = np.random.randint(0, width)
                spot_y = np.random.randint(0, height)
                spot_radius = np.random.randint(5, 15)

                spot_mask = (
                    (x_grid - spot_x)**2 + (y_grid - spot_y)**2
                ) <= spot_radius**2

                expr_data[g][spot_mask] = np.random.uniform(0.5, 1.5)

    else:  # random
        expr_data = np.random.rand(n_genes, width, height).astype(np.float32)

    # Save as zarr
    zarr_file = zarr_dir / f"bins_size_{resolution_um}.zarr.zip"
    store = zarr.storage.ZipStore(str(zarr_file), mode='w')
    z = zarr.create(
        store=store,
        shape=(n_genes, width, height),
        chunks=(1, width, height),
        dtype='float32'
    )
    z[:] = expr_data
    store.close()

    # Create haystack results (mark some genes as spatially variable)
    n_variable = int(n_genes * 0.7)  # 70% spatially variable
    haystack_data = {
        'gene': [f'Gene_{i:03d}' for i in range(n_genes)],
        'KLD': np.random.uniform(0, 2, n_genes),
        'pval': np.random.uniform(0, 1, n_genes),
        'logpval_adj': np.random.uniform(-5, 0, n_genes)
    }
    # Make top genes more significant
    haystack_data['logpval_adj'][:n_variable] = np.random.uniform(-5, -2.1, n_variable)

    haystack_df = pd.DataFrame(haystack_data)
    haystack_df.to_csv(sample_dir / "haystack_results.csv", index=False)

    # Create genes.csv
    with open(sample_dir / "genes.csv", 'w') as f:
        for i in range(n_genes):
            f.write(f"Gene_{i:03d}\n")

    return sample_dir


def test_encoder():
    """Test radial shell encoder."""
    print("\n" + "=" * 70)
    print("TEST 1: Radial Shell Encoder")
    print("=" * 70)

    encoder = RadialShellEncoder(n_shells=5)

    # Create test data
    n_bins = 100
    n_genes = 50

    # Circular patch
    angles = np.random.uniform(0, 2*np.pi, n_bins)
    radii = np.random.uniform(0, 50, n_bins)
    bin_coords = np.column_stack([
        radii * np.cos(angles),
        radii * np.sin(angles)
    ])

    gene_expression = np.random.rand(n_bins, n_genes).astype(np.float32)

    # Encode
    encoding = encoder.encode_patch(bin_coords, gene_expression)

    print(f"Input: {n_bins} bins, {n_genes} genes")
    print(f"Output: {encoding.shape[0]} dimensions")
    print(f"Expected: {n_genes * encoder.n_shells} dimensions")

    assert encoding.shape[0] == n_genes * encoder.n_shells, "Encoding dimension mismatch"
    assert encoding.dtype == np.float32, "Encoding dtype should be float32"

    print("✓ Encoder test passed")
    return True


def test_patch_generator():
    """Test patch generator."""
    print("\n" + "=" * 70)
    print("TEST 2: Patch Generator")
    print("=" * 70)

    generator = PatchGenerator(resolution_um=40)

    # Create sample bin coordinates
    width, height = 100, 100
    y_grid, x_grid = np.meshgrid(np.arange(height), np.arange(width), indexing='ij')
    bin_coords = np.column_stack([x_grid.flatten(), y_grid.flatten()])

    # Generate patches
    patches = generator.generate_patches(bin_coords, sample_id="test_sample")

    print(f"Generated {len(patches)} patches")
    print(f"Radii used: {generator.radii}")

    # Check patches
    assert len(patches) > 0, "No patches generated"
    for patch in patches[:3]:
        print(f"\nPatch: {patch['patch_id']}")
        print(f"  Radius: {patch['radius']}")
        print(f"  Center: {patch['center']}")
        print(f"  Bins: {patch['n_bins']}")

    print("\n✓ Patch generator test passed")
    return True


def test_database_build():
    """Test database building."""
    print("\n" + "=" * 70)
    print("TEST 3: Database Building")
    print("=" * 70)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Create synthetic samples
        print("Creating synthetic samples...")
        samples = []
        for i in range(3):
            sample_dir = create_synthetic_sample(
                temp_path / "samples",
                f"Sample_{i:02d}",
                resolution_um=40,
                n_genes=50,
                width=50,
                height=50,
                pattern='gradient' if i % 2 == 0 else 'spots'
            )
            samples.append(sample_dir)
            print(f"  Created: {sample_dir.name}")

        # Build database
        print("\nBuilding database...")
        db_dir = temp_path / "database"

        builder = MultiResolutionDatabaseBuilder(
            output_dir=db_dir,
            n_shells=5,
            pca_dims=32,  # Small for test
            use_variable_genes=False,  # Disable for synthetic test data
            batch_size=100
        )

        samples_by_resolution = {40: samples}
        builder.build_database(samples_by_resolution)

        # Check outputs
        res_dir = db_dir / "40um"
        assert res_dir.exists(), "Resolution directory not created"
        assert (res_dir / "patches_metadata.parquet").exists(), "Metadata not created"
        assert (res_dir / "pca_model_40um.pkl").exists(), "PCA model not created"

        # Check FAISS indices
        index_dir = res_dir / "faiss_indices"
        assert index_dir.exists(), "FAISS index directory not created"
        indices = list(index_dir.glob("faiss_r*.index"))
        assert len(indices) > 0, "No FAISS indices created"

        print(f"\n✓ Database build test passed")
        print(f"  Created {len(indices)} FAISS indices")

        return True


def test_search():
    """Test similarity search."""
    print("\n" + "=" * 70)
    print("TEST 4: Similarity Search")
    print("=" * 70)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Create samples
        print("Creating test database...")
        samples = []
        for i in range(5):
            sample_dir = create_synthetic_sample(
                temp_path / "samples",
                f"Sample_{i:02d}",
                resolution_um=40,
                n_genes=50,
                width=60,
                height=60,
                pattern='gradient'
            )
            samples.append(sample_dir)

        # Build database
        db_dir = temp_path / "database"
        builder = MultiResolutionDatabaseBuilder(
            output_dir=db_dir,
            n_shells=5,
            pca_dims=32,
            use_variable_genes=False,  # Disable for synthetic test data
            batch_size=100
        )
        builder.build_database({40: samples})

        # Search
        print("\nPerforming search...")
        searcher = MultiResolutionSpatialSearch(db_dir)

        # Search from coordinates
        query_sample = samples[0]
        results = searcher.search_from_coordinates(
            sample_path=query_sample,
            x_center=30,
            y_center=30,
            radius_physical=400,  # 10 bins at 40um
            resolution_um=40,
            k=10
        )

        print(f"\nFound {len(results)} results")
        assert len(results) > 0, "No search results"

        print("\nTop 3 results:")
        for idx, row in results.head(3).iterrows():
            print(f"  {idx+1}. {row['sample_id']}")
            print(f"     Similarity: {row['similarity']:.4f}")
            print(f"     Center: ({row['center_x']:.1f}, {row['center_y']:.1f})")

        # Check that query sample is in top results (should match itself)
        top_samples = results.head(5)['sample_id'].values
        assert query_sample.name in top_samples, "Query sample not in top results"

        print("\n✓ Search test passed")
        return True


def test_auto_radius_selection():
    """Test automatic radius selection."""
    print("\n" + "=" * 70)
    print("TEST 5: Auto Radius Selection")
    print("=" * 70)

    test_cases = [
        (100, 40),   # 100 bins at 40um
        (500, 40),   # 500 bins at 40um
        (50, 10),    # 50 bins at 10um
        (200, 160),  # 200 bins at 160um
    ]

    for n_bins, resolution in test_cases:
        radius = auto_select_radius(n_bins, resolution)
        physical_size = radius * resolution
        print(f"  {n_bins} bins at {resolution}μm → radius {radius} bins ({physical_size}μm)")

        # Check radius is reasonable
        assert radius > 0, "Radius must be positive"

    print("\n✓ Auto radius selection test passed")
    return True


def run_all_tests():
    """Run all tests."""
    logging.basicConfig(level=logging.WARNING)  # Suppress debug logs during tests

    print("=" * 70)
    print("RADIAL SHELL ENCODING SYSTEM - VALIDATION TESTS")
    print("=" * 70)

    tests = [
        ("Encoder", test_encoder),
        ("Patch Generator", test_patch_generator),
        ("Database Building", test_database_build),
        ("Similarity Search", test_search),
        ("Auto Radius Selection", test_auto_radius_selection),
    ]

    passed = 0
    failed = 0

    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"\n✗ {test_name} test FAILED: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Passed: {passed}/{len(tests)}")
    print(f"Failed: {failed}/{len(tests)}")

    if failed == 0:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(run_all_tests())
