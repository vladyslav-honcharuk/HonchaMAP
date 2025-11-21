#!/usr/bin/env python3
"""
Simple Radial Pattern Search - All in one script

Does 3 things:
1. Loads samples and extracts circular patches
2. Encodes each patch (gene expression in 5 concentric rings)
3. Finds similar patches

Usage:
    python simple_search.py --data-dir /path/to/samples
"""

import numpy as np
import zarr
from pathlib import Path
from scipy.spatial.distance import cosine
import argparse

# Configuration
RADIUS = 5  # patch radius in bins
N_SHELLS = 5  # number of concentric rings


# ============================================================================
# STEP 1: Load sample data
# ============================================================================

def load_sample(sample_path):
    """Load genes and expression from one sample."""
    # Read gene names
    with open(sample_path / "genes.csv") as f:
        genes = [line.strip() for line in f if line.strip()]

    # Load spatial expression: (genes, x, y)
    zarr_file = sample_path / "zarr" / "bins_size_160.zarr.zip"
    store = zarr.storage.ZipStore(str(zarr_file), mode='r')
    expression = np.asarray(zarr.open_array(store, mode='r'))
    store.close()

    return genes, expression


# ============================================================================
# STEP 2: Extract and encode patches
# ============================================================================

def extract_patch(expression, center_x, center_y, radius):
    """Get a circular patch of expression data."""
    n_genes, width, height = expression.shape

    # Find bins in circle
    bins = []
    for x in range(max(0, center_x - radius), min(width, center_x + radius + 1)):
        for y in range(max(0, center_y - radius), min(height, center_y + radius + 1)):
            dist = np.sqrt((x - center_x)**2 + (y - center_y)**2)
            if dist <= radius:
                bins.append((x, y))

    if not bins:
        return None

    # Extract expression for these bins
    x_coords = np.array([b[0] for b in bins])
    y_coords = np.array([b[1] for b in bins])
    patch_data = expression[:, x_coords, y_coords]  # (genes, n_bins)

    return x_coords, y_coords, patch_data


def encode_radial(x_coords, y_coords, patch_data, n_shells=5):
    """
    Encode patch into radial features.

    For each gene:
        - Divide patch into n_shells concentric rings
        - Calculate average expression in each ring
        - Result: [gene0_ring0, gene0_ring1, ..., gene0_ring4, gene1_ring0, ...]
    """
    n_genes = patch_data.shape[0]

    # Calculate distances from center
    cx, cy = x_coords.mean(), y_coords.mean()
    dists = np.sqrt((x_coords - cx)**2 + (y_coords - cy)**2)
    max_dist = dists.max()

    if max_dist == 0:
        return patch_data.mean(axis=1).repeat(n_shells)

    # Define ring boundaries
    ring_edges = np.linspace(0, max_dist, n_shells + 1)

    # Calculate average in each ring for each gene
    features = []
    for gene_idx in range(n_genes):
        for ring_idx in range(n_shells):
            # Which bins are in this ring?
            in_ring = (dists >= ring_edges[ring_idx]) & (dists < ring_edges[ring_idx + 1])

            if in_ring.sum() > 0:
                features.append(patch_data[gene_idx, in_ring].mean())
            else:
                features.append(0.0)

    return np.array(features, dtype=np.float32)


# ============================================================================
# STEP 3: Build database
# ============================================================================

def build_database(data_dir, radius):
    """Generate patches from all samples and encode them."""
    samples = [d for d in data_dir.iterdir() if d.is_dir() and (d / "zarr").exists()]
    print(f"Found {len(samples)} samples\n")

    # Get all genes (union across samples)
    all_genes = set()
    for s in samples:
        with open(s / "genes.csv") as f:
            all_genes.update(line.strip() for line in f if line.strip())

    global_genes = sorted(all_genes)
    print(f"Total genes: {len(global_genes)}\n")

    # Build database
    database = []

    for sample_dir in samples:
        print(f"Processing {sample_dir.name}...")
        genes, expression = load_sample(sample_dir)

        # Map to global genes (fill missing with zeros)
        n_genes, width, height = expression.shape
        global_expr = np.zeros((len(global_genes), width, height), dtype=expression.dtype)

        for i, gene in enumerate(genes):
            if gene in global_genes:
                idx = global_genes.index(gene)
                global_expr[idx] = expression[i]

        # Generate patches every 10 bins
        count = 0
        for cx in range(radius, width - radius, 10):
            for cy in range(radius, height - radius, 10):
                result = extract_patch(global_expr, cx, cy, radius)
                if result is None:
                    continue

                x_coords, y_coords, patch_data = result
                encoding = encode_radial(x_coords, y_coords, patch_data, N_SHELLS)

                database.append({
                    'sample': sample_dir.name,
                    'x': cx,
                    'y': cy,
                    'encoding': encoding
                })
                count += 1

        print(f"  → {count} patches\n")

    print(f"Total database: {len(database)} patches")
    return database, global_genes


# ============================================================================
# STEP 4: Search
# ============================================================================

def search(query_encoding, database, top_k=10):
    """Find most similar patches."""
    similarities = []

    for entry in database:
        sim = 1 - cosine(query_encoding, entry['encoding'])
        similarities.append(sim)

    # Get top K
    top_indices = np.argsort(similarities)[-top_k:][::-1]

    results = []
    for idx in top_indices:
        results.append({
            'sample': database[idx]['sample'],
            'x': database[idx]['x'],
            'y': database[idx]['y'],
            'similarity': similarities[idx]
        })

    return results


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description='Simple radial pattern search')
    parser.add_argument('--data-dir', required=True, help='Directory with samples')
    args = parser.parse_args()

    data_dir = Path(args.data_dir)

    print("="*60)
    print("SIMPLE RADIAL PATTERN SEARCH")
    print("="*60)
    print()

    # Build database
    print("Building database...")
    print("-"*60)
    database, genes = build_database(data_dir, RADIUS)
    print()

    # Create a query from first sample
    print("Creating query...")
    print("-"*60)
    samples = [d for d in data_dir.iterdir() if d.is_dir() and (d / "zarr").exists()]
    query_sample = samples[0]

    sample_genes, expression = load_sample(query_sample)
    n_genes, width, height = expression.shape

    # Find a valid query location (center of tissue with some expression)
    query_x = width // 2
    query_y = height // 2

    # Make sure we're within bounds
    query_x = max(RADIUS + 5, min(query_x, width - RADIUS - 5))
    query_y = max(RADIUS + 5, min(query_y, height - RADIUS - 5))

    print(f"Query sample: {query_sample.name}")
    print(f"Query location: ({query_x}, {query_y})")
    print(f"Query radius: {RADIUS} bins")
    print(f"Tissue size: {width} x {height} bins")

    # Map to global genes
    global_expr = np.zeros((len(genes), expression.shape[1], expression.shape[2]), dtype=expression.dtype)
    for i, gene in enumerate(sample_genes):
        if gene in genes:
            idx = genes.index(gene)
            global_expr[idx] = expression[i]

    # Extract query patch
    result = extract_patch(global_expr, query_x, query_y, RADIUS)
    if result is None:
        print("Failed to create query!")
        return

    x_coords, y_coords, patch_data = result
    query_encoding = encode_radial(x_coords, y_coords, patch_data, N_SHELLS)

    print(f"Query encoding: {len(query_encoding)} features ({len(genes)} genes × {N_SHELLS} rings)")
    print()

    # Search
    print("Searching...")
    print("-"*60)
    results = search(query_encoding, database, top_k=10)

    print("\nTop 10 matches:")
    for i, r in enumerate(results):
        print(f"{i+1:2d}. {r['sample']:12s} ({r['x']:3d}, {r['y']:3d})  similarity: {r['similarity']:.3f}")

    print()
    print("="*60)
    print("DONE!")
    print("="*60)


if __name__ == '__main__':
    main()
