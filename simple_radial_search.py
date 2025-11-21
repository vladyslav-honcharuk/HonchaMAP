#!/usr/bin/env python3
"""
Simple Radial Pattern Search - Find similar tissue regions based on gene expression

What it does:
1. Takes a circular patch of tissue
2. Divides it into concentric rings (like tree rings)
3. Calculates average gene expression in each ring
4. Finds other patches with similar patterns

Usage:
    python simple_radial_search.py --data-dir /path/to/samples --query-sample GSM7780155
"""

import numpy as np
import pandas as pd
import zarr
import json
from pathlib import Path
from scipy.spatial.distance import cosine
import matplotlib.pyplot as plt


# ============================================================================
# STEP 1: Load spatial gene expression data
# ============================================================================

def load_sample(sample_path, bin_size=160):
    """Load gene expression data from one sample."""
    sample_path = Path(sample_path)

    # Load genes
    with open(sample_path / "genes.csv") as f:
        genes = [line.strip() for line in f if line.strip()]

    # Load spatial expression data (genes × width × height)
    zarr_file = sample_path / "zarr" / f"bins_size_{bin_size}.zarr.zip"
    store = zarr.storage.ZipStore(str(zarr_file), mode='r')
    expression = zarr.open_array(store, mode='r')
    expression = np.asarray(expression)  # Convert to numpy: (genes, width, height)

    return genes, expression


# ============================================================================
# STEP 2: Extract circular patches and encode them
# ============================================================================

def make_circular_patch(expression, center_x, center_y, radius_bins):
    """Extract a circular patch from the expression data."""
    n_genes, width, height = expression.shape

    # Find all bins within the radius
    x_coords = []
    y_coords = []
    for x in range(max(0, center_x - radius_bins), min(width, center_x + radius_bins + 1)):
        for y in range(max(0, center_y - radius_bins), min(height, center_y + radius_bins + 1)):
            # Check if inside circle
            dist = np.sqrt((x - center_x)**2 + (y - center_y)**2)
            if dist <= radius_bins:
                x_coords.append(x)
                y_coords.append(y)

    if len(x_coords) == 0:
        return None, None, None

    # Extract expression for these bins: (genes, n_bins)
    patch_expression = expression[:, x_coords, y_coords]

    return np.array(x_coords), np.array(y_coords), patch_expression


def encode_patch_radial(x_coords, y_coords, patch_expression, n_shells=5):
    """
    Encode a patch by dividing it into concentric shells.

    Returns a vector: [gene0_shell0, gene0_shell1, ..., gene0_shell4, gene1_shell0, ...]
    """
    n_genes = patch_expression.shape[0]

    # Calculate center
    center_x = x_coords.mean()
    center_y = y_coords.mean()

    # Calculate distance from center for each bin
    distances = np.sqrt((x_coords - center_x)**2 + (y_coords - center_y)**2)
    max_dist = distances.max()

    if max_dist == 0:
        return patch_expression.mean(axis=1).repeat(n_shells)  # All bins at same location

    # Create shell boundaries
    shell_edges = np.linspace(0, max_dist, n_shells + 1)

    # For each gene, calculate average expression in each shell
    features = []
    for gene_idx in range(n_genes):
        gene_expr = patch_expression[gene_idx, :]

        for shell_idx in range(n_shells):
            # Find bins in this shell
            in_shell = (distances >= shell_edges[shell_idx]) & (distances < shell_edges[shell_idx + 1])

            if in_shell.sum() > 0:
                features.append(gene_expr[in_shell].mean())
            else:
                features.append(0.0)

    return np.array(features, dtype=np.float32)


# ============================================================================
# STEP 3: Build database of patches from all samples
# ============================================================================

def build_patch_database(data_dir, radius_bins=5, bin_size=160):
    """Generate patches from all samples and encode them."""
    data_dir = Path(data_dir)

    # Find all sample directories
    sample_dirs = [d for d in data_dir.iterdir() if d.is_dir() and (d / "zarr").exists()]
    print(f"Found {len(sample_dirs)} samples")

    # Collect all genes from all samples
    all_genes = set()
    for sample_dir in sample_dirs:
        with open(sample_dir / "genes.csv") as f:
            genes = [line.strip() for line in f if line.strip()]
            all_genes.update(genes)

    global_genes = sorted(list(all_genes))
    print(f"Total genes: {len(global_genes)}")

    # Generate patches from each sample
    database = []

    for sample_dir in sample_dirs:
        print(f"Processing {sample_dir.name}...")
        genes, expression = load_sample(sample_dir, bin_size)

        # Map sample genes to global genes
        gene_map = {g: i for i, g in enumerate(global_genes)}
        global_expression = np.zeros((len(global_genes), expression.shape[1], expression.shape[2]))
        for i, gene in enumerate(genes):
            if gene in gene_map:
                global_expression[gene_map[gene]] = expression[i]

        # Generate patches at regular intervals
        _, width, height = global_expression.shape
        step = radius_bins * 2  # Don't overlap too much

        for cx in range(radius_bins, width - radius_bins, step):
            for cy in range(radius_bins, height - radius_bins, step):
                x_coords, y_coords, patch_expr = make_circular_patch(global_expression, cx, cy, radius_bins)

                if patch_expr is None:
                    continue

                # Encode the patch
                encoding = encode_patch_radial(x_coords, y_coords, patch_expr)

                database.append({
                    'sample': sample_dir.name,
                    'center_x': cx,
                    'center_y': cy,
                    'radius': radius_bins,
                    'encoding': encoding
                })

        print(f"  Generated {len([d for d in database if d['sample'] == sample_dir.name])} patches")

    return database, global_genes


# ============================================================================
# STEP 4: Search for similar patches
# ============================================================================

def search_similar(query_encoding, database, top_k=10):
    """Find most similar patches using cosine similarity."""
    similarities = []

    for entry in database:
        # Cosine similarity (1 - cosine distance)
        sim = 1 - cosine(query_encoding, entry['encoding'])
        similarities.append(sim)

    # Sort by similarity (highest first)
    sorted_indices = np.argsort(similarities)[::-1][:top_k]

    results = []
    for idx in sorted_indices:
        results.append({
            'sample': database[idx]['sample'],
            'center_x': database[idx]['center_x'],
            'center_y': database[idx]['center_y'],
            'similarity': similarities[idx]
        })

    return results


# ============================================================================
# STEP 5: Visualize a single gene's radial pattern
# ============================================================================

def plot_gene_pattern(encoding, gene_idx, gene_name, n_shells=5):
    """Plot how one gene's expression changes from center to edge."""
    # Extract this gene's features (one per shell)
    gene_features = encoding[gene_idx * n_shells : (gene_idx + 1) * n_shells]

    plt.figure(figsize=(8, 5))
    plt.bar(range(n_shells), gene_features, color='steelblue', alpha=0.7)
    plt.xlabel('Shell (center → edge)')
    plt.ylabel('Average Expression')
    plt.title(f'{gene_name} - Radial Expression Pattern')
    plt.xticks(range(n_shells), [f'Shell {i+1}' for i in range(n_shells)])
    plt.grid(axis='y', alpha=0.3)

    return plt.gcf()


# ============================================================================
# MAIN SCRIPT
# ============================================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Simple radial pattern search')
    parser.add_argument('--data-dir', required=True, help='Directory with samples')
    parser.add_argument('--query-sample', required=True, help='Sample to query from')
    parser.add_argument('--query-x', type=int, default=50, help='Query center X')
    parser.add_argument('--query-y', type=int, default=50, help='Query center Y')
    parser.add_argument('--radius', type=int, default=5, help='Patch radius in bins')
    parser.add_argument('--gene', default='ACTA2', help='Gene to visualize')
    parser.add_argument('--top-k', type=int, default=10, help='Number of results')
    args = parser.parse_args()

    # Build database
    print("\n=== Building Database ===")
    database, global_genes = build_patch_database(args.data_dir, args.radius)
    print(f"Database: {len(database)} patches")

    # Create query
    print(f"\n=== Creating Query ===")
    print(f"Sample: {args.query_sample}")
    print(f"Location: ({args.query_x}, {args.query_y})")

    query_path = Path(args.data_dir) / args.query_sample
    genes, expression = load_sample(query_path)

    # Map to global genes
    gene_map = {g: i for i, g in enumerate(global_genes)}
    global_expression = np.zeros((len(global_genes), expression.shape[1], expression.shape[2]))
    for i, gene in enumerate(genes):
        if gene in gene_map:
            global_expression[gene_map[gene]] = expression[i]

    x_coords, y_coords, patch_expr = make_circular_patch(global_expression, args.query_x, args.query_y, args.radius)
    query_encoding = encode_patch_radial(x_coords, y_coords, patch_expr)

    # Search
    print(f"\n=== Searching ===")
    results = search_similar(query_encoding, database, args.top_k)

    print(f"\nTop {args.top_k} matches:")
    for i, result in enumerate(results):
        print(f"{i+1}. {result['sample']} at ({result['center_x']}, {result['center_y']}) - similarity: {result['similarity']:.3f}")

    # Visualize one gene
    if args.gene in global_genes:
        gene_idx = global_genes.index(args.gene)
        fig = plot_gene_pattern(query_encoding, gene_idx, args.gene)
        fig.savefig(f'{args.gene}_pattern.png', dpi=100, bbox_inches='tight')
        print(f"\nSaved visualization: {args.gene}_pattern.png")
    else:
        print(f"\nGene {args.gene} not found. Available: {global_genes[:10]}...")


if __name__ == '__main__':
    main()
