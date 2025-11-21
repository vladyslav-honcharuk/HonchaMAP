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
import matplotlib.pyplot as plt

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
# STEP 5: Visualize spatial data
# ============================================================================

def visualize_matches_spatial(data_dir, query_sample_name, query_x, query_y,
                               results, genes, gene_name='ACTA2', radius=5, n_show=4):
    """
    Visualize original spatial gene expression for query and top matches.

    Shows actual tissue bin locations (not circular patches) with gene expression as color.
    """
    # Find gene index
    if gene_name not in genes:
        print(f"Gene {gene_name} not in database. Using first gene: {genes[0]}")
        gene_name = genes[0]

    gene_idx = genes.index(gene_name)

    # Prepare subplots: query + top matches
    n_show = min(n_show, len(results) + 1)
    fig, axes = plt.subplots(1, n_show, figsize=(5*n_show, 5))
    if n_show == 1:
        axes = [axes]

    # Plot query
    print(f"\nCreating spatial visualization for gene: {gene_name}")
    print("-"*60)

    query_path = data_dir / query_sample_name
    sample_genes, expression = load_sample(query_path)

    # Map to global genes
    gene_map = {g: i for i, g in enumerate(genes)}
    global_expr = np.zeros((len(genes), expression.shape[1], expression.shape[2]))
    for i, gene in enumerate(sample_genes):
        if gene in gene_map:
            global_expr[gene_map[gene]] = expression[i]

    # Extract patch region (not just circular, but square around it for context)
    margin = radius + 3
    x_min = max(0, query_x - margin)
    x_max = min(global_expr.shape[1], query_x + margin)
    y_min = max(0, query_y - margin)
    y_max = min(global_expr.shape[2], query_y + margin)

    # Get all bins in region
    x_coords_all = []
    y_coords_all = []
    expression_vals = []

    for x in range(x_min, x_max):
        for y in range(y_min, y_max):
            x_coords_all.append(x)
            y_coords_all.append(y)
            expression_vals.append(global_expr[gene_idx, x, y])

    x_coords_all = np.array(x_coords_all)
    y_coords_all = np.array(y_coords_all)
    expression_vals = np.array(expression_vals)

    # Plot query
    ax = axes[0]
    scatter = ax.scatter(x_coords_all, y_coords_all, c=expression_vals,
                        cmap='YlOrRd', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)

    # Mark the patch center
    ax.scatter([query_x], [query_y], c='blue', s=200, marker='x', linewidths=3, label='Center')

    # Draw circle showing patch radius
    circle = plt.Circle((query_x, query_y), radius, color='blue', fill=False, linewidth=2, linestyle='--')
    ax.add_patch(circle)

    ax.set_title(f'QUERY\n{query_sample_name}\n({query_x}, {query_y})', fontsize=10, fontweight='bold')
    ax.set_xlabel('X coordinate (bins)')
    ax.set_ylabel('Y coordinate (bins)')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label=f'{gene_name} expression')

    print(f"Query: {query_sample_name} at ({query_x}, {query_y})")

    # Plot top matches
    for idx, result in enumerate(results[:n_show-1]):
        ax = axes[idx + 1]

        # Load match sample
        match_path = data_dir / result['sample']
        match_genes, match_expression = load_sample(match_path)

        # Map to global genes
        match_global_expr = np.zeros((len(genes), match_expression.shape[1], match_expression.shape[2]))
        for i, gene in enumerate(match_genes):
            if gene in gene_map:
                match_global_expr[gene_map[gene]] = match_expression[i]

        # Extract region around match
        match_x, match_y = result['x'], result['y']
        x_min = max(0, match_x - margin)
        x_max = min(match_global_expr.shape[1], match_x + margin)
        y_min = max(0, match_y - margin)
        y_max = min(match_global_expr.shape[2], match_y + margin)

        x_coords_match = []
        y_coords_match = []
        expression_vals_match = []

        for x in range(x_min, x_max):
            for y in range(y_min, y_max):
                x_coords_match.append(x)
                y_coords_match.append(y)
                expression_vals_match.append(match_global_expr[gene_idx, x, y])

        x_coords_match = np.array(x_coords_match)
        y_coords_match = np.array(y_coords_match)
        expression_vals_match = np.array(expression_vals_match)

        # Plot match
        scatter = ax.scatter(x_coords_match, y_coords_match, c=expression_vals_match,
                           cmap='YlOrRd', s=50, alpha=0.7, edgecolors='black', linewidth=0.5)

        # Mark center
        ax.scatter([match_x], [match_y], c='blue', s=200, marker='x', linewidths=3)

        # Draw circle
        circle = plt.Circle((match_x, match_y), radius, color='blue', fill=False, linewidth=2, linestyle='--')
        ax.add_patch(circle)

        ax.set_title(f'MATCH #{idx+1}\n{result["sample"]}\n({match_x}, {match_y})\nSimilarity: {result["similarity"]:.3f}',
                    fontsize=10, fontweight='bold')
        ax.set_xlabel('X coordinate (bins)')
        ax.set_ylabel('Y coordinate (bins)')
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        plt.colorbar(scatter, ax=ax, label=f'{gene_name} expression')

        print(f"Match #{idx+1}: {result['sample']} at ({match_x}, {match_y}) - similarity: {result['similarity']:.3f}")

    plt.suptitle(f'Spatial Expression of {gene_name}\n(Original tissue coordinates, blue circle = patch radius)',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    # Save
    output_file = f'spatial_matches_{gene_name}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved spatial visualization: {output_file}")

    return output_file


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

    # Show query encoding details
    print("\nQuery Encoding (first 3 genes):")
    print("-"*60)
    for gene_idx in range(min(3, len(genes))):
        gene_name = genes[gene_idx]
        gene_features = query_encoding[gene_idx * N_SHELLS : (gene_idx + 1) * N_SHELLS]
        print(f"{gene_name:10s}: [", end="")
        print(", ".join([f"{v:5.2f}" for v in gene_features]), end="")
        print("]  (ring 1→5: center→edge)")
    print()

    # Search
    print("Searching...")
    print("-"*60)
    results = search(query_encoding, database, top_k=10)

    print("\nTop 10 matches:")
    for i, r in enumerate(results):
        print(f"{i+1:2d}. {r['sample']:12s} ({r['x']:3d}, {r['y']:3d})  similarity: {r['similarity']:.3f}")

    # Show detailed comparison for top match (if not the query itself)
    if len(results) > 1:
        print("\n" + "="*60)
        print("DETAILED COMPARISON: Query vs Top Match")
        print("="*60)

        top_match = results[1] if results[0]['similarity'] > 0.99 else results[0]

        # Find the match in database
        match_encoding = None
        for entry in database:
            if (entry['sample'] == top_match['sample'] and
                entry['x'] == top_match['x'] and
                entry['y'] == top_match['y']):
                match_encoding = entry['encoding']
                break

        if match_encoding is not None:
            print(f"\nQuery:     {query_sample.name} at ({query_x}, {query_y})")
            print(f"Match:     {top_match['sample']} at ({top_match['x']}, {top_match['y']})")
            print(f"Similarity: {top_match['similarity']:.3f}\n")

            print("Gene expression comparison (first 5 genes, ring 1→5):")
            print("-"*80)
            print(f"{'Gene':<12s} | {'Query Rings':<30s} | {'Match Rings':<30s}")
            print("-"*80)

            for gene_idx in range(min(5, len(genes))):
                gene_name = genes[gene_idx]
                query_rings = query_encoding[gene_idx * N_SHELLS : (gene_idx + 1) * N_SHELLS]
                match_rings = match_encoding[gene_idx * N_SHELLS : (gene_idx + 1) * N_SHELLS]

                query_str = " ".join([f"{v:5.2f}" for v in query_rings])
                match_str = " ".join([f"{v:5.2f}" for v in match_rings])

                print(f"{gene_name:<12s} | {query_str:<30s} | {match_str:<30s}")

            print("-"*80)
            print("(Numbers show average gene expression in each ring)")
            print("(Ring 1 = center, Ring 5 = edge)")

    # Visualize spatial data for query and top matches
    print("\n" + "="*60)
    print("SPATIAL VISUALIZATION")
    print("="*60)

    # Use first gene that has non-zero expression in query
    viz_gene = None
    for gene_idx in range(min(20, len(genes))):
        gene_features = query_encoding[gene_idx * N_SHELLS : (gene_idx + 1) * N_SHELLS]
        if gene_features.sum() > 0.1:  # Has some expression
            viz_gene = genes[gene_idx]
            break

    if viz_gene is None:
        viz_gene = genes[0]  # Fallback to first gene

    visualize_matches_spatial(
        data_dir=data_dir,
        query_sample_name=query_sample.name,
        query_x=query_x,
        query_y=query_y,
        results=results,
        genes=genes,
        gene_name=viz_gene,
        radius=RADIUS,
        n_show=4  # Show query + top 3 matches
    )

    print()
    print("="*60)
    print("DONE!")
    print("="*60)


if __name__ == '__main__':
    main()
