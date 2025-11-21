#!/usr/bin/env python3
"""
Visualize Non-PCA Radial Shell Database

This script visualizes the database WITHOUT PCA compression,
showing individual gene patterns across radial shells.

Usage:
    python visualize_no_pca_database.py \
        --database spatial_database/no_pca_db \
        --resolution 160 \
        --gene-name TP53
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import json
import faiss


def load_embeddings_from_faiss(index_path: Path) -> np.ndarray:
    """Load embeddings directly from FAISS index."""
    index = faiss.read_index(str(index_path))

    # Extract vectors from the index
    n_vectors = index.ntotal
    d = index.d  # dimension

    # Reconstruct all vectors
    embeddings = np.zeros((n_vectors, d), dtype=np.float32)
    for i in range(n_vectors):
        embeddings[i] = index.reconstruct(i)

    return embeddings


def load_database_info(database_path: Path, resolution_um: int):
    """Load database configuration and metadata (no PCA version)."""
    res_dir = database_path / f"{resolution_um}um"

    if not res_dir.exists():
        raise FileNotFoundError(f"Resolution directory not found: {res_dir}")

    # Load configuration
    config_path = res_dir / "config.json"
    with open(config_path, 'r') as f:
        config = json.load(f)

    # Load metadata
    metadata_path = res_dir / "patches_metadata.parquet"
    metadata = pd.read_parquet(metadata_path)

    # Check if PCA was used
    use_pca = config.get('use_pca', True)

    print(f"Database Info for {resolution_um}μm resolution:")
    print(f"  Total genes: {len(config['global_genes'])}")
    print(f"  Variable genes: {'Yes' if config['use_variable_genes'] else 'No'}")
    if 'variable_gene_indices' in config:
        print(f"  Num variable genes: {len(config['variable_gene_indices'])}")
    print(f"  N shells: {config['n_shells']}")
    print(f"  PCA compression: {'Yes' if use_pca else 'No'}")
    if use_pca:
        print(f"  PCA dims: {config['pca_dims']}")
    else:
        expected_dims = len(config['global_genes']) * config['n_shells']
        print(f"  Full embedding dims: {expected_dims}")
    print(f"  Total patches: {len(metadata)}")
    print(f"  Radii available: {sorted(metadata['radius'].unique())}")
    print()

    return {
        'config': config,
        'metadata': metadata,
        'res_dir': res_dir
    }


def visualize_gene_heatmap(
    gene_name: str,
    gene_idx: int,
    n_shells: int,
    radii: list,
    embeddings_by_radius: dict,
    config: dict,
    max_patches: int = 100
):
    """
    Create comprehensive heatmaps for a single gene across different radii.

    Shows how the gene's expression varies across:
    - Different radial shells (center to edge)
    - Different patches (spatial locations)
    - Different patch sizes (radii)
    """
    n_radii = len(radii)
    fig, axes = plt.subplots(3, n_radii, figsize=(6*n_radii, 15))

    if n_radii == 1:
        axes = axes.reshape(3, 1)

    for col_idx, radius in enumerate(radii):
        embeddings = embeddings_by_radius[radius]
        n_patches = min(embeddings.shape[0], max_patches)

        # Extract features for this gene across all shells
        # Features are ordered as: [gene0_shell0, gene0_shell1, ..., gene1_shell0, ...]
        gene_features = []
        for shell_idx in range(n_shells):
            feature_idx = gene_idx * n_shells + shell_idx
            if feature_idx < embeddings.shape[1]:
                gene_features.append(embeddings[:n_patches, feature_idx])
            else:
                print(f"Warning: feature index {feature_idx} out of range for embeddings shape {embeddings.shape}")
                gene_features.append(np.zeros(n_patches))

        gene_features = np.array(gene_features)  # Shape: (n_shells, n_patches)

        # Top plot: Heatmap showing expression across shells and patches
        ax = axes[0, col_idx]
        im = ax.imshow(gene_features, aspect='auto', cmap='YlOrRd', interpolation='nearest')
        ax.set_ylabel('Shell (center→edge)', fontsize=11, fontweight='bold')
        ax.set_xlabel('Patch Index', fontsize=11)
        ax.set_title(f'Radius {radius} bins\n{n_patches} patches', fontsize=12, fontweight='bold')
        ax.set_yticks(range(n_shells))
        ax.set_yticklabels([f'{i+1}' for i in range(n_shells)])

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Expression', fontsize=10)

        # Middle plot: Average expression profile (center to edge)
        ax = axes[1, col_idx]
        avg_per_shell = np.mean(gene_features, axis=1)
        std_per_shell = np.std(gene_features, axis=1)

        shells = np.arange(n_shells)
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, n_shells))

        bars = ax.bar(shells, avg_per_shell, yerr=std_per_shell,
                      color=colors, alpha=0.8, capsize=5, edgecolor='black', linewidth=1.5)

        ax.set_xlabel('Shell Index (center→edge)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Mean Expression ± SD', fontsize=11, fontweight='bold')
        ax.set_title('Radial Expression Profile', fontsize=12, fontweight='bold')
        ax.set_xticks(shells)
        ax.set_xticklabels([f'{i+1}' for i in range(n_shells)])
        ax.grid(True, alpha=0.3, axis='y', linestyle='--')
        ax.set_axisbelow(True)

        # Add value labels on bars
        for i, (bar, val) in enumerate(zip(bars, avg_per_shell)):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.2f}',
                   ha='center', va='bottom', fontsize=9, fontweight='bold')

        # Bottom plot: Distribution of expression values per shell (violin plot)
        ax = axes[2, col_idx]

        # Prepare data for violin plot
        data_for_violin = []
        positions = []
        for shell_idx in range(n_shells):
            shell_data = gene_features[shell_idx, :]
            # Remove zeros for better visualization of actual expression
            non_zero_data = shell_data[shell_data > 0]
            if len(non_zero_data) > 0:
                data_for_violin.append(non_zero_data)
                positions.append(shell_idx)

        if len(data_for_violin) > 0:
            parts = ax.violinplot(data_for_violin, positions=positions,
                                 widths=0.7, showmeans=True, showextrema=True)

            # Color the violins
            for i, pc in enumerate(parts['bodies']):
                pc.set_facecolor(colors[positions[i]])
                pc.set_alpha(0.7)

            ax.set_xlabel('Shell Index (center→edge)', fontsize=11, fontweight='bold')
            ax.set_ylabel('Expression Distribution\n(non-zero values)', fontsize=11, fontweight='bold')
            ax.set_title('Expression Variability', fontsize=12, fontweight='bold')
            ax.set_xticks(range(n_shells))
            ax.set_xticklabels([f'{i+1}' for i in range(n_shells)])
            ax.grid(True, alpha=0.3, axis='y', linestyle='--')
            ax.set_axisbelow(True)
        else:
            ax.text(0.5, 0.5, 'No expression detected',
                   ha='center', va='center', transform=ax.transAxes,
                   fontsize=12, color='red')
            ax.set_xticks([])
            ax.set_yticks([])

    # Add main title
    fig.suptitle(f'{gene_name} - Radial Shell Expression Pattern',
                 fontsize=16, fontweight='bold', y=0.995)

    # Add interpretation guide
    fig.text(0.5, 0.02,
             'Interpretation: Top=Expression heatmap | Middle=Average radial profile | Bottom=Distribution per shell',
             ha='center', fontsize=10, style='italic', color='gray')

    plt.tight_layout(rect=[0, 0.03, 1, 0.99])

    return fig


def create_multi_gene_comparison(
    gene_names: list,
    embeddings_by_radius: dict,
    config: dict,
    radius: int = 5,
    max_patches: int = 50
):
    """Compare multiple genes side by side."""

    n_genes = len(gene_names)
    n_shells = config['n_shells']
    global_genes = config['global_genes']

    fig, axes = plt.subplots(n_genes, 2, figsize=(14, 4*n_genes))

    if n_genes == 1:
        axes = axes.reshape(1, 2)

    embeddings = embeddings_by_radius[radius]
    n_patches = min(embeddings.shape[0], max_patches)

    for gene_idx, gene_name in enumerate(gene_names):
        try:
            gene_id = global_genes.index(gene_name)
        except ValueError:
            print(f"Warning: Gene {gene_name} not found in global genes")
            continue

        # Extract gene features
        gene_features = []
        for shell_idx in range(n_shells):
            feature_idx = gene_id * n_shells + shell_idx
            if feature_idx < embeddings.shape[1]:
                gene_features.append(embeddings[:n_patches, feature_idx])

        gene_features = np.array(gene_features)

        # Left: Heatmap
        ax = axes[gene_idx, 0]
        im = ax.imshow(gene_features, aspect='auto', cmap='YlOrRd', interpolation='nearest')
        ax.set_ylabel(f'{gene_name}\nShell', fontsize=10, fontweight='bold')
        ax.set_xlabel('Patch Index', fontsize=10)
        ax.set_yticks(range(n_shells))
        ax.set_yticklabels([f'{i+1}' for i in range(n_shells)])
        plt.colorbar(im, ax=ax, label='Expression')

        # Right: Radial profile
        ax = axes[gene_idx, 1]
        avg_per_shell = np.mean(gene_features, axis=1)
        std_per_shell = np.std(gene_features, axis=1)

        shells = np.arange(n_shells)
        ax.bar(shells, avg_per_shell, yerr=std_per_shell,
               color='steelblue', alpha=0.7, capsize=5)
        ax.set_xlabel('Shell (center→edge)', fontsize=10, fontweight='bold')
        ax.set_ylabel('Expression', fontsize=10, fontweight='bold')
        ax.set_title(f'Radial Profile', fontsize=11, fontweight='bold')
        ax.set_xticks(shells)
        ax.set_xticklabels([f'{i+1}' for i in range(n_shells)])
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle(f'Gene Comparison - Radius {radius} bins',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    return fig


def main():
    parser = argparse.ArgumentParser(description='Visualize Non-PCA Radial Shell Database')
    parser.add_argument('--database', type=str, required=True,
                       help='Path to database directory')
    parser.add_argument('--resolution', type=int, default=160,
                       help='Resolution to visualize (default: 160)')
    parser.add_argument('--gene-name', type=str, default=None,
                       help='Gene name to visualize (default: first gene)')
    parser.add_argument('--compare-genes', type=str, nargs='+', default=None,
                       help='Multiple gene names to compare')
    parser.add_argument('--output-dir', type=str, default='database_plots',
                       help='Directory to save plots')
    parser.add_argument('--max-patches', type=int, default=100,
                       help='Maximum number of patches to show (default: 100)')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load database info
    print("Loading database information...")
    db_info = load_database_info(Path(args.database), args.resolution)

    config = db_info['config']
    metadata = db_info['metadata']
    res_dir = db_info['res_dir']

    # Check if PCA was used
    if config.get('use_pca', True):
        print("\nWARNING: This database uses PCA compression!")
        print("Individual gene information is not preserved.")
        print("Use --use-pca False when building the database to see gene-level patterns.")
        return

    # Get available radii
    radii = sorted(metadata['radius'].unique())
    print(f"Available radii: {radii}")

    # Load embeddings for each radius from FAISS indices
    print("\nLoading embeddings from FAISS indices...")
    embeddings_by_radius = {}
    faiss_dir = res_dir / "faiss_indices"

    for radius in radii:
        index_path = faiss_dir / f"faiss_r{radius}.index"
        if index_path.exists():
            embeddings = load_embeddings_from_faiss(index_path)
            embeddings_by_radius[radius] = embeddings
            print(f"  Radius {radius}: {embeddings.shape[0]} patches × {embeddings.shape[1]} features")
        else:
            print(f"  Warning: Index not found: {index_path}")

    if not embeddings_by_radius:
        print("\nNo FAISS indices found!")
        return

    # Single gene visualization
    if args.gene_name or not args.compare_genes:
        if args.gene_name:
            try:
                gene_idx = config['global_genes'].index(args.gene_name)
                gene_name = args.gene_name
            except ValueError:
                print(f"\nError: Gene '{args.gene_name}' not found!")
                print(f"Available genes: {config['global_genes'][:20]}...")
                return
        else:
            gene_idx = 0
            gene_name = config['global_genes'][0]

        print(f"\nCreating visualization for gene: {gene_name} (index {gene_idx})...")
        fig = visualize_gene_heatmap(
            gene_name=gene_name,
            gene_idx=gene_idx,
            n_shells=config['n_shells'],
            radii=radii,
            embeddings_by_radius=embeddings_by_radius,
            config=config,
            max_patches=args.max_patches
        )

        output_path = output_dir / f'{gene_name}_radial_pattern.png'
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close(fig)

    # Multi-gene comparison
    if args.compare_genes:
        print(f"\nCreating comparison for genes: {args.compare_genes}...")
        fig = create_multi_gene_comparison(
            gene_names=args.compare_genes,
            embeddings_by_radius=embeddings_by_radius,
            config=config,
            radius=radii[1] if len(radii) > 1 else radii[0],  # Use middle radius
            max_patches=args.max_patches
        )

        output_path = output_dir / f'gene_comparison_{"_".join(args.compare_genes)}.png'
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close(fig)

    print(f"\n✓ All visualizations saved to: {output_dir}")
    print("\nTo visualize a different gene:")
    print(f"  python {Path(__file__).name} --database {args.database} --resolution {args.resolution} --gene-name GENE_NAME")
    print("\nTo compare multiple genes:")
    print(f"  python {Path(__file__).name} --database {args.database} --resolution {args.resolution} --compare-genes GENE1 GENE2 GENE3")


if __name__ == '__main__':
    main()
