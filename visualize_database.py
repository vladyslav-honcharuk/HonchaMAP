#!/usr/bin/env python3
"""
Visualize Radial Shell Encoding Database

This script visualizes what's stored in the radial shell database:
- Shows the radial shell encoding structure
- Displays embeddings for different radii
- Creates heatmaps showing feature patterns

Usage:
    python visualize_database.py --database spatial_database --resolution 160
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import h5py
import json
from typing import Dict, List, Tuple


def load_database_info(database_path: Path, resolution_um: int) -> Dict:
    """Load database configuration and metadata."""
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

    # Load PCA models
    pca_path = res_dir / "pca_models.pkl"
    import pickle
    with open(pca_path, 'rb') as f:
        pca_models = pickle.load(f)

    print(f"Database Info for {resolution_um}μm resolution:")
    print(f"  Total genes: {len(config['global_genes'])}")
    print(f"  Variable genes: {'Yes' if config['use_variable_genes'] else 'No'}")
    if 'variable_gene_indices' in config:
        print(f"  Num variable genes: {len(config['variable_gene_indices'])}")
    print(f"  N shells: {config['n_shells']}")
    print(f"  PCA dims: {config['pca_dims']}")
    print(f"  Total patches: {len(metadata)}")
    print(f"  Radii available: {sorted(metadata['radius'].unique())}")
    print()

    return {
        'config': config,
        'metadata': metadata,
        'pca_models': pca_models,
        'res_dir': res_dir
    }


def load_raw_embeddings(res_dir: Path, radius: int) -> Tuple[np.ndarray, np.ndarray]:
    """Load raw embeddings for a specific radius."""
    h5_path = res_dir / f"raw_embeddings_r{radius}.h5"

    if not h5_path.exists():
        raise FileNotFoundError(f"Raw embeddings not found: {h5_path}")

    with h5py.File(h5_path, 'r') as f:
        embeddings = f['embeddings'][:]
        patch_ids = f['patch_ids'][:]

    return embeddings, patch_ids


def visualize_radial_shell_concept(n_shells: int = 5):
    """Visualize the radial shell encoding concept."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left plot: Show the shells
    ax = axes[0]
    center = np.array([0, 0])
    max_radius = 1.0

    # Draw concentric shells
    shell_boundaries = np.linspace(0, max_radius, n_shells + 1)
    colors = plt.cm.viridis(np.linspace(0, 1, n_shells))

    for i in range(n_shells):
        inner_r = shell_boundaries[i]
        outer_r = shell_boundaries[i + 1]

        # Draw shell as annulus
        theta = np.linspace(0, 2*np.pi, 100)

        # Outer circle
        x_outer = center[0] + outer_r * np.cos(theta)
        y_outer = center[1] + outer_r * np.sin(theta)

        if i == 0:
            # First shell is a filled circle
            ax.fill(x_outer, y_outer, color=colors[i], alpha=0.3, label=f'Shell {i+1}')
        else:
            # Other shells are annuli
            x_inner = center[0] + inner_r * np.cos(theta)
            y_inner = center[1] + inner_r * np.sin(theta)

            # Create polygon for annulus
            vertices = np.column_stack([
                np.concatenate([x_outer, x_inner[::-1]]),
                np.concatenate([y_outer, y_inner[::-1]])
            ])
            ax.fill(vertices[:, 0], vertices[:, 1], color=colors[i], alpha=0.3, label=f'Shell {i+1}')

    # Add some example bins
    np.random.seed(42)
    n_bins = 50
    angles = np.random.uniform(0, 2*np.pi, n_bins)
    radii = np.random.uniform(0, max_radius, n_bins)
    bin_x = center[0] + radii * np.cos(angles)
    bin_y = center[1] + radii * np.sin(angles)

    ax.scatter(bin_x, bin_y, c='black', s=30, alpha=0.5, zorder=10, label='Bins')
    ax.scatter(*center, c='red', s=100, marker='*', zorder=11, label='Centroid')

    ax.set_aspect('equal')
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xlabel('X (normalized)')
    ax.set_ylabel('Y (normalized)')
    ax.set_title(f'Radial Shell Encoding ({n_shells} shells)')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)

    # Right plot: Show encoding structure
    ax = axes[1]

    # Example: 10 genes, 5 shells
    n_genes_example = 10
    gene_names = [f'Gene{i+1}' for i in range(n_genes_example)]
    shell_names = [f'Shell{i+1}' for i in range(n_shells)]

    # Create example encoding matrix
    encoding_matrix = np.random.rand(n_genes_example, n_shells)

    # Create heatmap
    im = ax.imshow(encoding_matrix, aspect='auto', cmap='YlOrRd')

    # Set ticks
    ax.set_xticks(range(n_shells))
    ax.set_xticklabels(shell_names)
    ax.set_yticks(range(n_genes_example))
    ax.set_yticklabels(gene_names)

    ax.set_xlabel('Radial Shells')
    ax.set_ylabel('Genes')
    ax.set_title('Feature Encoding\n(Gene × Shell = Feature)')

    plt.colorbar(im, ax=ax, label='Expression')

    # Add grid
    for i in range(n_genes_example + 1):
        ax.axhline(i - 0.5, color='white', linewidth=0.5)
    for i in range(n_shells + 1):
        ax.axvline(i - 0.5, color='white', linewidth=0.5)

    plt.suptitle('Radial Shell Encoding: How It Works', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    return fig


def visualize_single_gene_encoding(
    gene_name: str,
    gene_idx: int,
    n_shells: int,
    radii: List[int],
    embeddings_by_radius: Dict[int, np.ndarray],
    config: Dict
):
    """Visualize how a single gene is encoded across different radii."""

    n_radii = len(radii)
    fig, axes = plt.subplots(2, n_radii, figsize=(5*n_radii, 10))

    if n_radii == 1:
        axes = axes.reshape(2, 1)

    for col_idx, radius in enumerate(radii):
        embeddings = embeddings_by_radius[radius]
        n_patches = min(embeddings.shape[0], 100)  # Show up to 100 patches

        # Extract features for this gene across all shells
        # Features are ordered as: [gene0_shell0, gene0_shell1, ..., gene1_shell0, ...]
        gene_features = []
        for shell_idx in range(n_shells):
            feature_idx = gene_idx * n_shells + shell_idx
            gene_features.append(embeddings[:n_patches, feature_idx])

        gene_features = np.array(gene_features)  # Shape: (n_shells, n_patches)

        # Top plot: Heatmap of this gene across shells and patches
        ax = axes[0, col_idx]
        im = ax.imshow(gene_features, aspect='auto', cmap='viridis', interpolation='nearest')
        ax.set_ylabel('Shell Index')
        ax.set_xlabel('Patch Index')
        ax.set_title(f'Radius {radius} bins\n{n_patches} patches')
        ax.set_yticks(range(n_shells))
        ax.set_yticklabels([f'Shell {i+1}' for i in range(n_shells)])
        plt.colorbar(im, ax=ax, label='Expression')

        # Bottom plot: Average expression per shell
        ax = axes[1, col_idx]
        avg_per_shell = np.mean(gene_features, axis=1)
        std_per_shell = np.std(gene_features, axis=1)

        shells = np.arange(n_shells)
        ax.bar(shells, avg_per_shell, yerr=std_per_shell,
               color='steelblue', alpha=0.7, capsize=5)
        ax.set_xlabel('Shell Index')
        ax.set_ylabel('Average Expression')
        ax.set_title(f'Mean ± SD per Shell')
        ax.set_xticks(shells)
        ax.set_xticklabels([f'{i+1}' for i in range(n_shells)])
        ax.grid(True, alpha=0.3, axis='y')

    # Add gene name to the figure
    fig.suptitle(f'Gene: {gene_name} - Radial Shell Encoding Across Different Radii',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    return fig


def visualize_feature_statistics(
    embeddings_by_radius: Dict[int, np.ndarray],
    config: Dict,
    n_genes_to_show: int = 20
):
    """Visualize statistics of features across different radii."""

    n_shells = config['n_shells']
    radii = sorted(embeddings_by_radius.keys())
    n_radii = len(radii)

    fig, axes = plt.subplots(2, n_radii, figsize=(5*n_radii, 10))

    if n_radii == 1:
        axes = axes.reshape(2, 1)

    for col_idx, radius in enumerate(radii):
        embeddings = embeddings_by_radius[radius]

        # Top plot: Feature distribution (first 100 features)
        ax = axes[0, col_idx]
        features_to_plot = embeddings[:, :100].T  # First 100 features across all patches

        # Create violin plot
        positions = range(min(20, features_to_plot.shape[0]))
        data_to_plot = [features_to_plot[i, :] for i in positions]

        parts = ax.violinplot(data_to_plot, positions=positions, widths=0.7,
                             showmeans=True, showextrema=True)

        ax.set_xlabel('Feature Index')
        ax.set_ylabel('Feature Value')
        ax.set_title(f'Radius {radius}: Feature Distribution\n(First 20 features)')
        ax.grid(True, alpha=0.3, axis='y')

        # Bottom plot: Average activation per gene
        ax = axes[1, col_idx]

        # Calculate average activation per gene (average across all shells)
        gene_activations = []
        for gene_idx in range(min(n_genes_to_show, embeddings.shape[1] // n_shells)):
            gene_features = embeddings[:, gene_idx*n_shells:(gene_idx+1)*n_shells]
            avg_activation = np.mean(gene_features)
            gene_activations.append(avg_activation)

        gene_indices = range(len(gene_activations))
        ax.bar(gene_indices, gene_activations, color='coral', alpha=0.7)
        ax.set_xlabel('Gene Index')
        ax.set_ylabel('Average Activation')
        ax.set_title(f'Mean Gene Activity\n(First {len(gene_activations)} genes)')
        ax.grid(True, alpha=0.3, axis='y')

        # Rotate x labels if many genes
        if len(gene_activations) > 10:
            ax.tick_params(axis='x', rotation=45)

    fig.suptitle('Feature Statistics Across Different Patch Radii',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()

    return fig


def create_database_summary_plot(db_info: Dict):
    """Create a comprehensive summary of the database."""

    metadata = db_info['metadata']
    config = db_info['config']

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Number of patches per radius
    ax = axes[0, 0]
    radius_counts = metadata['radius'].value_counts().sort_index()
    ax.bar(radius_counts.index, radius_counts.values, color='steelblue', alpha=0.7)
    ax.set_xlabel('Patch Radius (bins)')
    ax.set_ylabel('Number of Patches')
    ax.set_title('Patch Distribution by Radius')
    ax.grid(True, alpha=0.3, axis='y')

    # Add counts on bars
    for i, (radius, count) in enumerate(radius_counts.items()):
        ax.text(radius, count, f'{count:,}', ha='center', va='bottom')

    # Plot 2: Number of patches per sample
    ax = axes[0, 1]
    sample_counts = metadata['sample_id'].value_counts()
    ax.bar(range(len(sample_counts)), sample_counts.values, color='coral', alpha=0.7)
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Number of Patches')
    ax.set_title('Patch Distribution by Sample')
    ax.set_xticks(range(len(sample_counts)))
    ax.set_xticklabels([s[:12] for s in sample_counts.index], rotation=45, ha='right')
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 3: Patch size distribution
    ax = axes[1, 0]
    ax.hist(metadata['n_bins'], bins=30, color='green', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Number of Bins per Patch')
    ax.set_ylabel('Frequency')
    ax.set_title('Patch Size Distribution')
    ax.axvline(metadata['n_bins'].median(), color='red', linestyle='--',
               label=f'Median: {metadata["n_bins"].median():.0f}')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    # Plot 4: Database configuration info
    ax = axes[1, 1]
    ax.axis('off')

    info_text = f"""
DATABASE CONFIGURATION

Resolution: {config['resolution_um']} μm
Total Genes: {len(config['global_genes'])}
Variable Genes: {'Yes' if config['use_variable_genes'] else 'No'}
"""
    if 'variable_gene_indices' in config:
        info_text += f"Num Variable: {len(config['variable_gene_indices'])}\n"

    info_text += f"""
N Shells: {config['n_shells']}
PCA Dims: {config['pca_dims']}

Total Patches: {len(metadata):,}
Radii: {sorted(metadata['radius'].unique())}

Samples: {len(metadata['sample_id'].unique())}
"""

    ax.text(0.1, 0.5, info_text, fontsize=12, verticalalignment='center',
            fontfamily='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.suptitle('Radial Shell Database Summary', fontsize=14, fontweight='bold')
    plt.tight_layout()

    return fig


def main():
    parser = argparse.ArgumentParser(description='Visualize Radial Shell Database')
    parser.add_argument('--database', type=str, required=True,
                       help='Path to database directory')
    parser.add_argument('--resolution', type=int, default=160,
                       help='Resolution to visualize (default: 160)')
    parser.add_argument('--gene-name', type=str, default=None,
                       help='Gene name to visualize (default: first gene)')
    parser.add_argument('--output-dir', type=str, default='database_visualizations',
                       help='Directory to save plots')

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

    # Get available radii
    radii = sorted(metadata['radius'].unique())
    print(f"Available radii: {radii}")

    # Figure 1: Radial shell concept
    print("\nCreating radial shell concept visualization...")
    fig1 = visualize_radial_shell_concept(n_shells=config['n_shells'])
    fig1.savefig(output_dir / 'radial_shell_concept.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_dir / 'radial_shell_concept.png'}")

    # Figure 2: Database summary
    print("\nCreating database summary...")
    fig2 = create_database_summary_plot(db_info)
    fig2.savefig(output_dir / 'database_summary.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_dir / 'database_summary.png'}")

    # Load raw embeddings for each radius
    print("\nLoading raw embeddings...")
    embeddings_by_radius = {}
    for radius in radii:
        try:
            embeddings, patch_ids = load_raw_embeddings(res_dir, radius)
            embeddings_by_radius[radius] = embeddings
            print(f"  Radius {radius}: {embeddings.shape[0]} patches × {embeddings.shape[1]} features")
        except FileNotFoundError as e:
            print(f"  Warning: {e}")

    if not embeddings_by_radius:
        print("\nNo raw embeddings found. Did you build the database?")
        return

    # Figure 3: Feature statistics
    print("\nCreating feature statistics visualization...")
    fig3 = visualize_feature_statistics(embeddings_by_radius, config, n_genes_to_show=20)
    fig3.savefig(output_dir / 'feature_statistics.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_dir / 'feature_statistics.png'}")

    # Figure 4: Single gene encoding
    if args.gene_name:
        gene_idx = config['global_genes'].index(args.gene_name)
        gene_name = args.gene_name
    else:
        gene_idx = 0
        gene_name = config['global_genes'][0]

    print(f"\nCreating single gene visualization for: {gene_name} (index {gene_idx})...")
    fig4 = visualize_single_gene_encoding(
        gene_name=gene_name,
        gene_idx=gene_idx,
        n_shells=config['n_shells'],
        radii=radii,
        embeddings_by_radius=embeddings_by_radius,
        config=config
    )
    fig4.savefig(output_dir / f'gene_encoding_{gene_name}.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_dir / f'gene_encoding_{gene_name}.png'}")

    print(f"\n✓ All visualizations saved to: {output_dir}")
    print("\nGenerated plots:")
    print("  1. radial_shell_concept.png - How radial shell encoding works")
    print("  2. database_summary.png - Database statistics and configuration")
    print("  3. feature_statistics.png - Feature distributions across radii")
    print(f"  4. gene_encoding_{gene_name}.png - Single gene encoding details")

    # Show plots if in interactive mode
    try:
        plt.show()
    except:
        pass


if __name__ == '__main__':
    main()
