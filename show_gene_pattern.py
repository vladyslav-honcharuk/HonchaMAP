#!/usr/bin/env python3
"""
Show Gene Pattern - Super simple visualization of what's stored in the database

This shows you exactly what the database contains for one gene across shells.
"""

import numpy as np
import faiss
import json
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Show what is stored for one gene')
    parser.add_argument('--database', default='spatial_database/no_pca_db/160um')
    parser.add_argument('--gene', default='ACTA2')
    parser.add_argument('--radius', type=int, default=5, help='Patch radius (3, 5, or 8)')
    args = parser.parse_args()

    db_path = Path(args.database)

    # Load config to get gene list
    with open(db_path / 'config.json') as f:
        config = json.load(f)

    genes = config['global_genes']
    n_shells = config['n_shells']

    if args.gene not in genes:
        print(f"Gene {args.gene} not found!")
        print(f"Available genes: {genes[:20]}...")
        return

    gene_idx = genes.index(args.gene)
    print(f"Gene: {args.gene} (index {gene_idx})")
    print(f"Shells: {n_shells}")

    # Load FAISS index (contains all the encoded patches)
    index_file = db_path / 'faiss_indices' / f'faiss_r{args.radius}.index'
    index = faiss.read_index(str(index_file))

    print(f"\nDatabase has {index.ntotal} patches at radius {args.radius}")
    print(f"Each patch has {index.d} features total")
    print(f"  = {index.d // n_shells} genes × {n_shells} shells")

    # Extract all patches
    patches = np.zeros((index.ntotal, index.d), dtype=np.float32)
    for i in range(index.ntotal):
        patches[i] = index.reconstruct(i)

    # Extract features for this gene (5 features = 5 shells)
    gene_start = gene_idx * n_shells
    gene_end = gene_start + n_shells
    gene_features = patches[:, gene_start:gene_end]  # Shape: (n_patches, n_shells)

    print(f"\nExtracted {args.gene} features: {gene_features.shape[0]} patches × {gene_features.shape[1]} shells")

    # Create a simple plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Show first 50 patches as heatmap
    n_show = min(50, gene_features.shape[0])
    im = ax1.imshow(gene_features[:n_show].T, aspect='auto', cmap='YlOrRd')
    ax1.set_xlabel('Patch Number')
    ax1.set_ylabel('Shell (1=center, 5=edge)')
    ax1.set_yticks(range(n_shells))
    ax1.set_yticklabels([f'{i+1}' for i in range(n_shells)])
    ax1.set_title(f'{args.gene} in {n_show} patches\n(Yellow=high expression)')
    plt.colorbar(im, ax=ax1, label='Expression')

    # Right: Show average across all patches
    avg_per_shell = gene_features.mean(axis=0)
    std_per_shell = gene_features.std(axis=0)

    ax2.bar(range(n_shells), avg_per_shell, yerr=std_per_shell,
            color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd'][:n_shells],
            alpha=0.7, capsize=5)
    ax2.set_xlabel('Shell Number')
    ax2.set_ylabel('Average Expression')
    ax2.set_title(f'{args.gene} Average Pattern\n(averaged over all {gene_features.shape[0]} patches)')
    ax2.set_xticks(range(n_shells))
    ax2.set_xticklabels([f'{i+1}\n(center)' if i == 0 else f'{i+1}\n(edge)' if i == n_shells-1 else f'{i+1}'
                          for i in range(n_shells)])
    ax2.grid(axis='y', alpha=0.3)

    # Add text showing what this means
    fig.text(0.5, 0.02,
             'Left: Each column is one patch, each row is one shell | Right: Average shows if gene is higher in center or edge',
             ha='center', fontsize=10, style='italic')

    plt.suptitle(f'What is stored in the database for {args.gene}', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])

    output = f'{args.gene}_simple.png'
    plt.savefig(output, dpi=100, bbox_inches='tight')
    print(f"\nSaved: {output}")

    # Print interpretation
    print("\n" + "="*60)
    print("INTERPRETATION:")
    print("="*60)
    if avg_per_shell[0] > avg_per_shell[-1] * 1.2:
        print(f"→ {args.gene} is HIGHER in the CENTER of patches")
    elif avg_per_shell[-1] > avg_per_shell[0] * 1.2:
        print(f"→ {args.gene} is HIGHER at the EDGES of patches")
    else:
        print(f"→ {args.gene} has NO CLEAR radial pattern")

    print(f"\nShell values: {[f'{v:.2f}' for v in avg_per_shell]}")
    print("="*60)


if __name__ == '__main__':
    main()
