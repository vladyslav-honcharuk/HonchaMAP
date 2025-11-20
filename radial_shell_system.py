#!/usr/bin/env python3
"""
Radial Shell System - Main CLI interface for building and searching.

This is the main entry point for:
- Building multi-resolution databases
- Searching for similar spatial regions
- Managing database metadata

Author: XeniumMundus Project
"""

import argparse
import logging
from pathlib import Path
import sys
import json
from typing import Optional, List

from database_builder import (
    MultiResolutionDatabaseBuilder,
    discover_samples
)
from similarity_search import MultiResolutionSpatialSearch


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('radial_shell_system.log')
        ]
    )


def build_database(args):
    """Build multi-resolution database."""
    logger = logging.getLogger(__name__)
    logger.info("=" * 70)
    logger.info("BUILDING RADIAL SHELL ENCODING DATABASE")
    logger.info("=" * 70)

    # Discover samples
    logger.info(f"\nDiscovering samples in: {args.data_dir}")
    samples_by_resolution = discover_samples(
        Path(args.data_dir),
        min_resolutions=1
    )

    if not samples_by_resolution:
        logger.error("No samples found!")
        return 1

    # Filter to requested resolutions
    if args.resolutions:
        samples_by_resolution = {
            res: samples
            for res, samples in samples_by_resolution.items()
            if res in args.resolutions
        }

    if not samples_by_resolution:
        logger.error(f"No samples found for resolutions: {args.resolutions}")
        return 1

    # Build database
    builder = MultiResolutionDatabaseBuilder(
        output_dir=Path(args.output_dir),
        n_shells=args.n_shells,
        pca_dims=args.pca_dims,
        use_variable_genes=args.use_variable_genes,
        variable_gene_threshold=args.variable_gene_threshold,
        batch_size=args.batch_size
    )

    builder.build_database(samples_by_resolution)

    logger.info("\n" + "=" * 70)
    logger.info("DATABASE BUILD COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Output directory: {args.output_dir}")

    return 0


def search_database(args):
    """Search for similar regions."""
    logger = logging.getLogger(__name__)
    logger.info("=" * 70)
    logger.info("SEARCHING FOR SIMILAR SPATIAL REGIONS")
    logger.info("=" * 70)

    # Initialize search engine
    searcher = MultiResolutionSpatialSearch(Path(args.database))

    # Determine resolution
    resolution_um = args.resolution
    if resolution_um is None:
        # Use first available
        resolution_um = searcher.resolutions[0]
        logger.info(f"Auto-selected resolution: {resolution_um}μm")

    # Load sample path
    sample_path = Path(args.sample)
    if not sample_path.exists():
        logger.error(f"Sample not found: {sample_path}")
        return 1

    # Perform search based on input type
    if args.selection_file:
        # Search from selection bounds file
        with open(args.selection_file, 'r') as f:
            selection_bounds = json.load(f)

        results = searcher.search_from_selection(
            sample_path=sample_path,
            selection_bounds=selection_bounds,
            resolution_um=resolution_um,
            k=args.top_k,
            min_bins=args.min_bins,
            max_bins=args.max_bins,
            bin_size=args.bin_size
        )

    else:
        # Search from coordinates
        if args.x_center is None or args.y_center is None or args.radius is None:
            logger.error(
                "Must specify either --selection-file or "
                "(--x-center, --y-center, --radius)"
            )
            return 1

        results = searcher.search_from_coordinates(
            sample_path=sample_path,
            x_center=args.x_center,
            y_center=args.y_center,
            radius_physical=args.radius,
            resolution_um=resolution_um,
            k=args.top_k,
            min_bins=args.min_bins,
            max_bins=args.max_bins,
            bin_size=args.bin_size
        )

    # Display results
    if len(results) == 0:
        logger.warning("No results found")
        return 0

    logger.info("\n" + "=" * 70)
    logger.info(f"TOP {min(args.top_k, len(results))} RESULTS")
    logger.info("=" * 70)

    # Print formatted results
    for idx, row in results.head(args.top_k).iterrows():
        logger.info(
            f"\n{idx+1}. Sample: {row['sample_id']}"
        )
        logger.info(f"   Similarity: {row['similarity']:.4f}")
        logger.info(f"   Center: ({row['center_x']:.1f}, {row['center_y']:.1f})")
        logger.info(f"   Radius: {row['radius']} bins")
        logger.info(f"   Bins: {row['n_bins']}")

    # Save results if requested
    if args.output:
        output_path = Path(args.output)
        results.to_csv(output_path, index=False)
        logger.info(f"\nResults saved to: {output_path}")

    return 0


def info_database(args):
    """Display database information."""
    logger = logging.getLogger(__name__)
    db_path = Path(args.database)

    if not db_path.exists():
        logger.error(f"Database not found: {db_path}")
        return 1

    logger.info("=" * 70)
    logger.info(f"DATABASE INFO: {db_path}")
    logger.info("=" * 70)

    # Detect resolutions
    resolutions = []
    for item in db_path.iterdir():
        if item.is_dir() and item.name.endswith('um'):
            try:
                res = int(item.name[:-2])
                resolutions.append(res)
            except ValueError:
                continue

    resolutions.sort()

    if not resolutions:
        logger.warning("No valid resolutions found")
        return 1

    # Display info for each resolution
    for res in resolutions:
        logger.info(f"\n{res}μm Resolution:")
        res_dir = db_path / f"{res}um"

        # Metadata
        metadata_file = res_dir / "patches_metadata.parquet"
        if metadata_file.exists():
            import pandas as pd
            metadata = pd.read_parquet(metadata_file)
            logger.info(f"  Total patches: {len(metadata):,}")
            logger.info(f"  Samples: {metadata['sample_id'].nunique()}")

            # Stats by radius
            for radius in sorted(metadata['radius'].unique()):
                count = (metadata['radius'] == radius).sum()
                logger.info(f"    Radius {radius}: {count:,} patches")

        # PCA model
        pca_file = res_dir / f"pca_model_{res}um.pkl"
        if pca_file.exists():
            import pickle
            with open(pca_file, 'rb') as f:
                pca = pickle.load(f)
            variance = pca.explained_variance_ratio_.sum()
            logger.info(f"  PCA dimensions: {pca.n_components_}")
            logger.info(f"  Variance explained: {variance:.1%}")

        # Storage size
        total_size = sum(f.stat().st_size for f in res_dir.rglob('*') if f.is_file())
        logger.info(f"  Storage: {total_size / (1024**3):.2f} GB")

    return 0


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Radial Shell Encoding System for Spatial Transcriptomics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Build database from all samples
  python radial_shell_system.py build --data-dir /path/to/xenium/data --output-dir spatial_db

  # Build for specific resolutions
  python radial_shell_system.py build --data-dir /path/to/data --output-dir spatial_db --resolutions 20 40 80

  # Search from coordinates
  python radial_shell_system.py search --database spatial_db --sample /path/to/sample \\
    --x-center 5000 --y-center 5000 --radius 500 --resolution 40

  # Search from selection file
  python radial_shell_system.py search --database spatial_db --sample /path/to/sample \\
    --selection-file selection.json --resolution 40

  # Display database info
  python radial_shell_system.py info --database spatial_db
        """
    )

    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    subparsers = parser.add_subparsers(dest='command', help='Command to execute')

    # Build command
    build_parser = subparsers.add_parser('build', help='Build database')
    build_parser.add_argument('--data-dir', required=True, type=str,
                             help='Directory containing Xenium samples')
    build_parser.add_argument('--output-dir', required=True, type=str,
                             help='Output directory for database')
    build_parser.add_argument('--resolutions', type=int, nargs='+',
                             help='Resolutions to process (default: all available)')
    build_parser.add_argument('--n-shells', type=int, default=5,
                             help='Number of concentric shells (default: 5)')
    build_parser.add_argument('--pca-dims', type=int, default=256,
                             help='PCA compression dimensions (default: 256)')
    build_parser.add_argument('--use-variable-genes',
                             type=lambda x: x.lower() in ('true', '1', 'yes'),
                             default=True,
                             help='Use only spatially variable genes (default: True)')
    build_parser.add_argument('--variable-gene-threshold', type=float, default=-2.0,
                             help='Haystack logpval_adj threshold for variable genes (default: -2.0 for p_adj<0.01)')
    build_parser.add_argument('--batch-size', type=int, default=10000,
                             help='Batch size for PCA fitting (default: 10000)')

    # Search command
    search_parser = subparsers.add_parser('search', help='Search for similar regions')
    search_parser.add_argument('--database', required=True, type=str,
                              help='Path to database directory')
    search_parser.add_argument('--sample', required=True, type=str,
                              help='Path to query sample directory')
    search_parser.add_argument('--resolution', type=int,
                              help='Resolution to search (default: auto-detect)')
    search_parser.add_argument('--bin-size', type=int,
                              help='Bin size for zarr file (default: same as resolution)')

    # Search query options (mutually exclusive)
    query_group = search_parser.add_argument_group('query options')
    query_group.add_argument('--x-center', type=float,
                            help='X coordinate of query center (bin units)')
    query_group.add_argument('--y-center', type=float,
                            help='Y coordinate of query center (bin units)')
    query_group.add_argument('--radius', type=float,
                            help='Query radius (micrometers)')
    query_group.add_argument('--selection-file', type=str,
                            help='JSON file with selection bounds')

    # Search result options
    search_parser.add_argument('--top-k', type=int, default=100,
                              help='Number of results to return (default: 100)')
    search_parser.add_argument('--min-bins', type=int,
                              help='Minimum bins per result')
    search_parser.add_argument('--max-bins', type=int,
                              help='Maximum bins per result')
    search_parser.add_argument('--output', type=str,
                              help='Output file for results (CSV)')

    # Info command
    info_parser = subparsers.add_parser('info', help='Display database information')
    info_parser.add_argument('--database', required=True, type=str,
                            help='Path to database directory')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 1

    # Setup logging
    setup_logging(args.verbose)

    # Execute command
    if args.command == 'build':
        return build_database(args)
    elif args.command == 'search':
        return search_database(args)
    elif args.command == 'info':
        return info_database(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
