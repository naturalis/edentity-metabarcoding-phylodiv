# prune_tree.py
import util
import argparse
from pathlib import Path
import subprocess
import sys

def load_reference_ids(ref_ids_file, logger):
    """
    Load reference IDs from BLAST filtering results.
    
    Parameters:
    :param ref_ids_file: Path to file containing reference IDs
    :param logger: Logger instance
    
    Returns:
    :return: Set of reference IDs
    """
    try:
        with open(ref_ids_file) as f:
            ref_ids = set(line.strip() for line in f)
        logger.info(f"Loaded {len(ref_ids)} reference IDs")
        return ref_ids
    except Exception as e:
        logger.error(f"Error loading reference IDs: {str(e)}")
        raise

def initialize_database(tree_file, output_dir, logger):
    """
    Initialize sqlite database for the reference tree.
    
    Parameters:
    :param tree_file: Path to reference tree file
    :param output_dir: Output directory for database
    :param logger: Logger instance
    
    Returns:
    :return: Path to created database
    """
    try:
        db_path = output_dir / "reference_tree.db"
        
        # Remove existing database if present
        if db_path.exists():
            db_path.unlink()
        
        cmd = [
            "megatree-loader",
            "-i", str(tree_file),
            "-d", str(db_path)
        ]
        
        subprocess.run(cmd, check=True)
        logger.info(f"Initialized tree database at {db_path}")
        return db_path
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error initializing tree database: {str(e)}")
        raise

def extract_subtree(db_path, ref_ids_file, output_file, logger):
    """
    Extract subtree containing specified taxa using megatree-pruner.

    Parameters:
    :param db_path: Path to tree database
    :param ref_ids_file: Path to file containing reference IDs
    :param output_file: Path for output tree file
    :param logger: Logger instance
    """
    try:
        cmd = [
            "megatree-pruner",
            "-d", str(db_path),
            "-i", str(ref_ids_file)
        ]

        # Run command and redirect output to the specified file
        with open(output_file, 'w') as outfile:
            subprocess.run(cmd, stdout=outfile, check=True)
        
        logger.info(f"Extracted subtree to {output_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error extracting subtree: {str(e)}")
        raise


def main():
    parser = argparse.ArgumentParser(description='Extract subtree based on BLAST results.')
    parser.add_argument('-t', '--tree', required=True, help='Input reference tree file')
    parser.add_argument('-r', '--ref_ids', required=True, help='Reference IDs from BLAST filtering')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()
    
    logger = util.get_formatted_logger('prune_tree', args.verbosity)
    
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize tree database
        db_path = initialize_database(Path(args.tree), output_dir, logger)
        
        # Extract subtree
        output_tree = output_dir / "pruned_reference_tree.tre"
        extract_subtree(db_path, Path(args.ref_ids), output_tree, logger)
        
        logger.info("Tree pruning completed successfully")
        
    except Exception as e:
        logger.error(f"Error during tree pruning: {str(e)}")
        raise

if __name__ == "__main__":
    main()