# prune_tree.py
import util
import argparse
from pathlib import Path
import subprocess
import sys
from Bio import AlignIO

def extract_reference_ids_from_stockholm(stockholm_file, output_file, logger):
    """
    Extract reference sequence IDs from Stockholm alignment file.
    
    Parameters:
    :param stockholm_file: Path to Stockholm alignment file
    :param output_file: Path to write reference IDs
    :param logger: Logger instance
    
    Returns:
    :return: Set of reference IDs and path to output file
    """
    try:
        # Read Stockholm alignment
        logger.info(f"Reading Stockholm alignment from {stockholm_file}")
        alignment = AlignIO.read(stockholm_file, "stockholm")
        
        # Extract unique sequence IDs
        ref_ids = set(record.id for record in alignment)
        logger.info(f"Found {len(ref_ids)} unique reference IDs")
        
        # Write IDs to file
        with open(output_file, 'w') as f:
            for ref_id in sorted(ref_ids):
                f.write(f"{ref_id}\n")
        logger.info(f"Wrote reference IDs to {output_file}")
        
        return ref_ids, output_file
        
    except Exception as e:
        logger.error(f"Error extracting reference IDs: {str(e)}")
        raise

def load_reference_ids(ref_ids_file, logger):
    """
    Load reference IDs from file.
    
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
    parser = argparse.ArgumentParser(description='Extract subtree based on reference IDs from Stockholm alignment.')
    parser.add_argument('-t', '--tree', required=True, help='Input reference tree file')
    parser.add_argument('-a', '--alignment', required=True, help='Stockholm alignment file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()
    
    logger = util.get_formatted_logger('prune_tree', args.verbosity)
    
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract reference IDs from Stockholm file
        ref_ids_file = output_dir / "reference_ids.txt"
        _, ref_ids_file = extract_reference_ids_from_stockholm(
            Path(args.alignment),
            ref_ids_file,
            logger
        )
        
        # Initialize tree database
        db_path = initialize_database(Path(args.tree), output_dir, logger)
        
        # Extract subtree
        output_tree = output_dir / "pruned_reference_tree.tre"
        extract_subtree(db_path, ref_ids_file, output_tree, logger)
        
        # Clean up temporary files
        db_path.unlink()
        
        logger.info("Tree pruning completed successfully")
        
    except Exception as e:
        logger.error(f"Error during tree pruning: {str(e)}")
        raise

if __name__ == "__main__":
    main()