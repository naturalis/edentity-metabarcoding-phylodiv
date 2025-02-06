import util
import argparse
from pathlib import Path
import subprocess
import sys
import json
from ete3 import Tree

def check_placements(jplace_file, logger):
    """
    Check placement file for content and validity.
    
    Parameters:
    :param jplace_file: Path to jplace file
    :param logger: Logger instance
    
    Returns:
    :return: number of placements found
    """
    try:
        with open(jplace_file) as f:
            data = json.load(f)
            
        num_placements = len(data.get('placements', []))
        logger.info(f"Found {num_placements} placements")
            
        return num_placements
        
    except Exception as e:
        logger.error(f"Error checking placement file: {str(e)}")
        raise

def convert_to_newick(xml_file, output_file, logger):
    """
    Convert phyloXML tree to Newick format using ete3.
    
    Parameters:
    :param xml_file: Path to input XML tree file
    :param output_file: Path for output Newick file
    :param logger: Logger instance
    """
    try:
        # Read XML tree
        tree = Tree(str(xml_file), format=1)
        
        # Write Newick format
        tree.write(outfile=str(output_file), format=0)  # format=0 for Newick
        logger.info(f"Converted {xml_file} to Newick format: {output_file}")
        
    except Exception as e:
        logger.error(f"Error converting tree to Newick: {str(e)}")
        raise

def graft_samples(jplace_file, output_dir, logger):
    """
    Use guppy to graft placed samples onto the reference tree and convert to Newick format.
    
    Parameters:
    :param jplace_file: Path to jplace file from pplacer
    :param output_dir: Output directory path
    :param logger: Logger instance
    
    Returns:
    :return: Dictionary containing paths to generated output files
    """
    try:
        # Create output directory if it doesn't exist
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # First check the placements file
        num_placements = check_placements(jplace_file, logger)
        if num_placements == 0:
            raise ValueError("No placements found in jplace file")
        
        # Define output files
        output_paths = {
            'tog_tree_xml': output_dir / 'grafted_tree.xml',
            'tog_tree_nwk': output_dir / 'grafted_tree.nwk'
        }
        
        # Create tree with all placements as pendant edges
        logger.info("Creating grafted tree with placements...")
        cmd = [
            "guppy", "tog",
            "-o", str(output_paths['tog_tree_xml']),
            str(jplace_file)
        ]
        subprocess.run(cmd, check=True)
        
        # Convert to Newick
        convert_to_newick(output_paths['tog_tree_xml'], output_paths['tog_tree_nwk'], logger)
        
        return output_paths
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running guppy command: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error during tree grafting: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Graft placed samples onto reference tree using guppy.')
    parser.add_argument('-i', '--input', required=True, help='Input jplace file from pplacer')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('-v', '--verbosity', default='INFO', help='Logging verbosity level')
    args = parser.parse_args()
    
    # Setup logging
    logger = util.get_formatted_logger('graft_samples', args.verbosity)
    
    try:
        # Validate input file
        if not Path(args.input).exists():
            raise FileNotFoundError(f"Input file not found: {args.input}")
        
        # Graft samples and create tree
        output_files = graft_samples(
            Path(args.input),
            Path(args.output_dir),
            logger
        )
        
        # Log output file locations
        logger.info("Generated Newick format tree:")
        logger.info(f"- {output_files['tog_tree_nwk']}")
        
        logger.info("Sample grafting completed successfully")
        
    except Exception as e:
        logger.error(f"Error during sample grafting: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()