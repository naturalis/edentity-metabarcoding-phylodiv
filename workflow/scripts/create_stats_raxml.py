import argparse
import subprocess
from pathlib import Path
import logging
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import tempfile
import os
import re
import random
import sys

def setup_logger(verbosity):
    """Configure logging."""
    logging.basicConfig(
        level=verbosity,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger('raxml_stats')

def process_multiple_subsets(ref_tree, ref_align, n_samples, n_subsets, output_dir, threads, logger):
    """
    Process multiple subsets of the data and collect their RAxML stats.
    
    Parameters:
    :param ref_tree: Path to reference tree
    :param ref_align: Path to reference alignment
    :param n_samples: Number of samples per subset
    :param n_subsets: Number of subsets to process
    :param output_dir: Output directory
    :param threads: Number of threads for RAxML
    :param logger: Logger instance
    
    Returns:
    :return: List of paths to RAxML info files
    """
    raxml_info_files = []
    
    for i in range(n_subsets):
        subset_dir = Path(output_dir) / f"subset_{i+1}"
        subset_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Processing subset {i+1}/{n_subsets}")
        
        # Run subsampling
        subset_align, subset_tree = subsample_and_prune(
            ref_align,
            ref_tree,
            n_samples,
            subset_dir,
            logger
        )
        
        # Run RAxML
        raxml_info = run_raxml(
            subset_tree, 
            subset_align, 
            subset_dir, 
            threads, 
            logger
        )
        
        raxml_info_files.append(raxml_info)
        
    return raxml_info_files

def verify_files_exist(files_dict, logger):
    """Verify that all input files exist."""
    for name, path in files_dict.items():
        path = Path(path)
        if not path.exists():
            logger.error(f"{name} file not found: {path}")
            logger.error(f"Current working directory: {os.getcwd()}")
            raise FileNotFoundError(f"{name} file not found: {path}")
        else:
            logger.info(f"Found {name} file: {path}")
            logger.info(f"File size: {path.stat().st_size} bytes")

def verify_raxml_installed(logger):
    """Verify RAxML is installed and accessible."""
    try:
        result = subprocess.run(['raxmlHPC-PTHREADS', '-v'], 
                              capture_output=True, 
                              text=True)
        logger.info("RAxML is installed and accessible")
        logger.debug(f"RAxML version info: {result.stdout}")
    except FileNotFoundError:
        logger.error("RAxML (raxmlHPC-PTHREADS) not found in PATH")
        raise

def subsample_and_prune(stockholm_file, tree_file, n_samples, output_dir, logger):
    """
    Subsample sequences from Stockholm alignment and prune tree accordingly.
    """
    try:
        # Read full alignment
        logger.info(f"Reading alignment from {stockholm_file}")
        alignment = AlignIO.read(stockholm_file, "stockholm")
        total_seqs = len(alignment)
        
        logger.info(f"Total sequences in alignment: {total_seqs}")
        
        if n_samples >= total_seqs:
            logger.warning(f"Requested sample size {n_samples} >= total sequences {total_seqs}")
            return stockholm_file, tree_file
            
        # Randomly sample sequence indices and get their IDs
        sampled_indices = random.sample(range(total_seqs), n_samples)
        sampled_ids = [alignment[idx].id for idx in sorted(sampled_indices)]
        
        logger.info(f"Selected {len(sampled_ids)} sequences for subsampling")
        
        # Write sampled IDs to temporary file for tree pruning
        id_file = output_dir / "sampled_ids.txt"
        with open(id_file, 'w') as f:
            f.write('\n'.join(sampled_ids))
        logger.info(f"Wrote sampled IDs to {id_file}")
            
        # Initialize tree database
        db_path = output_dir / "reference_tree.db"
        if db_path.exists():
            db_path.unlink()
            
        # Initialize database
        cmd = [
            "megatree-loader",
            "-i", str(tree_file),
            "-d", str(db_path)
        ]
        logger.info(f"Running command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        logger.info("Initialized tree database")
        
        # Extract pruned subtree
        pruned_tree = output_dir / "pruned_subsampled_tree.tre"
        cmd = [
            "megatree-pruner",
            "-d", str(db_path),
            "-i", str(id_file)
        ]
        
        logger.info(f"Running command: {' '.join(cmd)}")
        with open(pruned_tree, 'w') as outfile:
            subprocess.run(cmd, stdout=outfile, check=True)
        logger.info(f"Extracted pruned subtree to {pruned_tree}")
        
        # Create subsampled alignment
        subsampled_align = output_dir / "subsampled_alignment.sto" 
        sampled_alignment = MultipleSeqAlignment([])
        for idx in sorted(sampled_indices):
            sampled_alignment.append(alignment[idx])
            
        AlignIO.write(sampled_alignment, subsampled_align, "stockholm")
        logger.info(f"Created subsampled alignment with {n_samples} sequences")
        
        # Clean up temporary files
        id_file.unlink()
        db_path.unlink()
        logger.info("Cleaned up temporary files")
        
        return subsampled_align, pruned_tree
        
    except Exception as e:
        logger.error(f"Error during subsampling and pruning: {str(e)}")
        raise

def convert_to_phylip(stockholm_file, output_dir, logger):
    """
    Convert Stockholm alignment to PHYLIP format
    """
    try:
        output_path = Path(output_dir) / "reference_alignment.phy"
        logger.info(f"Converting {stockholm_file} to PHYLIP format")
        
        # Read Stockholm and write PHYLIP
        AlignIO.convert(
            stockholm_file,
            "stockholm",
            str(output_path),
            "phylip-relaxed"  # Using relaxed format to handle longer sequence names
        )
        
        logger.info(f"Converted alignment saved to {output_path}")
        return output_path
        
    except Exception as e:
        logger.error(f"Error converting alignment format: {str(e)}")
        raise

def run_raxml(ref_tree, ref_align, output_dir, threads, logger):
    """Run RAxML to generate stats file"""
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert alignment if it's in Stockholm format
    if str(ref_align).endswith(('.sto', '.stockholm')):
        ref_align = convert_to_phylip(ref_align, output_dir, logger)
    
    # Create temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        # RAxML command
        cmd = [
            "raxmlHPC-PTHREADS",
            "-T", str(threads),
            "-f", "e",
            "-t", str(ref_tree),
            "-s", str(ref_align),
            "-m", "GTRGAMMA",
            "-n", "stats",
            "-w", temp_dir,
            "-p", "12345"
        ]
        
        try:
            logger.info(f"Running RAxML command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd, 
                check=True, 
                capture_output=True, 
                text=True
            )
            
            # Log RAxML output
            logger.debug("RAxML stdout:")
            logger.debug(result.stdout)
            if result.stderr:
                logger.warning("RAxML stderr:")
                logger.warning(result.stderr)
            
            # Copy RAxML output files to final destination
            raxml_info = Path(temp_dir) / "RAxML_info.stats"
            final_info = output_dir / "RAxML_info.stats"
            
            if raxml_info.exists():
                with open(raxml_info, 'r') as src, open(final_info, 'w') as dst:
                    dst.write(src.read())
                logger.info("RAxML completed successfully")
                return final_info
            else:
                raise FileNotFoundError("Expected RAxML output file not found")
                
        except subprocess.CalledProcessError as e:
            logger.error(f"RAxML failed with error:\nstdout: {e.stdout}\nstderr: {e.stderr}")
            raise
        except Exception as e:
            logger.error(f"Error running RAxML: {str(e)}")
            raise

def extract_parameters(raxml_info_file, logger):
    """Extract parameters from a single RAxML info file."""
    with open(raxml_info_file) as f:
        raxml_output = f.read()
    
    patterns = {
        'alpha': r'alpha:\s+([\d.]+)',
        'tree_length': r'Tree-Length:\s+([\d.]+)',
        'rates': (r'rate A <-> C:\s+([\d.]+)\s*'
                 r'rate A <-> G:\s+([\d.]+)\s*'
                 r'rate A <-> T:\s+([\d.]+)\s*'
                 r'rate C <-> G:\s+([\d.]+)\s*'
                 r'rate C <-> T:\s+([\d.]+)\s*'
                 r'rate G <-> T:\s+([\d.]+)'),
        'freqs': (r'freq pi\(A\):\s+([\d.]+)\s*'
                 r'freq pi\(C\):\s+([\d.]+)\s*'
                 r'freq pi\(G\):\s+([\d.]+)\s*'
                 r'freq pi\(T\):\s+([\d.]+)'),
        'likelihood': r'Final GAMMA\s+likelihood:\s+([-\d.]+)',
    }
    
    values = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, raxml_output)
        if match:
            values[key] = [float(x) for x in match.groups()]
        else:
            logger.warning(f"Could not find pattern for {key} in {raxml_info_file}")
            
    return values

def average_parameters(parameter_list, logger):
    """Average parameters from multiple runs."""
    if not parameter_list:
        raise ValueError("No parameters to average")
        
    # Initialize with first set of parameters
    avg_params = {}
    for key in parameter_list[0].keys():
        avg_params[key] = [0.0] * len(parameter_list[0][key])
        
    # Sum all parameters
    for params in parameter_list:
        for key in params:
            for i, value in enumerate(params[key]):
                avg_params[key][i] += value
                
    # Calculate averages
    n = len(parameter_list)
    for key in avg_params:
        avg_params[key] = [x/n for x in avg_params[key]]
        
    return avg_params

def format_raxml_output(raxml_info_files, template_file, output_file, logger):
    """Format averaged RAxML output according to template."""
    try:
        # Read template
        with open(template_file) as f:
            template = f.read()
            
        # Extract parameters from each file
        all_params = []
        for info_file in raxml_info_files:
            params = extract_parameters(info_file, logger)
            all_params.append(params)
            
        # Average parameters
        avg_params = average_parameters(all_params, logger)
        
        # Format template with averaged values
        formatted_output = template.format(
            alpha=avg_params['alpha'][0],
            tree_length=avg_params['tree_length'][0],
            rate_ac=avg_params['rates'][0],
            rate_ag=avg_params['rates'][1],
            rate_at=avg_params['rates'][2],
            rate_cg=avg_params['rates'][3],
            rate_ct=avg_params['rates'][4],
            rate_gt=avg_params['rates'][5],
            freq_a=avg_params['freqs'][0],
            freq_c=avg_params['freqs'][1],
            freq_g=avg_params['freqs'][2],
            freq_t=avg_params['freqs'][3],
            likelihood=avg_params['likelihood'][0]
        )
        
        # Write formatted output
        with open(output_file, 'w') as f:
            f.write(formatted_output)
            
        logger.info(f"Formatted averaged RAxML output written to {output_file}")
        
    except Exception as e:
        logger.error(f"Error formatting RAxML output: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Generate RAxML stats file')
    parser.add_argument('-t', '--tree', required=True, help='Reference tree file')
    parser.add_argument('-s', '--align', required=True, help='Reference alignment file')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument('--template', required=True, help='RAxML template file')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    parser.add_argument('--subsample', type=int, help='Number of sequences to subsample')
    parser.add_argument('--n-subsets', type=int, default=3, help='Number of subsets to process')
    parser.add_argument('-v', '--verbosity', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       help='Logging verbosity')
    args = parser.parse_args()
    
    logger = setup_logger(args.verbosity)
    
    try:
        logger.info("Starting RAxML stats generation...")
        
        # Verify inputs and RAxML installation
        verify_files_exist({
            'Tree': args.tree,
            'Alignment': args.align,
            'Template': args.template
        }, logger)
        verify_raxml_installed(logger)
        
        # Create output directory
        output_dir = Path(args.outdir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if args.subsample:
            # Process multiple subsets
            raxml_info_files = process_multiple_subsets(
                args.tree,
                args.align,
                args.subsample,
                args.n_subsets,
                output_dir,
                args.threads,
                logger
            )
        else:
            # Single run without subsampling
            raxml_info = run_raxml(
                args.tree,
                args.align,
                output_dir,
                args.threads,
                logger
            )
            raxml_info_files = [raxml_info]
        
        # Format output with averaged parameters
        output_file = output_dir / "RAxML_info_formatted.stats"
        format_raxml_output(raxml_info_files, args.template, output_file, logger)
        
        logger.info("RAxML stats generation completed successfully")
        
    except Exception as e:
        logger.error(f"Error during RAxML stats generation: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()