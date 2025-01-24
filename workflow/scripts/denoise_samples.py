import util
import argparse
from pathlib import Path
import subprocess
import yaml
import time
import re
from threading import Thread
import queue

def monitor_dada2_progress(process, logger, progress_queue):
    """Monitor DADA2 progress by parsing stdout"""
    current_sample = None
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
            
        line = line.decode().strip()
        if line:
            # Sample progress detection
            if "Processing sample" in line:
                current_sample = re.search(r"Processing sample '(.+)'", line).group(1)
                logger.info(f"Started processing sample: {current_sample}")
                
            # Error model progress
            elif "Learning error rates for forward reads:" in line:
                logger.info("Learning forward read error rates")
            elif "Learning error rates for reverse reads:" in line:
                logger.info("Learning reverse read error rates")
            elif "Conducting initial pass..." in line:
                logger.info("Initial error rate estimation pass")
            elif "Converging..." in line:
                logger.info("Converging error rates")


            elif "Inferring ASV sequences" in line:
                logger.info(f"Inferring sequences for sample: {current_sample}")
                
            # Merging progress    
            elif "Merging forward and reverse reads" in line:
                logger.info(f"Merging paired reads for sample: {current_sample}")
                
            # Chimera detection
            elif "Identifying chimeras" in line:
                logger.info("Running chimera detection")
                
            progress_queue.put(line)

def run_dada2_with_progress(cmd, logger):
    """Run DADA2 with progress monitoring"""
    progress_queue = queue.Queue()
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=False
    )
    
    # Start progress monitoring thread
    monitor_thread = Thread(
        target=monitor_dada2_progress,
        args=(process, logger, progress_queue)
    )
    monitor_thread.start()
    
    # Wait for completion
    process.wait()
    monitor_thread.join()
    
    # Check for errors
    if process.returncode != 0:
        error = process.stderr.read().decode()
        raise subprocess.CalledProcessError(process.returncode, cmd, error)
        
    return process.returncode

def build_dada2_command(input_path, outputs, params):
    """Build DADA2-specific command for sequence denoising."""
    cmd = [
        "qiime", "dada2", "denoise-paired",
        "--i-demultiplexed-seqs", str(input_path),
        "--o-table", str(outputs['table']),
        "--o-representative-sequences", str(outputs['rep_seqs']),
        "--o-denoising-stats", str(outputs['stats']),
        "--p-trunc-len-f", str(params['trunc_len']),
        "--p-trunc-len-r", str(params['trunc_len']),
        "--p-trim-left-f", str(params['trim_left']),
        "--p-trim-left-r", str(params['trim_left']),
        "--p-max-ee-f", str(params['max_ee']),
        "--p-max-ee-r", str(params['max_ee']),
        "--p-trunc-q", str(params['trunc_q']),
        "--p-chimera-method", str(params['chimera_method']),
        "--p-min-fold-parent-over-abundance", str(params['min_fold_parent_over_abundance']),
        "--p-n-threads", str(params['n_threads']),
        "--verbose"
    ]
    return cmd

def build_deblur_command(input_path, outputs, params):
    """
    Build Deblur-specific command for sequence denoising.
    
    Parameters:
    :param input_path: Input demultiplexed sequences path
    :param outputs: Dictionary of output paths for table, rep_seqs, and stats
    :param params: Deblur parameters dictionary
    
    Returns:
    :return: List of command arguments for subprocess
    """
    cmd = [
        "qiime", "deblur", "denoise-16S",
        "--i-demultiplexed-seqs", str(input_path),
        "--o-table", str(outputs['table']),
        "--o-representative-sequences", str(outputs['rep_seqs']),
        "--o-stats", str(outputs['stats']),
        "--p-trim-length", str(params['trim_length']),
        "--p-left-trim-len", str(params['left_trim_len']),
        "--p-min-reads", str(params['min_reads']),
        "--p-min-size", str(params['min_size']),
        "--p-jobs-to-start", str(params['jobs_to_start'])
    ]
    
    if params.get('hashed_feature_ids', False):
        cmd.append("--p-hashed-feature-ids")
    
    return cmd

def generate_visualizations(outputs, logger):
    """
    Generate QIIME2 visualizations for output artifacts.
    
    Parameters:
    :param outputs: Dictionary of output paths
    :param logger: Logger instance
    """
    for output_type, path in outputs.items():
        if output_type == 'stats':
            continue
            
        vis_path = str(path).replace('.qza', '.qzv')
        try:
            if output_type == 'table':
                subprocess.run(["qiime", "feature-table", "summarize",
                              "--i-table", str(path),
                              "--o-visualization", vis_path], check=True)
            elif output_type == 'rep_seqs':
                subprocess.run(["qiime", "feature-table", "tabulate-seqs",
                              "--i-data", str(path),
                              "--o-visualization", vis_path], check=True)
                
            logger.info(f"Generated visualization for {output_type}: {vis_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error generating visualization for {output_type}: {str(e)}")
            raise

def validate_dada2_params(params):
    """
    Validate DADA2 parameters.
    
    Parameters:
    :param params: Dictionary of DADA2 parameters
    
    Returns:
    :return: Validated parameters dictionary
    """
    required_params = {
        'trunc_len': int,
        'trim_left': int,
        'max_ee': (int, float),
        'trunc_q': int,
        'chimera_method': str,
        'min_fold_parent_over_abundance': (int, float),
        'n_threads': int
    }
    
    for param, param_type in required_params.items():
        if param not in params:
            raise KeyError(f"Missing required DADA2 parameter: {param}")
        if not isinstance(params[param], param_type):
            raise TypeError(f"Parameter {param} should be of type {param_type}")
    
    return params

def denoise(input_path, output_dir, method, params, logger):
    """Run denoising with specified method and progress monitoring"""
    if method not in ['dada2', 'deblur']:
        raise ValueError(f"Unsupported denoise method: {method}")
    
    if method == 'dada2':
        params = validate_dada2_params(params)
        
    outputs = {
        'table': output_dir / f'table_{method}.qza',
        'rep_seqs': output_dir / f'rep_seqs_{method}.qza',
        'stats': output_dir / f'stats_{method}.qza'
    }
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input artifact not found: {input_path}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Running {method.upper()} denoising")
    
    start_time = time.time()
    
    try:
        if method == 'dada2':
            cmd = build_dada2_command(input_path, outputs, params)
            run_dada2_with_progress(cmd, logger)
        else:
            cmd = build_deblur_command(input_path, outputs, params)
            subprocess.run(cmd, check=True)
            
        elapsed_time = time.time() - start_time
        logger.info(f"{method.upper()} denoising completed in {elapsed_time:.1f} seconds")
        
        generate_visualizations(outputs, logger)
        logger.info(f"All outputs saved to: {output_dir}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running {method.upper()}: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Denoise sequences using DADA2 or Deblur.')
    parser.add_argument('-i', '--input', required=True, help='Input demultiplexed sequences')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-c', '--config', required=True, help='Path to configuration YAML file')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    logger = util.get_formatted_logger('denoise_sequences', args.verbosity)
    
    try:
        with open(args.config) as f:
            config = yaml.safe_load(f)
        
        method = config['denoise_method'].lower()
        params = config[f'{method}_params']
        
        denoise(
            Path(args.input),
            Path(args.output_dir),
            method,
            params,
            logger
        )
        logger.info("Denoising completed successfully")
    except Exception as e:
        logger.error(f"Error during denoising: {str(e)}")
        raise