import util
import argparse
import pandas as pd
from pathlib import Path
import re
import subprocess


def parse_fastq_filename(filename):
    """
    Parse sample information from fastq filename.
    
    Parameters:
    :param filename: Input fastq filename to parse
    
    Returns:
    :return: tuple of (sample_id, read) or (None, None) if parsing fails 
    """
    base = Path(filename).stem.replace('.fastq', '').replace('.fq', '')
    
    if 'R1' in base:
        read = 'R1'
    elif 'R2' in base:
        read = 'R2' 
    else:
        return None, None

    sample_pattern = r'(.+?)_(?:R1|R2)'
    match = re.search(sample_pattern, base)
    if match:
        return match.group(1), read
        
    return None, None


def create_manifest(input_dir, output_file, logger):
    """
    Create a QIIME 2 manifest file from a directory of FastQ files.
    
    Parameters:
    :param input_dir: Directory containing fastq files (Path object)
    :param output_file: Output manifest file path (Path object)
    :param logger: Logger instance
    
    Returns:
    :return: Path to created manifest file
    
    Raises:
    :raises FileNotFoundError: If input directory doesn't exist
    """
    logger.info(f"Creating manifest file from directory: {input_dir}")
    
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    # Get all fastq.gz files
    fastq_files = sorted([f for f in input_dir.glob('*.fastq.gz')] + [f for f in input_dir.glob('*.fq.gz')])
    
    # Process files in pairs
    processed_samples = {}
    
    for filepath in fastq_files:
        sample_id, read = parse_fastq_filename(filepath)
        
        if sample_id is None:
            logger.warning(f"Could not parse sample info from filename: {filepath}")
            continue
            
        if sample_id not in processed_samples:
            processed_samples[sample_id] = {'R1': None, 'R2': None}
            
        processed_samples[sample_id][read] = str(filepath.absolute())
    
    # Create manifest data
    manifest_data = {
        'sample-id': [],
        'forward-absolute-filepath': [],
        'reverse-absolute-filepath': []
    }
    
    for sample_id, paths in processed_samples.items():
        if paths['R1'] and paths['R2']:  # Only include complete pairs
            manifest_data['sample-id'].append(sample_id)
            manifest_data['forward-absolute-filepath'].append(paths['R1'])
            manifest_data['reverse-absolute-filepath'].append(paths['R2'])
        else:
            logger.warning(f"Missing read for sample {sample_id}")
    
    # Create manifest DataFrame and save
    manifest_df = pd.DataFrame(manifest_data)
    manifest_df.to_csv(output_file, index=False, sep='\t')
    
    logger.info(f"Created manifest file with {len(manifest_data['sample-id'])} samples: {output_file}")
    validate_manifest(manifest_df, logger)
    
    return output_file


def validate_manifest(manifest_df, logger):
    """
    Validate the manifest file contents.
    
    Parameters:
    :param manifest_df: Manifest DataFrame to validate
    :param logger: Logger instance
    
    Logs warnings for any validation issues found.
    """
    # Check for duplicate sample IDs
    duplicates = manifest_df['sample-id'].duplicated()
    if duplicates.any():
        logger.error(f"Duplicate sample IDs found: {manifest_df['sample-id'][duplicates].tolist()}")
    
    # Verify file existence
    for column in ['forward-absolute-filepath', 'reverse-absolute-filepath']:
        missing_files = [f for f in manifest_df[column] if not Path(f).exists()]
        if missing_files:
            logger.error(f"Files not found for {column}: {missing_files}")


def import_sequences(manifest_file, output_dir, dataset_name, logger):
    """
    Import sequences into QIIME 2 using the manifest file.
    
    Parameters:
    :param manifest_file: Path to manifest file (Path object)
    :param output_dir: Output directory for QIIME2 artifacts (Path object)
    :param dataset_name: Name of the dataset
    :param logger: Logger instance
    
    Raises:
    :raises subprocess.CalledProcessError: If QIIME2 commands fail
    """
    # Define output paths
    output_artifact = output_dir / f'{dataset_name}_paired_end_demux.qza'
    vis_artifact = output_dir / f'{dataset_name}_paired_end_demux.qzv'
    
    try:
        # Import sequences
        logger.info("Importing sequences into QIIME 2")
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "SampleData[PairedEndSequencesWithQuality]",
            "--input-path", str(manifest_file),
            "--output-path", str(output_artifact),
            "--input-format", "PairedEndFastqManifestPhred33V2"
        ], check=True)
        
        # Generate visualization
        subprocess.run([
            "qiime", "demux", "summarize",
            "--i-data", str(output_artifact),
            "--o-visualization", str(vis_artifact)
        ], check=True)
        
        logger.info(f"Created QIIME 2 artifacts: {output_artifact} and {vis_artifact}")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during QIIME 2 import: {str(e)}")
        raise


if __name__ == "__main__":
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Prepare sequences for QIIME2 analysis.')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing fastq files')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for QIIME2 artifacts')
    parser.add_argument('-n', '--name', required=True, help='Dataset name')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logger
    logger = util.get_formatted_logger('prepare_sequences', args.verbosity)
    
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        manifest_file = output_dir / f'{args.name}_manifest.tsv'
        manifest_file = create_manifest(Path(args.input_dir), manifest_file, logger)
        import_sequences(manifest_file, output_dir, args.name, logger)
        
        logger.info("Sequence preparation completed successfully")
    except Exception as e:
        logger.error(f"Error during sequence preparation: {str(e)}")
        raise