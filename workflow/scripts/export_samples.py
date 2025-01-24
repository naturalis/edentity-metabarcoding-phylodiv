import argparse
from pathlib import Path
import subprocess
import logging

def export_sequences(input_qza, output_dir, logger):
    """
    Export sequences from QIIME2 artifact to FASTA format.
    
    Parameters:
    :param input_qza: Path to input .qza file
    :param output_dir: Path to output directory
    :param logger: Logger instance
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        "conda", "run", "-n", "qiime2-amplicon-2024.5",
        "qiime", "tools", "export",
        "--input-path", str(input_qza),
        "--output-path", str(output_dir)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        logger.info(f"Successfully exported sequences to {output_dir}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error exporting sequences: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Export sequences from QIIME2 artifact.')
    parser.add_argument('-i', '--input', required=True, help='Input QIIME2 artifact')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-v', '--verbosity', default='INFO', help='Log level (e.g. DEBUG)')
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=args.verbosity,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger('export_sequences')
    
    try:
        export_sequences(args.input, args.output_dir, logger)
    except Exception as e:
        logger.error(f"Error during sequence export: {str(e)}")
        raise