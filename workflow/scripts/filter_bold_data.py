import argparse
import pandas as pd
import os
import logging
from pathlib import Path

def setup_logging(log_file):
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def filter_bold_data(input_file, order, output_dir):
    """
    Filter BOLD database to include only entries from specified order.
    
    Args:
        input_file (str): Path to input BOLD TSV file
        order (str): Order to filter by (e.g., 'Lepidoptera')
        output_dir (str): Directory to save filtered data
    
    Returns:
        str: Path to output file
    """
    logger = setup_logging(os.path.join(output_dir, 'filter_bold.log'))
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"Reading BOLD database from {input_file}")
    logger.info(f"Filtering for order: {order}")
    
    # Read TSV file in chunks to handle large files
    chunk_size = 100000
    output_file = os.path.join(output_dir, f"{order.lower()}_entries.tsv")
    
    # Delete output file if it exists
    if os.path.exists(output_file):
        os.remove(output_file)
    
    # Process file in chunks
    chunks_processed = 0
    total_rows = 0
    filtered_rows = 0
    
    for chunk in pd.read_csv(input_file, sep='\t', chunksize=chunk_size):
        chunks_processed += 1
        total_rows += len(chunk)
        
        # Filter rows where order matches (case-insensitive)
        filtered_chunk = chunk[chunk['order'].str.lower() == order.lower()]
        filtered_rows += len(filtered_chunk)
        
        # Append to output file
        filtered_chunk.to_csv(
            output_file,
            sep='\t',
            mode='a',
            header=(chunks_processed == 1),
            index=False
        )
        
        if chunks_processed % 10 == 0:
            logger.info(f"Processed {chunks_processed} chunks ({total_rows:,} rows)")
    
    logger.info(f"Processing complete:")
    logger.info(f"Total rows processed: {total_rows:,}")
    logger.info(f"Rows matching order {order}: {filtered_rows:,}")
    logger.info(f"Output written to: {output_file}")
    
    return output_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter BOLD database by order')
    parser.add_argument('--input', required=True, help='Input BOLD TSV file')
    parser.add_argument('--order', required=True, help='Order to filter by')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()
    
    filter_bold_data(args.input, args.order, args.output_dir)
