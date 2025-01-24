import util
import argparse
import pandas as pd
from pathlib import Path


def filter_bold(order, infile, outfile, logger):
    """
    Filter BOLD data for a specific order and COI-5P marker.
    
    Parameters:
    :param order: Taxonomic order to filter for
    :param infile: Input BOLD TSV file path (Path object)
    :param outfile: Output filtered TSV file path (Path object) 
    :param logger: Logger instance for tracking execution
    """
    logger.info(f"Reading BOLD database from {infile}")
    logger.info(f"Filtering for order: {order} and marker: COI-5P")
    
    chunk_size = 100000
    usecols = ['processid', 'order', 'nuc', 'marker_code']
    
    if outfile.exists():
        outfile.unlink()
    
    order_lower = order.lower()
    seen_processids = set()
    stats = {
        'chunks_processed': 0,
        'total_rows': 0,
        'filtered_rows': 0,
        'duplicates': 0,
        'non_coi5p': 0
    }
    
    for chunk in pd.read_csv(infile, sep='\t', 
                           chunksize=chunk_size,
                           usecols=usecols,
                           on_bad_lines='skip', 
                           low_memory=False,
                           dtype={'order': 'string', 'marker_code': 'string'}):
        
        stats['chunks_processed'] += 1
        stats['total_rows'] += len(chunk)
        
        chunk = chunk.dropna(subset=['processid', 'marker_code', 'nuc'])
        filtered_chunk = chunk[
            (chunk['order'].str.lower().eq(order_lower)) &
            (chunk['marker_code'].str.contains('COI-5P', case=False, na=False))
        ]
        
        filtered_chunk = filtered_chunk[['processid', 'nuc']]
        
        stats['non_coi5p'] += len(chunk[chunk['order'].str.lower().eq(order_lower)]) - len(filtered_chunk)
        
        chunk_processids = set(filtered_chunk['processid'])
        new_duplicates = len(chunk_processids.intersection(seen_processids))
        stats['duplicates'] += new_duplicates
        seen_processids.update(chunk_processids)
        
        stats['filtered_rows'] += len(filtered_chunk)
        
        filtered_chunk.to_csv(
            outfile,
            sep='\t',
            mode='a',
            header=(stats['chunks_processed'] == 1),
            index=False
        )
        
        if stats['chunks_processed'] % 10 == 0:
            logger.info(f"Processed {stats['chunks_processed']} chunks ({stats['total_rows']:,} rows)")
    
    log_stats(stats, order, logger)


def log_stats(stats, order, logger):
    """
    Log final statistics of the filtering process.
    
    Parameters:
    :param stats: Dictionary containing processing statistics
    :param order: Taxonomic order that was filtered
    :param logger: Logger instance
    """
    logger.info("Processing complete:")
    logger.info(f"Total rows processed: {stats['total_rows']:,}")
    logger.info(f"Rows matching order {order} with COI-5P: {stats['filtered_rows']:,}")
    logger.info(f"Non-COI-5P entries excluded: {stats['non_coi5p']:,}")


if __name__ == "__main__":
    # Define command line arguments
    parser = argparse.ArgumentParser(description='Filter BOLD data for specific order and COI-5P marker.')
    parser.add_argument('-i', '--infile', required=True, help='Input BOLD TSV file')
    parser.add_argument('-o', '--outfile', required=True, help='Output filtered TSV file')
    parser.add_argument('-r', '--order', required=True, help='Taxonomic order to filter')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    # Configure logger
    logger = util.get_formatted_logger('filter_bold_data', args.verbosity)

    try:
        filter_bold(args.order, Path(args.infile), Path(args.outfile), logger)
        logger.info("BOLD data filtering completed successfully")
    except Exception as e:
        logger.error(f"Error during BOLD data filtering: {str(e)}")
        raise