import util
import argparse
import pandas as pd 
from pathlib import Path
import random
import re

def extract_processids(tree_file):
    """Extract processids from Newick tree file."""
    with open(tree_file) as f:
        tree_text = f.read()
    pattern = r'[,(]([A-Z0-9-]+):'
    processids = re.findall(pattern, tree_text)
    return set(processids)

def is_valid_sequence(sequence, required_length=500):
    """
    Check if sequence only contains uppercase A, C, T, G and is of required length.
    
    Parameters:
    :param sequence: Input DNA sequence string
    :param required_length: Required sequence length (default: 658 for COI-5P)
    
    Returns:
    :return: Boolean indicating if sequence is valid
    """
    if pd.isna(sequence):
        return False
        
    # Check sequence length
    if len(sequence) < required_length:
        return False
        
    # Check for valid characters (must be uppercase ACTG)
    return set(sequence).issubset({'A', 'C', 'T', 'G'})

def clean_sequence(sequence, required_length=658):
    """
    Clean sequence and check if valid.
    
    Parameters:
    :param sequence: Input DNA sequence string
    :param required_length: Required sequence length (default: 658 for COI-5P)
    
    Returns:
    :return: Cleaned sequence string if valid, None if invalid
    """
    if pd.isna(sequence):
        return None
    
    # Remove alignment dashes and any whitespace
    cleaned = sequence.replace('-', '').strip()
    
    # Check sequence length and valid characters
    if is_valid_sequence(cleaned, required_length):
        return cleaned
    return None

def write_fasta(sequences, output_file):
    """
    Write sequences to FASTA format.
    
    Parameters:
    :param sequences: DataFrame with processid and nuc columns
    :param output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as f:
        for _, row in sequences.iterrows():
            f.write(f">{row['processid']}\n{row['nuc']}\n")

def create_feature_table(regional_sequences, output_file):
    """
    Create QIIME2-style feature table for geographical samples.
    
    Parameters:
    :param regional_sequences: Dict of region -> DataFrame of sequences
    :param output_file: Path to output TSV file
    """
    # Initialize feature table with all sequence IDs as index
    all_sequences = pd.concat([df for df in regional_sequences.values()])
    feature_table = pd.DataFrame(0, 
                               index=all_sequences['processid'].unique(),
                               columns=['NA', 'EU', 'AS'])
    
    # Mark presence (1) for each sequence in its region
    for region, sequences in regional_sequences.items():
        feature_table.loc[sequences['processid'], region] = 1
    
    # Save feature table
    feature_table.to_csv(output_file, sep='\t')

def process_chunk_for_references(chunk, reference_ids, stats):
    """Process a chunk of data for reference sequences."""
    reference_matches = chunk[chunk['processid'].isin(reference_ids)].copy()
    
    # Clean and validate sequences
    reference_matches['nuc'] = reference_matches['nuc'].apply(clean_sequence)
    reference_matches = reference_matches.dropna(subset=['nuc'])
    
    # Update stats to track invalid sequences
    if 'invalid_sequences' not in stats:
        stats['invalid_sequences'] = 0
        stats['wrong_length'] = 0
        
    total_matches = len(chunk[chunk['processid'].isin(reference_ids)])
    valid_matches = len(reference_matches)
    stats['invalid_sequences'] += total_matches - valid_matches
    
    return reference_matches[['processid', 'nuc']]

def process_chunk_for_region(chunk, region, region_map, exclude_ids):
    """
    Process a chunk of data for geographical sampling.
    
    Parameters:
    :param chunk: DataFrame chunk to process
    :param region: Region identifier (NA, EU, AS)
    :param region_map: Dictionary mapping regions to countries
    :param exclude_ids: Set of processids to exclude (reference sequences)
    """
    # First filter by region and basic criteria
    filtered = chunk[
        (chunk['country/ocean'].isin(region_map[region])) &
        (chunk['order'].str.lower() == 'lepidoptera') &
        (chunk['marker_code'].str.contains('COI-5P', case=False, na=False))
    ].copy()
    
    # Remove sequences that are in the reference set
    filtered = filtered[~filtered['processid'].isin(exclude_ids)]
    
    # Clean and validate sequences
    filtered['nuc'] = filtered['nuc'].apply(clean_sequence)
    filtered = filtered.dropna(subset=['nuc'])
    
    # Select required columns
    return filtered[['processid', 'nuc']]

def filter_bold(tree_file, infile, output_dir, logger):
    """
    Filter BOLD data based on reference tree processids and create geographical samples.
    
    Parameters:
    :param tree_file: Path to reference tree file (.nwk format)
    :param infile: Input BOLD TSV file path
    :param output_dir: Output directory for all files
    :param logger: Logger instance
    """
    logger.info(f"Reading reference tree from {tree_file}")
    reference_ids = extract_processids(tree_file)
    logger.info(f"Extracted {len(reference_ids)} processids from reference tree")
    
    logger.info(f"Reading BOLD database from {infile}")
    
    chunk_size = 100000
    usecols = ['processid', 'order', 'nuc', 'marker_code', 'country/ocean']
    
    # Create output directories
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_dir = output_dir / "sample_sequences"
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    stats = {
        'chunks_processed': 0,
        'total_rows': 0,
        'reference_matches': 0,
        'invalid_sequences': 0
    }

    # First pass: collect reference sequences
    reference_sequences = pd.DataFrame(columns=['processid', 'nuc'])
    
    logger.info("Processing reference sequences...")
    for chunk in pd.read_csv(infile, sep='\t', usecols=usecols, chunksize=chunk_size, 
                           on_bad_lines='skip', low_memory=False,
                           dtype={'order': 'string', 'marker_code': 'string'}):
        stats['chunks_processed'] += 1
        stats['total_rows'] += len(chunk)
        
        reference_matches = process_chunk_for_references(chunk, reference_ids, stats)
        reference_sequences = pd.concat([reference_sequences, reference_matches])
        
        if stats['chunks_processed'] % 10 == 0:
            logger.info(f"Processed {stats['chunks_processed']} chunks ({stats['total_rows']:,} rows)")
    
    # Save reference sequences as FASTA
    reference_fasta = sample_dir / "reference_sequences.fasta"
    write_fasta(reference_sequences, reference_fasta)
    stats['reference_matches'] = len(reference_sequences)
    stats['invalid_sequences'] = len(reference_ids) - len(reference_sequences)
    
    # Second pass: collect geographical samples
    region_map = {
        'NA': ['Costa Rica', 'United States', 'Canada', 'Mexico', 'Panama', 'Trinidad and Tobago', 
               'Guatemala', 'Honduras', 'Dominican Republic', 'Cuba', 'Bahamas'],
        'EU': ['Germany', 'Italy', 'France', 'Norway', 'Austria', 'Spain', 'United Kingdom', 
               'Netherlands', 'Portugal', 'Finland', 'Sweden', 'Croatia', 'Latvia', 'Slovakia', 
               'Switzerland', 'Romania', 'Czech Republic', 'Denmark', 'Slovenia', 'Estonia', 
               'Lithuania', 'Hungary', 'Poland', 'Belgium', 'Ireland', 'North Macedonia'],
        'AS': ['China', 'Thailand', 'Indonesia', 'Malaysia', 'Pakistan', 'Philippines', 'Vietnam', 
               'India', 'South Korea', 'Taiwan', 'Turkey', 'Iran', 'Japan']
    }
    
    regional_sequences = {region: pd.DataFrame(columns=['processid', 'nuc']) for region in region_map}
    
    logger.info("Processing geographical samples...")
    for chunk in pd.read_csv(infile, sep='\t', usecols=usecols, chunksize=chunk_size,
                           on_bad_lines='skip', low_memory=False,
                           dtype={'order': 'string', 'marker_code': 'string'}):
        for region in region_map:
            matches = process_chunk_for_region(chunk, region, region_map, reference_sequences['processid'])
            regional_sequences[region] = pd.concat([regional_sequences[region], matches])
    
    # Sample sequences and save combined FASTA
    sample_size = 300
    sampled_sequences = {}
    for region in region_map:
        df = regional_sequences[region]
        if len(df) < sample_size:
            logger.warning(f"Only {len(df)} sequences available for {region} (requested {sample_size})")
            sampled_sequences[region] = df
        else:
            sampled_sequences[region] = df.sample(n=sample_size, random_state=42)
    
    # Write combined FASTA
    combined_fasta = sample_dir / "geographical_samples.fasta"
    combined_sequences = pd.concat([df for df in sampled_sequences.values()])
    write_fasta(combined_sequences, combined_fasta)
    
    # Create and save feature table
    feature_table = sample_dir / "feature_table.tsv"
    create_feature_table(sampled_sequences, feature_table)
    
    logger.info(f"Saved combined samples to: {combined_fasta}")
    logger.info(f"Saved feature table to: {feature_table}")
    
    log_stats(stats, logger)

def log_stats(stats, logger):
    """Log final statistics of the filtering process."""
    logger.info("Processing complete:")
    logger.info(f"Total rows processed: {stats['total_rows']:,}")
    logger.info(f"Reference sequences found: {stats['reference_matches']:,}")
    logger.info(f"Invalid sequences: {stats['invalid_sequences']:,}")
    logger.info("  - Common issues:")
    logger.info(f"    * Wrong sequence length (not 658bp)")
    logger.info(f"    * Invalid characters (non ACTG or lowercase)")
    logger.info(f"Invalid/missing sequences in reference set: {stats['invalid_sequences']:,}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter BOLD data based on reference tree and create geographical samples.')
    parser.add_argument('-t', '--tree', required=True, help='Input reference tree file (.nwk format)')
    parser.add_argument('-i', '--infile', required=True, help='Input BOLD TSV file')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    logger = util.get_formatted_logger('filter_bold_data', args.verbosity)

    try:
        filter_bold(Path(args.tree), Path(args.infile), Path(args.outdir), logger)
        logger.info("BOLD data filtering completed successfully")
    except Exception as e:
        logger.error(f"Error during BOLD data filtering: {str(e)}")
        raise