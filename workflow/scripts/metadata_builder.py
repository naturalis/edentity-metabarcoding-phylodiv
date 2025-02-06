import pandas as pd
from pathlib import Path
import argparse
import qiime2
import logging

def setup_logger():
    """Set up logger for debugging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger('metadata_builder')

def create_qiime2_files(feature_table_path, output_dir, logger):
    """
    Create QIIME 2 metadata file and properly formatted feature table.
    
    Parameters:
    feature_table_path: Path to input feature table TSV
    output_dir: Directory to write output files
    logger: Logger instance
    """
    # Convert to Path objects
    feature_table_path = Path(feature_table_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read original feature table
    ft = pd.read_csv(feature_table_path, sep='\t', index_col=0)
    logger.info(f"Read original feature table with shape: {ft.shape}")
    
    # Features are the sample IDs, samples are the geographical regions
    # No need to transpose as it's already in correct format
    feature_table_path = output_dir / 'feature-table.tsv'
    ft.to_csv(feature_table_path, sep='\t')
    logger.info(f"Created feature table with {len(ft.index)} features and {len(ft.columns)} samples")
    
    # Create metadata for the samples (NA, EU, AS)
    sample_metadata = pd.DataFrame(index=ft.columns)
    sample_metadata.index.name = '#SampleID'
    
    # Add descriptive information for each sample
    sample_metadata['sample-type'] = 'geographical-region'
    sample_metadata['description'] = sample_metadata.index.map({
        'NA': 'North America',
        'EU': 'Europe',
        'AS': 'Asia'
    })
    
    # Write sample metadata file
    metadata_path = output_dir / 'sample-metadata.tsv'
    with open(metadata_path, 'w') as f:
        f.write(f"#SampleID\tsample-type\tdescription\n")
        f.write("#q2:types\tcategorical\tcategorical\n")
        sample_metadata.to_csv(f, sep='\t', header=False)
        
    logger.info(f"Created metadata file at: {metadata_path}")
    logger.info("Next steps:")
    logger.info("1. Import the feature table into QIIME 2:")
    logger.info(f"   qiime tools import --input-path {feature_table_path} \\\n" + 
               "      --type 'FeatureTable[Frequency]' \\\n" +
               "      --input-format BIOMV210Format \\\n" +
               f"      --output-path {output_dir}/feature-table.qza")

if __name__ == "__main__":
    logger = setup_logger()
    
    parser = argparse.ArgumentParser(description='Create QIIME 2 files')
    parser.add_argument('feature_table', help='Path to feature table TSV file')
    parser.add_argument('output_dir', help='Directory to write output files')
    
    args = parser.parse_args()
    
    try:
        create_qiime2_files(args.feature_table, args.output_dir, logger)
    except Exception as e:
        logger.error(f"Error creating files: {str(e)}")
        raise