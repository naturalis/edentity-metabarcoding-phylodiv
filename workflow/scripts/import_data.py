import util
import argparse
from pathlib import Path
import subprocess
import pandas as pd
import biom

def tsv_to_biom(tsv_path, biom_path, logger):
    """
    Convert TSV feature table to BIOM format.
    
    Parameters:
    :param tsv_path: Path to input TSV file
    :param biom_path: Path to output BIOM file
    :param logger: Logger instance
    """
    try:
        # Read TSV file
        data = pd.read_csv(tsv_path, sep='\t', index_col=0)
        
        # Convert to biom table
        table = biom.table.Table(
            data=data.values,
            observation_ids=data.index,
            sample_ids=data.columns
        )
        
        # Write to file
        with biom.util.biom_open(str(biom_path), 'w') as f:
            table.to_hdf5(f, generated_by="feature-table-import")
            
        logger.info(f"Converted TSV to BIOM: {biom_path}")
        
    except Exception as e:
        logger.error(f"Error converting TSV to BIOM: {str(e)}")
        raise

def import_biom_to_qiime2(biom_path, output_path, logger):
    """
    Import BIOM file into QIIME2.
    
    Parameters:
    :param biom_path: Path to input BIOM file
    :param output_path: Path to output QIIME2 artifact
    :param logger: Logger instance
    """
    try:
        cmd = [
            "qiime", "tools", "import",
            "--input-path", str(biom_path),
            "--type", "FeatureTable[Frequency]",
            "--input-format", "BIOMV210Format",
            "--output-path", str(output_path)
        ]
        
        subprocess.run(cmd, check=True)
        logger.info(f"Imported BIOM to QIIME2: {output_path}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error importing BIOM to QIIME2: {str(e)}")
        raise

def import_tree_to_qiime2(tree_path, output_path, is_rooted, logger):
    """
    Import Newick tree into QIIME2.
    
    Parameters:
    :param tree_path: Path to input Newick tree file
    :param output_path: Path to output QIIME2 artifact
    :param is_rooted: Whether tree is rooted or not
    :param logger: Logger instance
    """
    try:
        tree_type = "Phylogeny[Rooted]" if is_rooted else "Phylogeny[Unrooted]"
        
        cmd = [
            "qiime", "tools", "import",
            "--input-path", str(tree_path),
            "--type", tree_type,
            "--output-path", str(output_path)
        ]
        
        subprocess.run(cmd, check=True)
        logger.info(f"Imported tree to QIIME2: {output_path}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error importing tree to QIIME2: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Import TSV feature table and tree into QIIME2')
    parser.add_argument('-i', '--input_tsv', required=True, help='Input TSV feature table')
    parser.add_argument('-t', '--input_tree', required=True, help='Input Newick tree file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('--rooted', action='store_true', help='Specify if tree is rooted')
    parser.add_argument('-v', '--verbosity', default='INFO', help='Log level')
    args = parser.parse_args()
    
    # Setup logging
    logger = util.get_formatted_logger('import_qiime2', args.verbosity)
    
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Convert TSV to BIOM
        biom_path = output_dir / "feature-table.biom"
        tsv_to_biom(Path(args.input_tsv), biom_path, logger)
        
        # Import BIOM to QIIME2
        feature_table_qza = output_dir / "feature-table.qza"
        import_biom_to_qiime2(biom_path, feature_table_qza, logger)
        
        # Import tree to QIIME2
        tree_qza = output_dir / "tree.qza"
        import_tree_to_qiime2(Path(args.input_tree), tree_qza, args.rooted, logger)
        
        logger.info("Import to QIIME2 completed successfully")
        
    except Exception as e:
        logger.error(f"Error during import: {str(e)}")
        raise

if __name__ == "__main__":
    main()