# blast_filter.py
import util
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import subprocess
import numpy as np

def create_blast_db(reference_fasta, output_path, logger):
    """
    Create BLAST database from FASTA file.
    
    Parameters:
    :param reference_fasta: Path to reference FASTA file
    :param output_path: Path for output directory
    :param logger: Logger instance
    
    Returns:
    :return: Path to reference FASTA file
    """
    try:
        # Create BLAST database directly from FASTA
        cmd = [
            "makeblastdb",
            "-in", str(reference_fasta),
            "-dbtype", "nucl",
            "-out", str(output_path / "reference_db")
        ]
        subprocess.run(cmd, check=True)
        logger.info(f"Created BLAST database at {output_path}")
        
        return reference_fasta
        
    except Exception as e:
        logger.error(f"Error creating BLAST database: {str(e)}")
        raise

def run_blast(query_fasta, db_path, output_file, threads, logger):
    """
    Run BLAST against reference database.
    
    Parameters:
    :param query_fasta: Path to query FASTA file
    :param db_path: Path to BLAST database
    :param output_file: Path to output file
    :param threads: Number of threads to use
    :param logger: Logger instance
    """
    try:
        cmd = [
            "blastn",
            "-query", str(query_fasta),
            "-db", str(db_path / "reference_db"),
            "-out", str(output_file),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-num_threads", str(threads)
        ]
        subprocess.run(cmd, check=True)
        logger.info(f"BLAST search completed, results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during BLAST search: {str(e)}")
        raise

def filter_blast_results(blast_results, max_evalue, logger):
    """
    Filter BLAST results based on E-value threshold.
    
    Parameters:
    :param blast_results: Path to BLAST results file
    :param max_evalue: Maximum E-value threshold
    :param logger: Logger instance
    
    Returns:
    :return: DataFrame with filtered results
    """
    try:
        # Read BLAST results
        cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        df = pd.read_csv(blast_results, sep='\t', names=cols)
        
        # Filter based on E-value
        filtered_df = df[df['evalue'] <= max_evalue]
        #CHANGE BEST_HITS TO filtered_df

        # Keep only best hit per query (lowest E-value)
        best_hits = filtered_df.sort_values('evalue').groupby('qseqid').head(5)
        
        return best_hits
        
    except Exception as e:
        logger.error(f"Error filtering results: {str(e)}")
        raise

def write_filtered_sequences(query_fasta, filtered_results, output_fasta, logger):
    """
    Write filtered sequences to new FASTA file based on BLAST results.
    
    Parameters:
    :param query_fasta: Path to query sequences FASTA file
    :param filtered_results: DataFrame of filtered BLAST results
    :param output_fasta: Path to output filtered FASTA file
    :param logger: Logger instance
    """
    try:
        # Get unique query sequence IDs from filtered results
        query_ids = set(filtered_results['qseqid'].unique())
        
        # Write matching sequences to new file
        count = 0
        total_seqs = 0
        written_ids = set()  # Track which IDs we've written
        
        with open(output_fasta, 'w') as out_handle:
            for record in SeqIO.parse(query_fasta, "fasta"):
                total_seqs += 1
                if record.id in query_ids:
                    if record.id not in written_ids:  # Only write each ID once
                        out_handle.write(f">{record.id}\n{str(record.seq)}\n")
                        written_ids.add(record.id)
                        count += 1
        
        logger.info(f"Found {len(query_ids)} unique query IDs with BLAST hits")
        logger.info(f"Processed {total_seqs} total input sequences")
        logger.info(f"Wrote {count} filtered sequences to {output_fasta}")
        
        if count == 0:
            raise ValueError("No sequences passed BLAST filtering criteria")
            
    except Exception as e:
        logger.error(f"Error writing filtered sequences: {str(e)}")
        raise

def extract_reference_sequences(ref_fasta, ref_ids, output_fasta, logger):
    """
    Extract reference sequences from reference FASTA for matched IDs.
    
    Parameters:
    :param ref_fasta: Path to reference FASTA file
    :param ref_ids: Set or list of reference IDs to extract
    :param output_fasta: Path to output FASTA file
    :param logger: Logger instance
    """
    try:
        # Convert ref_ids to set for faster lookup
        ref_ids = set(ref_ids)
        
        # Write matching sequences to new file
        count = 0
        written_ids = set()  # Track which IDs we've written
        with open(output_fasta, 'w') as out_handle:
            for record in SeqIO.parse(ref_fasta, "fasta"):
                if record.id in ref_ids:
                    if record.id not in written_ids:  # Only write each ID once
                        out_handle.write(f">{record.id}\n{str(record.seq)}\n")
                        written_ids.add(record.id)
                        count += 1
        
        logger.info(f"Found {len(ref_ids)} unique reference IDs")
        logger.info(f"Extracted {count} reference sequences to {output_fasta}")
        
        if count == 0:
            raise ValueError("No matching reference sequences found")
            
    except Exception as e:
        logger.error(f"Error extracting reference sequences: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='BLAST and filter sequences.')
    parser.add_argument('-i', '--input', required=True, help='Input sequences in FASTA format')
    parser.add_argument('-r', '--reference', required=True, help='Reference sequences TSV file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('--max_evalue', type=float, default=1e-10, help='Maximum E-value threshold')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads for BLAST')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()

    logger = util.get_formatted_logger('blast_filter', args.verbosity)
    
    try:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create BLAST database and get reference FASTA path
        logger.info(f"Creating BLAST database from reference TSV")
        ref_fasta = create_blast_db(Path(args.reference), output_dir, logger)
        
        # Run BLAST
        blast_output = output_dir / "blast_results.txt"
        logger.info(f"Running BLAST search against database")
        run_blast(Path(args.input), output_dir, blast_output, args.threads, logger)
        
        # Filter results
        filtered_results = filter_blast_results(
            blast_output,
            args.max_evalue,
            logger
        )
        
        # Write reference IDs
        ref_ids = filtered_results['sseqid'].unique()
        ref_range_file = output_dir / "reference_ids.txt"
        with open(ref_range_file, 'w') as f:
            f.write('\n'.join(ref_ids))
        
        # Write filtered sequences
        filtered_seqs = output_dir / "filtered_sequences.fasta"
        write_filtered_sequences(Path(args.input), filtered_results, filtered_seqs, logger)
        
        # Extract matching reference sequences
        reference_subset = output_dir / "reference_subset.fasta"
        extract_reference_sequences(ref_fasta, ref_ids, reference_subset, logger)
        
        logger.info("Processing completed successfully")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise

if __name__ == "__main__":
    main()