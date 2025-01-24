import util
import argparse
from pathlib import Path
import subprocess
import tempfile
from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
import os

def align_sequences_hmm(sequences, hmm_model, output_file, logger, output_format='fasta'):
    """
    Align sequences using HMMER's hmmalign against a provided HMM model.
    
    Parameters:
    :param sequences: Path to input sequences FASTA file
    :param hmm_model: Path to HMM model file
    :param output_file: Path to output aligned file
    :param logger: Logger instance 
    :param output_format: Format for output ('fasta' or 'afa')
    """
    try:
        # Create temporary files for intermediate steps
        with tempfile.NamedTemporaryFile(mode='w+') as temp_align:
            # Run hmmalign
            cmd = [
                'hmmalign',
                '--trim',  # Trim alignment to match model
                '--outformat', 'afa',  # Always use aligned FASTA for hmmalign output
                '-o', temp_align.name,
                str(hmm_model),
                str(sequences)
            ]
            subprocess.run(cmd, check=True)
            
            # Read aligned sequences
            alignment = read_alignment(temp_align.name, "fasta")
            
            # Write to output file based on requested format
            if output_format.lower() == 'afa':
                # For AFA format, preserve alignment characters
                with open(output_file, 'w') as f:
                    for record in alignment:
                        f.write(f">{record.id}\n{str(record.seq)}\n")
            else:
                # For regular FASTA format, remove gaps
                with open(output_file, 'w') as f:
                    for record in alignment:
                        sequence = str(record.seq).replace('-', '')
                        f.write(f">{record.id}\n{sequence}\n")
                    
        logger.info(f"Successfully aligned sequences to {output_file} in {output_format} format")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during alignment with hmmalign: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error processing alignment: {str(e)}")
        raise

def check_reverse_complement(sequence, hmm_model, logger):
    """
    Check if sequence should be reverse complemented based on HMM alignment score.
    
    Parameters:
    :param sequence: SeqIO sequence record
    :param hmm_model: Path to HMM model file
    :param logger: Logger instance
    
    Returns:
    :return: SeqIO record in correct orientation
    """
    try:
        # Score original orientation
        with tempfile.NamedTemporaryFile(mode='w+') as temp_fasta, \
             tempfile.NamedTemporaryFile(mode='w+') as temp_stockholm:
            
            # Write sequence to temp file
            SeqIO.write(sequence, temp_fasta.name, 'fasta')
            
            # Align and get score
            subprocess.run([
                'hmmalign',
                '--trim',
                '-o', temp_stockholm.name,
                str(hmm_model),
                temp_fasta.name
            ], check=True)
            
            forward_alignment = read_alignment(temp_stockholm.name, "stockholm")
            forward_score = sum(1 for c in forward_alignment.column_annotations['posterior_probability'] if c == '*')
        
        # Score reverse complement
        sequence.seq = sequence.seq.reverse_complement()
        with tempfile.NamedTemporaryFile(mode='w+') as temp_fasta, \
             tempfile.NamedTemporaryFile(mode='w+') as temp_stockholm:
            
            SeqIO.write(sequence, temp_fasta.name, 'fasta')
            subprocess.run([
                'hmmalign',
                '--trim',
                '-o', temp_stockholm.name,
                str(hmm_model),
                temp_fasta.name
            ], check=True)
            
            reverse_alignment = read_alignment(temp_stockholm.name, "stockholm")
            reverse_score = sum(1 for c in reverse_alignment.column_annotations['posterior_probability'] if c == '*')
        
        # Return sequence in best orientation
        if forward_score >= reverse_score:
            sequence.seq = sequence.seq.reverse_complement()  # Revert back
            logger.debug(f"Keeping forward orientation for {sequence.id}")
        else:
            logger.debug(f"Using reverse complement for {sequence.id}")
            
        return sequence
        
    except Exception as e:
        logger.error(f"Error checking sequence orientation: {str(e)}")
        raise

def process_sequences(input_fasta, hmm_model, output_file, logger, output_format='fasta', check_orientation=True):
    """
    Process sequences through HMM alignment pipeline.
    
    Parameters:
    :param input_fasta: Path to input FASTA file
    :param hmm_model: Path to HMM model file
    :param output_file: Path to output aligned file
    :param logger: Logger instance
    :param output_format: Format for output ('fasta' or 'afa')
    :param check_orientation: Whether to check sequence orientation
    """
    try:
        if check_orientation:
            # Process sequences and check orientation
            temp_fasta = output_file.parent / f"{output_file.stem}_oriented.fasta"
            with open(temp_fasta, 'w') as f:
                for record in SeqIO.parse(input_fasta, "fasta"):
                    oriented_record = check_reverse_complement(record, hmm_model, logger)
                    f.write(f">{oriented_record.id}\n{str(oriented_record.seq)}\n")
            
            # Align oriented sequences
            align_sequences_hmm(temp_fasta, hmm_model, output_file, logger, output_format)
            
            # Clean up
            os.remove(temp_fasta)
        else:
            # Direct alignment without orientation check
            align_sequences_hmm(input_fasta, hmm_model, output_file, logger, output_format)
            
    except Exception as e:
        logger.error(f"Error processing sequences: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Align sequences using HMM model.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-m', '--model', required=True, help='HMM model file')
    parser.add_argument('-o', '--output', required=True, help='Output aligned file')
    parser.add_argument('-f', '--format', choices=['fasta', 'afa'], default='fasta',
                      help='Output format (fasta or afa)')
    parser.add_argument('--skip-orientation', action='store_true', help='Skip orientation check')
    parser.add_argument('-v', '--verbosity', required=True, help='Log level (e.g. DEBUG)')
    args = parser.parse_args()
    
    # Setup logging
    logger = util.get_formatted_logger('sequence_alignment', args.verbosity)
    
    try:
        # Convert paths
        input_file = Path(args.input)
        hmm_model = Path(args.model)
        output_file = Path(args.output)
        
        # Create output directory if needed
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Process sequences
        process_sequences(
            input_file,
            hmm_model,
            output_file,
            logger,
            args.format,
            not args.skip_orientation
        )
        
        logger.info("Sequence alignment completed successfully")
        
    except Exception as e:
        logger.error(f"Error during sequence alignment: {str(e)}")
        raise

if __name__ == "__main__":
    main()