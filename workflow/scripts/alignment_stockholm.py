import util
import argparse
from pathlib import Path
import subprocess
import tempfile
from Bio import SeqIO
from Bio.AlignIO import read as read_alignment
import os

def align_score(sequence, hmm_model, logger):
    """
    Align sequence and calculate score using hmmalign's posterior probability output.
    Uses Stockholm format for scoring purposes.
    
    Parameters:
    :param sequence: SeqIO sequence record
    :param hmm_model: Path to HMM model file
    :param logger: Logger instance
    
    Returns:
    :return: alignment score (float)
    """
    try:
        with tempfile.NamedTemporaryFile(mode='w+') as temp_fasta, \
             tempfile.NamedTemporaryFile(mode='w+') as temp_stockholm:
            
            # Write sequence to temp file
            SeqIO.write(sequence, temp_fasta.name, 'fasta')
            
            # Run hmmalign in Stockholm format for scoring
            subprocess.run([
                'hmmalign',
                '--trim',
                '-o', temp_stockholm.name,
                str(hmm_model),
                temp_fasta.name
            ], check=True)
            
            # Read alignment and get posterior probabilities
            alignment = read_alignment(temp_stockholm.name, "stockholm")
            quality_string = alignment.column_annotations['posterior_probability']
            
            # Calculate weighted score using more sophisticated scoring
            # '.' = 0.0, '*' = 10.0, numbers weighted by value
            count = 0
            dot_count = quality_string.count('.')
            star_count = quality_string.count('*')
            
            # Weight stars and dots
            count += dot_count * 0  # No confidence
            count += star_count * 10  # Maximum confidence
            
            # Add numerical values for medium confidence positions
            digit_sum = sum(int(char) for char in quality_string if char.isdigit())
            count += digit_sum
            
            # Calculate average score
            average_score = count / len(quality_string)
            
            return average_score
            
    except Exception as e:
        logger.error(f"Error calculating alignment score: {str(e)}")
        raise

def align_sequences_hmm(sequences, hmm_model, output_file, logger, output_format='fasta'):
    """
    Align sequences using HMMER's hmmalign against a provided HMM model.
    
    Parameters:
    :param sequences: Path to input sequences FASTA file
    :param hmm_model: Path to HMM model file
    :param output_file: Path to output aligned file
    :param logger: Logger instance 
    :param output_format: Format for output ('fasta', 'afa', or 'stockholm')
    """
    try:
        # Always use Stockholm format for hmmalign output initially
        with tempfile.NamedTemporaryFile(mode='w+') as temp_align:
            # Run hmmalign
            cmd = [
                'hmmalign',
                '--trim',  # Trim alignment to match model
                '-o', temp_align.name,
                str(hmm_model),
                str(sequences)
            ]
            subprocess.run(cmd, check=True)
            
            # Read alignment in Stockholm format
            alignment = read_alignment(temp_align.name, "stockholm")
            
            # Write output based on requested format
            if output_format.lower() == 'stockholm':
                # For Stockholm format, preserve original output
                with open(output_file, 'w') as f:
                    f.write(temp_align.read())
            elif output_format.lower() == 'afa':
                # For aligned FASTA format
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

def check_reverse_complement(sequence, hmm_model, logger, min_score_threshold=3.0):
    """
    Check if sequence should be reverse complemented based on HMM alignment score.
    Also filters out sequences with poor alignment in both directions.
    
    Parameters:
    :param sequence: SeqIO sequence record
    :param hmm_model: Path to HMM model file
    :param logger: Logger instance
    :param min_score_threshold: Minimum acceptable alignment score (default: 3.0)
    
    Returns:
    :return: Tuple of (SeqIO record in best orientation or None, best_score)
    """
    try:
        # Score original orientation
        forward_score = align_score(sequence, hmm_model, logger)
        logger.debug(f'Forward alignment score {forward_score} for {sequence.id}')
        
        # Score reverse complement
        sequence.seq = sequence.seq.reverse_complement()
        reverse_score = align_score(sequence, hmm_model, logger)
        logger.debug(f'Reverse complemented alignment score {reverse_score} for {sequence.id}')
        
        # Get best score and corresponding sequence orientation
        best_score = max(forward_score, reverse_score)
        
        # Check if both scores are below threshold
        if best_score < min_score_threshold:
            logger.warning(f"Sequence {sequence.id} filtered out due to poor alignment scores "
                         f"(forward: {forward_score:.2f}, reverse: {reverse_score:.2f})")
            return None, best_score
            
        # Return sequence in best orientation
        if forward_score >= reverse_score:
            sequence.seq = sequence.seq.reverse_complement()  # Revert back
            logger.debug(f"Keeping forward orientation for {sequence.id}")
        else:
            logger.debug(f"Using reverse complement for {sequence.id}")
            
        return sequence, best_score
            
    except Exception as e:
        logger.error(f"Error checking sequence orientation: {str(e)}")
        raise

def process_sequences(input_fasta, hmm_model, output_file, logger, output_format='fasta', 
                     check_orientation=True, min_score_threshold=3.0):
    """
    Process sequences through HMM alignment pipeline with quality filtering.
    
    Parameters:
    :param input_fasta: Path to input FASTA file
    :param hmm_model: Path to HMM model file
    :param output_file: Path to output aligned file
    :param logger: Logger instance
    :param output_format: Format for output ('fasta', 'afa', or 'stockholm')
    :param check_orientation: Whether to check sequence orientation
    :param min_score_threshold: Minimum acceptable alignment score
    """
    try:
        if check_orientation:
            # Process sequences and check orientation
            reverse_complement_count = 0
            filtered_count = 0
            total_sequences = 0
            scores = []
            retained_sequences = 0
            
            # Write oriented sequences
            temp_fasta = output_file.parent / f"{output_file.stem}_oriented.fasta"
            with open(temp_fasta, 'w') as f:
                for record in SeqIO.parse(input_fasta, "fasta"):
                    total_sequences += 1
                    original_seq = str(record.seq)
                    
                    # Check orientation and filter low quality sequences
                    result = check_reverse_complement(
                        record, 
                        hmm_model, 
                        logger,
                        min_score_threshold
                    )
                    
                    # Skip if sequence was filtered out
                    if result[0] is None:
                        filtered_count += 1
                        continue
                        
                    oriented_record, score = result
                    scores.append(score)
                    
                    if str(oriented_record.seq) != original_seq:
                        reverse_complement_count += 1
                        
                    # Write oriented sequence
                    f.write(f">{oriented_record.id}\n{str(oriented_record.seq)}\n")
                    retained_sequences += 1
            
            # Log statistics
            logger.info(f'Total sequences processed: {total_sequences}')
            logger.info(f'Sequences filtered due to poor alignment: {filtered_count} '
                       f'({filtered_count/total_sequences*100:.1f}%)')
            logger.info(f'Sequences retained: {retained_sequences} '
                       f'({retained_sequences/total_sequences*100:.1f}%)')
            logger.info(f'Sequences reverse complemented: {reverse_complement_count} '
                       f'({reverse_complement_count/retained_sequences*100:.1f}% of retained)')
            
            if scores:
                logger.info(f'Alignment score statistics:')
                logger.info(f'  Mean score: {sum(scores)/len(scores):.2f}')
                logger.info(f'  Min score: {min(scores):.2f}')
                logger.info(f'  Max score: {max(scores):.2f}')
            
            # Only proceed with alignment if we have sequences left
            if retained_sequences > 0:
                # Align oriented sequences
                align_sequences_hmm(temp_fasta, hmm_model, output_file, logger, output_format)
                # Clean up
                os.remove(temp_fasta)
            else:
                logger.error("No sequences passed quality filtering")
                raise ValueError("All sequences were filtered out due to poor alignment scores")
            
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
    parser.add_argument('-f', '--format', choices=['fasta', 'afa', 'stockholm'], default='fasta',
                      help='Output format (fasta, afa, or stockholm)')
    parser.add_argument('--skip-orientation', action='store_true', help='Skip orientation check')
    parser.add_argument('--min-score', type=float, default=3.0,
                      help='Minimum alignment score threshold (default: 3.0)')
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
            not args.skip_orientation,
            args.min_score
        )
        
        logger.info("Sequence alignment completed successfully")
        
    except Exception as e:
        logger.error(f"Error during sequence alignment: {str(e)}")
        raise

if __name__ == "__main__":
    main()