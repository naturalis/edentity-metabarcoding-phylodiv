import util
import argparse
from pathlib import Path
import subprocess
import sys
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from pathlib import Path

def build_hmm_profile(ref_align, hmm_profile, logger):
    """Build HMM profile with explicit RF marking."""
    cmd = [
        "hmmbuild",
        "--dna",
        "--hand",  # Use RF annotation if present
        str(hmm_profile),
        str(ref_align)
    ]
    subprocess.run(cmd, check=True)
    logger.info(f"Built HMM profile: {hmm_profile}")

def align_sequences(query_seqs, hmm_profile, output_align, logger):
    """Align sequences to HMM profile."""
    cmd = [
        "hmmalign",
        "--outformat", "Stockholm",
        "--trim",
        "--dna",
        str(hmm_profile),
        str(query_seqs)
    ]
    subprocess.run(cmd, check=True, stdout=open(output_align, 'w'))
    logger.info(f"Aligned sequences to profile: {output_align}")

def merge_alignments_fallback(ref_align_path, query_align_path, output_path, logger):
    """
    Merge alignments using BioPython as fallback method when esl-alimerge fails.
    """
    # Read both alignments
    ref_align = AlignIO.read(ref_align_path, "stockholm")
    query_align = AlignIO.read(query_align_path, "stockholm")
    
    # Create new alignment with all sequences
    merged = MultipleSeqAlignment([])
    
    # Add reference sequences
    for record in ref_align:
        merged.append(record)
        
    # Add query sequences
    for record in query_align:
        merged.append(record)
    
    # Ensure all sequences are same length by padding if needed
    max_length = max(len(ref_align[0]), len(query_align[0]))
    for record in merged:
        if len(record.seq) < max_length:
            record.seq = record.seq + '-' * (max_length - len(record.seq))
    
    # Write merged alignment
    AlignIO.write(merged, output_path, "stockholm")
    logger.info(f"Merged alignments using fallback method: {output_path}")
    return output_path

def merge_alignments(ref_align, query_align, output_path, logger):
    """
    Try to merge alignments using esl-alimerge first, fall back to Python method if it fails.
    """
    try:
        # First attempt with esl-alimerge
        cmd = [
            "esl-alimerge",
            "--outformat", "Stockholm",
            "-o", str(output_path),
            str(ref_align),
            str(query_align)
        ]
        subprocess.run(cmd, check=True)
        logger.info(f"Merged alignments using esl-alimerge: {output_path}")
        
    except subprocess.CalledProcessError:
        logger.warning("esl-alimerge failed, trying fallback merge method...")
        merge_alignments_fallback(ref_align, query_align, output_path, logger)

def prepare_alignments(ref_align, query_seqs, output_dir, logger):
    """
    Prepare and merge alignments for pplacer.
    """
    try:
        # Setup paths
        hmm_profile = output_dir / "ref.hmm"
        query_aligned = output_dir / "query_aligned.sto"
        combined_align = output_dir / "combined_alignment.sto"
        
        # Build HMM profile with RF annotation
        build_hmm_profile(ref_align, hmm_profile, logger)
        
        # Align query sequences
        align_sequences(query_seqs, hmm_profile, query_aligned, logger)
        
        # Merge alignments with fallback
        merge_alignments(ref_align, query_aligned, combined_align, logger)
        
        return combined_align
        
    except Exception as e:
        logger.error(f"Error preparing alignments: {str(e)}")
        raise

def preprocess_tree(tree_path, output_dir, logger):
    """
    Preprocess tree to fix common issues.
    """
    try:
        from Bio import Phylo
        import io
        
        # Read tree
        tree = Phylo.read(str(tree_path), 'newick')
        
        # Add small branch lengths where missing
        for clade in tree.find_clades():
            if clade.branch_length == 0 or clade.branch_length is None:
                clade.branch_length = 0.00001
                
        # Ensure tree is properly rooted
        tree.root_at_midpoint()
        
        # Write processed tree
        processed_tree = output_dir / "processed_reference_tree.tre"
        Phylo.write(tree, str(processed_tree), 'newick')
        
        logger.info(f"Preprocessed tree written to: {processed_tree}")
        return processed_tree
        
    except Exception as e:
        logger.error(f"Error preprocessing tree: {str(e)}")
        raise

def run_pplacer_with_stats(ref_align, query_seqs, ref_tree, stats_file, output_dir, threads, logger):
    """
    Run pplacer using an existing stats file.
    Optimized for large reference trees with memory constraints.
    """
    try:
        # Create output directory
        output_dir = Path(output_dir).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup paths
        output_jplace = output_dir / "placements.jplace"
        mmap_file = output_dir / "pplacer.mmap"
        
        # Convert to absolute paths
        ref_tree = Path(ref_tree).resolve()
        ref_align = Path(ref_align).resolve()
        query_seqs = Path(query_seqs).resolve()
        stats_file = Path(stats_file).resolve()
        
        # Setup paths and prepare alignments
        combined_align = prepare_alignments(ref_align, query_seqs, output_dir, logger)

        # In run_pplacer_with_stats:
        # Before running pplacer
        processed_tree = preprocess_tree(ref_tree, output_dir, logger)
        # Use processed_tree instead of ref_tree in pplacer command
        
        # Run pplacer with memory-optimized parameters
        logger.info("Running pplacer...")
        cmd = [
            "pplacer",
            "-t", str(processed_tree),           # Reference tree
            "-s", str(stats_file),         # RAxML stats file
            "--groups", "10",              # Split into more groups for lower memory usage
            "--max-strikes", "3",          # Reduced for speed and memory
            "--strike-box", "3.0",         # Default strike box size
            "--max-pitches", "20",         # Reduced number of pitches
            "-j", str(threads),            # Number of threads
            "--keep-at-most", "5",         # Reduced number of placements to keep
            "--keep-factor", "0.01",       # Default keep factor
            "--mmap-file", str(mmap_file), # Use memory mapping
            "-o", str(output_jplace),      # Output file
            str(combined_align)            # Input alignment
        ]
        
        logger.debug(f"Running command: {' '.join(cmd)}")
        
        # Run with increased process limits if possible
        try:
            import resource
            # Try to increase memory limit
            resource.setrlimit(resource.RLIMIT_AS, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        except (ImportError, resource.error):
            logger.warning("Could not adjust process resource limits")
        
        subprocess.run(cmd, check=True)
        
        logger.info(f"Placements written to: {output_jplace}")
        return output_jplace
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running command: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error during placement: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Place query sequences using pplacer with existing stats file.')
    parser.add_argument('-r', '--ref-align', required=True, help='Reference alignment in Stockholm format')
    parser.add_argument('-q', '--query', required=True, help='Query sequences in FASTA format')
    parser.add_argument('-t', '--tree', required=True, help='Reference tree in Newick format')
    parser.add_argument('-s', '--stats', required=True, help='Existing RAxML stats file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('-j', '--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('-v', '--verbosity', default='INFO', help='Logging verbosity level')
    args = parser.parse_args()

    logger = util.get_formatted_logger('pplacer_direct', args.verbosity)
    
    try:
        # Validate input files
        for path_arg in [args.ref_align, args.query, args.tree, args.stats]:
            if not Path(path_arg).exists():
                raise FileNotFoundError(f"Input file not found: {path_arg}")
        
        # Run pplacer with memory management
        output_jplace = run_pplacer_with_stats(
            Path(args.ref_align),
            Path(args.query), 
            Path(args.tree),
            Path(args.stats),
            Path(args.output_dir),
            args.threads,
            logger
        )
        
        logger.info("Phylogenetic placement completed successfully")
        
    except Exception as e:
        logger.error(f"Error during phylogenetic placement: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()