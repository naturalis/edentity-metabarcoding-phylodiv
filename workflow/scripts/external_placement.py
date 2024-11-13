import argparse
import subprocess
import os
import util


def run_sepp(query_sequences, reference_alignment, reference_tree, output_dir):
    """
    Run SEPP for phylogenetic placement.

    :param query_sequences: Path to query sequences file
    :param reference_alignment: Path to reference alignment file
    :param reference_tree: Path to reference tree file
    :param output_dir: Directory to save output files
    """
    logger = util.setup_logging('external_placement.log')
    logger.info("Starting SEPP phylogenetic placement")

    util.ensure_dir(output_dir)
    output_file = os.path.join(output_dir, "placement.json")

    command = [
        "run_sepp.py",
        "-s", query_sequences,
        "-a", reference_alignment,
        "-t", reference_tree,
        "-o", output_file
    ]

    logger.info(f"Running SEPP command: {' '.join(command)}")
    subprocess.run(command, check=True)


def convert_to_newick(placement_json, output_newick):
    """
    Convert SEPP JSON output to Newick format.

    :param placement_json: Path to SEPP JSON output
    :param output_newick: Path to save Newick tree
    """
    logger = util.setup_logging('external_placement.log')
    logger.info("Converting SEPP output to Newick format")

    command = [
        "guppy", "tog",
        "-o", output_newick,
        placement_json
    ]

    logger.info(f"Running guppy command: {' '.join(command)}")
    subprocess.run(command, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run external phylogenetic placement')
    parser.add_argument('--query-seqs', required=True, help='Path to query sequences file')
    parser.add_argument('--ref-alignment', required=True, help='Path to reference alignment file')
    parser.add_argument('--ref-tree', required=True, help='Path to reference tree file')
    parser.add_argument('--output-dir', required=True, help='Directory to save output files')
    args = parser.parse_args()

    util.validate_inputs([args.query_seqs, args.ref_alignment, args.ref_tree])
    run_sepp(args.query_seqs, args.ref_alignment, args.ref_tree, args.output_dir)

    placement_json = os.path.join(args.output_dir, "placement.json")
    output_newick = os.path.join(args.output_dir, "placement_tree.nwk")
    convert_to_newick(placement_json, output_newick)