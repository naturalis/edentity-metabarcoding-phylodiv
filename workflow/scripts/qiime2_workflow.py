import argparse
from qiime2.plugins import demux, dada2, feature_classifier, phylogeny
import util


def run_qiime2_workflow(input_sequences, sample_metadata, classifier, output_dir):
    """
    Run the QIIME 2 workflow including denoising, taxonomic classification, and tree building.

    :param input_sequences: Path to input sequences QIIME 2 artifact
    :param sample_metadata: Path to sample metadata file
    :param classifier: Path to trained classifier
    :param output_dir: Directory to save output artifacts
    """
    logger = util.setup_logging('qiime2_workflow.log')
    logger.info("Starting QIIME 2 workflow")

    util.ensure_dir(output_dir)

    # Import sequences
    logger.info("Loading input sequences")
    demuxed_seqs = util.load_qiime2_artifact(input_sequences)

    # Denoise with DADA2
    logger.info("Running DADA2 denoising")
    dada2_result = dada2.methods.denoise_single(
        demultiplexed_seqs=demuxed_seqs,
        trunc_len=150,
        trim_left=0
    )

    # Taxonomic classification
    logger.info("Performing taxonomic classification")
    classifier_artifact = util.load_qiime2_artifact(classifier)
    taxonomy_result = feature_classifier.methods.classify_sklearn(
        reads=dada2_result.representative_sequences,
        classifier=classifier_artifact
    )

    # Phylogenetic tree construction
    logger.info("Constructing phylogenetic tree")
    mafft_result = phylogeny.methods.align_to_tree_mafft_fasttree(
        sequences=dada2_result.representative_sequences
    )

    # Save outputs
    logger.info("Saving output artifacts")
    util.save_qiime2_artifact(dada2_result.table, f"{output_dir}/feature_table.qza")
    util.save_qiime2_artifact(dada2_result.representative_sequences, f"{output_dir}/rep_seqs.qza")
    util.save_qiime2_artifact(taxonomy_result.classification, f"{output_dir}/taxonomy.qza")
    util.save_qiime2_artifact(mafft_result.alignment, f"{output_dir}/aligned_rep_seqs.qza")
    util.save_qiime2_artifact(mafft_result.tree, f"{output_dir}/tree.qza")

    # Generate feature table summary
    util.summarize_feature_table(f"{output_dir}/feature_table.qza", sample_metadata, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run QIIME 2 workflow')
    parser.add_argument('--input-seqs', required=True, help='Path to input sequences QIIME 2 artifact')
    parser.add_argument('--sample-metadata', required=True, help='Path to sample metadata file')
    parser.add_argument('--classifier', required=True, help='Path to trained classifier')
    parser.add_argument('--output-dir', required=True, help='Directory to save output artifacts')
    args = parser.parse_args()

    util.validate_inputs([args.input_seqs, args.sample_metadata, args.classifier])
    run_qiime2_workflow(args.input_seqs, args.sample_metadata, args.classifier, args.output_dir)