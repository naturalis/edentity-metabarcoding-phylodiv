import argparse
from qiime2.plugins import feature_classifier
import util


def prepare_reference_set(reference_sequences, reference_taxonomy, output_classifier):
    """
    Prepare the reference sequence set and train the classifier.

    :param reference_sequences: Path to reference sequences QIIME 2 artifact
    :param reference_taxonomy: Path to reference taxonomy QIIME 2 artifact
    :param output_classifier: Path to output trained classifier
    """
    logger = util.setup_logging('prepare_reference_set.log')
    logger.info("Starting reference set preparation")

    # Load reference sequences and taxonomy
    ref_seqs = util.load_qiime2_artifact(reference_sequences)
    ref_tax = util.load_qiime2_artifact(reference_taxonomy)

    # Train the classifier
    logger.info("Training the classifier")
    classifier = feature_classifier.methods.fit_classifier_naive_bayes(
        reference_reads=ref_seqs,
        reference_taxonomy=ref_tax
    )

    # Save the trained classifier
    logger.info(f"Saving the trained classifier to {output_classifier}")
    util.save_qiime2_artifact(classifier.classifier, output_classifier)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare reference set and train classifier')
    parser.add_argument('--ref-seqs', required=True, help='Path to reference sequences QIIME 2 artifact')
    parser.add_argument('--ref-tax', required=True, help='Path to reference taxonomy QIIME 2 artifact')
    parser.add_argument('--output', required=True, help='Path to output trained classifier')
    args = parser.parse_args()

    util.validate_inputs([args.ref_seqs, args.ref_tax])
    prepare_reference_set(args.ref_seqs, args.ref_tax, args.output)