import argparse
from qiime2.plugins import diversity
import util


def run_diversity_analysis(table, tree, sample_metadata, metrics, sampling_depth, output_dir):
    """
    Run diversity analysis using QIIME 2.

    :param table: Path to feature table QIIME 2 artifact
    :param tree: Path to phylogenetic tree QIIME 2 artifact
    :param sample_metadata: Path to sample metadata file
    :param metrics: List of diversity metrics to compute
    :param sampling_depth: Rarefaction depth
    :param output_dir: Directory to save output artifacts
    """
    logger = util.setup_logging('diversity_analysis.log')
    logger.info("Starting diversity analysis")

    util.ensure_dir(output_dir)

    # Load artifacts
    feature_table = util.load_qiime2_artifact(table)
    phylogeny = util.load_qiime2_artifact(tree)

    # Run core diversity metrics
    logger.info(f"Running core diversity metrics with sampling depth {sampling_depth}")
    diversity_results = diversity.pipelines.core_metrics_phylogenetic(
        table=feature_table,
        phylogeny=phylogeny,
        sampling_depth=sampling_depth,
        metadata=sample_metadata
    )

    # Save diversity results
    for name, artifact in diversity_results.artifacts.items():
        output_path = f"{output_dir}/{name}.qza"
        logger.info(f"Saving diversity artifact: {output_path}")
        util.save_qiime2_artifact(artifact, output_path)

    # Compute additional alpha diversity metrics if specified
    alpha_metrics = set(metrics) - set(['faith_pd', 'evenness', 'shannon'])
    if alpha_metrics:
        logger.info(f"Computing additional alpha diversity metrics: {alpha_metrics}")
        alpha_results = diversity.methods.alpha_diversity(
            table=feature_table,
            metric=list(alpha_metrics)
        )
        for metric, artifact in alpha_results.items():
            output_path = f"{output_dir}/alpha_{metric}.qza"
            logger.info(f"Saving alpha diversity artifact: {output_path}")
            util.save_qiime2_artifact(artifact, output_path)

    # Compute additional beta diversity metrics if specified
    beta_metrics = set(metrics) - set(['unweighted_unifrac', 'weighted_unifrac', 'jaccard', 'bray_curtis'])
    if beta_metrics:
        logger.info(f"Computing additional beta diversity metrics: {beta_metrics}")
        beta_results = diversity.methods.beta_diversity(
            table=feature_table,
            metric=list(beta_metrics),
            phylogeny=phylogeny
        )
        for metric, artifact in beta_results.items():
            output_path = f"{output_dir}/beta_{metric}.qza"
            logger.info(f"Saving beta diversity artifact: {output_path}")
            util.save_qiime2_artifact(artifact, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run diversity analysis')
    parser.add_argument('--table', required=True, help='Path to feature table QIIME 2 artifact')
    parser.add_argument('--tree', required=True, help='Path to phylogenetic tree QIIME 2 artifact')
    parser.add_argument('--sample-metadata', required=True, help='Path to sample metadata file')
    parser.add_argument('--metrics', nargs='+',
                        default=['faith_pd', 'shannon', 'observed_features', 'weighted_unifrac', 'unweighted_unifrac'],
                        help='Diversity metrics to compute')
    parser.add_argument('--sampling-depth', type=int, required=True, help='Rarefaction depth')
    parser.add_argument('--output-dir', required=True, help='Directory to save output artifacts')
    args = parser.parse_args()

    util.validate_inputs([args.table, args.tree, args.sample_metadata])
    run_diversity_analysis(args.table, args.tree, args.sample_metadata, args.metrics, args.sampling_depth,
                           args.output_dir)