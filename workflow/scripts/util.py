import os
import logging
from qiime2 import Artifact
from qiime2.plugins import feature_table, metadata


def setup_logging(log_file, log_level='INFO'):
    """
    Set up logging for the scripts.

    :param log_file: Path to the log file
    :param log_level: Logging level (default: INFO)
    :return: Configured logger object
    """
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        filename=log_file,
        filemode='w'
    )
    logger = logging.getLogger(__name__)
    return logger


def ensure_dir(directory):
    """
    Make sure the specified directory exists, create it if it doesn't.

    :param directory: Path to the directory
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def load_qiime2_artifact(filepath):
    """
    Load a QIIME 2 artifact from a file.

    :param filepath: Path to the QIIME 2 artifact
    :return: Loaded QIIME 2 artifact
    """
    return Artifact.load(filepath)


def save_qiime2_artifact(artifact, filepath):
    """
    Save a QIIME 2 artifact to a file.

    :param artifact: QIIME 2 artifact to save
    :param filepath: Path to save the artifact
    """
    artifact.save(filepath)


def merge_feature_tables(table_filepaths, output_filepath):
    """
    Merge multiple QIIME 2 feature tables into a single table.

    :param table_filepaths: List of paths to feature table artifacts
    :param output_filepath: Path to save the merged feature table
    """
    tables = [load_qiime2_artifact(filepath) for filepath in table_filepaths]
    merged_table = feature_table.methods.merge(tables=tables)
    save_qiime2_artifact(merged_table.merged_table, output_filepath)


def filter_features(table_filepath, min_frequency, min_samples, output_filepath):
    """
    Filter features from a QIIME 2 feature table based on frequency and prevalence.

    :param table_filepath: Path to the input feature table artifact
    :param min_frequency: Minimum total frequency for a feature to be retained
    :param min_samples: Minimum number of samples a feature must be observed in to be retained
    :param output_filepath: Path to save the filtered feature table
    """
    table = load_qiime2_artifact(table_filepath)
    filtered_table = feature_table.methods.filter_features(
        table=table,
        min_frequency=min_frequency,
        min_samples=min_samples
    )
    save_qiime2_artifact(filtered_table.filtered_table, output_filepath)


def summarize_feature_table(table_filepath, sample_metadata_filepath, output_dir):
    """
    Generate a summary of a QIIME 2 feature table.

    :param table_filepath: Path to the feature table artifact
    :param sample_metadata_filepath: Path to the sample metadata file
    :param output_dir: Directory to save the summary
    """
    table = load_qiime2_artifact(table_filepath)
    sample_metadata = metadata.Metadata.load(sample_metadata_filepath)
    summary = feature_table.visualizers.summarize(
        table=table,
        sample_metadata=sample_metadata
    )
    summary.visualization.save(os.path.join(output_dir, 'feature_table_summary.qzv'))


def validate_inputs(filepaths):
    """
    Validate that all input files exist.

    :param filepaths: List of file paths to validate
    :raises FileNotFoundError: If any file is not found
    """
    for filepath in filepaths:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Input file not found: {filepath}")


def get_sampling_depth(feature_table_summary):
    """
    Determine an appropriate sampling depth for rarefaction based on the feature table summary.
    This is a placeholder function and should be implemented based on your specific criteria.

    :param feature_table_summary: Path to the feature table summary
    :return: Suggested sampling depth
    """
    # This is a placeholder implementation
    # In a real scenario, you might parse the summary file and apply some logic
    # to determine an appropriate sampling depth
    return 10000  # Example return value

# Add more utility functions as needed for your workflow