# QIIME 2 Phylogenetic Placement Workflow

![workflow](https://github.com/yourusername/qiime2-phylogenetic-placement/actions/workflows/python-package-conda.yml/badge.svg)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

![Logo](doc/logo.png)

## Overview

This repository contains a workflow for phylogenetic placement and diversity analysis using QIIME 2 and external tools. The workflow integrates QIIME 2's standard analysis pipeline with custom steps for using external reference phylogenies, such as those created by the BACTRIA pipeline.

Key features:
- Integration of external reference phylogenies into QIIME 2 workflow
- Custom sequence placement using SEPP or similar tools
- Calculation of phylogenetic diversity metrics based on placements
- Scalable analysis for large metabarcoding datasets

## Installation

This workflow uses Conda for dependency management. To set up the environment:

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/qiime2-phylogenetic-placement.git
   cd qiime2-phylogenetic-placement
   ```

2. Create and activate the Conda environment:
   ```
   conda env create -f workflow/envs/environment.yml
   conda activate qiime2-phylo-placement
   ```

## Usage

The workflow is managed using Snakemake. To run the entire pipeline:

```
snakemake --use-conda --cores all
```

For more detailed usage instructions, see the [documentation](workflow/documentation.md).

## Configuration

Adjust the `config/config.yaml` file to set parameters for your analysis, including:
- Input data locations
- Reference phylogeny file
- QIIME 2 parameters
- Diversity metric calculations

## Results

Output files will be generated in the `results/` directory, including:
- QIIME 2 artifacts
- Phylogenetic placements
- Diversity metric calculations

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for more details.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE.md](LICENSE.md) file for details.

