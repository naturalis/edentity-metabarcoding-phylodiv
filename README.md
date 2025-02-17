# Edentity Metabarcoding Phylodiversity

## Overview

This repository contains a workflow for phylogenetic placement and diversity analysis using QIIME 2 and external tools. The workflow integrates QIIME 2's standard analysis pipeline with custom steps for using external reference phylogenies, such as those created by the BACTRIA pipeline.

Key features:
- Integration of external reference phylogenies into QIIME 2 workflow
- Custom sequence placement
- Calculation of phylogenetic diversity metrics based on placements
- Scalable analysis for large metabarcoding datasets

## Installation

This workflow uses Conda for dependency management. To set up the environment:

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/qiime2-phylogenetic-placement.git
   cd qiime2-phylogenetic-placement
   ```

2. At the moment environments are installed per script used:
   ```
   conda env create -f workflow/envs/script.yaml
   conda activate script
   ```

## Usage

The workflow is managed using Snakemake. To run the entire pipeline:

NOT YET IMPLEMENTED

All scripts will need to be executed manually, for all commands and input files both Inputs_Outputs_Structure.txt & Workflow_commands_internship.txt can be used to find the necessary commands

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

