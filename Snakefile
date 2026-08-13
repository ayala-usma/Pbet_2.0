###############################################################################
# File: Snakefile
# Version: 1.3
#
# Purpose:
#   Main workflow entry point for the Pbet2.0 Snakemake pipeline.
#   Defines the workflow structure, loads configuration files,
#   and includes all rule modules.
#
# Part of:
#   Pbet2.0
#
# Author:
#   Aurelia Ayala Usma (ayala.usma@gmail.com)
#
# Created:
#   2026-07-20
#
# Last Updated:
#   2026-08-13
#
# License:
#   GNU GPL-3.0
#
# Dependencies:
#   - Mamba or Conda
#   - Snakemake
#   - For the dependencies of each subpipeline, please check the /workflow/envs
#     directory
#  
###############################################################################


# -----------------------------------------------------------------------------
#### Reading the general config file

configfile: "config/project_config.yaml"


# -----------------------------------------------------------------------------
#### Load subpipelines

include: "workflow/rules/download_ref.smk"
include: "workflow/rules/download_reads.smk"
include: "workflow/rules/read_preprocessing.smk"


# -----------------------------------------------------------------------------
#### Progress messages

onsuccess:
    print("\n--- WORKFLOW FINISHED! ---\n")


onerror:
    print("\n--- :o An error occurred! ---\n")


# -----------------------------------------------------------------------------
#### Target rule

rule all:
    default_target: True
    input:
        rules.subpipeline_download_ref.output,
        rules.subpipeline_download_reads.output,
        rules.subpipeline_read_preprocessing.output