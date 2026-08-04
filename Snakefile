###############################################################################
# File: Snakefile
# Version: 1.0
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
#   2026-08-03
#
# License:
#   GNU GPL-3.0
#
# Dependencies:
#   - Snakemake
#   - Please check the conda environment specifications in the envs/ directory
#
###############################################################################


# -----------------------------------------------------------------------------
#### Reading the general config file

configfile: "config/project_config.yaml"


# -----------------------------------------------------------------------------
#### Load rules

include: "workflow/rules/download_ref.smk"
include: "workflow/rules/download_reads.smk"


# -----------------------------------------------------------------------------
#### Progress messages

onstart:
    print("\n--- WORKFLOW STARTED ---\n")


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
        rules.subpipeline_download_reads.output