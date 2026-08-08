###############################################################################
# File: read_preprocessing.smk
# Version: 0.1
#
# Purpose:
# Preprocessing and quality filtering
# long and short reads for the sequencing project of P. betacei.
#
# Part of:
#   Pbet2.0
#
# Author:
#   Aurelia Ayala Usma (ayala.usma@gmail.com)
#
# Created:
#   2026-08-08
#
# Last Updated:
#   2026-08-08
#
# License:
#   GNU GPL-3.0
#
# Dependencies:
#   - Conda environment specifications at envs/downloads.yaml
#
###############################################################################


# -----------------------------------------------------------------------------
#### Reading the accession list with polars

import polars as pl

long_reads = lookup(within=config, dpath="get_reads/lr")
short_reads = lookup(within=config, dpath="get_reads/sr")

LR_ACCESSIONS = pl.read_csv(long_reads).select('acc').to_series().to_list()
SR_ACCESSIONS = pl.read_csv(short_reads).select('acc').to_series().to_list()


# -----------------------------------------------------------------------------
#### Declaration of the rules


rule generate_hifi_zmw: 
    params:
        out_dir = "results/data/lr"
    output: 
        "results/data/lr/{long_read}"
    log:
        "logs/download_reads/lr/download_{long_read}.log"
    conda:
        "../envs/raw_data_processing.yaml"
    message:
        "Downloading the PacBio long reads for the P. betacei \
        P8084 WGS project hosted in an AWS s3 bucket."
    shell:
        "TO DO \
        ccs m54167_171009_163752.subreads.bam m54167_171009_163752.ccs.bam --log-level INFO --all --all-kinetics --subread-fallback --log-file m54167_171009_163752.ccs.log -j 8
        extracthifi m54167_171010_024731.ccs.bam m54167_171010_024731.hifi.bam -j 8 &

        "


# -----------------------------------------------------------------------------
#### Overall subpipeline rule

rule subpipeline_read_preprocessing:
    input:
        "TO DO"
    output:
        "logs/read_preprocessing/preprocessing.log"
    shell:
        "TO DO"