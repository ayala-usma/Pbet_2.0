###############################################################################
# File: download_reads.smk
# Version: 1.0
#
# Purpose:
#   Download the long and short reads for the sequencing project of P. betacei
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
#   2026-07-20
#
# License:
#   GNU GPL-3.0
#
# Dependencies:
#   - Please check the conda environment specifications at envs/downloads.yaml
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


rule complete_downloads:
    input:
        


rule download_lr:
    output: 
        expand("results/data/lr/{lr_acc}.txt", lr_acc=LR_ACCESSIONS)
    log:
        "logs/download_lr.log"
    conda:
        "../envs/downloads.yaml"
    message:
        "Downloading the PacBio long reads for the P. betacei \
        P8084 WGS project."
    shell:
        "touch {output}"


rule download_sr:
    output: 
        expand("results/data/sr/{sr_acc}.txt", sr_acc=SR_ACCESSIONS)
    log:
        "logs/download_sr.log"
    conda:
        "../envs/downloads.yaml"
    message:
        "Downloading the Illumina raw short reads for the P. betacei \
        P8084 WGS project."
    shell:
        "touch {output}"
