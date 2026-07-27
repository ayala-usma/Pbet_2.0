###############################################################################
# File: download_ref.smk
# Version: 1.0
#
# Purpose:
#   Download the reference genomes of P. betacei and P. infestans mitochondrion
#   using entrez-direct tools
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
#   2026-07-22
#
# License:
#   GNU GPL-3.0
#
# Dependencies:
#   - Please check the conda environment specifications at envs/downloads.yaml
#
###############################################################################

# -----------------------------------------------------------------------------
#### Declaration of the rules


rule download_assembly:
    params:
        acc = lookup(within = config,
                        dpath = "get_references/P8084_assembly")
    output:
        ref_gen = "results/ref/geno_P8084.fa"
    log:
        "logs/download_P8084_assembly.log"
    conda:
        "../envs/downloads.yaml"
    message:
        "Downloading the reference assembly for P. betacei P8084"
        
    shell:
        "elink -db assembly -id {params.acc} -target nuccore | \
        efetch -format fasta 2> {log} 1> {output.ref_gen}"


rule download_mitochondria:
    params:
        acc = lookup(within = config,
                    dpath = "get_references/Pinf_mitochondrion")
    output:
        ref_mito = "results/ref/mito_Pinf.fa"
    log:
        "logs/download_Pinf_mito.log"
    conda:
        "../envs/downloads.yaml"
    message:
        "Downloading the reference mitochondrial genome for P. infestans"
    shell:
        "efetch -db nuccore -id {params.acc} -format fasta \
        2>{log} 1>{output.ref_mito}"
