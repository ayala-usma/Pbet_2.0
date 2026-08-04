###############################################################################
# File: download_reads.smk
# Version: 2.0
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
#   2026-08-04
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

# This rules presumes that the user is logged in to the AWS CLI. Please do that
# before running the pipeline.

rule download_lr: 
    params:
        out_dir = "results/data/lr"
    output: 
        "results/data/lr/{long_read}"
    log:
        "logs/download_reads/lr/download_{long_read}.log"
    conda:
        "../envs/downloads.yaml"
    message:
        "Downloading the PacBio long reads for the P. betacei \
        P8084 WGS project hosted in an AWS s3 bucket."
    shell:
        "export AWS_MAX_ATTEMPTS=10 &\
        aws s3 cp \
        s3://pbet-annotation/Pbet_SQII_subreads/{wildcards.long_read} \
        {params.out_dir}/{wildcards.long_read} --quiet &> {log}"


rule prefetch_sr:
    params:
        out_dir = "results/data/sr"
    output: 
        "results/data/sr/{short_read}/{short_read}.sra"
    log:
        "logs/download_reads/sr/download_{short_read}/{short_read}.log"
    conda:
        "../envs/downloads.yaml"
    message:
        "SRA prefetching the Illumina raw short reads for the P. betacei \
        P8084 WGS project."
    shell:
        "prefetch -v -L info --max-size 1t -O {params.out_dir} \
        {wildcards.short_read} &> {log}"
        
rule fasterq_sr:
    params:
        out_dir = "results/data/sr"
    input:
        rules.prefetch_sr.output
    output: 
        "results/data/sr/{short_read}/{short_read}.fastq"
    log:
        "logs/download_reads/sr/download_{short_read}/{short_read}.log"
    conda:
        "../envs/downloads.yaml"
    threads:
        workflow.cores * 0.8
    message:
        "Extraction of the Illumina raw short reads for the P. betacei \
        P8084 WGS project from the SRA prefetch file."
    shell:
        "fasterq-dump -e {threads} --details --skip-technical --split-spot \
        --outdir {params.out_dir}/{wildcards.short_read} \
        {params.out_dir}/{wildcards.short_read} &>> {log}"

rule compressing_sr:
    params:
        out_dir = "results/data/sr"
    input:
        reads = "results/data/sr/{short_read}/{short_read}.fastq",
        sra = rules.prefetch_sr.output
    output: 
        "results/data/sr/{short_read}/{short_read}.fq.zd"
    log:
        "logs/download_reads/sr/download_{short_read}/{short_read}.log"
    conda:
        "../envs/downloads.yaml"
    threads:
        workflow.cores * 0.8
    message:
        "Compressing the Illumina raw short reads for the P. betacei \
        P8084 WGS project."
    shell:
        "zstd --adapt --rm -f -v -T{threads} -o {input.reads}.tmp \
        {input.reads} &>> {log} & rm {input.sra} & \
        mv {input.reads}.tmp {output}"


# -----------------------------------------------------------------------------
#### Overall subpipeline rule

rule subpipeline_download_reads:
    input:
        expand("results/data/lr/{lr_acc}", lr_acc=LR_ACCESSIONS),
        expand("results/data/sr/{sr_acc}/{sr_acc}.fq.zd",
                sr_acc=SR_ACCESSIONS)
    output:
        "logs/download_reads/process_download_reads.log"
    shell:
        "echo '***Reads downloaded successfully!***' > {output}"