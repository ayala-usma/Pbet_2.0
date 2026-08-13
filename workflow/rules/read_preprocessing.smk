###############################################################################
# File: read_preprocessing.smk
# Version: 0.3
#
# Purpose:
# Preprocessing and quality filtering long and short reads for the sequencing 
# project of P. betacei.
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
#   2026-08-13
#
# License:
#   GNU GPL-3.0
#
# Dependencies:
#   - Conda environment specifications at envs/raw_data_preprocessing.yaml
#
###############################################################################

# -----------------------------------------------------------------------------
#### Reading the input file paths

import polars as pl

long_reads = lookup(within=config, dpath="get_reads/lr")
short_reads = lookup(within=config, dpath="get_reads/sr")

LR_PATHS = pl.read_csv(long_reads).select('acc').to_series().to_list()
LR_PATHS = [
    'results/data/lr/' + element for element in LR_PATHS
    ]

SR_PATHS = pl.read_csv(short_reads).select('acc').to_series().to_list()
SR_PATHS = [
    'results/data/sr/' + element + '/' + element for element in SR_PATHS
    ]


# -----------------------------------------------------------------------------
#### Declaration of the rules

rule lr_ccs: 
    params:
        snr = 4.5
    input:
        "results/data/lr/{lr_path}.subreads.bam"
    output: 
        ccs_bam = "results/data/lr/{lr_path}.ccs.bam",
        index= "results/data/lr/{lr_path}.ccs.bam.pbi",
        hifi_bam = "results/data/lr/{lr_path}.hifi.bam",
        met = "results/data/lr/{lr_path}.ccs.zmw_metrics.json.gz"
    log:
        "logs/read_preprocessing/lr/ccs_{lr_path}.log"
    conda:
        "../envs/raw_data_preprocessing.yaml"
    threads:
        workflow.cores * 0.9
    message:
        "Performing the CCS pipeline from PacBio with custom SNR threshold"
    shell:
        """
        ccs --min-snr {params.snr} --log-level INFO --all --all-kinetics \
            --subread-fallback --log-file {log} -j {threads} \
            {input} {output.ccs_bam}

        extracthifi {output.ccs_bam} {output.hifi_bam} -j {threads}
        """

rule ccs_filtering:
    params:
        min_len = 2500
    input:
        hifi_bam = rules.lr_ccs.output.hifi_bam,
        ccs_bam =  rules.lr_ccs.output.ccs_bam,
        met = rules.lr_ccs.output.met
    output:
        "results/data/lr/{lr_path}.filtered_reads.bam"
    log:
        "logs/read_preprocessing/lr/filtering_{lr_path}.log"
    conda:
        "../envs/raw_data_preprocessing.yaml"
    threads:
        workflow.cores * 0.9
    message:
        "Performing the quality filtering of the CCS bam file according to \
        minimum read length and SNR threshold."
    shell:
        """
        samtools view {input.hifi_bam} | \
            awk '{{print substr($NF, 6)}}' > {input.hifi_bam}.excl.mask

        zcat {input.met} | \
        jq '.[][] | select(.insert_size >= {params.min_len} and
            (.status == "LOW_PASS_SHORTCUT" or
             .status == "SUCCESS")) .zmw' | \
             sed -E 's/\".*\/([0-9]+)\"/\\1/g' > {input.met}.incl.mask

        zmwfilter -j {threads} --log-level INFO --log-file {log} \
        --exclude {input.hifi_bam}.excl.mask {input.ccs_bam} {input.hifi_bam}.excl.temp.bam

        pbindex -j {threads} {input.hifi_bam}.excl.temp.bam

        zmwfilter -j {threads} --log-level INFO --log-file {log} \
        --include {input.met}.incl.mask {input.hifi_bam}.excl.temp.bam {output}

        pbindex -j {threads} {output}

        rm {input.hifi_bam}.excl.mask {input.hifi_bam}.excl.temp.bam
        rm {input.ccs_bam} {input.ccs_bam}.pbi {input.met}.incl.mask
        rm {input.hifi_bam}.excl.temp.bam.pbi
        """
        
# rule sr_PE_fastp:
#     params:
#         min_qual = 30,
#         min_len = 50,
#     input:
#     output:
#     log:
#         "logs/read_preprocessing/sr/filtering_{sr_path}.log"
#     conda:
#         "../envs/raw_data_preprocessing.yaml"
#     threads:
#         workflow.cores * 0.9
#     message:
#         "Performing the quality filtering of the paired-end short reads \
#         using the fastp with base trimming using minimum base quality and \
#         keeping a minimum read length"
#     shell:
#         """
#         TO DO
#         zstd -dc SRR14352918/SRR14352918.fq.zst | \
#         fastp --interleaved_in --stdin -q {params.min_qual} \
#         -o SRR14352918_R1.fq.gz -O SRR14352918_R2.fq.gz -h SRR14352918.html \
#         --cut_front --cut_tail -l {params.min_len} -M {params.min_qual}
#         -j SRR14352918.json -w {threads}
#         """

# rule sr_SE_fastp:
#     params:
#         min_qual = 30,
#         min_len = 50,
#     input:
#     output:
#     log:
#         "logs/read_preprocessing/sr/filtering_{sr_path}.log"
#     conda:
#         "../envs/raw_data_preprocessing.yaml"
#     threads:
#         workflow.cores * 0.9
#     message:
#         "Performing the quality filtering of the paired-end short reads \
#         using the fastp with base trimming using minimum base quality and \
#         keeping a minimum read length"
#     shell:
#         """
#         TO DO
#         zstd -dc SRR14352918/SRR14352918.fq.zst | \
#         fastp --interleaved_in --stdin -q {params.min_qual} \
#         -o SRR14352918_R1.fq.gz -O SRR14352918_R2.fq.gz -h SRR14352918.html \
#         --cut_front --cut_tail -l {params.min_len} -M {params.min_qual}
#         -j SRR14352918.json -w {threads}
#         """

# -----------------------------------------------------------------------------
#### Overall subpipeline rule

rule subpipeline_read_preprocessing:
    input:
        expand("{lr_full_path}.filtered_reads.bam", lr_full_path=LR_PATHS)
    output:
        "logs/read_preprocessing/preprocessing.log"
    shell:
        "echo '***Reads successfully preprocessed!***' > {output}"