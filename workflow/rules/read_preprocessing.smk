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
#   2026-08-15
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
        snr = 4
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
        min_len = 2000
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
        --exclude {input.hifi_bam}.excl.mask {input.ccs_bam} \
        {input.hifi_bam}.excl.temp.bam

        pbindex -j {threads} {input.hifi_bam}.excl.temp.bam

        zmwfilter -j {threads} --log-level INFO --log-file {log} \
        --include {input.met}.incl.mask {input.hifi_bam}.excl.temp.bam {output}

        pbindex -j {threads} {output}

        rm {input.hifi_bam}.excl.mask {input.hifi_bam}.excl.temp.bam
        rm {input.ccs_bam} {input.ccs_bam}.pbi {input.met}.incl.mask
        rm {input.hifi_bam}.excl.temp.bam.pbi
        """

rule download_adapter_seq:
    params:
        adapter_url = "https://raw.githubusercontent.com/BioInfoTools/BBMap/refs/heads/master/resources/adapters.fa"
    output:
        "results/ref/BBMap_adapters.fa"
    log:
        "logs/read_preprocessing/sr/adapters.log"
    conda:
        "../envs/raw_data_preprocessing.yaml"
    message:
        "Downloading BBMap fasta file with the Illumina adapters used in \
        multiple sequencing applications."
    shell:
        """
        curl -o {output} {params.adapter_url}
        """

rule sr_fastp_filtering:
    params:
        min_qual = 30,
        min_len = 80,
        max_unq = 5
    input:
        raw_rds = "results/data/sr/{sr_accID}/{sr_accID}.fq.zst",
        adapters = rules.download_adapter_seq.output
    output:
        clean_rds = "results/data/sr/{sr_accID}/{sr_accID}_clean.fq.gz",
        html_rep = "results/data/sr/{sr_accID}/{sr_accID}.html",
        json_rep = "results/data/sr/{sr_accID}/{sr_accID}.json"
    log:
        "logs/read_preprocessing/sr/filtering_{sr_accID}.log"
    conda:
        "../envs/raw_data_preprocessing.yaml"
    threads:
        workflow.cores * 0.9
    message:
        "Performing the quality filtering of the Illumina WGS short reads \
        using fastp with base trimming using minimum base quality and \
        keeping a minimum read length"
    shell:
        """
        if [ "{wildcards.sr_accID}" == "SRR14352918" ]; then
            zstd -dc {input.raw_rds} | \
            fastp --interleaved_in --stdin -o {output.clean_rds} \
            -O "results/data/sr/{wildcards.sr_accID}/{wildcards.sr_accID}_clean_R2.fq.gz" \
            -q {params.min_qual} -h {output.html_rep} --cut_front \
            --cut_tail -l {params.min_len} -M {params.min_qual} \
            -u {params.max_unq} -j {output.json_rep} --detect_adapter_for_pe \
            --adapter_fasta {input.adapters} -w {threads}

        else
            zstd -dc {input.raw_rds} | \
            fastp --stdin -q {params.min_qual} -o {output.clean_rds} \
            -h {output.html_rep} --cut_front --cut_tail -l {params.min_len} \
            -M {params.min_qual} -u {params.max_unq} -j {output.json_rep} \
            --adapter_fasta {input.adapters} -w {threads}
        fi
        """

# -----------------------------------------------------------------------------
#### Overall subpipeline rule

rule subpipeline_read_preprocessing:
    input:
        expand("{lr_full_path}.filtered_reads.bam", lr_full_path=LR_PATHS),
        "results/ref/BBMap_adapters.fa",
        expand("{sr_full_path}_clean.fq.gz", sr_full_path=SR_PATHS)
    output:
        "logs/read_preprocessing/preprocessing.log"
    shell:
        "echo '***Reads successfully preprocessed!***' > {output}"