###############################################################################
# File: read_preprocessing.smk
# Version: 0.4
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
#   2026-11-12
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

# From Config Files
long_reads = lookup(within=config, dpath="get_reads/lr")
short_reads = lookup(within=config, dpath="get_reads/sr")

def get_accessions(accession_file: str):
    return pl.read_csv(accession_file).select('acc').to_series().to_list()

LR_PATHS = ['results/data/lr/' + str(i) for i in get_accessions(long_reads)]
SR_PATHS = ['results/data/sr/' + str(i) + '/' + str(i) 
            for i in get_accessions(short_reads)]


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
        --include {input.met}.incl.mask {input.hifi_bam}.excl.temp.bam \
        {output}

        pbindex -j {threads} {output}

        rm {input.hifi_bam}.excl.mask {input.hifi_bam}.excl.temp.bam
        rm {input.ccs_bam} {input.ccs_bam}.pbi {input.met}.incl.mask
        rm {input.hifi_bam}.excl.temp.bam.pbi
        """

rule sr_bbduk_decontamination:
    params:
        compression_level = 6,
        maq = 20,
        minlen = 80
    input: 
        raw_rds = "results/data/sr/{sr_accID}/{sr_accID}.fq.zst"
    output:
        bb_clean = "results/data/sr/{sr_accID}/{sr_accID}_bbduk.fq.gz",
        bb_contamination = "results/data/sr/{sr_accID}/{sr_accID}_bbduk_contam.fq.gz",
        contamination_stats="results/data/sr/{sr_accID}/{sr_accID}_bbduk.stats",
    log:
        "logs/read_preprocessing/sr/bbduk_{sr_accID}.log"
    conda:
        "../envs/raw_data_preprocessing.yaml"
    threads:
        workflow.cores * 0.9
    message:
        "Performing Illumina adapter trimming and \
        sequencing contaminants removal in short reads using BBDuk."
    shell:
        """
        if [ "{wildcards.sr_accID}" == "SRR14352918" ]; then
            zstd -dc {input.raw_rds} | \
            bbduk.sh in=stdin.fq out={output.bb_clean} interleaved=t \
            outm={output.bb_contamination} \
            zl={params.compression_level} \
            ref=adapters,phix,pjet,mtst,kapa,artifacts \
            ktrim=r k=27 mink=11 hdist=1 qtrim=rl trimq=10 ftm=5 \
            minlen={params.minlen} maq={params.maq} \
            tpe tbo stats={output.contamination_stats} threads={threads}

        else
            zstd -dc {input.raw_rds} | \
            bbduk.sh in=stdin.fq out={output.bb_clean} interleaved=f \
            outm={output.bb_contamination} \
            zl={params.compression_level} \
            ref=adapters,phix,pjet,mtst,kapa,artifacts \
            ktrim=r k=27 mink=11 hdist=1 qtrim=rl trimq=10 ftm=5 \
            minlen={params.minlen} maq={params.maq} \
            tpe tbo stats={output.contamination_stats} threads={threads}
        fi
        """

rule sr_fastp_filtering:
    params:
        min_qual = 30,
        min_len = 80,
        max_unq = 5
    input:
        bb_clean = rules.sr_bbduk_decontamination.output.bb_clean,
    output:
        html_rep = "results/data/sr/{sr_accID}/{sr_accID}.html",
        json_rep = "results/data/sr/{sr_accID}/{sr_accID}.json"
    log:
        "logs/read_preprocessing/sr/fastp_{sr_accID}.log"
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
            fastp --interleaved_in -i {input.bb_clean} \
            -o "results/data/sr/{wildcards.sr_accID}/{wildcards.sr_accID}_bb_fastp_R1.fq.gz" \
            -O "results/data/sr/{wildcards.sr_accID}/{wildcards.sr_accID}_bb_fastp_R2.fq.gz" \
            -q {params.min_qual} -h {output.html_rep} --cut_front \
            --cut_tail -l {params.min_len} -M {params.min_qual} \
            -u {params.max_unq} -j {output.json_rep} -w {threads}

        else
            fastp -i {input.bb_clean} \
            -o "results/data/sr/{wildcards.sr_accID}/{wildcards.sr_accID}_bb_fastp.fq.gz" \
            -q {params.min_qual} -h {output.html_rep} --cut_front --cut_tail \
            -l {params.min_len} -M {params.min_qual} -u {params.max_unq} \
            -j {output.json_rep} -w {threads}
        fi

        rm {input.bb_clean}
        """

# -----------------------------------------------------------------------------
#### Overall subpipeline rule

rule subpipeline_read_preprocessing:
    input:
        expand("{lr_full_path}.filtered_reads.bam", lr_full_path=LR_PATHS),
        expand("{sr_full_path}.html", sr_full_path=SR_PATHS)
    output:
        "logs/read_preprocessing/preprocessing.log"
    shell:
        "echo '***Reads successfully preprocessed!***' > {output}"