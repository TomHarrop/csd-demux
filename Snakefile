#!/usr/bin/env python3

import os
import pathlib2


#############
# FUNCTIONS #
#############

def resolve_path(x):
    return str(pathlib2.Path(x).resolve())


###########
# GLOBALS #
###########

barcodes = 'data/EXP-PBC096_barcodes.fasta'
read_dir = 'data/raw'

# containers
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'
last_container = 'shub://TomHarrop/singularity-containers:last_973'
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'
pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
guppy_container = 'shub://TomHarrop/singularity-containers:guppy_2.3.5'

########
# MAIN #
########

# find raw fastq files
raw_fastq_files = []
for dirpath, dirnames, filenames in os.walk(read_dir):
    for filename in filenames:
        if filename.endswith('fastq'):
            raw_fastq_files.append(os.path.join(dirpath, filename))

#########
# RULES #
#########

# rule map_to_barcodes:
#     input:
#         filtered_reads = 'output/010_raw/filtered_reads.fastq',
#         barcodes = barcodes
#     output:
#         'output/020_mapped/aln.sam'
#     log:
#         'output/000_logs/020_map_to_barcodes.log'
#     threads:
#         50
#     singularity:
#         minimap_container
#     shell:
#         'minimap2 '
#         '-t {threads} '
#         '-ax map-ont '
#         '-N 1 '
#         '{input.barcodes} '
#         '{input.filtered_reads} '
#         '> {output} '
#         '2> {log}'

rule manual_demultiplex_guppy_results:
    input:
        guppy_results = 'output/035_guppy-manual-demux/barcoding_summary_filtered.txt',
        filtered_reads = 'output/010_raw/filtered_reads.fastq'
    output:
        expand('output/035_guppy-manual-demux/BC{bc}.fastq',
               bc=[f'{i:02}' for i in range(1, 97)])
    params:
        outdir = 'output/035_guppy-manual-demux'
    threads:
        50
    script:
        'src/manual_demultiplex_guppy_results.py'

rule filter_guppy_results:
    input:
        guppy_results = 'output/030_guppy-barcoder/barcoding_summary.txt',
    output:
        guppy_results = 'output/035_guppy-manual-demux/barcoding_summary_filtered.txt'
    log:
        'output/000_logs/035_filter_guppy_results.log'
    threads:
        1
    singularity:
        r_container
    script:
        'src/filter_guppy_results.R'

rule guppy_demux:
    input:
        read_dir
    output:
        'output/030_guppy-barcoder/barcoding_summary.txt'
    log:
        'output/000_logs/030_guppy-barcoder.log'
    threads:
        50
    singularity:
        guppy_container
    shell:
        'guppy_barcoder '
        '-t {threads} '
        '-i {input} '
        '-s {output} '
        '--barcode_kits "EXP-PBC096" '
        '&> {log}'

rule convert_to_tab:
    input:
        'output/020_mapped/aln.maf'
    output:
        'output/020_mapped/aln.tab'
    singularity:
        last_container
    shell:
        'maf-convert tab {input} > {output}'

rule map_to_barcodes:
    input:
        filtered_reads = 'output/010_raw/filtered_reads.fastq',
        db = 'output/020_mapped/barcodes.fasta'
    output:
        'output/020_mapped/aln.maf'
    log:
        'output/000_logs/020_map_to_barcodes.log'
    threads:
        50
    singularity:
        last_container
    shell:
        'lastal '
        '-P{threads} '
        '-Q 1 '
        '-N 1 '
        '{input.db} '
        '{input.filtered_reads} '
        '> {output} '
        '2> {log}'

rule lastdb:
    input:
        resolve_path(barcodes)
    output:
        expand('output/020_mapped/barcodes.fasta{suffix}',
               suffix=['', '.bck', '.des', '.prj', '.sds', '.ssp',
                       '.suf', '.tis'])
    params:
        wd = 'output/020_mapped'
    singularity:
        last_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'cp {input} ./barcodes.fasta ; '
        'lastdb barcodes.fasta barcodes.fasta'


rule filter_reads:
    input:
        summary = 'data/raw/sequencing_summary.txt',
        all_reads = 'output/010_raw/all_reads.fastq'
    output:
        to_keep = temp('output/010_raw/to_keep.txt'),
        filtered_reads = 'output/010_raw/filtered_reads.fastq'
    log:
        'output/000_logs/010_filter_reads.log'
    singularity:
        bbduk_container
    shell:
        'grep "TRUE" {input.summary} | cut -f 2 > {output.to_keep} '
        '; '
        'filterbyname.sh '
        'in={input.all_reads} '
        'names={output.to_keep} '
        'include=t '
        'out={output.filtered_reads} '
        '&> {log}'

rule combine_reads:
    input:
        raw_fastq_files
    output:
        'output/010_raw/all_reads.fastq'
    shell:
        'cat {input} > {output}'


