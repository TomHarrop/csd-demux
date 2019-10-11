#!/usr/bin/env python3

import os
import pathlib2
import multiprocessing


#############
# FUNCTIONS #
#############

# function for checking which samples were demuxed
def aggregate_input(wildcards):
    print(wildcards)
    checkpoint_output = checkpoints.manual_demultiplex_guppy_results.get(
        **wildcards).output['outdir']
    return expand('output/050_mapped/BC{bc}_sorted.bam',
                  bc=glob_wildcards(
                    os.path.join(checkpoint_output, "barcode{bc}.fastq")).bc)


def resolve_path(x):
    return str(pathlib2.Path(x).resolve())


###########
# GLOBALS #
###########

barcodes = 'data/EXP-PBC096_barcodes.fasta'
read_dir = 'data/raw2'
honeybee_ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
csd_chrom = "NC_037640.1"

# containers
r_container = 'shub://TomHarrop/singularity-containers:r_3.5.1'
last_container = 'shub://TomHarrop/singularity-containers:last_973'
minimap_container = 'shub://TomHarrop/singularity-containers:minimap2_2.11r797'
pigz_container = 'shub://TomHarrop/singularity-containers:pigz_2.4.0'
bwa_container = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
guppy_container = 'shub://TomHarrop/singularity-containers:guppy_2.3.7'
sambamba_container = 'shub://TomHarrop/singularity-containers:sambamba_0.6.8'
freebayes_container = 'shub://TomHarrop/singularity-containers:freebayes_1.2.0'
vcflib_container = 'shub://TomHarrop/singularity-containers:vcflib_1.0.0-rc2'
minionqc_container = 'shub://TomHarrop/singularity-containers:minionqc_1.4.1'
qcat_container = 'shub://TomHarrop/singularity-containers:qcat_1.0.1'

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

rule target:
    input:
        'output/015_minionqc',
        'output/060_freebayes/long_haplotypes.vcf'

rule extract_csd_chrom:
    input:
        vcf = 'output/060_freebayes/variants.vcf'
    output:
        'output/060_freebayes/variants_csd_chrom.vcf'
    shell:
        'grep "^#" {input.vcf} > {output} ; '
        'grep {csd_chrom} {input.vcf} >> {output}'


rule freebayes_lr:
    input:
        bam = aggregate_input,
        vcf = 'output/060_freebayes/variants_filtered.vcf',
        fa = honeybee_ref
    output:
        bg = 'output/060_freebayes/variants_filtered.vcf.gz',
        index = 'output/060_freebayes/variants_filtered.vcf.gz.tbi',
        vcf = 'output/060_freebayes/long_haplotypes.vcf'
    params:
        haplotype_length = '10'
    log:
        'output/000_logs/060_freebayes/freebayes.log'
    singularity:
        freebayes_container
    shell:
        'bgzip --force --stdout {input.vcf} > {output.bg} ; '
        'tabix -p vcf {output.bg} ; '
        'freebayes '
        '-f {input.fa} '
        '--haplotype-length {params.haplotype_length} '
        '--haplotype-basis-alleles {output.bg} '
        '{input.bam} '
        '> {output.vcf} '
        '2> {log}'


rule filter_vcf:
    input:
        'output/060_freebayes/variants.vcf'
    output:
        'output/060_freebayes/variants_filtered.vcf'
    params:
        filter = "QUAL > 20"
    log:
        'output/000_logs/060_freebayes/freebayes_filter.log'
    singularity:
        vcflib_container
    shell:
        'vcffilter -f "{params.filter}" {input} > {output} 2> {log}'


rule freebayes:
    input:
        bam = aggregate_input,
        fa = honeybee_ref
    output:
        vcf = 'output/060_freebayes/variants.vcf'
    log:
        'output/000_logs/060_freebayes/freebayes.log'
    singularity:
        freebayes_container
    shell:
        'freebayes '
        '--region NC_037640.1:11771679-11781139 '
        '-f {input.fa} '
        '{input.bam} '
        '> {output} '
        '2> {log}'


rule sort:
    input:
        'output/050_mapped/BC{bc}.sam'
    output:
        bam = 'output/050_mapped/BC{bc}_sorted.bam'
    log:
        'output/000_logs/050_mapped/{bc}_sort.log'
    threads:
        1
    singularity:
        sambamba_container
    shell:
        'sambamba view '
        '{input} '
        '-f "bam" '
        '--sam-input '
        '-l 0 '
        '2> {log} '
        '| '
        'sambamba sort '
        '-o {output.bam} '
        '-l 9 '
        '-t {threads} '
        '/dev/stdin '
        '2>> {log} '
        '; '
        'sambamba index '
        '{output.bam} '
        '2>> {log}'

rule map_to_csd:
    input:
        fq = 'output/035_guppy-manual-demux/demuxed/barcode{bc}.fastq',
        ref = 'output/010_raw/honeybee_ref.mmi'
    output:
        temp('output/050_mapped/BC{bc}.sam')
    params:
        rg = '\'@RG\\tID:BC{bc}\\tSM:BC{bc}\''
    log:
        'output/000_logs/050_mapped/{bc}.log'
    threads:
        1
    singularity:
        minimap_container
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax map-ont '
        '-N 1 '
        '-R {params.rg} '
        '{input.ref} '
        '{input.fq} '
        '> {output} '
        '2> {log}'

rule prepare_ref:
    input:
        honeybee_ref
    output:
        'output/010_raw/honeybee_ref.mmi'
    log:
        'output/000_logs/010_raw/prepare_ref.log'
    threads:
        3
    singularity:
        minimap_container
    shell:
        'minimap2 '
        '-x map-ont '
        '-d {output} '
        '{input} '
        '2> {log}'


# everything from here on is testing
rule demultiplex_stats:
    input:
        'output/035_guppy-manual-demux/demuxed/barcode{bc}.fastq'
    output:
        'output/040_stats/BC{bc}_readlength.txt'
    log:
        'output/000_logs/040_stats/BC{bc}_readlength.log'
    threads:
        1
    singularity:
        bbduk_container
    shell:
        'readlength.sh '
        'in={input} '
        'out={output} '
        'binsize=50 '
        '2> {log}'

checkpoint manual_demultiplex_guppy_results:
    input:
        guppy_results = 'output/035_guppy-manual-demux/barcoding_summary_filtered.txt',
        filtered_reads = 'output/010_raw/filtered_reads.fastq'
    output:
        outdir = directory('output/035_guppy-manual-demux/demuxed')
    threads:
        multiprocessing.cpu_count()
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
    params:
        outdir = 'output/030_guppy-barcoder'
    threads:
        multiprocessing.cpu_count()
    singularity:
        guppy_container
    shell:
        'guppy_barcoder '
        '-t {threads} '
        '-i {input} '
        '-s {params.outdir} '
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
        multiprocessing.cpu_count()
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

rule qcat:
    input:
        all_reads = 'output/010_raw/all_reads.fastq'
    output:
        bc_dir = directory('output/037_qcat/demuxed'),
        tsv = 'output/037_qcat/barcodes.tsv'
    log:
        'output/000_logs/037_qcat.log'
    threads:
        multiprocessing.cpu_count()
        # 1
    singularity:
        qcat_container
    shell:
        'qcat '
        '--fastq {input.all_reads} '
        '--barcode_dir {output.bc_dir} '
        '--detect-middle '
        '--threads {threads} '
        # '--tsv '
        '--trim '
        '--kit PBC096 '
        '--guppy '
        # '--epi2me '
        '> {output.tsv} '
        '2> {log}'

rule combine_reads:
    input:
        raw_fastq_files
    output:
        'output/010_raw/all_reads.fastq'
    shell:
        'cat {input} > {output}'


rule minionqc:
    input:
        str(pathlib2.Path(read_dir, 'sequencing_summary.txt'))
    output:
        directory('output/015_minionqc')
    threads:
        1
    priority:
        1
    singularity:
        minionqc_container
    log:
        'output/000_logs/015_minionqc.log'
    shell:
        'MinIONQC.R '
        '--input={input} '
        '--outputdirectory={output} '
        '&> {log}'

