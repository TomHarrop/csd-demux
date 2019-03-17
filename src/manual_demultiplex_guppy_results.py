#!/usr/bin/env python3

from Bio import SeqIO
from functools import partial
from multiprocessing import Pool
import pandas
import pathlib2

#############
# FUNCTIONS #
#############

# manual trimming with exact match
def manual_trim_by_location(
        seq_id,
        seq_index,
        guppy_filtered):
    # extract the seq from seq_index and the row from guppy
    my_seq = seq_index[seq_id]
    my_row = guppy_filtered.loc[guppy_filtered['read_id'] == seq_id]
    # extract the found barcode sequences
    my_f_bc = my_row['barcode_front_foundseq'].iloc[0]
    my_r_bc = my_row['barcode_rear_foundseq'].iloc[0]
    # find where the barcodes are in the read
    my_seqstart = my_seq.seq.find(my_f_bc) + len(my_f_bc)
    my_seqend = my_seq.seq.find(my_r_bc) - 1
    # trim
    my_trimmed_seq = my_seq[my_seqstart:my_seqend]
    return(my_trimmed_seq)


def write_fastq_from_row(
        row,
        seq_index,
        guppy_filtered,
        outdir):
    my_filename = '{0}/{1}.fastq'.format(outdir, row['barcode_front_bc'])
    print('Demultiplexing {0} reads to {1}'.format(
        row['barcode_front_bc'],
        my_filename))
    my_readlist = row[0]
    with open(my_filename, 'wt') as f:
        for rec in my_readlist:
            if rec in seq_index.keys():
                my_seq = manual_trim_by_location(
                    rec,
                    seq_index,
                    guppy_filtered)
                SeqIO.write(my_seq, f, 'fastq')


###############
# GLOBALS     #
# (fix later) #
###############

# reads = 'test2/reads.fastq'
# reads = 'output/010_raw/filtered_reads.fastq'
# guppy_summary = 'output/030_guppy-barcoder/barcoding_summary.txt'
# threads = 8

reads = snakemake.input['filtered_reads']
guppy_summary = snakemake.input['guppy_results']
outdir = snakemake.output['outdir']
threads = snakemake.threads

########
# MAIN #
########


def main():
    # check for the output directory
    with pathlib2.Path(outdir) as p:
        if not p.is_dir():
            p.mkdir(parents=True)

    # index the fastq. 
    # have to read the whole thing into memory if we want to serialize it.
    seq_index = SeqIO.to_dict(SeqIO.parse(reads, 'fastq'))

    # read the guppy results
    guppy_filtered = pandas.read_csv(guppy_summary, sep='\s+')

    # group by barcode
    matched_bc_grouped = guppy_filtered.groupby('barcode_front_bc')

    # extract the read IDs in each group
    reads_by_bc = matched_bc_grouped.apply(
        lambda group: list(group['read_id']))

    # run write_fastq_from_row in parallel for each row
    rows = (y for x, y in reads_by_bc.reset_index().iterrows())
    with Pool(threads) as pool:
        pool.map(partial(
            write_fastq_from_row,
            seq_index=seq_index,
            guppy_filtered=guppy_filtered,
            outdir=outdir), rows)


if __name__ == '__main__':
    main()
