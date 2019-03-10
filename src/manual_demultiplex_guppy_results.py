#!/usr/bin/env python3

from Bio import SeqIO
from functools import partial
from multiprocessing import Pool
import pandas

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
        guppy_filtered):
    print('Working on {}'.format(row['barcode_front_bc']))
    my_filename = 'test2/{}.fastq'.format(row['barcode_front_bc'])
    my_readlist = row[0]
    # my_seq = []
    # for rec in my_readlist:
    #     my_seq.append(manual_trim_by_location(
    #         rec,
    #         seq_index,
    #         guppy_filtered))
    # return({row['barcode_front_bc']: my_seq})
    # SeqIO.write(my_seq, my_filename, 'fastq')
    with open(my_filename, 'wt') as f:
        for rec in my_readlist:
            my_seq = manual_trim_by_location(
                rec,
                seq_index,
                guppy_filtered)
            SeqIO.write(my_seq, f, 'fastq')


###############
# GLOBALS     #
# (fix later) #
###############

reads = 'output/010_raw/filtered_reads.fastq'
# reads = 'test2/reads.fastq'
guppy_summary = 'output/030_guppy-barcoder/barcoding_summary.txt'
threads = 8

########
# MAIN #
########


def main():
    # index the fastq. 
    # have to read the whole thing into memory if we want to serialize it.
    seq_index = SeqIO.to_dict(SeqIO.parse(reads, 'fastq'))
    all_ids = list(seq_index.keys())

    # read the guppy results
    guppy = pandas.read_csv(guppy_summary, sep='\s+')

    # only work on reads that are in the q-filtered fastq
    guppy_filtered = guppy.loc[guppy['read_id'].isin(all_ids)]

    # which barcode was recognised at each end of read
    guppy_filtered['barcode_front_bc'] = guppy_filtered.apply(
        lambda row: row['barcode_front_id'].split("_")[0], axis=1)
    guppy_filtered['barcode_rear_bc'] = guppy_filtered.apply(
        lambda row: row['barcode_rear_id'].split("_")[0], axis=1)

    # get a subset of reads where front barcode and rear barcode match
    matched_bc = guppy_filtered.loc[
        guppy_filtered['barcode_front_bc'] ==
        guppy_filtered['barcode_rear_bc']]

    # group by barcode
    matched_bc_grouped = matched_bc.groupby('barcode_front_bc')

    # extract the read IDs in each group
    reads_by_bc = matched_bc_grouped.apply(
        lambda group: list(group['read_id']))

    # write the reads per barcode
    # reads_by_bc.reset_index().apply(write_fastq_from_row, axis=1)

    # run write_fastq_from_row in parallel
    # joblib_results = joblib.Parallel(
    #     n_jobs=threads,
    #     prefer='threads')(
    #         joblib.delayed(write_fastq_from_row)(y)
    #         for x, y in reads_by_bc.reset_index().iterrows())

    rows = (y for x, y in reads_by_bc.reset_index().iterrows())
    with Pool(threads) as pool:
        all_res = pool.map(partial(
            write_fastq_from_row,
            seq_index=seq_index,
            guppy_filtered=guppy_filtered), [x for x in rows])


if __name__ == '__main__':
    main()
