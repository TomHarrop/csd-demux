#!/usr/bin/env python3

from Bio import SeqIO
import pandas

#############
# FUNCTIONS #
#############

# manual trimming with exact match
def manual_trim_by_location(seq_id):
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


def write_fastq_from_row(row):
    my_filename = 'test2/{}.fastq'.format(row['barcode_front_bc'])
    my_readlist = row[0]
    with open(my_filename, 'wt') as f:
        for rec in my_readlist:
            my_seq = manual_trim_by_location(rec)
            SeqIO.write(my_seq, f, 'fastq')


###############
# GLOBALS     #
# (fix later) #
###############

reads = 'output/010_raw/filtered_reads.fastq'
guppy_summary = 'output/030_guppy-barcoder/barcoding_summary.txt'

########
# MAIN #
########

# index the fastq
seq_index = SeqIO.index(reads, 'fastq')
all_ids = list(seq_index.keys())

# read the guppy results
guppy = pandas.read_csv(guppy_summary, sep='\s+')

# only work on reads that are in the q-filtered fastq
guppy_filtered = guppy.loc[guppy['read_id'].isin(all_ids)]

# which barcode was recognised at each end of read
guppy_filtered['barcode_front_bc'] = guppy_filtered.apply(lambda row: row['barcode_front_id'].split("_")[0], axis=1)
guppy_filtered['barcode_rear_bc'] = guppy_filtered.apply(lambda row: row['barcode_rear_id'].split("_")[0], axis=1)

# get a subset of reads where front barcode and rear barcode match
matched_bc = guppy_filtered.loc[
    guppy_filtered['barcode_front_bc'] == guppy_filtered['barcode_rear_bc']]

# group by barcode
matched_bc_grouped = matched_bc.groupby('barcode_front_bc')

# extract the read IDs in each group
reads_by_bc = matched_bc_grouped.apply(lambda group: list(group['read_id']))

# write the reads per barcode
reads_by_bc.reset_index().apply(write_fastq_from_row, axis=1)


