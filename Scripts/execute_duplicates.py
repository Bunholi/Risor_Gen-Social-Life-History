import os
from datetime import datetime

# HISTORICAL PREPROCESSING STEP
# This Python 2.7 wrapper documents how PCR duplicates were removed from one
# original pooled read pair using ParseDBR_ddRAD (commit 8fc0edc).
# Users starting from the sample-specific FASTQ files in SRA BioProject PRJNA1427305 should NOT run this script and
# should begin with the de novo STACKS assembly documented in denovo_map.py.

# USER CONFIGURATION
# For each original pooled read pair, update only the ParseFastQ.py location
# and the input/output paths below. These paths are examples of the required
# directory structure and are not analytical parameters.
filename = '../ParseDBR_ddRAD/ParseFastQ.py'
read_file_1 = '../seq_batch1/POOLED_READ1.fq.gz'
read_file_2 = '../seq_batch1/POOLED_READ2.fq.gz'
read1_out = '../seq_batch1/POOLED_READ1_NODP.fq'
read2_out = '../seq_batch1/POOLED_READ2_NODP.fq'
dropped_reads = '../seq_batch1/dropped_reads.txt'

# PUBLISHED LIBRARY SETTINGS
# The original pools were processed separately with their corresponding DBR
# index and adapter-length values.
# The active values below document the DBR01 configuration.
index_seq = 'ATCACG'
adpt_len = '0'
r2_enzyme = 'CG'

log_name = datetime.now().strftime("%y%m%d_%H%M")

command = "python {} " \
        "-r {} " \
        "-R {} " \
        "-i {} " \
        "-e {} " \
        "-n {} " \
        "-N {} " \
        "--drop {} " \
        "-Z " \
        "-l {} " \
        "> {}.log".format(filename, read_file_1, read_file_2, index_seq, r2_enzyme, read1_out, read2_out, dropped_reads, adpt_len, log_name)

os.system(command)

print("Finished run")
