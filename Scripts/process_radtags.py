import os
from datetime import datetime

filename = 'process_radtags'
input_type = 'fastq'

# HISTORICAL PREPROCESSING STEP
# This wrapper documents how the original pooled batch-1 reads were
# demultiplexed and quality filtered. Users starting from the sample-specific paired FASTQ files in
# SRA BioProject PRJNA1427305 should NOT run this script; they should begin
# with the denovo STACKS assembly documented in denovo_map.py.

# USER CONFIGURATION: for the original processing, these paths were changed
# to identify the matching pooled read pair, output directory, and barcode
# file.
read_file_1 = '../seq_batch1/POOLED_READ1_NODP.fq'
read_file_2 = '../seq_batch1/POOLED_READ2_NODP.fq'
output = '../seq_batch1_demultiplexed'

# The P1 and P8 barcode combinations were processed separately. Set `barcode`
# to the file corresponding to the pooled read pair being processed, and run
# the script once for each applicable combination.
barcode = '../data/sequence_processing_inputs/barcode_1st_P1.txt'
# barcode = '../data/sequence_processing_inputs/barcode_1st_P8.txt'  # use for the matching P8 reads

# PUBLISHED PARAMETERS: MluCI and MspI were used for the published analysis.
enzyme_1 = 'mluCI'
enzyme_2 = 'mspI'

#log_name = datetime.now().strftime("%y%m%d_%H%M")

command = "{} " \
        "-i {} " \
        "-1 {} " \
        "-2 {} " \
        "-o {} " \
        "-b {} " \
        "-c " \
        "-q " \
        "-r " \
        "--disable-rad-check " \
        "--paired " \
        "--inline_index " \
        "--renz_1 {} " \
        "--renz_2 {} ".format(filename, input_type, read_file_1, read_file_2, output, barcode, enzyme_1, enzyme_2)

os.system(command)

print("Finished run")
