import os

filename = 'denovo_map.pl'

# PUBLISHED PARAMETERS: these are the final STACKS values used for the
# manuscript and should not be changed for reproduction. The value of m was
# fixed at 5; M and n were selected through the documented optimization.
m = '5'
M_ustacks = '6'
n_cstacks = '7'

# USER CONFIGURATION: update only these paths and, if needed, the number of
# threads. "samples" must point to the per-sample files produced by
# process_radtags (deposited under NCBI BioProject [PRJNA1427305])
samples = '../seq_batch1_demultiplexed'
popmap = '../data/sequence_processing_inputs/popmap_batch1.txt'
output = '../denovo.m5M6n7'
threads = '10'

command = "{} " \
        "-m {} " \
        "-M {} " \
        "-n {} " \
        "--samples {} " \
        "--popmap {} " \
        "-o {} " \
        "-T {} " \
        "--paired ".format(
            filename, m, M_ustacks, n_cstacks,
            samples, popmap, output, threads
        )

print(command)
os.system(command)

print("Finished run")
