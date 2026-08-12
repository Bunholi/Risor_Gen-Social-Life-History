import os
from datetime import datetime

filename = 'populations'

# PRELIMINARY RUN USED TO GENERATE THE LOCUS WHITELIST
# Before the final runs configured below, populations was run once with
# popmap_zero.txt, which assigns every sample to a single population. This
# preliminary run did not use --whitelist or -r because its purpose was to
# produce populations.sumstats.tsv. The commands in
# Risor_Stacks_VCFtools_script.sh use that file to retain loci containing
# 1–10 SNPs and generate whitelist_loci.tsv.
#
# Configuration used for that preliminary run:
#   path_stacks = '../denovo.m5M6n7'
#   popmap = '../data/sequence_processing_inputs/popmap_zero.txt'
#   output = '../pop_zero'
#
# For the preliminary run, omit the `--whitelist {}` and `-r {}` portions of
# the command below. After whitelist_loci.tsv has been generated, restore
# those command portions and use the final-run configuration below.

# USER CONFIGURATION: update these paths for the run being reproduced.
# The settings below are for the final population-structure run (r = 0.70).
path_stacks = '../denovo.m5M6n7'
popmap = '../data/sequence_processing_inputs/popmap_batch1.txt'
output = '../denovo.m5M6n7/pop_wl_r07'
whitelist = '../denovo.m5M6n7/whitelist_loci.tsv'

# PUBLISHED PARAMETER: use 0.70 for the population-structure dataset.
# For the separate relatedness run, use 0.95, the relatedness population map
# (`../data/sequence_processing_inputs/popmap_batch1_relatedness.txt`, excluding the nine high-F
# individuals), and output name pan_populations.
min_samples = '0.70'

# Thread count is machine-specific and does not change the analysis.
threads = '8'

# Before running, confirm that path_stacks, popmap, output, whitelist, and
# min_samples describe the intended population-structure or relatedness run.

#log_name = datetime.now().strftime("%y%m%d_%H%M")

command = f"{filename} " \
        f"-P {path_stacks} " \
        f"-M {popmap} " \
        f"--whitelist {whitelist} " \
        f"-r {min_samples} " \
        f"-t {threads} " \
        f"-O {output} " \
        f"--vcf " \
        f"--structure " \
        f"--plink "

print(command)
os.system(command)

print("Finished run")
