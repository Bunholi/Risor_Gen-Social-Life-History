#!/bin/bash

# ddRAD Risor workflow: public SRA reads to final filtered VCFs

# This file contains the workflow needed to reproduce the genomic
# processing from the sample-specific paired FASTQ files deposited under NCBI
# BioProject PRJNA1427305. Historical pooled-read preprocessing is documented
# separately in the README, execute_duplicates.py, and process_radtags.py.
#
# USER CONFIGURATION
# Users may change only local paths, filenames, output directories, conda
# environment names, and thread counts to match their computing environment.
# These changes do not alter the analysis.
#
# PUBLISHED ANALYTICAL PARAMETERS
# The following are the exact settings used for the manuscript and should not
# be changed when reproducing the published results:
#   STACKS assembly: m=5, M=6, n=7
#   locus presence: r=0.70 (population structure), r=0.95 (relatedness)
#   minimum allele count: MAC=5 (population structure), MAC=16 (relatedness)
#   excess-heterozygosity threshold: P<1e-4
#   LD pruning: r2=0.2, window=300 bp, retain one site per window


## 1 - Run the final de novo STACKS assembly

# Download the sample-specific paired FASTQs from PRJNA1427305 and place them
# in the sample directory specified at the top of denovo_map.py. The supplied
# popmap_batch1.txt defines the exact samples retained in the published
# assembly. Confirm only the sample, popmap, output, and thread settings in
# denovo_map.py; its published m=5, M=6, and n=7 values should not be changed.

# Activate a local environment containing STACKS 2.65.
# Example from the authors’ server:
conda activate stacks_2.65
python denovo_map.py

### Check assembly output

### Sequence coverage from USTACKS

stacks-dist-extract denovo_map.log cov_per_sample |
  grep -v '^#' > cov_per_sample_m5M6n7.tsv

### Sequence coverage from GSTACKS and other important results

cat gstacks.log | grep -A 3 'Genotyped' # Overall coverage

stacks-dist-extract gstacks.log.distribs effective_coverages_per_sample |
  grep -v '^#' > effective_cov_per_sample_m5M6n7.tsv # Per-sample coverage


## 2 - Run populations and generate the locus whitelist

# Refer to populations_popgen.py.
# Three populations runs were used:
#   1. Preliminary whitelist run: popmap_zero.txt, without --whitelist or -r;
#      populations.sumstats.tsv was used to generate whitelist_loci.tsv.
#   2. Population-structure run: popmap_batch1.txt, whitelist_loci.tsv, r=0.70;
#      output populations.snps.vcf.
#   3. Relatedness run: popmap_batch1_relatedness.txt, whitelist_loci.tsv,
#      r=0.95; output pan_populations.snps.vcf.

### First run

# Run without additional filtering parameters and use a popmap with only
# one population (popmap_zero.txt). Configure the preliminary run exactly as
# documented at the top of populations_popgen.py, omit --whitelist and -r,
# and run the wrapper to generate populations.sumstats.tsv:
python populations_popgen.py

### Build whitelist

cat populations.sumstats.tsv |
  grep -v "^#" |
  cut -f 1 |
  uniq -c |
  sort |
  awk '$1 >= 1 && $1 <= 10' |
  sort -n > whitelist.tsv

awk '{sub(/.* /, "", $0); print}' whitelist.tsv > whitelist_loci.tsv

### Second run

# Run populations with the whitelist, the full population map (popmap_batch1.txt), and r=0.70.
# This full-sample run was evaluated first. The separate relatedness dataset
# was constructed later, after individual inbreeding coefficients revealed
# nine individuals with unusually high F values (see below).

# Set the configuration block in populations_popgen.py for the population-
# structure run (r=0.70), then run:
python populations_popgen.py

## 3 - VCFtools filtering: full population-structure dataset

# Activate a local environment containing the reported VCFtools and BCFtools versions.
# The environment name below is specific to the authors’ server.
conda activate vcftools

vcftools --vcf populations.snps.vcf #version 0.1.14

### Retain biallelic loci

vcftools \
  --vcf populations.snps.vcf \
  --max-alleles 2 \
  --min-alleles 2 \
  --recode \
  --out pop_m5_r07_wl_biallelic

### Apply minimum allele count filters

# MAC = 5 final full population-structure dataset

vcftools \
  --vcf pop_m5_r07_wl_biallelic.recode.vcf \
  --mac 5 \
  --recode \
  --out pop_m5_r07_wl_biallelic_mac5


### Hardy-Weinberg equilibrium filtering

# Evaluate the MAC 5 population-structure dataset

vcftools \
  --vcf pop_m5_r07_wl_biallelic_mac5.recode.vcf \
  --hardy \
  --out pop_m5_r07_hwe_mac5

# Excess-heterozygosity outliers:
# p < 1e-4 and observed heterozygosity > expected heterozygosity

awk 'NR>1 && $NF != "nan" && $NF < 1e-4 {print $1"\t"$2}' \
  pop_m5_r07_hwe_mac5.hwe |
  sort -u > pop_m5_r07_hwe_excess_mac5_p1e4.sites


### Remove excess-heterozygosity loci

# MAC 5 final full population-structure dataset using p < 1e-4

vcftools \
  --vcf pop_m5_r07_wl_biallelic_mac5.recode.vcf \
  --exclude-positions pop_m5_r07_hwe_excess_mac5_p1e4.sites \
  --recode \
  --recode-INFO-all \
  --out pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4

### Check homozygosity and inbreeding coefficients

vcftools \
  --vcf pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4.recode.vcf \
  --het \
  --out post_mac5


### Identify individuals with high inbreeding coefficients

# Inspection of post_mac5.het showed that the following nine samples had
# F > 0.9: golden_2, golden_3, golden_21, golden_64, grunt_127, grunt_159,
# grunt_181, long_209, and long_219.
#
# This result led to the decision to exclude these samples from the
# relatedness analysis. The decision was therefore made after the first
# full-sample populations run and the homozygosity/inbreeding assessment,
# not before the initial STACKS analysis.
#
# A separate reduced population map, popmap_batch1_relatedness.txt, was made
# by removing these nine rows from the full population map. All retained
# population assignments remained unchanged. No separate orange_clade.txt
# input is required to reproduce either final VCF.

# Check depth

vcftools \
  --vcf pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4.recode.vcf \
  --depth \
  --out depth


### LD pruning

# MAC 5 final full population-structure dataset

bcftools +prune \
  pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4.recode.vcf \
  -Ov \
  -m r2=0.2 \
  -w 300 \
  -n 1 \
  -o pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4_LDpruned.vcf

vcftools \
  --vcf pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4_LDpruned.vcf #final VCF file for population-structure analysis


## 4 - VCFtools filtering: relatedness dataset

# After the high-F values were identified in the full-sample dataset above,
# the nine affected individuals were removed from the dataset for the relatedness analysis,
# which is now represented by popmap_batch1_relatedness.txt. This reduced population map was supplied to
# STACKS populations, generating pan_populations.snps.vcf.
# Loci were required to be present in at least 95% of retained individuals.

# Use the published relatedness configuration documented in
# populations_popgen.py: r=0.95, popmap=popmap_batch1_relatedness.txt, and
# output=pan_populations. These analytical settings should not be changed;
# only their local paths may need adjustment. Then run:
python populations_popgen.py

### Retain biallelic loci and apply the conservative MAC = 16 filter

# MAC = 16 was selected for the final relatedness analysis as a conservative
# allele-count filter to reduce the number of rare alleles and minimize the potential for false-positive relatedness estimates.

vcftools \
  --vcf pan_populations.snps.vcf \
  --max-alleles 2 \
  --min-alleles 2 \
  --mac 16 \
  --recode \
  --out rel_m5_r095_wl_biallelic_mac16

### Hardy-Weinberg equilibrium test

vcftools \
  --vcf rel_m5_r095_wl_biallelic_mac16.recode.vcf \
  --hardy \
  --out rel_m5_r095_hwe_mac16

# Excess-heterozygosity outliers:
# p < 1e-4 and observed heterozygosity > expected heterozygosity

awk 'NR>1 && $NF != "nan" && $NF < 1e-4 {print $1"\t"$2}' \
  rel_m5_r095_hwe_mac16.hwe |
  sort -u > rel_m5_r095_hwe_excess_mac16_p1e4.sites


### Remove excess-heterozygosity loci

vcftools \
  --vcf rel_m5_r095_wl_biallelic_mac16.recode.vcf \
  --exclude-positions rel_m5_r095_hwe_excess_mac16_p1e4.sites \
  --recode \
  --recode-INFO-all \
  --out rel_m5_r095_wl_biallelic_mac16.noHWEexcess_p1e4


### Check homozygosity and inbreeding coefficients

vcftools \
  --vcf rel_m5_r095_wl_biallelic_mac16.noHWEexcess_p1e4.recode.vcf \
  --het \
  --out rel_post


### LD pruning

bcftools +prune \
  rel_m5_r095_wl_biallelic_mac16.noHWEexcess_p1e4.recode.vcf \
  -Ov \
  -m r2=0.2 \
  -w 300 \
  -n 1 \
  -o rel_m5_r095_wl_biallelic_mac16.noHWEexcess_p1e4.LDpruned.vcf

vcftools \
  --vcf rel_m5_r095_wl_biallelic_mac16.noHWEexcess_p1e4.LDpruned.vcf #final VCF file for relatedness analysis
