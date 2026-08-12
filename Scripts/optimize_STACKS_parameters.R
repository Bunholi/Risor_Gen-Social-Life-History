# STACKS parameter evaluation for the Risor ruber ddRAD dataset
#
# Purpose
# -------
# This script documents the comparisons used to select the STACKS parameters
# M (within-individual mismatch distance) and n (between-individual catalog
# mismatch distance). The minimum stack depth m was fixed at 5 as a conservative setting recommended
# for relatedness inference. The final analysis used m = 5, M = 6, and n = 7.
# Those values are already set in denovo_map.py; this parameter-selection
# script is not required to rerun the final assembly.
#
# Inputs
# ------
# The VCF files below are the unfiltered populations.snps.vcf outputs from
# historical preliminary STACKS runs used to compare M = 1–8 and n = 6–8.
# These are the candidate files retained from the original analysis. 

library(vcfR)
library(RADstackshelpR)

# USER CONFIGURATION
# Change only these directory paths to match the local directory structure.
# No working directory is set; all input and output paths are built below.
input_dir <- file.path("..", "parameter_optimization_vcfs")
output_dir <- file.path("..", "outputs", "parameter_optimization")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_vcf <- function(filename) file.path(input_dir, filename)
output_file <- function(filename) file.path(output_dir, filename)

# PUBLISHED PARAMETER COMPARISONS
# Keep these filenames and parameter ranges unchanged when documenting the
# comparisons used for the manuscript.

#### 1. Evaluate big M: M = 1, 2, 3, 4, 5, 6, 7, and 8 ####

M.out <- optimize_bigM(
  M1 = input_vcf("bigM1.vcf"),
  M2 = input_vcf("bigM2.vcf"),
  M3 = input_vcf("bigM3.vcf"),
  M4 = input_vcf("bigM4.vcf"),
  M5 = input_vcf("bigM5.vcf"),
  M6 = input_vcf("bigM6.vcf"),
  M7 = input_vcf("bigM7.vcf"),
  M8 = input_vcf("bigM8.vcf")
)

png(output_file("vis_loci_bigM.png"), width = 1600, height = 1200, res = 200)
print(vis_loci(output = M.out, stacks_param = "M"))
dev.off()

#### 2. Evaluate candidate n values: n = 6, 7, and 8 ####

# The RADstackshelpR function requires three arguments named relative to M.
# The filenames supplied here reproduce the historical R analysis exactly:
# n6.vcf, n7.vcf, and n8.vcf. The final selected value was n = 7. 
n.out <- optimize_n(
  nequalsMminus1 = input_vcf("n6.vcf"),
  nequalsM = input_vcf("n7.vcf"),
  nequalsMplus1 = input_vcf("n8.vcf")
)

png(output_file("vis_loci_n.png"), width = 1600, height = 1200, res = 200)
print(vis_loci(output = n.out, stacks_param = "n"))
dev.off()
