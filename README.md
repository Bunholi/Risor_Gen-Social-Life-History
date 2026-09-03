# Risor_Social-Life-History

This README file was updated on 2026-09-03 by Ingrid Vasconcellos Bunholi.

This archive contains the processed data, metadata, bioinformatics workflows, and R scripts used to reproduce the analyses, tables, and figures presented in Bunholi et al. (2026), *Ecological constraints of habitat specialization shape the social organization and life-history of a sponge-dwelling goby* - The American Naturalist.

## GENERAL INFORMATION

**Title of dataset:** Risor_Gen-Social-Life-History

**Brief description of the study:** This study examines how specialization on sponge habitats shapes habitat occupancy, social organization, relatedness, and life-history traits in the sponge-dwelling goby *Risor ruber*. Fish and host sponges were sampled at four patch reefs near Carrie Bow Cay, Belize, in July 2023.

**Name:** Ingrid Vasconcellos Bunholi
**ORCID:** 0000-0001-5489-276X 
**Email:** ingrid.bunholi@utexas.edu
**Institution:** The University of Texas at Austin Marine Science Institute
**Person responsible for code:** Ingrid Vasconcellos Bunholi
**Person responsible for data:** Ingrid Vasconcellos Bunholi

## CONTENTS OF THIS ARCHIVE

The project includes:

1. **Analysis data:** Input files are organized by workflow under `data/social_life_history/` and `data/relatedness_popgen/`.
2. **Sequence-processing inputs:** Population maps and historical barcode files are provided separately in `data/sequence_processing_inputs/`.
3. **Bioinformatics workflows:** Scripts and documented commands used to process the raw ddRAD sequence data with STACKS, VCFtools, and associated software.
4. **Statistical scripts:** R Markdown files used to perform the analyses and reproduce the tables, figures, and summary statistics presented in the manuscript.

Sequence reads are not included in this archive because of their file size. Sample-specific paired-end FASTQ files are deposited in the NCBI Sequence Read Archive under BioProject [PRJNA1427305](https://www.ncbi.nlm.nih.gov/sra/PRJNA1427305).

## DATA AND FILE OVERVIEW

The archive supports three independent starting points:

1. **Social and life-history analyses:** Run `Social_LifeHistory_Risor_0812.Rmd` using only the files in `data/social_life_history/`.
2. **Processed genomic analyses:** Run `Relatedness_PopGen_Risor_0812.Rmd` using the files in `data/relatedness_popgen/`.
3. **Sequence processing:** Start from the sample-specific FASTQ files in PRJNA1427305 to recreate the processed VCF files before following entry point 2.

```text
data/
├── social_life_history/         # inputs for Social_LifeHistory_Risor_0812.Rmd
├── relatedness_popgen/          # inputs for Relatedness_PopGen_Risor_0812.Rmd
└── sequence_processing_inputs/  # population maps and historical barcodes
```

### Raw sequence data

The original sequencing output consisted of pooled paired-end ddRAD files. PCR duplicates were removed using ParseDBR_ddRAD, and the reads were subsequently demultiplexed and quality filtered using `process_radtags` in STACKS. The scripts and project-specific inputs used for these upstream steps are retained to document the complete original workflow.

The files deposited under NCBI BioProject [PRJNA1427305](https://www.ncbi.nlm.nih.gov/sra/PRJNA1427305) are organized as separate paired-end FASTQ files for each sample. For example, sample `grunt_137` is represented by `grunt_137.1.fq.gz` and `grunt_137.2.fq.gz`. Therefore, users starting from the public SRA files do not need to repeat barcode-based demultiplexing. The public sequence-data workflow begins with these sample-specific FASTQ files at the de novo STACKS assembly step described below.

### 1. Sequence-data entry point (optional)

#### Inputs and scripts for reproduction from deposited sequences

* `popmap_zero.txt`, `popmap_batch1.txt`, and `popmap_batch1_relatedness.txt`  
  Tab-separated sample-to-population files provided in `data/sequence_processing_inputs/`. `popmap_zero.txt` assigns all samples to one population for the preliminary `populations` run used to generate the locus whitelist. `popmap_batch1.txt` contains the complete sample set for the population-structure analysis. `popmap_batch1_relatedness.txt` contains the reduced sample set for the relatedness analysis.

* `denovo_map.py`  
  Python wrapper used to run `denovo_map.pl`. Users edit only the sample-directory, population-map, output, and thread variables labeled `USER CONFIGURATION`. The published parameters (`m = 5`, `M = 6`, and `n = 7`) are already specified and should not be changed for reproduction.

* `populations_popgen.py`  
  Python wrapper used to run STACKS `populations`. It documents the preliminary whitelist run, the population-structure run (`r = 0.70`), and the relatedness run (`r = 0.95`, `popmap_batch1_relatedness.txt`, and `pan_populations` output).

* `Risor_Stacks_VCFtools_script.sh`  
  Command-line workflow for reproducing the genomic processing from the sample-specific FASTQ files deposited under PRJNA1427305. It documents the final de novo STACKS assembly, locus-whitelist generation, population-structure and relatedness `populations` runs, VCF filtering, Hardy–Weinberg filtering, and linkage-disequilibrium pruning.

#### Additional files retained to document original preprocessing and parameter exploration

The following files document steps performed before the publicly deposited sample-specific FASTQ files and final STACKS settings were available. They are included for documentation but are not required to reproduce the published results from PRJNA1427305 or from the processed VCF files.

* `execute_duplicates.py`  
  Project-specific Python 2.7 wrapper used to remove PCR duplicates with ParseDBR_ddRAD ([Eljensen/ParseDBR_ddRAD](https://github.com/Eljensen/ParseDBR_ddRAD), commit `8fc0edc`). This script is included to document the index sequence, R2 enzyme overhang, adapter-length setting, and output files used for the original pooled sequence data.

* `process_radtags.py`  
  Python wrapper used to demultiplex and quality filter the original pooled batch-1 reads with `process_radtags` in STACKS 2.65.

* `barcode_1st_P1.txt` and `barcode_1st_P8.txt`  
  Barcode files provided in `data/sequence_processing_inputs/` and retained in their original separate form. They document the demultiplexing of the original pooled reads and are not needed when starting from the sample-specific SRA FASTQs.

* `optimize_STACKS_parameters.R`  
  R script using RADstackshelpR and vcfR to document the historical comparison of candidate values `M = 1–8` and `n = 6–8`. The minimum stack depth was fixed at `m = 5` for the published analysis. The final analysis used `m = 5`, `M = 6`, and `n = 7`, as recorded in `denovo_map.py`. The resulting `vis_loci_bigM.png` and `vis_loci_n.png` files are provided in `outputs/parameter_optimization/`. The large preliminary VCF inputs are not included because this exploratory step is not required to rerun the final pipeline.

#### Sequence-data workflow

The complete historical workflow and the public SRA starting point are:

```text
Original pooled batch-1 FASTQ files
        |
        v
FastQC
        |
        v
ParseDBR_ddRAD using execute_duplicates.py
        |
        v
process_radtags.py using the original barcode files
        |
        v
Sample-specific paired-end FASTQ files
        |
        |<----- PUBLIC REPRODUCTION STARTS HERE:
        |       download per-sample FASTQs from PRJNA1427305
        v
De novo STACKS assembly using denovo_map.py
with the optimized parameters m = 5, M = 6, and n = 7
        |
        v
STACKS populations and whitelist generation
        |
        v
VCFtools and BCFtools filtering documented in
Risor_Stacks_VCFtools_script.sh
        |
        v
Final filtered population-structure and relatedness VCFs
        |
        v
Continue with "2. Processed genomic data and analysis" below
```

### 2. Processed genomic data and analysis

Users who do not need to repeat the sequence-processing and variant-filtering pipeline can begin with the following analysis-ready VCF and associated tabular files in `data/relatedness_popgen/`.

* `pop_m5_r07_wl_biallelic_mac5.noHWEexcess_p1e4_LDpruned.vcf`  
  Filtered variant-calling file used for population-genomic analyses. The dataset contains 42 biological samples and two technical replicates. The filename indicates: STACKS minimum depth `m = 5`; loci present in at least 70% of samples (`r07`); whitelist applied (`wl`); biallelic loci retained; minimum allele count 5 (`mac5`); excess-heterozygosity loci at *P* < 1 × 10^-4 removed; and linkage-disequilibrium pruning applied.

* `rel_m5_r095_wl_biallelic_mac16.noHWEexcess_p1e4.LDpruned.vcf`  
  Filtered variant-calling file used for relatedness analyses. The dataset contains 33 individuals. The filename indicates: STACKS minimum depth `m = 5`; loci present in at least 95% of samples (`r095`); whitelist applied (`wl`); biallelic loci retained; minimum allele count 16 (`mac16`); excess-heterozygosity loci at *P* < 1 × 10^-4 removed; and linkage-disequilibrium pruning applied.

These processed VCF files are the starting inputs for the downstream analyses presented in `Relatedness_PopGen_Risor_0812.Rmd`.

#### Associated files and variable definitions

Identifiers, categorical codes, lineage assignments, and relatedness estimators are unitless. Missing values are represented by blank cells, `NA`, or `na`.

* `Risor_metadata_geno_subset_with_all_variables.csv`: individual metadata for 42 biological samples and two `_d` technical replicates.  
  `id`, `id_num`: sample and biological-individual identifiers; `site`, `Reef`: reef; `lat_dd`, `lon_dd`: decimal-degree coordinates; `depth`, `Depth_m`, `Depth`: depth (m); `sponge_id`, `Sponge_id`: host identifier; `S_width_cm`, `Sponge_width_cm`: sponge width (cm); `Log_S_width`: natural log of sponge width; `S_species`, `S_colour`: sponge species and color; `weight`: body mass (g); `sl`: standard length (mm); `tl`, `TL_mm`: total length (mm); `Group.size`, `Group_size`: number of fish in the sponge; `Rank`: size rank (1 = largest); `R1_TL`–`R14_TL`: length (mm) at each rank; `AVG.TL_mm`: group mean length (mm); `mt_lineage`, `nc_lineage`: lineage assignments; `is_male`: 1 = male, 0 = not classified as male; `egg_count`: egg number; `eggs`: 1 = present, 0 = absent; `age`: estimated age (days).

* `Risor_ruber_2023_Pairwise_Distances.csv`: geographic distances for ordered individual pairs.  
  `IN_ind_id`, `NEAR_ind_id`: individual identifiers; `IN_sponge_id`, `NEAR_sponge_id`: host identifiers; `ind_id_Pair`, `sponge_id_Pair`: pair identifiers; `Distance_m`: distance between host sponges (m).

* `Risor_rel_pan_1664SNPs_withDist.csv`: relatedness and distance data for 528 dyads among 33 individuals.  
  `Dyad`: row identifier; `Individual1`, `Individual2`, `PopID1&2`: individuals and population pair; `TrioML`, `Wang`, `LynchLi`, `LynchRd`, `Ritland`, `QuelletGt`, `DyadML`: unitless COANCESTRY relatedness estimates (negative estimates are valid); `ind1_num`, `ind2_num`, `lo`, `hi`, `canon_pair`: identifiers used to standardize dyad order; `Distance_m`: distance (m); remaining `IN_`, `NEAR_`, and pair columns are matching identifiers inherited from the distance table.

Run `Relatedness_PopGen_Risor_0812.Rmd` from `Scripts/` to reproduce the population-genomic and relatedness analyses. It reads all required inputs from `../data/relatedness_popgen/` and writes generated files to `../outputs/`.

### 3. Social and life-history data and analysis (no genomic data required)

Both required files are provided in `data/social_life_history/`.

* `social_data_all_complete.csv`: individual social and life-history data for 223 fish.  
  `ind_id`: fish identifier; `Sponge_id`: host identifier; `Reef`: reef; `Sponge_width_cm`: sponge width (cm); `Depth`: depth (m); `S_species`, `S_colour`: sponge species and color; `Group_size`: number of fish in the sponge; `TL_mm`: total length (mm); `Rank`: size rank (1 = largest); `AVG.TL_mm`: group mean length (mm); `mt_lineage`, `nc_lineage`: lineage assignments; `is_male`: 1 = male, 0 = not classified as male; `egg_count`: egg number; `eggs`: 1 = present, 0 = absent; `age`: estimated age (days).

* `sponge_occupancy_df.csv`: sponge-level occupancy data.  
  `Sponge_ID`: sponge identifier; `Max_width_cm`: maximum width (cm); `Reef`: sampling reef; `S_occ`: 1 = occupied, 0 = unoccupied; `S_unocc`: 1 = unoccupied, 0 = occupied; `R_ruber`: `Present` or `Absent`.

Run `Social_LifeHistory_Risor_0812.Rmd` from `Scripts/` using only these two CSV files from `../data/social_life_history/`.

Rendered HTML reports are included in `outputs/` and end with `sessionInfo()`, which records the R and package versions used.

## METHODOLOGICAL INFORMATION

### Field collection

* **Location and date:** Four patch reefs near Carrie Bow Cay, Belize, sampled during July 2023 at depths shallower than 10 m.

* **Host sponges:** *Ircinia campana*, *I. felix*, *I. strobilina*, and *I. ruetzleri*, including sponges occupied and unoccupied by *Risor ruber*.

* **Field workflow:** Each sponge was photographed with a scale, measured, assigned a depth, and temporarily tagged to prevent resampling. All tags were removed after sampling. Resident fishes were anesthetized by applying a 5:1 ethanol:clove-oil solution around the sponge excurrent openings and were collected using SCUBA.

* **Sample storage:** Fin clips were preserved in 95% ethanol.

* **DNA extraction:** Genomic DNA was extracted using the QIAGEN DNeasy Blood and Tissue Kit following the manufacturer's protocol.

* **Sequencing:** ddRAD libraries were sequenced on an Illumina NovaSeq platform using paired-end 2 × 150-bp sequencing.

### Bioinformatics processing

Raw reads were processed with FastQC, ParseDBR_ddRAD, and STACKS. The final STACKS parameters were m = 5, M = 6, and n = 7. Separate datasets were produced for population-genomic and relatedness analyses, followed by whitelist filtering, retention of biallelic loci, analysis-specific minimum allele-count filtering, removal of loci showing excess heterozygosity, and linkage-disequilibrium pruning. Exact commands, parameters, and filtering order are provided in `Risor_Stacks_VCFtools_script.sh`, and the resulting VCF filenames are explained under “Processed genomic data and analysis.”

For reproduction, users should edit only variables labeled `USER CONFIGURATION` at the beginning of the Python scripts. These variables specify local paths, output locations, population-map files, and thread counts. Variables labeled `PUBLISHED PARAMETERS` contain the analytical settings used for the manuscript and should not be changed.

The relatedness VCF was converted by `Relatedness_PopGen_Risor_0812.Rmd` to the allele-coded COANCESTRY input deposited as `outputs/new_coancestry_input_filtered_dataset_related.txt`. COANCESTRY was then run outside R to calculate the TrioML, Wang, Lynch–Li, Lynch–Ritland, Ritland, Queller–Goodnight, and DyadML estimators. Its dyad-level results were combined with geographic distances in `data/relatedness_popgen/Risor_rel_pan_1664SNPs_withDist.csv`. Correlations among estimators were evaluated in the R Markdown workflow, and the Wang estimator was used for the reported downstream models.

### Software versions

* **FastQC:** 0.12.1
* **ParseDBR_ddRAD:** GitHub commit `8fc0edc` from [Eljensen/ParseDBR_ddRAD](https://github.com/Eljensen/ParseDBR_ddRAD), executed with Python 2.7
* **STACKS:** 2.65
* **VCFtools:** 0.1.14
* **BCFtools:** 1.21
* **COANCESTRY:** 1.0.2.0
* **Python used for the STACKS wrappers:** 3.7
* **R:** 4.5

STACKS, VCFtools, and BCFtools must be installed and available on the command line before running the sequence-processing workflow. Conda environment names shown in `Risor_Stacks_VCFtools_script.sh` correspond to the authors' server and may differ on another system. Users may change environment names and thread counts without changing the analysis.

The R version, operating system, and R-package versions used for the verified rerun are reported at the end of each included HTML report. The required packages are also declared in the setup section of each corresponding R Markdown file.

## ADDITIONAL INFORMATION

The shell file documents the public sequence-processing workflow beginning with the sample-specific FASTQ files in PRJNA1427305. It is not a software-installation script or a fully automated, one-command pipeline. Machine-specific paths and thread counts are identified in the configuration sections of the associated Python scripts.

Generated intermediate files, including STACKS working directories and the large preliminary VCFs used during exploratory parameter evaluation, are not duplicated in the archive because they are not required to rerun the final pipeline with the selected published parameters. The two parameter-selection figures are retained in `outputs/parameter_optimization/` for documentation.

For complete descriptions of field collection, genomic processing, filtering decisions, and statistical analyses, see the Materials and Methods section of the associated manuscript.
