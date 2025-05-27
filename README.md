# SeqSorcerer_May27_2025


Run this code in your terminal and it will prompt you with instructions that will guide you through an RNA seq pipeline for paired end reads. Note: depending on your annotation files/desired counting methods you will have to modify the featurecounts parameters in the run_feature_counts function. 


**AI Summary:** 

This script automates end-to-end processing of paired-end RNA-seq data, handling everything from raw FASTQ files to gene count tables. It performs the following steps:

**Dependency Checks & Installation:**

Verifies presence of required tools: hisat2, trim_galore, samtools, and featureCounts (via Subread).

Installs missing tools using conda.

**Version Reporting:**

Prints versions of key tools (hisat2, trim_galore, featureCounts).

**Input Collection:**

Prompts the user for input FASTQ folder, output directory, reference genome FASTA, and GTF/GFF file.

**Annotation Preparation:**

If a .gff file is provided, it's converted to .gtf using gffread.

Verifies that sequence headers in the FASTA and GTF files are compatible.

**Reference Indexing:**

Builds a HISAT2 index from the reference genome FASTA.

**Trimming & Quality Control:**

Trims adapters and low-quality bases using trim_galore with paired-end settings and FastQC.

Outputs are organized into per-sample subdirectories.

**Alignment:**

Aligns trimmed reads to the reference genome using hisat2 with strand-specific flags.

Converts SAM to BAM and sorts alignments by name using samtools.

**Quantification:**

Runs featureCounts on the sorted BAM files to generate gene-level count tables.

Uses stranded mode (-s 2) and counts CDS regions by transcript_id.
