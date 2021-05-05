# HBVIF
HBV Integration Finder (HBVIF) pipeline.

This bioinformatics pipeline was designed to directly extract the HBV-DNA integrated sequences from genomic breakpoints of the host genome (Human Genome).
This pipeline takes as input directly the raw reads, performs the trimming (customizable/skippable), and processes them generating several tab-separated files that can be easily visualized by researchers.
The pipeline is based on the previous installation of the following software:
1) fastqc (v0.11.8)
2) trimmomatic (v.0.39)
3) bwa (>= 0.7.17-r1188)
4) samtools (v. 1.9)
5) picard tools (MarkDuplicates v. 2.20.1-SNAPSHOT)
6) bedtools (v.2.28.0)
7) cap3 (VersionDate: 02/10/15)
8) cd-hit (v.4.8.1)
9) blast (v.2.9.0-2)

Please be sure that all the listed software are installed on your machine before running the pipeline.

This tool was developed by Domenico Giosa (dgiosa@unime.it) at the University of Messina, in collaboration with Sequentia Biotech SL (https://www.sequentiabiotech.com/).
