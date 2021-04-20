# HBVIF
HBV Integration Finder (HBVIF) pipeline.

This bioinformatics-pipeline was designed to directly extract the HBV-DNA integrated sequences from genomic breakpoints of the host genome (Human Genome).
This pipeline takes as input directly the raw reads, performs the trimming (customizable/skippable), and process them generating several tab-separated files than can be easly see by reserchers.
The pipeline is based on the previous installation of the following software:
1) trimmomatic (v.0.39)
2) fastqc
3) bwa (v.)
4) samtools (v.)
5) picard tools (v.)
6) bedtools (v.)
7) cap3
8) cd-hit (v.)
9) blast (v.)

Please be sure that all the listed softwares are installed on your machine before run the pipeline.
