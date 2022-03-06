#!/bin/bash

makeblastdb -in data/16S_rRNA_library_filtered_edit.fasta -dbtype nucl
blastn -query data/clean/asv_seq.fasta -db /mnt/nuuk/2021/HRGM/16S_rRNA/16S_rRNA_library_filtered_edit.fasta -perc_identity 97 -qcov_hsp_perc 95 -outfmt 6 > analysis/v1/files/asvs_to_HRGM.m8
