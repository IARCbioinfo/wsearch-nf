#!/usr/bin/env bash

all_SF=$1
all_labelled=$2
chimera=$3
silva=$4
threads=$5

mkdir 8A.otu
mkdir 8B.asv
mkdir Taxonomy


#dereplication
vsearch -fastx_uniques $all_SF -fastaout all_SF_DR.fasta -sizeout

#OTU clustering
cd 8A.otu
echo "\n ||| OTU clustering ||| \n"

#clusterisation (!= UPARSE)
vsearch -cluster_size ../all_SF_DR.fasta -id 0.97 -centroids clustered.fasta -relabel OTU
#remove chimeras
vsearch --uchime3_denovo clustered.fasta --sizein --sizeout --nonchimeras c_unchimered.fasta
vsearch --uchime_ref c_unchimered.fasta --sizein --sizeout --nonchimeras uparse_otus.fasta --db ../$chimera -threads $threads
#OTU table creation
vsearch -usearch_global ../$all_labelled -db uparse_otus.fasta -strand both -id 0.97 -otutabout uparse_otu_tab.txt -biomout uparse_otu_biom.biom -threads $threads


#ASV denoising
cd ../8B.asv
echo "\n ||| ASV denoising ||| \n"

#denoising (UNOISE3)
vsearch -cluster_unoise ../all_SF_DR.fasta --sizein -sizeout --centroids denoised.fasta -relabel ASV  -threads $threads
#remove chimeras
vsearch --uchime3_denovo denoised.fasta --sizein --sizeout --nonchimeras d_unchimered.fasta  -threads $threads
vsearch --uchime_ref d_unchimered.fasta --sizein --sizeout --nonchimeras unoise_zotus.fasta --db ../$chimera  -threads $threads
#ASV table creation
vsearch -usearch_global ../$all_labelled -db unoise_zotus.fasta -strand both -id 0.97 -otutabout unoise_zotus_tab.txt -biomout unoise_zotus_biom.biom  -threads $threads


#Taxonomy
cd ../Taxonomy
echo "\n ||| Taxonomy ||| \n"

#assign tax to OTU sequences
vsearch -sintax ../8A.otu/uparse_otus.fasta -db ../$silva -tabbedout sintax_tax_otu.sintax -strand both -sintax_cutoff 0.8
convert_tax.py sintax_tax_otu.sintax sintax_newTax_otu.txt #transform tax_table in order to input in phyloseq after

#assign tax to ASV sequences
vsearch -sintax ../8B.asv/unoise_zotus.fasta -db ../$silva -tabbedout sintax_tax_asv.sintax -strand both -sintax_cutoff 0.8
convert_tax.py sintax_tax_asv.sintax sintax_newTax_asv.txt





