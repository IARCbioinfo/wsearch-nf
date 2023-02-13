#!/usr/bin/env Rscript

# install.packages("dplyr")     # To manipulate dataframes
# install.packages("readxl")    # To read Excel files into R
# install.packages("ggplot2")   # for high quality graphics
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# BiocManager::install("phyloseq")
# BiocManager::install("microbiomeMarker")
# 
# install.packages("dplyr")     # To manipulate dataframes
# install.packages("readxl")    # To read Excel files into R
# install.packages("ggplot2")   # for high quality graphics
# 
# 
# library("phyloseq")
# library("ggplot2")
# library(tidyverse)
# library(readxl)
args = commandArgs(trailingOnly=TRUE)

library("dplyr")
library("readxl")
library("ggplot2")
library("phyloseq")
library("tidyverse")

if (!(file.exists(args[4]))){
  stop(" /!\ Five arguments must be supplied : \n 1- Otu table \n 2- Otu taxonomy table \n 3- Asv table \n 4- Asv taxonomy table \n")
}

otu_tab = read_tsv(args[1])
print(otu_tab)
otu_tax = read_tsv(args[2])
print(otu_tax)
asv_tab = read_tsv(args[3])
print(asv_tab)
asv_tax = read_tsv(args[4])
print(asv_tax)


#open tables from USEARCH

#output : ./Analyse/
dir.create("Analysis")

#############################################################################################################
# Phyloseq
#############################################################################################################

# asv_tab = read_tsv("/data/Microbiome/work/bertrandp/nextflow/Microbiome-nf/All-samples_WSEARCH/9.Denoising-Taxonomy/9B.asv/unoise_zotus_tab.txt")
# otu_tab = read_tsv("/data/Microbiome/work/bertrandp/nextflow/Microbiome-nf/All-samples_WSEARCH/9.Denoising-Taxonomy/9A.otu/uparse_otu_tab.txt")
# 
# asv_tax = read_tsv("/data/Microbiome/work/bertrandp/nextflow/Microbiome-nf/All-samples_WSEARCH/9.Denoising-Taxonomy/Taxonomy/sure_sintax_newTax_asv.txt")
# otu_tax = read_tsv("/data/Microbiome/work/bertrandp/nextflow/Microbiome-nf/All-samples_WSEARCH/9.Denoising-Taxonomy/Taxonomy/sure_sintax_newTax_otu.txt")

#row names = otu column
otu_tab = otu_tab %>%
  tibble::column_to_rownames("#OTU ID")
asv_tab = asv_tab %>%
  tibble::column_to_rownames("#OTU ID")

otu_tax = otu_tax %>%
  tibble::column_to_rownames("#OTU")
asv_tax = asv_tax %>%
  tibble::column_to_rownames("#OTU")


#transform into matrix
otu_tab = as.matrix(otu_tab)
asv_tab = as.matrix(asv_tab)

otu_tax = as.matrix(otu_tax)
asv_tax = as.matrix(asv_tax)


rownames(otu_tax)<-gsub(";size=[0-9]*","",rownames(otu_tax)) #discard the common suffix in all the samples names
rownames(asv_tax)<-gsub(";size=[0-9]*","",rownames(asv_tax)) #discard the common suffix in all the samples names

colnames(otu_tab)<-gsub("_S[0-9]*_L001__contig_All-samples_run1-2","",colnames(otu_tab)) #discard the common suffix in all the samples names
colnames(asv_tab)<-gsub("_S[0-9]*_L001__contig_All-samples_run1-2","",colnames(asv_tab)) #discard the common suffix in all the samples names



#Creation of the phyloseq objects
OTU = otu_table(otu_tab, taxa_are_rows = T)
colnames(OTU) = sub("_S[^_]*", "", colnames(OTU)) #remove Usearch tag after the sample name
OTU_TAX = tax_table(otu_tax)

otu_ps = phyloseq(OTU,OTU_TAX)


ASV = otu_table(asv_tab, taxa_are_rows = T)
colnames(ASV) = sub("_S[^_]*", "", colnames(ASV)) #remove Usearch tag after the sample name
ASV_TAX = tax_table(asv_tax)

asv_ps = phyloseq(ASV,ASV_TAX)

######################################################################

## ANALYSE OTU ##

#Normalisation
otu_total = median(sample_sums(otu_ps))
standf = function(OTU) OTU/sum(OTU)
otu_ps_norm = transform_sample_counts(otu_ps, standf)


#Alpha Diversity

otu_alpha_div=estimate_richness(otu_ps,measures=c("InvSimpson", "Shannon")) #Simpson and Shannon not reliant on singletons
saveRDS(otu_alpha_div, "Analysis/usearch_otu_alphaDiv_AllSamples2_run1-2.rds")

png(file="Analysis/usearch_otu_alpha.png")
plot_richness(otu_ps_norm, measures=c("InvSimpson", "Shannon"), x = "samples", color="samples", title="Alpha diversity otu")
dev.off()

#Ordination
otu_ps.ord_bray <- ordinate(otu_ps, "NMDS", "bray")

otu_beta_div = distance(otu_ps, method = "bray")
saveRDS(otu_beta_div, "Analysis/usearch_otu_betaDiv_AllSamples2_run1-2.rds")

# plot the distance between samples
png(file="Analysis/usearch_otu_plotOrd_samples.png")
plot_ordination(otu_ps, otu_ps.ord_bray,
                              title="Beta diversity otu") + geom_point(size=3) # Maybe interesting ?
dev.off()

# #Normalize abundances
otu_top <- names(sort(taxa_sums(otu_ps), decreasing=TRUE)) #how many OTU do we want to keep
otu_ps.top <- transform_sample_counts(otu_ps, function(OTU) OTU/sum(OTU)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
otu_ps_norm <- prune_taxa(otu_top, otu_ps.top)

#plots representing the Phylum composition
png(file="Analysis/usearch_otu_compo.png", width=2000, height=1000)
otu_phylumGlommed = tax_glom(otu_ps_norm, "PHYLUM")
plot_bar(otu_phylumGlommed, fill="PHYLUM", title = "Phylum Composition of each sample")
dev.off()


###################################################################################


## ANALYSE ASV ##

#Normalisation
asv_total = median(sample_sums(asv_ps))
asv_standf = function(ASV) ASV/sum(ASV)
asv_ps_norm = transform_sample_counts(asv_ps, asv_standf)


#Alpha Diversity

asv_alpha_div=estimate_richness(asv_ps,measures=c("InvSimpson", "Shannon")) #Simpson and Shannon not reliant on singletons
saveRDS(asv_alpha_div, "Analysis/usearch_asv_alphaDiv_AllSamples2_run1-2.rds")

png(file="Analysis/usearch_asv_alpha.png")
plot_richness(asv_ps_norm, measures=c("InvSimpson", "Shannon"), x = "samples", color="samples")
dev.off()

#Ordination
asv_ps.ord_bray <- ordinate(asv_ps, "NMDS", "bray")

asv_beta_div = distance(asv_ps, method = "bray")
saveRDS(asv_beta_div, "Analysis/usearch_asv_betaDiv_AllSamples2_run1-2.rds")

# plot the distance between samples
png(file="Analysis/usearch_asv_plotOrd_samples.png")
plot_ordination(asv_ps, asv_ps.ord_bray, 
                title="Beta diversity asv") + geom_point(size=3) 
dev.off()

# #Normalize abundances
asv_top <- names(sort(taxa_sums(asv_ps), decreasing=TRUE)) #how many OTU do we want to keep
asv_ps.top <- transform_sample_counts(asv_ps, function(ASV) ASV/sum(ASV)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
asv_ps_norm <- prune_taxa(asv_top, asv_ps.top)


#plots representing the Phylum composition
png(file="Analysis/usearch_asv_compo.png", width=2000, height=1000)
asv_phylumGlommed = tax_glom(asv_ps_norm, "PHYLUM")
plot_bar(asv_phylumGlommed, fill="PHYLUM", title = "Phylum composition of each samples")
dev.off()


