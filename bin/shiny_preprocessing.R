#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("dplyr")
library("readxl")
library("ggplot2")
library("phyloseq")
library("tidyverse")
library(plotly)
library(gapminder)
library(pals)
library(psadd)
library(vegan)


otu_tab = as.data.frame(read_tsv(args[1]))
head(otu_tab)

otu_tax = as.data.frame(read_tsv(args[2],na = c("", "NA")))
head(otu_tax)

asv_tab = as.data.frame(read_tsv(args[3]))
head(asv_tab)

asv_tax = as.data.frame(read_tsv(args[4], na = c("", "NA")))
head(asv_tax)

samples_xl = as.data.frame(read_excel(args[5]))
head(samples_xl)

experiment = as.character(args[6])
print(experiment)


# asv_tab = read_tsv("/data/Microbiome_16S_CRC/work/bertrandp/93AS34deleted_wsearch_r12/8.Denoising-Taxonomy/8B.asv/unoise_zotus_tab.txt")
# otu_tab = read_tsv("/data/Microbiome_16S_CRC/work/bertrandp/93AS34deleted_wsearch_r12/8.Denoising-Taxonomy/8A.otu/uparse_otu_tab.txt")
# 
# asv_tax = read_tsv("/data/Microbiome_16S_CRC/work/bertrandp/93AS34deleted_wsearch_r12/8.Denoising-Taxonomy/Taxonomy/sure_sintax_newTax_asv.txt", na = c("", "NA"))
# otu_tax = read_tsv("/data/Microbiome_16S_CRC/work/bertrandp/93AS34deleted_wsearch_r12/8.Denoising-Taxonomy/Taxonomy/sure_sintax_newTax_otu.txt", na = c("", "NA"))
# 
# samples_xl = read_excel("/data/Microbiome_16S_CRC/work/bertrandp/data/metadata_merged.xlsx")


samples <- samples_xl %>% 
  tibble::column_to_rownames("Sample")
as.factor(samples$Stage)

samples['sample_id'] <- row.names(samples)
samples['id_run'] = paste(samples$sample_id,samples$Run,sep="_run")
samples['patient_stage'] = paste(samples$Participant,samples$Stage,sep="_")
samples['stage_treatment'] = paste(samples$Stage,samples$Treatment,sep="_")
samples['stage_hpv'] = paste(samples$Stage,samples$HPV16_other,sep="_")
samples['participant_stage'] = paste(samples$Participant,samples$Stage,sep="_")

print("inputs preprocessing")

#row names = otu column
otu_tab = otu_tab %>%
  tibble::column_to_rownames("#OTU ID")
asv_tab = asv_tab %>%
  tibble::column_to_rownames("#OTU ID")

otu_tax = otu_tax %>%
  tibble::column_to_rownames("#OTU")
asv_tax = asv_tax %>%
  tibble::column_to_rownames("#OTU")

otu_tax["FULL"] = paste(otu_tax$PHYLUM,otu_tax$CLASS,otu_tax$ORDER,otu_tax$FAMILY,otu_tax$GENUS,otu_tax$SPECIE,sep="_")
asv_tax["FULL"] = paste(asv_tax$PHYLUM,asv_tax$CLASS,asv_tax$ORDER,asv_tax$FAMILY,asv_tax$GENUS,asv_tax$SPECIE,sep="_")

#transform into matrix
otu_tab = as.matrix(otu_tab)
asv_tab = as.matrix(asv_tab)

otu_tax = as.matrix(otu_tax)
asv_tax = as.matrix(asv_tax)


rownames(otu_tax)<-gsub(";size=[0-9]*","",rownames(otu_tax)) #discard the common suffix in all the samples names
rownames(asv_tax)<-gsub(";size=[0-9]*","",rownames(asv_tax)) #discard the common suffix in all the samples names

colnames(otu_tab)<-gsub("_S[0-9]*_L001*","",colnames(otu_tab)) #discard the common suffix in all the samples names
colnames(asv_tab)<-gsub("_S[0-9]*_L001*","",colnames(asv_tab)) #discard the common suffix in all the samples names


#Creation of the phyloseq objects
print("phyloseq objects")

OTU = otu_table(otu_tab, taxa_are_rows = T)
colnames(OTU) = sub("_S[^_]*", "", colnames(OTU)) #remove Usearch tag after the sample name
OTU_TAX = tax_table(otu_tax)

otu_ps = phyloseq(otu_table(OTU),tax_table(OTU_TAX),sample_data(samples))


ASV = otu_table(asv_tab, taxa_are_rows = T)
colnames(ASV) = sub("_S[^_]*", "", colnames(ASV)) #remove Usearch tag after the sample name
ASV_TAX = tax_table(asv_tax)

asv_ps = phyloseq(otu_table(ASV),tax_table(ASV_TAX),sample_data(samples))


#Normalisation
print("normalisation")

otu_top <- names(sort(taxa_sums(otu_ps), decreasing=TRUE)) #how many OTU do we want to keep
otu_ps.top <- transform_sample_counts(otu_ps, function(OTU) OTU/sum(OTU)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
otu_ps_norm <- prune_taxa(otu_top, otu_ps.top)


asv_top <- names(sort(taxa_sums(asv_ps), decreasing=TRUE)) #how many OTU do we want to keep
asv_ps.top <- transform_sample_counts(asv_ps, function(OTU) OTU/sum(OTU)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
asv_ps_norm <- prune_taxa(asv_top, asv_ps.top)


# Creation de otu_df 
# otu_df est le dataframe contenant toutes les informations de l'objet phyloseq otu_ps 
print("otu_df creation")

otu_df <- psmelt(otu_ps)

otu_df = otu_df %>% group_by(Stage) %>% mutate(AbundNormStage = Abundance/sum(Abundance)) #normalize abund by stage

otu_df = otu_df %>% group_by(sample_id) %>% mutate(AbundNormSample = Abundance/sum(Abundance))


alpha_shannon = data.frame(diversity(t(otu_tab), index="shannon")) 
alpha_invsimpson = data.frame(diversity(t(otu_tab), index="invsimpson"))

rownames(alpha_shannon)<-gsub("_S[0-9]*","",rownames(alpha_shannon))
alpha_shannon["sample_id"] = row.names(alpha_shannon)
colnames(alpha_shannon) = c("alpha_shannon", "sample_id")

rownames(alpha_invsimpson)<-gsub("_S[0-9]*","",rownames(alpha_invsimpson))
alpha_invsimpson["sample_id"] = row.names(alpha_invsimpson)
colnames(alpha_invsimpson) = c("alpha_invsimpson", "sample_id")

if(length(grep("alpha_shannon",names(otu_df),value=TRUE)) == 0) { #si la l'alpha n'est pas deja calculé
  otu_df = full_join(alpha_shannon, otu_df, by="sample_id")
}

if(length(grep("alpha_invsimpson",names(otu_df),value=TRUE)) == 0) {
  otu_df = full_join(alpha_invsimpson, otu_df, by="sample_id")
}

# Creation de asv_df 
# asv_df est le dataframe contenant toutes les informations de l'objet phyloseq asv_ps 
print("asv_df creation")


asv_df <- psmelt(asv_ps)

asv_df = asv_df %>% group_by(Stage) %>% mutate(AbundNormStage = Abundance/sum(Abundance)) #normalize abund by stage

asv_df = asv_df %>% group_by(sample_id) %>% mutate(AbundNormSample = Abundance/sum(Abundance))


alpha_shannon = data.frame(diversity(t(asv_tab), index="shannon")) 
alpha_invsimpson = data.frame(diversity(t(asv_tab), index="invsimpson"))

rownames(alpha_shannon)<-gsub("_S[0-9]*","",rownames(alpha_shannon))
alpha_shannon["sample_id"] = row.names(alpha_shannon)
colnames(alpha_shannon) = c("alpha_shannon", "sample_id")

rownames(alpha_invsimpson)<-gsub("_S[0-9]*","",rownames(alpha_invsimpson))
alpha_invsimpson["sample_id"] = row.names(alpha_invsimpson)
colnames(alpha_invsimpson) = c("alpha_invsimpson", "sample_id")

if(length(grep("alpha_shannon",names(asv_df),value=TRUE)) == 0) { #si la l'alpha n'est pas deja calculé
  asv_df = full_join(alpha_shannon, asv_df, by="sample_id")
}

if(length(grep("alpha_invsimpson",names(asv_df),value=TRUE)) == 0) {
  asv_df = full_join(alpha_invsimpson, asv_df, by="sample_id")
}

print("saveRDS")
saveRDS(asv_ps, file = sprintf("./%s_asv_ps.Rds",experiment))
saveRDS(asv_df, file = sprintf("./%s_asv_df.Rds",experiment))
saveRDS(otu_df, file = sprintf("./%s_otu_df.Rds",experiment))
saveRDS(otu_ps, file = sprintf("./%s_otu_ps.Rds",experiment))
saveRDS(samples, file = sprintf("./%s_samples.Rds",experiment))
