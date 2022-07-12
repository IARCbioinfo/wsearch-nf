---
title: "Wsearch"
author: "Pierre BERTRAND"
date: '2022-07-07'
output: 
  html_document :
    fig_width: 15
    fig_height: 8
    fig_caption: true
---

# Wsearch

1. Processing microbiome analysis with Usearch and Vsearch
2. Plots analysis with Phyloseq
3. Other analysis to suit the study & example

### Wsearch pipeline

**This pipeline will use Usearch and Vsearch in order to analyse the microbiome sequenced from Gene Marker 16S RNA. **
It is using **Usearch 32bit** (free-version) for firsts step because of the better results. Then, steps which need high computing performances are performed with **Vsearch**, which is an open-source software created in order to replace the paid-version (64bit) of Usearch.
The default parameters (see 1.2 Inputs and 1.3 Parameters) were optimized for analyse paired reads of 300bp from Colorectal Cancer patients, but can be applied on many studies.

![Wsearch_pipeline_steps_options](Microbiome_CRC_results/Wsearch_pipeline_steps_options.png)


The goal of this pipeline is to compute the alpha and beta-diversity of the dataset and assign taxonomies in order to indentify the composition of each sample. Then, it will plot thoses measurement and composition to have an overview of the samples analyse, and guide decision for futurs analysis. 

### Prerequisites

Usearch v11 : https://www.drive5.com/usearch/download.html  
Vsearch v2.21.1_linux_x86_64 : https://github.com/torognes/vsearch/releases

```bash
install.packages("dplyr")
install.packages("readxl")
install.packages("ggplot2")
install.packages("tidyverse")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("phyloseq")
```

## 1. Processing microbiome analysis with Usearch and Vsearch

### 1.1 Command line 

```bash
nextflow Wsearch.nf --input_folder "fastq_folder" --output_folder "output_folder" --experiment "Microbiome_CRC"
```

### 1.2 Inputs

 | Mandatory parameters   | Example value       |     Description                                                          |
  |:-----------------------|:--------------------|:-------------------------------------------------------------------------|
  | --input_folder         | "/data/Microbiome/work/bertrandp/deziped-fastq_run1-2"| Path to the folder containing .fastq files (forward and reverse) of each samples  |
  | --output_folder                 | "Microbiome_CRC_WSEARCH"              | Folder where the pipeline will perform and output                                         |
  | --experiment                 | "Microbiome_CRC" | Name of the experiment   |

 
* *input_folder* is the path to the folder containing all .fastq files, this folder must include : 
  * 2 files per samples : one forward and one reverse
  * Names of the reads must be default name from Illumina Miseq sequencing : SampleName_R1_001.fastq (forward) and SampleName_R2_001.fastq (reverse)
    For instance : T1-1_S67_L001_R1_001.fastq (forward) and T1-1_S67_L001_R2_001.fastq (reverse)


* *output_folder* is the path to the folder where the pipeline will create the finals files, at the end of the run, this folder will contain :
  * 1.reads_summary : For each file, print some statistics
  * 2.merge_reads : For each samples, merge the forward and the reverse file in order to create one contig file for each sample. Add also a .log file for each contig with some statistics on the new contig.
  * 3.filtered_contigs : For each sample, create one file containing reads passed filter's parameters.
  * 4.labeled_data : For each sample, add the name of the sample in the header of the reads.
  * 5.dereplicated_data : For each sample, keep only one copie of each replicated reads. Add a label 'size' in the header of each uniq read in order to take into account the abundance later.
  * 6.low-abund_data : For each sample, keep only reads with a size=1 (singletons). It means reads not duplicated, usually, those reads are erroneous and will be remove in the next step.
  * 7.Singleton-free : Remove singleton (reads in the 6.low-abund) from the dereplicated reads (reads in the 5.dereplicated)
  * 8.Denoising-Taxonomy : merge into one file all reads from all samples (all_labelled.fasta) and from data Singleton-free (all_SF.fasta). 
    * 8A.otu : Fasta sequence of OTU (uparse_otus.fasta), OTU table with the abundance of OTU in each samples (uparse_otu_tab.txt, uparse_otu_biom.biom)
    * 8B.asv : Fasta sequence of ASV (unoise_zotus.fasta), ASV table with the abundance of ASV in each samples (unoise_zotus_tab.txt, unoise_zotus_biom.biom)
    * Taxonomy : Taxonomy of ASV and OTU in the sintrax format (.sintax), and for each denoising method (OTU & ASV) a tabbed taxonomy table : sure[...]_newTax_[...].txt contain only taxonomies assigned with a score higher than the cutoff (default = 0.8), in oposition to notSure[...]_newTax_[...].txt which contain all taxonomies found, even under the cutoff score. Only the "sure[...]_newTax_[...].txt" will be use in the next step.
  * Analysis : For each denoising method compute the alpha diversity (Shannon and InverseSimpson index) for each samples and the beta diversity between samples (Bray-curtis distance). .rds file are saved R object, they can be load in R in order to perform statistical test (see Example). Several .png are also created displaying the alpha and beta diversity, but also the Phylum composition of each sample.


* The *experiment* name will be write in the file name of outputted files.

### 1.3 Parameters

  | Optional parameters    | Example value   |     Description                                                          |
  |:-----------------------|:----------------|:-------------------------------------------------------------------------|
  | --maxdiff            | "40"         | Maximum number of missmatch between the forward and the reverse read, in the overlap region. Can be decreased if reads have great quality [default=40] |
  | --overlap              | "50"       | Minimum length (base pairs) of the overlap between the forward and the reverse reads. Couple of reads with a smaller overlap are removed [default=50]                                 |
  | --identity              | "70"      | Percentage of identity minimum between the overlap region of the forward and reverse reads. Can be increased with reads have good quality. [default=70]                                 |
  | --min_contig_length               | "50"     | Minimum length of the contig, after merge of the forward and reverse read [default=50] |
  | --max_contig_length   | "550"            | Maximum length of the contig, after merge of the forward and reverse read [default=550] |
  | --max_ee              | "1"             | Maximum of expected errors allowed in each contigs. Expected errors are linked to the quality score from the fastq, and permit to estimate the number of errors in the contig. Can be increased if reads have a bad quality. https://drive5.com/usearch/manual8.0/expected_errors.html      |
  | --threads                 | "20"            | Number of threads pipeline can use. Used only for steps need lot of computing power (Clustering, Denoising, Taxonomy, Analysis)                   |
  
  
Thoses default parameters come from a benchmarking of severals option within Usearch. 
**The main logic of this pipeline is to use flexible parameters in the merging step and to be strict during the filtering of the contigs.** This way, most of reads will be able to merge, and then improve the quality of the overlap region. Because, during the merging step, the algorithm will choose the nucleotide with the best quality for each position along the overlap region.
Then, after the merging, we can use strict parameters in order to keep only contigs with a good quality, and remove incorrect contigs.

## 2. Analysis output

Those plots and tables are automatically generated by the pipeline, you can find them in the "Analysis" folder.

![usearch_asv Alpha Diversity](Microbiome_CRC_results/usearch_asv_alpha.png)





![usearch_asv Phylum composition](Microbiome_CRC_results/usearch_asv_compo.png)




![usearch_asv Beta diversity](Microbiome_CRC_results/usearch_asv_plotOrd_samples.png)




![usearch_asv Beta diversity & taxonomy](Microbiome_CRC_results/usearch_asv_plotOrd_samples_tax.png)




  
  Alpha-diversity table :
  
  | Samples | Shannon | Inverse Simpson |
  |:--------|:--------------|:------------|
  | N10 | 7.500938  | 602.80860  |
  | N11 | 7.887776  | 1205.19233  |
  | N12   |       8.020262| 1129.59413|
|N13     |     7.408143 | 724.69644|
|N14      |    7.629771 | 641.82615|

  Beta-diversity table :

||N10    |   N11    |   N12     |  N13    |
|:--------|:--------------|:------------|:------------|:------------|
N11    |      0.8843791   ||||                                                         
N12     |     0.9203999 | 0.6940893 |||                                                  
N13     |     0.9117488 | 0.8499791 | 0.7356067 ||                                       
N14     |     0.9013869 |0.7143146| 0.5632045| 0.7815333  |


## 3.Other analysis to suit the study

In this part, we will introduce a script in order to improve plots and the analysis of our data.
The main idea is to add a metadata table to plot based on one or several variables.

### Metadata spreadsheet

It will contain all data about our sample. In the case of the colorectal cancer analysis, we had the antibiotic intake, the group the 'Stage' belong to (NC, Before therapy, After therapy), etc.
This spreadsheet must be in excel format to be input in this script.

Here, an overview of how the data must be organized in the metadata file :

 | Sample | Stage | Antibiotic |
  |:--------|:--------------|:------------|
  | N10 | NC  | No  |
  | N11 | NC  | No  |
  | T1-1   |       Before| No|
|T5-1     |     Before | Yes |
|T5-2     |    After | No|

**Note that the 'Sample' column must not be changed.**

The Metadata spreadsheet used in the example below is available in the folder.

### Example | Analysis of the Colorectal Cancer microbiome before and after therapy 

#### Prerequires

```{r, results='hide', message=FALSE, warning=FALSE}
library("dplyr")
library("readxl")
library("ggplot2")
library("phyloseq")
library("tidyverse")
```

#### Inputs

Firstly, we need to recover OTU/ASV table and Taxonomies table outputed by Wsearch
```{r, results='hide', warning=FALSE, message=FALSE}
asv_tab = read_tsv("Microbiome_CRC_results/unoise_zotus_tab.txt")
otu_tab = read_tsv("Microbiome_CRC_results/uparse_otu_tab.txt")
asv_tax = read_tsv("Microbiome_CRC_results/sure_sintax_newTax_asv.txt")
otu_tax = read_tsv("Microbiome_CRC_results/sure_sintax_newTax_otu.txt")

alpha_asv = readRDS("Microbiome_CRC_results/usearch_asv_alphaDiv_AllSamples2_run1-2.rds")
alpha_otu = readRDS("Microbiome_CRC_results/usearch_otu_alphaDiv_AllSamples2_run1-2.rds")
```
Secondly, to have great graphs and perform analysis at the best, we need a metadata table. It should show data about our samples. For instance, in which group the sample became, if the patient took antibiotics, or her diet.

```{r, results='hide', warning=FALSE, message=FALSE}
samples_xl <- read_excel("Metadata_table.xlsx","Metadata")
```

After loading table, we need to set rownames of each dataset
```{r, results='hide', warning=FALSE}
samples <- samples_xl %>% 
  tibble::column_to_rownames("Sample")
as.factor(samples$Stage)

#row names = otu column
otu_tab = otu_tab %>%
  tibble::column_to_rownames("#OTU ID")
asv_tab = asv_tab %>%
  tibble::column_to_rownames("#OTU ID")

otu_tax = otu_tax %>%
  tibble::column_to_rownames("#OTU")
asv_tax = asv_tax %>%
  tibble::column_to_rownames("#OTU")
```

and transform dataset into matrix
```{r, results='hide', warning=FALSE}
#transform into matrix
otu_tab = as.matrix(otu_tab)
asv_tab = as.matrix(asv_tab)

otu_tax = as.matrix(otu_tax)
asv_tax = as.matrix(asv_tax)
```

Then, we want to keep only the ID (Name before the "_S" followed by number) of the sample as samples name.  
For example : T10_S79_L001 into T10
```{r, results='hide', warning=FALSE}
rownames(otu_tax)<-gsub(";size=[0-9]*","",rownames(otu_tax)) #discard the common suffix in all the samples names
rownames(asv_tax)<-gsub(";size=[0-9]*","",rownames(asv_tax)) #discard the common suffix in all the samples names

colnames(otu_tab)<-gsub("_S[0-9]*_L001__contig_All-samples_run1-2","",colnames(otu_tab)) #discard the common suffix in all the samples names
colnames(asv_tab)<-gsub("_S[0-9]*_L001__contig_All-samples_run1-2","",colnames(asv_tab)) #discard the common suffix in all the samples names
```

Finally, we can create phyloseq objects
```{r, results='hide', warning=FALSE}
#Creation of the phyloseq objects
OTU = otu_table(otu_tab, taxa_are_rows = T)
colnames(OTU) = sub("_S[^_]*", "", colnames(OTU)) #remove Usearch tag after the sample name
OTU_TAX = tax_table(otu_tax)

otu_ps = phyloseq(OTU,OTU_TAX,sample_data(samples))

ASV = otu_table(asv_tab, taxa_are_rows = T)
colnames(ASV) = sub("_S[^_]*", "", colnames(OTU)) #remove Usearch tag after the sample name
ASV_TAX = tax_table(asv_tax)

asv_ps = phyloseq(ASV,ASV_TAX,sample_data(samples))
```

#### OTU ANALYSIS

We will start with the analysis of the OTU table and taxonomy.

Firstly, we want to normalize abundances into percentage
```{r, warning=FALSE}
#Normalisation
otu_total = median(sample_sums(otu_ps))
standf = function(OTU) OTU/sum(OTU)
otu_ps_norm = transform_sample_counts(otu_ps, standf)
```

Calculation of the alpha diversity (Inverse-Simpson and Shannon index), and plot according to the group
```{r, results='hide', warning=FALSE}
#Alpha Diversity, group by "Stage" from the metadata file
plot_richness(otu_ps_norm, measures=c("InvSimpson", "Shannon"), x = "Stage", color="Stage") #Simpson and Shannon index not reliant to singletons
```

Calculation of the beta diversity (Bray-Curtis distance), and display
```{r, results='hide', warning=FALSE}
#Ordination
otu_ps.ord_bray <- ordinate(otu_ps, "NMDS", "bray")

# plot the distance between samples, color by "Stage" and display patient took antibiotics
plot_ordination(otu_ps, otu_ps.ord_bray, type="samples", color="Stage", 
                                  shape="Antibiotic", title="Samples") + geom_point(size=3) 

#plot distances between samples and taxonomies as well
sample_data(otu_ps)['sample_id'] <- row.names(sample_data(otu_ps)) 

plot_ordination(otu_ps, otu_ps.ord_bray, type="split", color="CLASS", shape="Stage", title="biplot", label = "sample_id") + geom_point(size=3) #interesting
```

To continue, we need to transform data into factor
```{r, warning=FALSE}
#transform data to factor
otu_df <- as.data.frame(lapply(sample_data(otu_ps),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(otu_df) <- sample_names(otu_ps)
sample_data(otu_ps) <- sample_data(otu_df)
```

To do specific plot, we need to merge samples according to the "Stage", and re normalize dataset
```{r, warning=FALSE}
#merge samples per 'stages' (NC, before, after)
otu_ps_fraction <- merge_samples(otu_ps, "Stage", fun=mean)

# #Normalize abundances
otu_top <- names(sort(taxa_sums(otu_ps), decreasing=TRUE)) #how many OTU do we want to keep
otu_ps.top <- transform_sample_counts(otu_ps, function(OTU) OTU/sum(OTU)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
otu_ps_norm <- prune_taxa(otu_top, otu_ps.top)

otu_top_f <- names(sort(taxa_sums(otu_ps_fraction), decreasing=TRUE)) #how many OTU do we want to keep 
otu_ps.top_f <- transform_sample_counts(otu_ps_fraction, function(OTU) OTU/sum(OTU)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
otu_ps_fract_norm <- prune_taxa(otu_top_f, otu_ps.top_f)
```

Now, we want to plot the composition of each samples, and separate them according to the "Stage"
```{r, warning=FALSE}
#plots representing the Phylum composition
otu_phylumGlommed = tax_glom(otu_ps_norm, "PHYLUM")
plot_bar(otu_phylumGlommed, fill="PHYLUM", title = "Plot_ab" ,facet_grid=~Stage)
```

Here, we will plot the mean composition in each "Stage"
```{r, warning=FALSE}
plot_bar(otu_ps_fract_norm, title="usearch_otu_stack_compo" ,fill = "PHYLUM") + 
  geom_bar(aes(color=PHYLUM, fill=PHYLUM), stat="identity", position="stack")
```

#### ASV ANALYSIS

We will perform the same script from the ASV table
```{r, results='hide', warning=FALSE, fig.show='hide'}
#Normalisation
asv_total = median(sample_sums(asv_ps))
asv_standf = function(ASV) ASV/sum(ASV)
asv_ps_norm = transform_sample_counts(asv_ps, asv_standf)

#Alpha Diversity
plot_richness(asv_ps_norm, measures=c("InvSimpson", "Shannon"), x = "Stage", color="Stage") #Simpson and Shannon index not reliant to singletons

#Ordination
asv_ps.ord_bray <- ordinate(asv_ps, "NMDS", "bray")

asv_beta_div = distance(asv_ps, method = "bray")

# plot the distance between samples
plot_ordination(asv_ps, asv_ps.ord_bray, type="samples", color="Stage", 
                shape="Antibiotic", title="Samples") + geom_point(size=3) 

#plot distances between samples and taxonomy of ASV
sample_data(asv_ps)['sample_id'] <- row.names(sample_data(asv_ps)) 

plot_ordination(asv_ps, asv_ps.ord_bray, type="split", color="CLASS", shape="Stage", title="biplot", label = "sample_id") + geom_point(size=3)

#transform data to factor
asv_df <- as.data.frame(lapply(sample_data(asv_ps),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(asv_df) <- sample_names(asv_ps)
sample_data(asv_ps) <- sample_data(asv_df)

#merge samples per 'stages' (NC, before, after)
asv_ps_fraction <- merge_samples(asv_ps, "Stage", fun=mean)

# #Normalize abundances
asv_top <- names(sort(taxa_sums(asv_ps), decreasing=TRUE)) #how many OTU do we want to keep
asv_ps.top <- transform_sample_counts(asv_ps, function(ASV) ASV/sum(ASV)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
asv_ps_norm <- prune_taxa(asv_top, asv_ps.top)

asv_top_f <- names(sort(taxa_sums(asv_ps_fraction), decreasing=TRUE)) #how many OTU do we want to keep 
asv_ps.top_f <- transform_sample_counts(asv_ps_fraction, function(ASV) ASV/sum(ASV)) # transform counts into pourcentages -> for each samples : abundance/sum(abundances)
asv_ps_fract_norm <- prune_taxa(asv_top_f, asv_ps.top_f)


#plots representing the Phylum composition
asv_phylumGlommed = tax_glom(asv_ps_norm, "PHYLUM")
plot_bar(asv_phylumGlommed, fill="PHYLUM", title = "Plot_ab" ,facet_grid=~Stage)

plot_bar(asv_ps_fract_norm, title="usearch_asv_stack_compo" ,fill = "PHYLUM") + 
  geom_bar(aes(color=PHYLUM, fill=PHYLUM), stat="identity", position="stack")
```

#### Statistical tests

The goal of this part is to test a significant diffenrence between groups. Which could represent changes in the Microbiome due to the Colorectal cancer and therapy. 

```{r setup, results='hide', message=FALSE, warning=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(rstatix)
library(ggpubr)
```

##### Loading


```{r}
#replace rownames by a columns
alpha_asv <- tibble::rownames_to_column(alpha_asv, "sample")
alpha_otu <- tibble::rownames_to_column(alpha_otu, "sample")

#add the group for each samples: samples names starting with 'N' belong to 'NC' group, samples names ending by '-1' belong to the 'Before' group, and those ending with '-2' to the 'After' group
asv<-alpha_asv %>%
  mutate( group=case_when( grepl("N",sample) ~ "NC", grepl("-1", sample) ~ "Before", grepl("-2",sample) ~ "After"  )) %>%
  mutate(group=factor(group,levels = c("NC","Before","After")))


otu<-alpha_otu %>%
  mutate( group=case_when( grepl("N",sample) ~ "NC", grepl("-1",sample) ~ "Before", grepl("-2",sample) ~ "After"  )) %>%
  mutate(group=factor(group,levels = c("NC","Before","After")))
```

##### Normality test

Before to use any staistical test, we have to check the normality of the dataset
```{r, fig.show='hide'}
asv %>% group_by(group) %>% shapiro_test(Shannon)
ggqqplot(asv, x = "Shannon", facet.by = "group")

otu %>% group_by(group) %>% shapiro_test(Shannon)
ggqqplot(otu, x = "Shannon", facet.by = "group")
```

If the p.value is lower than the alpha risk (0.05), then the distribution isn't normal

##### Compare multiple Variance

We also compare variances
```{r}
asv %>% levene_test(Shannon ~ group)
otu %>% levene_test(Shannon ~ group)
```

Population variances are not equal if p < 0.05.

##### Identify aberrant values

We looking for the aberrant values as well
```{r}
asv %>% identify_outliers(Shannon)
otu %>% identify_outliers(Shannon)
```


##### Wilcox test

In our case, we have 3 groups (NC, Before and After), some of them has few samples (Before and After has less than 30 samples).
Some groups do not follow the normal distribution.

So we used a wilcoxon test :

```{r, results='hide', message=FALSE, warning=FALSE}
stat.test = asv %>% wilcox_test(Shannon~group, paired=FALSE, p.adjust.method = "fdr") %>% add_xy_position(x = "group")
stat.test 
ggboxplot(asv, x="group", y="Shannon", add=c("jitter","mean_sd"), fill="group", title = "ASV Alpha diversity changes between groups") +
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01)


stat.test = otu %>% wilcox_test(Shannon~group, paired=FALSE, p.adjust.method = "fdr") %>% add_xy_position(x = "group")
stat.test 
ggboxplot(otu, x="group", y="Shannon", add=c("jitter","mean_sd"), fill="group", title = "OTU Alpha diversity changes between groups") +
  stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01)
```

With the ASV method and the Shannon index, we see a significant difference between each group. It means there are changes in the microbiome due to the colorectal cancer and therapies.

