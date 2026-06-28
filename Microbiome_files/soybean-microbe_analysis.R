library(phyloseq)
library(biomformat)
library(BiocManager)
library(file2meco)
library(MicrobiomeStat)
library(WGCNA)
library(ggtree)
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)
library(microeco)
library(ape)
library(plyr)
library(magrittr)
library(tidygraph)
library(ggcor)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(vegan)
library(devtools)
library(MicEco)
library(psych)
library(igraph)
library(ggpubr)
library(Hmisc)
library(minpack.lm)
library(stats4)
library(EcoSimR) # load EcoSimR library
library(devEMF)
library(NST)
library(meconetcomp)
library(magrittr)
library(igraph)
library(Biostrings)
devtools::install_github("kassambara/ggcorrplot")

## Importing files

soy_16s_biom = import_biom("Microbiome_files/Input files/Soybean_16s.biom")

metadata = import_qiime_sample_data("Microbiome_files/Input files/Metadata_16S_Soybean.txt")

tree = read_tree("Microbiome_files/Input files/unrooted_tree.nwk")

rep_fasta = readDNAStringSet("Microbiome_files/Input files/mspb_seq.fasta", format = "fasta")

soy_16s_biom = merge_phyloseq(soy_16s_biom, metadata, rep_fasta, tree)

#rename columns in taxonomy table
colnames(tax_table(soy_16s_biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

colnames(tax_table(soy_16s_biom))

soy_16s_mecodataset <- phyloseq2meco(soy_16s_biom)

soy_16s_mecodataset$tax_table

soy_16s_mecodataset

phyis <- readxl::read_xlsx("Microbiome_files/Input files/Soybean_env.xlsx")
phyis <- as.data.frame(phyis)

rownames(phyis) = phyis[, 1]

phyis = phyis[ ,-1]

# add_data is used to add the environmental data
env_soy <- trans_env$new(dataset = soy_16s_mecodataset, add_data = phyis)

?trans_norm

#---------------------------------------------------------------------------------------------------------------------
# Rarefaction #

# OTU/ASV table from microeco object
otu <- t(soy_16s_mecodataset$otu_table)   # vegan wants samples as rows, ASVs as columns
depth <- rowSums(otu)

summary(depth)
sort(depth)
hist(depth, breaks = 30, main = "Sequencing depth", xlab = "Reads per sample")

quantile(depth, c(0, 0.05, 0.10, 0.25, 0.50))

#---------------------------------------------------------------------------------------------------------------------
# Rarefaction #
