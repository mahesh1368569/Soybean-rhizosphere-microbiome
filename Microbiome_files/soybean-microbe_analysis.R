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

rownames(phyis) = phyis[, 1]

phyis = phyis[ ,-1]

phyis_physiology <- phyis[, c(1:12)]

phyis_plant <- phyis[, c(13:18)]

phyis_soil = phyis[, c(19:21)]


# add_data is used to add the environmental data
env_sujan <- trans_env$new(dataset = meco_dataset, add_data = phyis)

env_sujan_phyis <- trans_env$new(dataset = meco_dataset, add_data = phyis_physiology)

env_sujan_plant <- trans_env$new(dataset = meco_dataset, add_data = phyis_plant)

env_sujan_soil <- trans_env$new(dataset = meco_dataset, add_data = phyis_soil)

