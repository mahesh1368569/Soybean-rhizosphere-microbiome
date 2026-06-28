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
library(ggdendro)
install.packages("ggdendro")

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

soy_16s_mecodataset$otu_table <- round(soy_16s_mecodataset$otu_table)
storage.mode(soy_16s_mecodataset$otu_table) <- "integer"

phyis <- readxl::read_xlsx("Microbiome_files/Input files/Soybean_env.xlsx")
phyis <- as.data.frame(phyis)

rownames(phyis) = phyis[, 1]

phyis = phyis[ ,-1]

# add_data is used to add the environmental data
env_soy <- trans_env$new(dataset = soy_16s_mecodataset, add_data = phyis)

##Check integer ###
otu <- soy_16s_mecodataset$otu_table

# microeco OTU table is usually ASVs as rows, samples as columns
sum(otu != round(otu), na.rm = TRUE)
colSums(otu)[1:10]

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
# Alpha #
soy_16s_mecodataset$cal_abund()
soy_16s_mecodataset$cal_alphadiv()
soy_16s_mecodataset$cal_betadiv()

alpha <- trans_alpha$new(dataset = soy_16s_mecodataset, group = "Nitrogen")
alpha_S <- trans_alpha$new(dataset = soy_16s_mecodataset, group = "Sulfur")

alpha_NS <- trans_alpha$new(dataset = soy_16s_mecodataset, group = "Nitrogen", by_group = "Sulfur")

alpha_NS$cal_diff(method = "KW")


soy_N_chao1 = alpha$plot_alpha(measure = "Chao1", y_increase = 0.3)+ 
  geom_boxplot(aes(fill = Nitrogen),  position = position_dodge(0.8)) + 
  geom_jitter(aes(fill = Nitrogen), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Shannon diversity index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.position = "top",# Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(0.8, "cm"), # Adjust legend key size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

soy_N_chao1

ggsave("soy_N_chao1.pdf", plot =soy_N_chao1, device = "pdf", width = 10, height = 8, units = "in", dpi = 1000)


soy_N_shannon = alpha$plot_alpha(measure = "Shannon", y_increase = 0.3) + 
  geom_boxplot(aes(fill = Nitrogen),  position = position_dodge(0.8)) + 
  geom_jitter(aes(fill = Nitrogen), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Shannon diversity index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.position = "top",# Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(0.8, "cm"), # Adjust legend key size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

soy_N_shannon

soy_S_shannon = alpha_S$plot_alpha(measure = "Shannon", y_increase = 0.3) + 
  geom_boxplot(aes(fill = Sulfur),  position = position_dodge(0.8)) + 
  geom_jitter(aes(fill = Sulfur), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), alpha = 0.5) +
  labs(title = "",
       x = "",
       y = "Shannon diversity index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.position = "top",# Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(0.8, "cm"), # Adjust legend key size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

soy_S_shannon

soy_NS_shannon = alpha_NS$plot_alpha(measure = "Shannon", y_increase = 0.3) +
  geom_boxplot(aes(fill = Nitrogen),  position = position_dodge(0.8)) +
  geom_point(aes(fill = Nitrogen), position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8)) +
  labs(title = "",
       x = "",
       y = "Shannon diversity index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.position = "top",# Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(0.8, "cm"), # Adjust legend key size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))


soy_NS_chao1 = alpha_NS$plot_alpha(measure = "Chao1", y_increase = 0.3) +
  geom_boxplot(aes(fill = Nitrogen),  position = position_dodge(0.8)) +
  geom_point(aes(fill = Nitrogen), position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8)) +
  labs(title = "",
       x = "",
       y = "Chao1 index") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.position = "top",# Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(0.8, "cm"), # Adjust legend key size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

soy_NS_chao1

ggsave("Microbiome_files/Output/soy_NS_chao1.pdf", plot = soy_NS_chao1, width = 9, 
       height = 7, units = "in", dpi = 1000)
ggsave("Microbiome_files/Output/soy_NS_shannon.pdf", plot = soy_NS_shannon, width = 9, 
       height = 7, units = "in", dpi = 1000)

#---------------------------------------------------------------------------------------------------------------------
# Beta #

soy_beta <- trans_beta$new(dataset = soy_16s_mecodataset, group = "Nitrogen", measure = "bray")

soy_beta$cal_ordination(method = "NMDS")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
soy_beta_plot = soy_beta$plot_ordination(plot_color = "Nitrogen", plot_shape = "Sulfur", plot_type = c("point", "ellipse"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.position = "right",# Adjust legend title size
        legend.text = element_text(size = 12),  # Adjust legend text size
        legend.key.size = unit(0.8, "cm"), # Adjust legend key size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave("Microbiome_files/Output/soy_beta_plot.pdf", plot = soy_beta_plot, width =8, 
       height = 5, units = "in", dpi = 1000)

#---------------------------------------------------------------------------------------------------------------------
# Clustering #

# extract a part of data
tmp <- clone(soy_16s_mecodataset)
tmp$tidy_dataset()

t1 <- trans_beta$new(
  dataset = tmp,
  group = "Nitrogen",
  measure = "bray"
)

soy_N_cluster <- t1$plot_clustering(
  group = "Nitrogen",
  replace_name = c("Nitrogen")
)

soy_N_cluster

t2 <- trans_beta$new(
  dataset = tmp,
  group = "Sulfur",
  measure = "bray"
)

soy_S_cluster <- t1$plot_clustering(
  group = "Sulfur",
  replace_name = c("Sulfur")
)

soy_S_cluster

#---------------------------------------------------------------------------------------------------------------------
# Clustering #