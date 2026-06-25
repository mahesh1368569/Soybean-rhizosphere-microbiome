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

## Importing files

soy_16s_biom = import_biom("Soybean_16smod.biom")

soy_16s_biom1 = import_biom("Soybean_16s1.biom")

metadata = import_qiime_sample_data("Metadata_16S_Soybean.txt")

tree = read_tree("unrooted_tree.nwk")

rep_fasta = readDNAStringSet("mspb_seq.fasta", format = "fasta")

soy_16s_biom1 = merge_phyloseq(soy_16s_biom1, metadata, rep_fasta, tree)

#rename columns in taxonomy table
colnames(tax_table(soy_16s_biom1)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

colnames(tax_table(soy_16s_biom))


soy_16s_mecodataset <- phyloseq2meco(soy_16s_biom)

soy_16s_mecodataset1 <- phyloseq2meco(soy_16s_biom1)

soy_16s_mecodataset$tax_table

soy_16s_mecodataset

## Composition analysis

soy_16s_mecodataset$cal_abund()

abun = trans_abund$new(dataset = soy_16s_mecodataset, taxrank = "Phylum", ntaxa = 15)

# The groupmean parameter can be used to obtain the group-mean barplot.
abun_trt_N <- trans_abund$new(dataset = soy_16s_mecodataset, taxrank = "Phylum", groupmean = "Nitrogen")

abun_trt_N

abundance_N = abun$plot_box(group = "Nitrogen", xtext_angle = 30) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Increase x-axis text size
        axis.text.y = element_text(size = 18), # Increase y-axis text size
        axis.title.x = element_text(size = 20), # Increase x-axis label size
        axis.title.y = element_text(size = 20), # Increase y-axis label size
        strip.text = element_text(size = 18),
        legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

ggsave("abundance_N.pdf", plot = abundance_N, device = "pdf", width = 16, height = 10, units = "in", dpi = 1000)


## Ven diagram

# merge samples as one community for each group
venn_1 <- soy_16s_mecodataset$merge_samples("Nitrogen")
# tmp is a new microtable object

# create trans_venn objec
Fig_venn<- trans_venn$new(venn_1, ratio = "seqratio")

Fig_venn$plot_venn()

# Alpha diversity metrics

alpha <- trans_alpha$new(dataset = soy_16s_mecodataset) # Alpha diversity metrics

alpha <- trans_alpha$new(dataset = soy_16s_mecodataset, group = "Nitrogen")

alpha1 <- trans_alpha$new(dataset = soy_16s_mecodataset1, group = "Nitrogen")

#calculating significance tests

alpha$cal_diff(method = "KW")

alpha1$cal_diff(method = "KW")

soy_N_chao1 = alpha1$plot_alpha(measure = "Chao1", y_increase = 0.3)+ 
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

ggsave("soy_N_shannon.pdf", plot =soy_N_shannon, device = "pdf", width = 10, height = 8, units = "in", dpi = 1000)



#calculating significance tests

alpha$cal_diff(method = "KW")

alpha$plot_alpha(measure = "Chao1", y_increase = 0.3)

alpha$plot_alpha(measure = "Shannon", y_increase = 0.3)

## Betadiversity metrics

soy_16s_mecodataset$cal_betadiv()

beta_N <- trans_beta$new(dataset = soy_16s_mecodataset, group = "Nitrogen", measure = "bray")

beta_N$cal_ordination(method = "PCoA")

# Custom palette with 8 distinct colors
custom_colors = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")



sulfur_gradient = c("#1f77b4", "#ff7f0e", "#9467bd", "#d62728")



# plot the PCoA result with confidence ellipse
N_beta = beta_N$plot_ordination(plot_color = "Nitrogen", plot_type = c("point", "ellipse"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Increase x-axis text size
        axis.text.y = element_text(size = 18), # Increase y-axis text size
        axis.title.x = element_text(size = 20), # Increase x-axis label size
        axis.title.y = element_text(size = 20), # Increase y-axis label size
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

N_beta

ggsave("MSPB_N_beta.pdf", plot =N_beta, device = "pdf", width = 14, height = 8, units = "in", dpi = 1000)

#sulfur

beta_S <- trans_beta$new(dataset = soy_16s_mecodataset, group = "Sulfur", measure = "bray")
beta_S$cal_ordination(method = "PCoA")

# plot the PCoA result with confidence ellipse
S_beta = beta_S$plot_ordination(plot_color = "Sulfur", color_values = sulfur_gradient,  plot_type = c("point", "ellipse"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Increase x-axis text size
        axis.text.y = element_text(size = 18), # Increase y-axis text size
        axis.title.x = element_text(size = 20), # Increase x-axis label size
        axis.title.y = element_text(size = 20), # Increase y-axis label size
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

S_beta

ggsave("MSPB_S_beta.pdf", plot =S_beta, device = "pdf", width = 14, height = 8, units = "in", dpi = 1000)

#-------------------------------Network analysis in microeco--------------------------------#
# 1. First constructing two correlation networks for treatment groups (Control and Cover crop)
soybean_network = list()

no_N = clone(soy_16s_mecodataset)

## change sample_table directly

no_N$sample_table %<>% subset(Nitrogen == "0_Nitrogen")

# trim all files in the objectNitrogen# trim all files in the object
no_N$tidy_dataset()

# use filter_thres parameter to filter the feature with low relative abundance
no_N_network <- trans_network$new(dataset = no_N, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold

no_N_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# put the network into the list
soybean_network$No_Nitrogen <- no_N_network

# Cover crop samples

high_N = clone(soy_16s_mecodataset)

## change sample_table directly

high_N$sample_table %<>% subset(Nitrogen == "200_Nitrogen")

# trim all files in the object
high_N$tidy_dataset()

# use filter_thres parameter to filter the feature with low relative abundance
high_N_network <- trans_network$new(dataset = high_N, cor_method = "spearman", filter_thres = 0.0005)
# COR_p_thres represents the p value threshold
# COR_cut denotes the correlation coefficient threshold

high_N_network$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# put the network into the list
soybean_network$High_Nitrogen <- high_N_network

#1.1 Network modularity for all networks

soybean_network %<>% cal_module(undirected_method = "cluster_fast_greedy")

soybean_network[["No_Nitrogen"]]$save_network(filepath = "No-Nitrogen.gexf")

soybean_network[["High_Nitrogen"]]$save_network(filepath = "High-Nitrogen.gexf")

#1.2 Network topological attributes for all networks

network_atr_soy = cal_network_attr(soybean_network)

write.csv(network_atr_soy, "network_attributes_soybean.csv")

#1.3 Node and edge properties extraction for all networks

soybean_network %<>% get_node_table(node_roles = TRUE) %>% get_edge_table

#1.4 compare nodes across networks

# obtain the node distributions by searching the res_node_table in the object
node_dist <- node_comp(soybean_network, property = "name")

# obtain nodes intersection
node_intersection <- trans_venn$new(node_dist, ratio = "numratio")

node_intersection

node_intersection_plot <- node_intersection$plot_venn(fill_color = FALSE)

node_intersection_plot

ggsave("soil_amp_node_overlap.pdf", g1, width = 7, height = 6)

# calculate jaccard distance to reflect the overall differences of networks
node_dist$cal_betadiv(method = "jaccard")

node_dist$beta_diversity$jaccard

#1.5 compare edges across networks

# get the edge distributions across networks
edge_dist <- edge_comp(soybean_network)

# obtain edges intersection
edge_intersection <- trans_venn$new(edge_dist, ratio = "numratio")

edge_interection_plot <- edge_intersection$plot_venn(fill_color = FALSE)

edge_interection_plot 

# calculate jaccard distance
edge_dist$cal_betadiv(method = "jaccard")

edge_dist$beta_diversity$jaccard

#1.6 Extract overlapped edges of netoworks to a new network

# first obtain edges distribution and intersection
edge_dist

edge_dist_plot = trans_venn$new(edge_dist)

edge_dist_plot

# convert intersection result to a microtable object
intersection = edge_dist_plot$trans_comm()

# extract the intersection of all the two networks (cover crop and control)
# please use colnames(tmp2$otu_table) to find the required name
Intersec_all <- subset_network(soybean_network, venn = intersection, name = "no_N&high_N")

# Intersec_all is a trans_network object
# for example, save Intersec_all as gexf format
Intersec_all$save_network("Intersec_all.gexf")

#1.7 Compare phylogeneitc distances of paired nodes in edges

# filter useless features to speed up the calculation
node_names <- unique(unlist(lapply(cc_ms_network, function(x){colnames(x$data_abund)})))

filter_cc_ms<- microeco::clone(CCMS_mecodataset)

filter_cc_ms$otu_table <- filter_cc_ms$otu_table[node_names, ]

filter_cc_ms$tidy_dataset()

# obtain phylogenetic distance matrix
phylogenetic_distance_ccms <- as.matrix(cophenetic(filter_cc_ms$phylo_tree))

# use both the positive and negative labels
ccedge.compare <- edge_node_distance$new(network_list = cc_ms_network, dis_matrix = phylogenetic_distance_ccms, label = c("+", "-"))

ccedge.compare$cal_diff(method = "anova")

# visualization

ccedge.compare_plot <- ccedge.compare$plot( add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")

ccedge.compare_plot

ggsave("soil_amp_phylo_distance.pdf", g1, width = 7, height = 6)

# show different modules with at least 10 nodes and positive edges
modules <- edge_node_distance$new(network_list = cc_ms_network, dis_matrix = phylogenetic_distance_ccms, 
                                  label = "+", with_module = TRUE, module_thres = 10)
modules$cal_diff(method = "anova")

modules_plot <- modules$plot(add_sig = TRUE, add_sig_text_size = 5) + ylab("Phylogenetic distance")

modules_plot 
ggsave("soil_amp_phylo_distance_modules.pdf", g1, width = 8, height = 6)


#1.8 Compare node sources of edges across networks
soil_cc_network_edgetax <- edge_tax_comp(cc_ms_network, taxrank = "Phylum", label = "+", rel = TRUE)

# filter the features with small number
soil_cc_network_edgetax <- soil_cc_network_edgetax[apply(soil_cc_network_edgetax, 1, mean) > 0.01, ]

# visualization
g1 <- pheatmap::pheatmap(soil_cc_network_edgetax, display_numbers = TRUE)

ggsave("soil_amp_edge_tax_comp.pdf", g1, width = 7, height = 7)

#1.9 compare topological properties of sub-networks
# calculate global properties of all sub-networks
topological <- subnet_property(cc_ms_network)

# then prepare the data for the correlation analysis
# use sample names (second column) as rownames
rownames(topological) <- topological[, 2]

# delete first two columns (network name and sample name)
topological <- topological[, -c(1:2)]


# generate correlation heatmap
g1 <- tmp1$plot_cor()

ggsave("soil_amp_subnet_property.pdf", g1, width = 11, height = 5)

#1.10 Robustness

?robustness

robustness_cc <- robustness$new(cc_ms_network, remove_strategy = c("edge_rand", "edge_strong", "node_rand", "node_degree_high"), 
                                remove_ratio = seq(0, 0.99, 0.1), measure = c("Eff", "Eigen", "Pcr"), run = 10)

View(robustness_cc$res_table)

View(robustness_cc$res_summary)

robustness_cc$plot(linewidth = 1)

robustness_df = robustness_cc[["res_summary"]]

# Assuming your dataframe is named df
ggplot(data = robustness_df, aes(x = factor(remove_ratio), y = Mean, fill = Network)) +
  geom_boxplot() +
  facet_grid(measure ~ .) +  # Faceting by measure on rows
  labs(title = "Boxplot of Network Robustness by Network Type",
       x = "Removal Ratio",
       y = "Robustness (Mean)") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")  # Using a different palette for better distinction between Network types


#1.11 Vulnerability of nodes

vul_table_cc <- vulnerability(cc_ms_network)

View(vul_table_cc)

############################## Functional analysis ##########################################

t2 <- trans_func$new(soy_16s_mecodataset1)

t2$cal_spe_func(prok_database = "FAPROTAX")

t2$res_spe_func

###### bar plot ########

t2$cal_spe_func_perc(abundance_weighted = F)

t2$trans_spe_func_perc()

faprotax_df_t <- as.data.frame(t2$res_spe_func_perc_trans)

#functional_long <- faprotax_df_t %>%
#rownames_to_column(var = "Function") %>%  # Convert rownames (functions) into a column
#pivot_longer(cols = -Function, names_to = "Sample_ID", values_to = "Relative_Abundance")

meta_func = import_qiime_sample_data("meta-function.txt")

# Step 4: Merge Functional Data with Metadata
faprotax_df_t <- faprotax_df_t  %>%
  left_join(meta_func, by = "sampname")  # Ensure metadata has a column named "Sample"

str(faprotax_df_t)

faprotax_df_t$Nitrogen <- as.factor(faprotax_df_t$Nitrogen)

faprotax_df_t$Sulfur <- as.factor(faprotax_df_t$Sulfur)

faprotax_df_t$group

plot_N <- ggplot(faprotax_df_t %>% filter(group == "N-cycle"), 
                 aes(x = variable, y = value, fill = Nitrogen)) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  theme_minimal() +
  labs(title = "N-Cycle Functional Pathways", x = "Pathway", y = "Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top")

plot_N

ggsave("plot_N.pdf", plot_N, width = 7, height = 6, dpi = 1000)

plot_C <- ggplot(faprotax_df_t %>% filter(group == "C-cycle"), 
                 aes(x = variable, y = value, fill = Nitrogen)) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  theme_minimal() +
  labs(title = "C-Cycle Functional Pathways", x = "Pathway", y = "Relative Abundance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.position = "top")

plot_C

ggsave("plot_C.pdf", plot_C, width = 7, height = 6, dpi = 1000)