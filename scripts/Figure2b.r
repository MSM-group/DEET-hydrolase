# Load packages
# Install packages if not already installed
library(tidyverse)
library(treedataverse)
library(RColorBrewer)
library(readxl)
library(ggnewscale)
library(ggstar)
library(ggtree)

# Set color palette and define functions
color_palette <- read_csv("data/taxonomy/taxonomic_phylum_color_palette_original.csv") %>%
  bind_rows(data.frame(hexcode = "gray80", phylum = ""))

rename_seq_lab <- function(name){
  dplyr::case_when(str_detect(name, "k121")|str_detect(name, "DEET")|str_detect(name, "Cocaine") ~ name,
                   str_detect(name, "WP") ~ str_remove(name, "\\."),
                   TRUE ~ word(name, sep = "_", start = 1, end = 1))
}
generate_tip_labs <- function(lab, p_lab, deet, methylbenozoate, fournpb, fournpt){
  dplyr::case_when(lab == "DEET_hydrolase" ~ "DEET hydrolase",
                   lab == "Cocaine_esterase" ~ "Cocaine esterase",
                   deet == TRUE|methylbenozoate == TRUE|fournpb == TRUE|fournpt == TRUE ~ p_lab,
                   TRUE ~ "")
}
add_phylum <- function(name, phy){
  case_when(name %in% c("DEET_hydrolase", "1MJ5", "3HI4", "6QHV", "k121_72495_6")|str_detect(name, "WP") ~ "Proteobacteria", # naming updated to Pseudomonadota
            name %in% c("Cocaine_esterase", "k121_1133021_48") ~ "Actinobacteriota", # naming updated to Actinomycetota
            name == "k121_530580_16" ~ "Gemmatimonadota",
            TRUE ~ phy)
}
get_star_color <- function(cmpnd){
  dplyr::case_when(cmpnd == TRUE ~ "yes",
                   TRUE ~ "no")
}
get_deet_star <- function(deet, x_tip){
  case_when(deet == TRUE ~ x_tip + 0.3)
}
get_methylbenzoate_star_x <- function(deet, methylbenzoate, x_tip){
  case_when(methylbenzoate == TRUE ~ case_when(deet == TRUE ~ x_tip + 0.6,
                                               TRUE ~ x_tip + 0.3),
            TRUE ~ NA)
}
get_4np_butyrate_star_x <- function(deet, methylbenzoate, fournpb, x_tip){
  case_when(fournpb == TRUE ~ case_when(deet == TRUE & methylbenzoate == TRUE ~ x_tip + 0.9,
                                        xor(deet == TRUE, methylbenzoate == TRUE) ~ x_tip + 0.6,
                                        TRUE ~ x_tip + 0.3),
            TRUE ~ NA)
  
}
get_4np_trimethylacetate_star_x <- function(deet, methylbenozoate, fournpb, fournpt, x_tip){
  case_when(fournpt == TRUE ~ case_when(deet == TRUE & methylbenozoate == TRUE & fournpb == TRUE ~ x_tip + 1.2,
                                        (deet == TRUE & methylbenozoate == TRUE & fournpb == FALSE)|(deet == TRUE & methylbenozoate == FALSE & fournpb == TRUE)|(deet == FALSE & methylbenozoate == TRUE & fournpb == TRUE) ~ x_tip + 0.9,
                                        (deet == TRUE & methylbenozoate == FALSE & fournpb == FALSE)|(deet == FALSE & methylbenozoate == FALSE & fournpb == TRUE)|(deet == FALSE & methylbenozoate == TRUE & fournpb == FALSE) ~ x_tip + 0.6,
                                        TRUE ~ x_tip + 0.3),
            TRUE ~ NA)
  
}

# Read in the tree labels
p_labels <- readxl::read_excel("data/taxonomy/table_s2_original.xlsx", skip = 1) %>%
  dplyr::select(id, name) %>%
  dplyr::mutate(p_name = id,
                seqid = stringr::word(name, sep = "_", start = 1, end = 3)) %>%
  dplyr::select(seqid, p_name) %>%
  dplyr::add_row(seqid = "DEET_hydrolase", p_name = "DH")

assay_res <- read_rds("data/enzyme_activity_data/assay_tree_data.rds")

# Read in the taxonomy data
metadat <- read_excel("data/taxonomy/table_s2_original.xlsx", skip = 1) %>%
  mutate(seqid = word(name, sep = "_", start = 1, end = 3)) %>%
  select(seqid, phylum)
metadat$phylum[metadat$phylum == "Verrucomicrobia"] <- "Verrucomicrobiota"

# Sequence-based tree
seq_based_tree <- ape::read.tree("data/phylogeny/DH_tree_bootstrap_100.treefile") %>%
  treeio::as_tibble() %>%
  dplyr::mutate(label = rename_seq_lab(label)) %>%
  left_join(., metadat, by = c("label" = "seqid")) %>%
  left_join(., p_labels, by = c("label" = "seqid")) %>%
  left_join(., assay_res, by = c("p_name" = "id")) %>%
  dplyr::mutate(phylum = add_phylum(label, phylum)) %>%
  dplyr::mutate(lab = generate_tip_labs(label, p_name, diethyltulouamide, methylbenzoate, x4np_butyrate, x4np_trimethylacetate)) %>%
  treeio::as.treedata()
seq_based_tree_plot <- ggtree::ggtree(seq_based_tree)

# Structure-based tree
rmsd <- readr::read_rds("data/phylogeny/new_struc_dist.rds")
hc <- nj(rmsd)
tree <- hc %>% as.treedata()

struc_based_tree_no_outgroup <- ggtree::ggtree(tree)
struc_based_tree_no_outgroup
struc_based_tree_no_outgroup$data$label
struc_based_tree <- ggtree(ape::root(phy = tree, outgroup = c("6QHV", "1MJ5", "3HI4"))) #edgelabel = T)) #resolve.root = TRUE))

struc_based_tree$data$label <- c(struc_based_tree_no_outgroup$data$label)
struc_based_tree$data


dat_seq <- seq_based_tree_plot$data %>%
  dplyr::mutate(color_deet_star = get_star_color(diethyltulouamide),
                color_methylbenzoate_star = get_star_color(methylbenzoate),
                color_4np_butyrate_star = get_star_color(x4np_butyrate),
                color_4np_trimethylacetate_star = get_star_color(x4np_trimethylacetate),
                x_deet_star = get_deet_star(diethyltulouamide, x),
                x_methylbenzoate_star = get_methylbenzoate_star_x(diethyltulouamide, methylbenzoate, x),
                x_4np_butyrate_star = get_4np_butyrate_star_x(diethyltulouamide, methylbenzoate, x4np_butyrate, x),
                x_4np_trimethylacetate_star = get_4np_trimethylacetate_star_x(diethyltulouamide, methylbenzoate, x4np_butyrate, x4np_trimethylacetate, x),
                x_tip_lab_with_star = case_when(lab == "Cocaine esterase" ~ x + 0.3,
                                                TRUE ~ pmax(x_deet_star, x_methylbenzoate_star, x_4np_butyrate_star, x_4np_trimethylacetate_star, na.rm = TRUE)+0.3))

dat_struc <- struc_based_tree_no_outgroup$data %>%
  mutate(x = max(x) - x + max(dat_seq$x) + 5) %>%
  left_join(., metadat, by = c("label" = "seqid")) %>%
  left_join(., p_labels, by = c("label" = "seqid")) %>%
  left_join(., assay_res, by = c("p_name" = "id")) %>%
  dplyr::mutate(phylum = add_phylum(label, phylum),
         lab = generate_tip_labs(label, p_name, diethyltulouamide, methylbenzoate, x4np_butyrate, x4np_trimethylacetate))
dat_combined <- bind_rows(dplyr::select(dat_seq, !ends_with("star")), dat_struc) %>%
  tidyr::drop_na(label) %>%
  filter(parent != node)

struc_based_tree_plot_2 <- ggtree::ggtree(dat_struc, color = "gray")

combined_tree <- struc_based_tree_plot_2 +
  geom_nodelab() +
  geom_tree(dat_seq, color = "gray") +
  geom_line(aes(x = x, y = y, group = label, color = phylum), data = filter(dat_combined, isTip == TRUE), linewidth = 0.5, alpha = 0.5) +
  geom_tiplab(aes(label = lab), size = 2.5, hjust = 1, offset = -0.1) +
  geom_tiplab(data = dat_seq, aes(x = x_tip_lab_with_star, label = lab), size = 2.5) +
  geom_point(data = drop_na(dat_combined, phylum), aes(x = x, y = y, color = phylum), size = 1, shape = 15) +
  scale_color_manual(name = "Phylum", values = pull(color_palette, hexcode)) +
  ggnewscale::new_scale_color() +
  geom_nodepoint(data = filter(dat_seq, isTip == FALSE & parent != node), aes(color = as.numeric(label)), size = 1) +
  scale_color_gradient(low = "gray80", high = "black", name = "Bootstrap support (%)") +
  geom_star(data = dat_seq, aes(x = x_deet_star, y = y, fill = color_deet_star), starstroke = 0.1) +
  scale_fill_manual(values = c("red", "white"), breaks = c("yes", "no"), guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_star(data = dat_seq, aes(x = x_methylbenzoate_star, y = y, fill = color_methylbenzoate_star), starstroke = 0.1) +
  scale_fill_manual(values = c("blue", "white"), breaks = c("yes", "no"), guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_star(data = dat_seq, aes(x = x_4np_butyrate_star, y = y, fill = color_4np_butyrate_star), starstroke = 0.1) +
  scale_fill_manual(values = c("orange", "white"), breaks = c("yes", "no"), guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_star(data = dat_seq, aes(x = x_4np_trimethylacetate_star, y = y, fill = color_4np_trimethylacetate_star), starstroke = 0.1) +
  scale_fill_manual(values = c("darkgreen", "white"), breaks = c("yes", "no"), guide = "none")


ggplot2::ggsave(combined_tree,
                filename = "output/Figure2b.png",
                dpi = 320,
                device = "png",
                units = "in",
                width = 12,
                height = 5)


