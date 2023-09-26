# Install packages
pacman::p_load("tidyverse", "Biostrings", "viridis", "ggpubr",
               "pals", "colorspace", "readxl", "pheatmap", "ggpmisc")

# Read in the datasets
megahit1 <- read_csv("data/DEET_hydrolase_BLAST_results/DEET_hydrolase_megahit_1_GG_readmerge.csv")
megahit2 <- read_csv("data/DEET_hydrolase_BLAST_results/DEET_hydrolase_megahit_2_non_GG_readmerge.csv")
readcounts <- bind_rows(megahit1, megahit2)

# Filter for samples of interest
percwat <- readcounts %>%
  dplyr::select(Experiment, Sample_name, Sample_type, Time, gene_prop) %>%
  dplyr::filter(grepl("BT_WW|BT_CB|V1|WW.._R._DNA_D.", Sample_name)) %>% # excludes ultrafiltration
  dplyr::mutate(expt_treatment = case_when(grepl("CB", Sample_name) ~ 0, # Stream water is 0% WW
                                           grepl("WW[0-9][0-9]", Sample_name) ~ as.numeric(gsub("WW", "", word(Sample_name, sep = "_", 1))),
                                           TRUE ~ 100)) %>%
  dplyr::mutate(expt = ifelse(grepl("D[0-9][0-9]", Sample_name), "Expt2", "Expt1")) %>%
  dplyr::filter(!is.na(Experiment)) %>%
  dplyr::filter(!duplicated(.))

# Plot
pl2 <- ggplot(percwat, aes(expt_treatment, gene_prop, color = factor(expt), fill = factor(expt))) +
  geom_point(aes(color = factor(expt))) + 
  scale_color_manual(values = c("gray60", "black")) +
  scale_fill_manual(values = c("gray80", "black")) +
  geom_smooth(method="lm", level = 0.95) +
  stat_poly_line() +
  stat_poly_eq() +
  theme_pubr() +
  theme(legend.title = element_blank()) +
  xlab("% Treated WW") +
  ylab("Normalized counts of DH \n homolog hits in metagenomes") +
  scale_x_continuous(limits = c(0, 100), breaks=c(0, 10, 30, 80, 100))
pl2

ggsave("output/Figure1a.png", pl2, width = 6, height = 6)
