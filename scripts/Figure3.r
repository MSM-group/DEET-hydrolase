# Install packages
pacman::p_load("tidyverse", "ggpubr", "ggdist", "readxl", "Biostrings")

# Read in the data
fils <- list.files("data/enzyme_activity_data/", full.names = T, pattern = "calculated_slopes.csv")
fils

rawdat <- tibble(filename = fils) %>%
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(cols = c(file_contents)) %>%
  dplyr::mutate(cmpnd =  paste0(word(filename, 6, sep = "_"), "_", word(filename, 7, sep = "_"))) %>%
  dplyr::mutate(enzyme = word(org, 1, sep = "_")) %>%
  dplyr::filter(enzyme != "NA") %>%
  dplyr::mutate(biological_replicate = gsub("bio", "", word(org, 2, sep = "_"))) %>%
  dplyr::mutate(technical_replicate = word(org, 3, sep = "_")) %>%
  dplyr::filter(!grepl("20230309", filename))

# Read in taxonomy information
tax <- read_csv("data/taxonomy/kaiju_taxonomy_align.csv") %>%
  dplyr::mutate(enzyme = word(name, sep = "Batch353_", -1))

# Write to file
maprdat_raw <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) %>%
  left_join(., tax, by = "enzyme")
table(maprdat_raw$phylum)

# Fix the naming
maprdat_raw$phylum[maprdat_raw$enzyme == "DH"] <- "Proteobacteria"
maprdat_raw$phylum[maprdat_raw$phylum == "Verrucomicrobia"] <- "Verrucomicrobiota"
maprdat_raw$phylum[is.na(maprdat_raw$phylum)] <- " "
maprdat_raw$phylum[maprdat_raw$enzyme == "EV"] <- " Negative control"
maprdat_raw$enzyme[maprdat_raw$enzyme == "EV"] <- "NEG"

# Draw line on graph for average activity
thresh <- summary(maprdat_raw$log_slope)[4] # mean = 2.135  

tokeep <- maprdat_raw %>%
  dplyr::group_by(enzyme, cmpnd) %>%
  mutate(flag = ifelse(log_slope <= thresh, 1, 0)) %>%
  summarise(any_flagged = max(flag)) %>%
  dplyr::filter(any_flagged == 0) %>%
  dplyr::pull(enzyme) %>%
  unique()

maprdat_long <- maprdat_raw %>%
  dplyr::filter(enzyme %in% c(tokeep, "NEG"))

# Order by activity
testr <- maprdat_long[order(maprdat_long$log_slope, decreasing = T),]
alph <- unique(testr$enzyme)
alph

# Reorder the factors
maprdat_long$enzyme <- factor(maprdat_long$enzyme,
                              levels = alph, ordered = TRUE)
maprdat_long$phylum <- as.factor(maprdat_long$phylum)
maprdat_long$biological_replicate <- as.factor(maprdat_long$biological_replicate)
table(maprdat_long$enzyme, maprdat_long$cmpnd)

pal2 <- read_csv("data/taxonomy/taxonomic_phylum_color_palette.csv") %>%
  dplyr::filter(phylum %in% levels(maprdat_long$phylum))

pal <- c("gray60", pal2$hexcode)

maprdat_long <- maprdat_long %>%
  dplyr::filter(!is.na(enzyme))

p3 <- ggplot(data = maprdat_long) + 
  facet_grid(cols = vars(cmpnd)) +
  geom_boxplot(aes(y = enzyme, x = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_point(aes(y = enzyme, x = log_slope, color = phylum, shape = biological_replicate), alpha = 0.5, size = 2) + 
  labs_pubr() +
  theme_pubr() +
  xlab("Enzyme activity") +
  ylab("Enzyme") +
  scale_color_manual(name='Phylum', values = pal) +
  theme(axis.title.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.title = element_text(),
    legend.position = "left") +
  scale_x_continuous(limits = c(1,4), breaks = c(1,1.5,2, 2.5, 3, 3.5, 4)) +
  scale_y_discrete(position = "right", limits = rev(levels(maprdat_long$enzyme))) +
  geom_vline(aes(xintercept = 2.135), linetype = "dashed", color = "gray60")
p3

ggsave( "output/Figure3.png", p3, width = 8, height = 5)

