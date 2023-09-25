# Install packages
pacman::p_load("tidyverse", "ggpubr", "ggdist", "readxl", "Biostrings")

# Read in the data
fils <- list.files("data/enzyme_activity_data/", full.names = T, pattern = "calculated_slopes.csv")

fils <- fils[!grepl(paste0(c("DEV", "SHINS1", "Thermostability", "pH", "°C"), collapse = "|"), fils)]
fils

rawdat <- tibble(filename = fils) %>%
  mutate(file_contents = map(filename,          # read files into
                             ~ read_csv(file.path(.))) # a new data column
  ) %>%
  unnest(cols = c(file_contents)) %>%
  dplyr::mutate(cmpnd =  paste0(word(filename, 4, sep = "_"), "_", word(filename, 5, sep = "_"))) %>%
  dplyr::mutate(enzyme = word(org, 1, sep = "_")) %>%
  dplyr::filter(enzyme != "NA") %>%
  dplyr::mutate(biological_replicate = gsub("bio", "", word(org, 2, sep = "_"))) %>%
  dplyr::mutate(technical_replicate = word(org, 3, sep = "_")) %>%
  dplyr::filter(!grepl("20230309", filename))
length(unique(rawdat$enzyme))

# Read in taxonomy information
tax <- read_csv("data/taxonomy/kaiju_taxonomy_align.csv") %>%
  dplyr::mutate(enzyme = word(name, sep = "Batch353_", -1))

# Write to file
maprdat_raw <- rawdat %>%
  dplyr::mutate(hr_slope = max_slope * 10 * 60) %>%
  dplyr::mutate(log_slope = log10(hr_slope)) %>%
  left_join(., tax, by = "enzyme")
table(maprdat_raw$phylum)
table(maprdat_raw$family)
maprdat_raw$phylum[maprdat_raw$enzyme == "DH"] <- "Proteobacteria"
maprdat_raw$phylum[maprdat_raw$phylum == "Verrucomicrobia"] <- "Verrucomicrobiota"
maprdat_raw$phylum[is.na(maprdat_raw$phylum)] <- " "
#maprdat_raw$phylum[maprdat_raw$enzyme == "EV"] <- "EV"


# Set threshold for activity
thresh <- 2.135  #summary(maprdat_raw$log_slope)
summary(maprdat_raw$log_slope)

tokeep <- maprdat_raw %>%
  dplyr::group_by(enzyme, cmpnd) %>%
  mutate(flag = ifelse(log_slope <= thresh, 1, 0)) %>%
  summarise(any_flagged = max(flag)) %>%
  dplyr::filter(any_flagged == 0) %>%
  dplyr::pull(enzyme) %>%
  unique()


maprdat_long <- maprdat_raw %>%
  dplyr::filter(enzyme %in% c(tokeep, "EV"))

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
# Remove after the empty vector
testr <- maprdat_long[order(desc(maprdat_long$mean_slope)),]
testr
alph <- unique(testr$enzyme)
alph
b2 <- alph[1:grep("NEG", alph)]
b2
maprdat_long$enzyme <- factor(maprdat_long$enzyme,
                              levels = b2, ordered = TRUE)
levels(maprdat_long$phylum)

pal2 <- read_csv("data/taxonomy/taxonomic_phylum_color_palette.csv") %>%
  dplyr::filter(phylum %in% levels(maprdat_long$phylum))

pal <- c("gray60", pal2$hexcode)
pal

maprdat_long <- maprdat_long %>%
  dplyr::filter(!is.na(enzyme))


pdf("output/draft_fig_combined_slopes_trimmed.pdf")
p3 <- ggplot(data = maprdat_long) + 
  facet_grid(cols = vars(cmpnd)) +
  geom_boxplot(aes(y = enzyme, x = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_point(aes(y = enzyme, x = log_slope, color = phylum, shape = biological_replicate), alpha = 0.5, size = 2) + #position=position_jitter(width=.3, height=0.0))
  #   geom_jitter(aes(y = enzyme, x = log_slope, color = phylum, shape = biological_replicate), alpha = 0.5, size = 2, position=position_jitter(width=0.0, height=0.0)) +
  labs_pubr() +
  theme_pubr() +
  xlab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  ylab("Enzyme") +
  scale_color_manual(name='Phylum', values = pal) +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
    axis.title.y = element_blank(),
    #axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
    legend.title = element_text(),
    legend.position = "left") +
  #  axis.text.x = element_text(angle = -75)) +
  scale_x_continuous(limits = c(1,4), breaks = c(1,1.5,2, 2.5, 3, 3.5, 4)) +
  scale_y_discrete(position = "right", limits = rev(levels(maprdat_long$enzyme)))+
  geom_vline(aes(xintercept = 2.135), linetype = "dashed", color = "gray60")
p3
dev.off()
rev(levels(maprdat_long$enzyme))
levels(maprdat_long$enzyme)
ggsave( "output/first_plot_idea.png", p3, width = 8, height = 5)



pdf("output/4NP_butyrate_draft_fig_combined_slopes.pdf")
ggplot(subset(maprdat_long, cmpnd %in% c("4NP_butyrate"))) + 
  #geom_violin(aes(x = enzyme, y = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_boxplot(aes(x = enzyme, y = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_jitter(aes(x = enzyme, y = log_slope, color = as.factor(filename), shape= as.factor(biological_replicate)), alpha = 0.5, size = 2, position=position_jitter(width=.3, height=0.0)) +
  # geom_point(aes(x = cmpnd, y = log_slope_mean_slope),
  #            fill = "firebrick", color = "black", shape = 23, size = 3) +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Enzyme") +
  guides(shape=guide_legend("Biological replicate")) +
  scale_color_manual(values = c("firebrick", "cornflowerblue", "forestgreen", "pink", "orange")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "right")
dev.off()

pdf("output/4NP_trimethylacetate_draft_fig_combined_slopes.pdf")
ggplot(subset(maprdat_long, cmpnd %in% c("4NP_trimethylacetate"))) + 
  #geom_violin(aes(x = enzyme, y = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_boxplot(aes(x = enzyme, y = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_jitter(aes(x = enzyme, y = log_slope, color = as.factor(filename), shape= as.factor(biological_replicate)), alpha = 0.5, size = 2, position=position_jitter(width=.3, height=0.0)) +
  # geom_point(aes(x = cmpnd, y = log_slope_mean_slope),
  #            fill = "firebrick", color = "black", shape = 23, size = 3) +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Enzyme") +
  guides(shape=guide_legend("Biological replicate")) +
  scale_color_manual(values = c("firebrick", "cornflowerblue", "forestgreen", "pink", "orange")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
#legend.title = element_blank()) 
dev.off()


ggplot(data = maprdat_long) + 
  ggdist::stat_halfeye(aes(x = enzyme, y = log_slope)) + #color = "firebrick") +
  #geom_boxplot(aes(x = enzyme, y = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_jitter(aes(x = enzyme, y = log_slope, shape = as.factor(biological_replicate)), alpha = 0.5, size = 2, position=position_jitter(width=.3, height=0.0)) +
  # geom_point(aes(x = cmpnd, y = log_slope_mean_slope),
  #            fill = "firebrick", color = "black", shape = 23, size = 3) +
  labs_pubr() +
  theme_pubr() +
  ylab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  xlab("Enzyme") +
  # scale_color_manual(values = c("firebrick", "cornflowerblue", "forestgreen")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.title = element_blank()) 

pdf("output/4NP_butyrate_distribution_plot.pdf") 
ggplot(data = maprdat_long) + 
  ggdist::stat_halfeye(aes(y = enzyme, x = log_slope)) + #color = "firebrick") +
  #geom_boxplot(aes(x = enzyme, y = log_slope), outlier.color = NA, colour = "gray60", fill = "gray80", alpha = 0.5) +
  geom_jitter(aes(y = enzyme, x = log_slope, shape = as.factor(biological_replicate)), alpha = 0.5, size = 2, position=position_jitter(width=.3, height=0.0)) +
  # geom_point(aes(x = cmpnd, y = log_slope_mean_slope),
  #            fill = "firebrick", color = "black", shape = 23, size = 3) +
  labs_pubr() +
  theme_pubr() +
  xlab("Enzyme activity log(nmol pNP/ OD 1/ hr)") +
  ylab("Enzyme") +
  
  # scale_color_manual(values = c("firebrick", "cornflowerblue", "forestgreen")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.6, hjust = 0.8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
# legend.title = element_blank()) 
dev.off()


