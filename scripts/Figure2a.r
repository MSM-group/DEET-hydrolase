# Install packages
pacman::p_load("tidyverse", "readxl", "DECIPHER", 
               "ggplot2", "ggpubr")

# Load data
taxdat <- read_excel("data/taxonomy/table_s2_corrected.xlsx")

sumdat <- taxdat %>%
  group_by(phylum) %>%
  add_count(phylum) %>%
  dplyr::mutate(perc = n/nrow(.) * 100) %>%
  dplyr::distinct(phylum, .keep_all = T) %>%
  dplyr::mutate(role = "")

# Set color palette
pal1 <- read_csv("data/taxonomy/taxonomic_phylum_color_palette_updated.csv")

# Make barplot of taxonomic distribution
pl1 <- ggplot(sumdat, aes(x = role, y = perc, fill = phylum)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(labels = pal1$phylum, values = pal1$hexcode) +
  theme_pubr() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(t = 20,  # Top margin
                             r = 20,  # Right margin
                             b = 20,  # Bottom margin
                             l = 20), # Left margin
        legend.text = element_text(size = 12)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) 
pl1
ggsave(pl1, filename = "output/Figure2a.png", width = 4, height = 4)
