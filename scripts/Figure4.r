# Install packages
pacman::p_load("tidyverse", "ggpubr", "ggdist", "readxl", "Biostrings")

# Read in the data
fils <- list.files("data/HPLC_data/", full.names = T, pattern = ".csv")
fils

rawdat <- tibble(filename = fils) %>%
  mutate(file_contents = map(filename,         
                             ~ read_csv(file.path(.)))) %>%
  unnest(cols = c(file_contents)) %>%
  dplyr::mutate(biorep =  word(filename, -3, sep = "_")) %>%
  dplyr::mutate(enzyme = word(variable, 1, sep = "_")) %>%
  dplyr::mutate(compound = word(variable, 2, sep = "_")) %>%
  mutate(time_hrs = as.numeric(gsub("h", "", time))) %>%
  dplyr::mutate(idvar = paste0(enzyme, "_", time_hrs, "_", compound)) %>%
  dplyr::filter(sample != "boiled_enzyme") %>%
  dplyr::filter(!grepl("DH_24|DH_48", idvar))
  
summdat <- rawdat %>%
  dplyr::group_by(idvar) %>%
  dplyr::mutate(mean = mean(value, na.rm = T)) %>%
  dplyr::mutate(sd = sd(value, na.rm = T)) %>%
  dplyr::mutate(substrate = case_when(compound %in% c("DEET", "3MB") ~ "DEET",
                                      compound %in% c("MB", "B") ~ "MB")) %>%
  ungroup() %>%
  mutate(compound = factor(compound, levels = c("DEET", "3MB", "MB", "B")))


pal3 <- c("maroon", "red", "darkblue", "cornflowerblue")

plt2 <- ggplot2::ggplot(summdat, aes(x = time_hrs, y = mean, shape = enzyme, color = compound)) +
  geom_point() +
  geom_errorbar(aes(ymax=mean + sd, ymin = mean - sd, width=0)) + 
  geom_line() +
  ggplot2::facet_grid(cols = vars(enzyme), scales = 'free') + 
  ggpubr::theme_pubr() +
  theme(strip.text.x = element_blank()) +
  ggplot2::labs(x = "Time (h)", y = "Concentration (µM)", 
                shape = "Enzyme", color = "Compound") +
  ggplot2::scale_color_manual(values = pal3) 
plt2
ggsave( "output/Figure4.png", plt2, width = 6.5, height = 2.5)
