library(dplyr)
library(ggplot2)

getwd()
setwd("../../Downloads/Computational Results and codes-20251102T212810Z-1-001/Computational Results and codes")

#For RNA RMSD
folder_path <- "./1. RNA_RMSD/"
files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)
merged_df <- do.call(dplyr::bind_rows, lapply(files, readr::read_csv))


merged_df %>%
  select(file, rmsd_A) %>%
  mutate(pair = case_when(
    grepl("at4g03250", file) ~ "at4g03250_ath-miR166e",
    grepl("at2g27550", file) ~ "atc_ath_miR398a",
    grepl("med21at4g04780", file) ~ "med21_ath_miR161.2",
    grepl("pranav_new", file) ~ "negative_control",
    grepl("at5g20230bcbp_ath_mir398a", file) ~ "BCBP_noncanonical",
    grepl("at2g34710_ath_mir166e", file) ~ "positive_control",
    TRUE ~ "unknown"
  )) %>%
  mutate(pair = factor(pair,
                       levels = c("positive_control",
                                  "at4g03250_ath-miR166e",
                                  "med21_ath_miR161.2",
                                  "negative_control",
                                  "BCBP_noncanonical",
                                  "atc_ath_miR398a"))) %>%
  filter(!grepl("fold_2024_12_02_ago10_at5g20230bcbp_ath_mir398a_holygrail_trial_6_model_[0-4]", file)) %>%
  filter(rmsd_A > 0.001, pair != "unknown") %>%
  group_by(pair) %>%
  summarise(
    mean_rmsd = mean(rmsd_A),
    sd_rmsd = sd(rmsd_A),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = pair, y = mean_rmsd, fill = pair)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_rmsd - sd_rmsd,
                    ymax = mean_rmsd + sd_rmsd),
                width = 0.2, linewidth = 0.8) +
  #geom_text(aes(label = paste0("n = ", n), y = mean_rmsd + sd_rmsd + 0.4), size = 6) +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 1)) +
  scale_fill_brewer(palette = "Blues") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(x = "Complex Pair",
    y = "Mean RMSD (Å)",
    title = "Mean RMSD ± SD ~ RNA duplex")

#For Protein RMSD
folder_path_protein <- "./2. Domain_RMSD/"
files_protein <- list.files(path = folder_path_protein, pattern = "\\.csv$", full.names = TRUE)
merged_df_protein <- do.call(dplyr::bind_rows, lapply(files_protein, readr::read_csv))

merged_df_protein %>%
  select(af_file, domain, rmsd) %>% 
  mutate(pair = case_when(
    grepl("at4g03250", af_file) ~ "at4g03250_ath-miR166e",
    grepl("at2g27550", af_file) ~ "atc_ath_miR398a",
    grepl("med21at4g04780", af_file) ~ "med21_ath_miR161.2",
    grepl("pranav_new", af_file) ~ "negative_control",
    grepl("at5g20230bcbp_ath_mir398a", af_file) ~ "BCBP_noncanonical",
    grepl("at2g34710_ath_mir166e", af_file) ~ "positive_control",
    TRUE ~ "unknown"
  )) %>%
  mutate(pair = factor(pair,
                       levels = c("positive_control",
                                  "at4g03250_ath-miR166e",
                                  "med21_ath_miR161.2",
                                  "negative_control",
                                  "BCBP_noncanonical",
                                  "atc_ath_miR398a"))) %>%
  filter(!grepl("fold_2024_12_02_ago10_at5g20230bcbp_ath_mir398a_holygrail_trial_6_model_[0-4]", af_file)) %>%
  filter(!domain %in% c("PAZ", "Piwi_N_1")) %>%
  filter(pair != "unknown") %>%
  group_by(pair, domain) %>%
  summarise(
    mean_rmsd = mean(rmsd),
    sd_rmsd = sd(rmsd),
    n = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = pair, y = mean_rmsd, fill = domain)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_rmsd - sd_rmsd,
                    ymax = mean_rmsd + sd_rmsd),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8) +
  #geom_text(aes(label = paste0("n = ", n), y = mean_rmsd + sd_rmsd + 0.3), position = position_dodge(width = 0.8), size = 5) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  theme_classic(base_size = 20) +
  theme(legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(
    x = "Complex Pair",
    y = "Mean RMSD (Å)",
    fill = "AGO10 Domains",
    title = "Mean RMSD ± SD ~ AGO10 Protein")

###For COM analysis
folder_path_com <- "./3. RNA(COM)_DOM(COM)/"
files_com <- list.files(path = folder_path_com, pattern = "\\.csv$", full.names = TRUE)
merged_df_com <- do.call(dplyr::bind_rows, lapply(files_com, readr::read_csv))

merged_df_com %>%
  select(file, com_distance_A, domain_label) %>% 
  mutate(pair = case_when(
    grepl("at4g03250", file) ~ "at4g03250_ath-miR166e",
    grepl("at2g27550", file) ~ "atc_ath_miR398a",
    grepl("med21at4g04780", file) ~ "med21_ath_miR161.2",
    grepl("pranav_new", file) ~ "negative_control",
    grepl("at5g20230bcbp_ath_mir398a", file) ~ "BCBP_noncanonical",
    grepl("at2g34710_ath_mir166e", file) ~ "positive_control",
    TRUE ~ "unknown"
  )) %>%
  mutate(pair = factor(pair,
                       levels = c("positive_control",
                                  "at4g03250_ath-miR166e",
                                  "med21_ath_miR161.2",
                                  "negative_control",
                                  "BCBP_noncanonical",
                                  "atc_ath_miR398a"))) %>%
  filter(!grepl("fold_2024_12_02_ago10_at5g20230bcbp_ath_mir398a_holygrail_trial_6_model_[0-4]", file)) %>%
  filter(!domain_label %in% c("PAZ", "Piwi_N_1")) %>%
  filter(pair != "unknown") %>%
  group_by(pair, domain_label) %>%
  summarise(
    mean_com_distanceA = mean(com_distance_A),
    sd_com_distanceA = sd(com_distance_A),
    n = n(),
    .groups = "drop"
  ) %>%
    ggplot(aes(x = pair, y = mean_com_distanceA, fill = domain_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_com_distanceA - sd_com_distanceA,
                    ymax = mean_com_distanceA + sd_com_distanceA),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.8) +
  #geom_text(aes(label = paste0("n = ", n), y = mean_rmsd + sd_rmsd + 0.3), position = position_dodge(width = 0.8), size = 5) +
  scale_y_continuous(limits = c(0, 27), breaks = seq(0, 27, by = 3)) +
  theme_classic(base_size = 20) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  labs(
    x = "Complex Pair",
    y = "Mean Distance (Å)",
    fill = "AGO10 Domains",
    title = "Distance btw AGO10 domains and active target sites")
