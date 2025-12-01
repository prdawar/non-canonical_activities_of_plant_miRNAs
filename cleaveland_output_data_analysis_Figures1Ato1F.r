library(dplyr)
library(tidyr)
library(tidyverse)
library(patchwork)
library(gridExtra)

getwd()
setwd("set your working directory")

degradome_all <- read.delim("./degradome_results_all.txt", 
                            header = TRUE)

degradome_all_clean <- degradome_all %>%
  select(Query.and.Transcript, contains("Degradome")) %>%
  select(
    -contains("Concateted"),
    -contains("Pseudodegradome")
  ) %>%
  mutate(unique_id = make.unique(Query.and.Transcript)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>% 
  filter(signif_count > 0)
#write.csv(degradome_all_clean, file = "degradome_all_clean_11_17_25.csv")

known_targets <- read.delim("./known_miRNA_targets.txt", header = TRUE)
known_targets <- known_targets %>%
  mutate(Query.and.Transcript = paste(known_targets$Interacting.miRNA, known_targets$Target, sep = "&")) %>% 
  select(Query.and.Transcript, Target, Interacting.miRNA)

degradome_all_clean <- degradome_all_clean %>%
  mutate(Query.and.Transcript = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                                     replacement = "\\1&\\2",
                                     Query.and.Transcript))

degradome_all_clean_meta <- degradome_all %>%
  select(Query.and.Transcript, MFEratio, AllenScore, Structure, Sequence, contains("Degradome")) %>%
  select(
    -contains("Concateted"),
    -contains("Pseudodegradome")
  ) %>%
  mutate(unique_id = make.unique(Query.and.Transcript)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.1, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>%
  filter(signif_count > 0) %>%
  mutate(Query.and.Transcript = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                                     replacement = "\\1&\\2",
                                     Query.and.Transcript)) %>%
  rownames_to_column(var = "ID") %>%
  select(Query.and.Transcript, MFEratio, AllenScore, Structure, Sequence)

degradome_all_clean_meta_new <- degradome_all_clean_meta %>%
  group_by(Query.and.Transcript) %>%
  filter(MFEratio == max(MFEratio)) %>%
  filter(AllenScore == min(AllenScore)) %>%
  ungroup() %>% 
  unique()

degradome_all_clean_canonical <- degradome_all_clean %>%
  filter(apply(degradome_all_clean, 1, function(row) {
    any(sapply(known_targets$Query.and.Transcript, function(pattern) {
      grepl(pattern, paste(row, collapse = " "))
    }))
  }))
#write.csv(degradome_all_clean_meta, "degradome_all_clean_meta_11_17_25.csv")
#write.csv(degradome_all_clean_meta_new, "degradome_all_clean_meta_new_11_17_25.csv")

final_list_canonical <- degradome_all_clean_canonical %>%
  rownames_to_column(var = "ID") %>%
  select(Query.and.Transcript, signif_count) %>%
  group_by(Query.and.Transcript) %>%
  slice_max(signif_count, n = 1, with_ties = FALSE) %>% 
  ungroup() %>%
  inner_join(y = degradome_all_clean_meta_new, by = "Query.and.Transcript")

#Figure 1A
final_list_canonical %>%
  select(Query.and.Transcript, signif_count, MFEratio, AllenScore) %>%
  mutate(signif_group = if_else(signif_count >= 2, "TRUE", "FALSE"),
         mfe_group = if_else(MFEratio >= 0.73, "TRUE", "FALSE"),
         allen_group = if_else(AllenScore <= 4, "TRUE", "FALSE")) %>%
  pivot_longer(cols = c(signif_group, mfe_group, allen_group),
               names_to = "parameter", values_to = "status") %>%
  mutate(parameter = recode(parameter,
                            signif_group = "Degradome Evidence",
                            mfe_group = "MFEratio",
                            allen_group = "AllenScore")) %>%
  dplyr::count(parameter, status, name = "n") %>%
ggplot(aes(x = status, y = n, fill = parameter)) +
  geom_bar(stat = "identity", 
           position = position_dodge(width = 0.62), width = 0.6) +
  labs(x = "", y = "Count (n)") +
  scale_y_continuous(limits = c(0, 175), breaks = seq(0,175, 25)) +
  scale_fill_manual(values = c("Degradome Evidence" = "#AFEEEE",
                               "MFEratio" = "#DDA0DD",
                               "AllenScore" = "#9ACD32")) +
  theme_minimal(base_size = 20)
#ggsave("known_targets_metrics_11_17_25.png", width = 8, height = 6, bg = "white")
#write.csv(final_list_canonical, "final_list_canonical_11_17_25.csv")

degradome_all_clean_noncanonical <- degradome_all %>%
  select(Query.and.Transcript, contains("Degradome")) %>%
  select(
    -contains("Concateted"),
    -contains("Pseudodegradome")
  ) %>%
  mutate(unique_id = make.unique(Query.and.Transcript)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>%
  filter(signif_count > 0) %>%
  mutate(Query.and.Transcript = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                                     replacement = "\\1&\\2",
                                     Query.and.Transcript)) %>%
  filter(apply(., 1, function(row) {
    all(sapply(known_targets$Query.and.Transcript, function(pattern) {
      !grepl(pattern, paste(row, collapse = " "))
    }))
  }))

final_list_noncanonical <- degradome_all_clean_noncanonical %>%
  rownames_to_column(var = "ID") %>%
  select(Query.and.Transcript, signif_count) %>%
  group_by(Query.and.Transcript) %>%
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup()

final_list_noncanonical <- final_list_noncanonical %>%
  filter(signif_count > 1)

degradome_all_clean_noncanonical_meta <- degradome_all %>%
  select(Query.and.Transcript, MFEratio, AllenScore, Structure, Sequence, contains("Degradome")) %>%
  select(
    -contains("Concateted"),
    -contains("Pseudodegradome")
  ) %>%
  mutate(unique_id = make.unique(Query.and.Transcript)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.1, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>%
  filter(signif_count > 0) %>%
  mutate(Query.and.Transcript = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                                     replacement = "\\1&\\2",
                                     Query.and.Transcript)) %>%
  filter(apply(., 1, function(row) {
    all(sapply(known_targets$Query.and.Transcript, function(pattern) {
      !grepl(pattern, paste(row, collapse = " "))
    }))
  })) %>%
  rownames_to_column(var = "ID") %>%
  select(Query.and.Transcript, MFEratio, AllenScore, Structure, Sequence)

degradome_all_clean_noncanonical_meta <- degradome_all_clean_noncanonical_meta %>%
  group_by(Query.and.Transcript) %>%
  filter(MFEratio == max(MFEratio)) %>%
  filter(AllenScore == min(AllenScore)) %>%
  ungroup() %>% 
  unique()

final_list_noncanonical <- final_list_noncanonical %>% 
  inner_join(degradome_all_clean_noncanonical_meta, by = "Query.and.Transcript")
#write.csv(final_list_noncanonical, file = "final_list_noncanonical_11_17_25.csv")

degradome_all_clean_noncanonical_cc_1 <- degradome_all %>%
  select(Query.and.Transcript, contains("Degradome")) %>%
  select(
    -contains("Concateted"),
    -contains("Pseudodegradome")
  ) %>%
  mutate(unique_id = make.unique(Query.and.Transcript)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 1, .names = NULL),
    na.rm = TRUE
  )) %>%
  filter(signif_count > 0) %>%
  mutate(Query.and.Transcript = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                                     replacement = "\\1&\\2",
                                     Query.and.Transcript)) %>%
  filter(apply(., 1, function(row) {
    all(sapply(known_targets$Query.and.Transcript, function(pattern) {
      !grepl(pattern, paste(row, collapse = " "))
    }))
  }))

final_list_noncanonical_cc_1 <- degradome_all_clean_noncanonical_cc_1 %>%
  rownames_to_column(var = "ID") %>%
  select(Query.and.Transcript, signif_count) %>%
  group_by(Query.and.Transcript) %>%
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup()

final_list_noncanonical_cc_1 <- final_list_noncanonical_cc_1 %>%
  filter(signif_count > 1)

final_list_noncanonical_cc_1 <- final_list_noncanonical_cc_1 %>% 
  inner_join(degradome_all_clean_noncanonical_meta, by = "Query.and.Transcript")

degradome_all_clean_noncanonical_cc_0 <- degradome_all %>%
  select(Query.and.Transcript, contains("Degradome")) %>%
  select(
    -contains("Concateted"),
    -contains("Pseudodegradome")
  ) %>%
  mutate(unique_id = make.unique(Query.and.Transcript)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x < 1, .names = NULL),
    na.rm = TRUE
  )) %>%
  filter(signif_count > 0) %>%
  mutate(Query.and.Transcript = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                                     replacement = "\\1&\\2",
                                     Query.and.Transcript)) %>%
  filter(apply(., 1, function(row) {
    all(sapply(known_targets$Query.and.Transcript, function(pattern) {
      !grepl(pattern, paste(row, collapse = " "))
    }))
  }))

final_list_noncanonical_cc_0 <- degradome_all_clean_noncanonical_cc_0 %>%
  rownames_to_column(var = "ID") %>%
  select(Query.and.Transcript, signif_count) %>%
  group_by(Query.and.Transcript) %>%
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup()

final_list_noncanonical_cc_0 <- final_list_noncanonical_cc_0 %>%
  filter(signif_count > 1)

final_list_noncanonical_cc_0 <- final_list_noncanonical_cc_0 %>% 
  inner_join(degradome_all_clean_noncanonical_meta, by = "Query.and.Transcript")

#Figure 1B
final_list_canonical %>%
  mutate(nc = MFEratio < 0.73 | AllenScore > 4) %>%
  ggplot(aes(x = MFEratio, y = AllenScore)) +
  geom_point(aes(shape = nc), alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0,15,3)) +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6,1,0.1)) +
  #geom_vline(xintercept = 0.73, linetype = "dashed", color = "grey", linewidth = 0.7) +
  #geom_hline(yintercept = 4, linetype = "dashed", color = "grey", linewidth = 0.7) +
  theme_minimal(base_size = 20)
#ggsave("relationship_Allenscore_MFEratio_known_11_17_25.png", width = 8, height = 6, bg = "white")

#Figure 1C
final_list_noncanonical %>%
  filter(signif_count > 1) %>%
  mutate(nc = MFEratio < 0.73 | AllenScore > 4) %>%
  ggplot(aes(x = MFEratio, y = AllenScore)) +
  geom_point(aes(shape = nc), alpha = 0.5, size = 2)  +
  geom_smooth(method = "loess", se = TRUE, color = "red") +
  scale_y_continuous(limits = c(0, 30), breaks = seq(0,30,5)) +
  scale_x_continuous(limits = c(0.6, 1), breaks = seq(0.6,1,0.1)) +
  #geom_vline(xintercept = 0.73, linetype = "dashed", color = "darkgrey", linewidth = 0.7) +
  #geom_hline(yintercept = 4, linetype = "dashed", color = "darkgrey", linewidth = 0.7) +
  theme_minimal(base_size = 20)
#ggsave("relationship_Allenscore_MFEratio_novel_11_17_25.png", width = 8, height = 6, bg = "white")

final_list_noncanonical %>%
  filter(signif_count > 1) %>%
  select(Query.and.Transcript, signif_count, MFEratio, AllenScore) %>%
  mutate(signif_group = if_else(signif_count >= 2, "TRUE", "FALSE"),
         mfe_group = if_else(MFEratio >= 0.73, "TRUE", "FALSE"),
         allen_group = if_else(AllenScore <= 4, "TRUE", "FALSE")) %>%
  pivot_longer(cols = c(signif_group, mfe_group, allen_group),
               names_to = "parameter", values_to = "status") %>%
  mutate(parameter = recode(parameter,
                            signif_group = "Degradome Evidence",
                            mfe_group = "MFEratio",
                            allen_group = "AllenScore")) %>% 
  dplyr::count(status, parameter)

T_plots_0.6_to_0.65 <- read.delim("./T_plots_MFE_0.6_0.65.txt", 
                                  header = TRUE)
T_plots_0.6_to_0.65_meta <-  T_plots_0.6_to_0.65 %>%
  select(Criteria, MFEratio, AllenScore, Structure, Sequence) %>%
  mutate(Criteria = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                         replacement = "\\1&\\2",
                         Criteria)) %>%
  group_by(Criteria) %>%
  filter(MFEratio == max(MFEratio)) %>%
  filter(AllenScore == min(AllenScore)) %>%
  ungroup() %>% 
  unique()
#write.csv(T_plots_0.6_to_0.65_meta, "T_plots_0.6_to_0.65_meta_11_17_25.csv")

#Figure 1E
T_plots_0.6_to_0.65 %>%
  select(Criteria, contains("Degradome")) %>%
  select(
    -contains("cat_degradome"),
    -contains("FP_pseudodegradome_new")
  ) %>%
  mutate(unique_id = make.unique(Criteria)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>% 
  filter(signif_count > 1) %>% 
  mutate(Criteria = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                         replacement = "\\1&\\2",
                         Criteria)) %>%
  group_by(Criteria) %>% 
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  inner_join(T_plots_0.6_to_0.65_meta, by = "Criteria") %>%
  select(Criteria, MFEratio, AllenScore, signif_count) %>%
  #rename(Query.and.Transcript = Criteria) %>%
  filter(signif_count > 1) %>%
  mutate(bin = cut(MFEratio,
                   breaks = seq(0.60, 0.65, by = 0.01),
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = c("0.60–0.61", 
                              "0.61–0.62", 
                              "0.62–0.63", 
                              "0.63–0.64", 
                              "0.64–0.65"))) %>%
  group_by(bin) %>%
  summarise(total_count = n()) %>%
  ggplot(aes(x = bin, y = total_count)) +
  geom_col(fill = "#AFEEEE") +
  labs(x = "MFE ratio",
       y = "Count (n)") +
  scale_y_continuous(limits = c(0, 30), 
                     breaks = seq(0, 30, by = 5)) +
  theme_classic(base_size = 20)
#ggsave("MFEratio_distribution_0.6-0.65_11_17_25.png", width = 8, height = 6, bg = "white")

T_plots_0.6_to_0.65_sig <- T_plots_0.6_to_0.65 %>%
  select(Criteria, contains("Degradome"), -number_of_degradome_hits) %>%
  select(
    -contains("cat_degradome"),
    -contains("FP_pseudodegradome_new")
  ) %>%
  mutate(unique_id = make.unique(Criteria)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>% 
  filter(signif_count > 1) %>%
  mutate(Criteria = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                         replacement = "\\1&\\2",
                         Criteria)) %>%
  select(Criteria, signif_count) %>% 
  group_by(Criteria) %>%
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  inner_join(T_plots_0.6_to_0.65_meta, by = "Criteria")
#write.csv(T_plots_0.6_to_0.65_sig, "T_plots_0.6_to_0.65_sig_11_17_25.csv")

#Figure 1D
T_plots_0.6_to_0.65 %>%
  select(Criteria, contains("Degradome")) %>%
  mutate(unique_id = make.unique(Criteria)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>% 
  filter(signif_count > 1) %>% 
  mutate(Criteria = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                         replacement = "\\1&\\2",
                         Criteria)) %>%
  group_by(Criteria) %>% 
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  inner_join(T_plots_0.6_to_0.65_meta, by = "Criteria") %>%
  select(Criteria, MFEratio, AllenScore, signif_count) %>%
  #rename(Query.and.Transcript = Criteria) %>%
  bind_rows(final_list_noncanonical) %>%
  bind_rows(final_list_canonical) %>%
  filter(signif_count > 1) %>% 
  mutate(bin = cut(MFEratio,
                   breaks = seq(0.6, 1, by = 0.05),
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = c("0.60–0.65", 
                              "0.65–0.70", 
                              "0.70–0.75", 
                              "0.75–0.80", 
                              "0.80–0.85",
                              "0.85-0.90",
                              "0.90-0.95",
                              "0.95-1.00"))) %>%
  group_by(bin) %>%
  summarise(total_count = n()) %>%
  ggplot(aes(x = bin, y = total_count)) +
  geom_col(fill = "#87CEEB") +
  labs(x = "MFE ratio",
       y = "Count (n)") +
  scale_y_continuous(limits = c(0, 150), 
                     breaks = seq(0, 150, by = 30)) +
  theme_classic(base_size = 20)
#ggsave("MFEratio_distribution_all_11_17_25.png", width = 10, height = 6, bg = "white")

#Figure 1F
T_plots_0.6_to_0.65 %>%
  select(Criteria, contains("Degradome")) %>%
  mutate(unique_id = make.unique(Criteria)) %>%
  column_to_rownames("unique_id") %>%
  mutate(signif_count = rowSums(
    across(names(.)[grepl("DegradomePval", names(.))], ~ .x < 0.05, .names = NULL) &
      across(names(.)[grepl("DegradomeCategory", names(.))], ~ .x <= 2, .names = NULL),
    na.rm = TRUE
  )) %>% 
  filter(signif_count > 1) %>% 
  rownames_to_column(var = "unique_id") %>%
  mutate(Criteria = gsub(pattern = "(ath-miR\\d+)(?:\\.\\d+)?[a-z]*-?\\d*p?&([A-Z0-9]+)\\.\\d+", 
                         replacement = "\\1&\\2",
                         Criteria)) %>%
  group_by(Criteria) %>%
  slice_max(signif_count, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  inner_join(T_plots_0.6_to_0.65_meta, by = "Criteria") %>% 
  select(Criteria, MFEratio, AllenScore, signif_count) %>% 
  #rename(Query.and.Transcript = Criteria) %>% 
  bind_rows(final_list_noncanonical) %>%
  bind_rows(final_list_canonical) %>%
  filter(signif_count > 1) %>%
  mutate(bin = cut(AllenScore,
                   breaks = seq(0, 30, by = 3),
                   include.lowest = TRUE,
                   right = FALSE,
                   labels = c("0-3",
                              "3-6",
                              "6-9",
                              "9-12",
                              "12-15",
                              "15-18",
                              "18-21",
                              "21-24",
                              "24-27",
                              "27-30"))) %>%
  group_by(bin) %>%
  summarise(total_count = n()) %>%
  ggplot(aes(x = bin, y = total_count)) +
  geom_col(fill = "#FFDAB9") +
  labs(x = "Allen Score",
       y = "Count (n)") +
  scale_y_continuous(limits = c(0, 200), 
                     breaks = seq(0, 200, by = 40)) +
  theme_classic(base_size = 20)
#ggsave("Allenscore_distribution_11_17_25.png", width = 8, height = 6, bg = "white")
