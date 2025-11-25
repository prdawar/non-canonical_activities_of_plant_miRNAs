library(dplyr)
library(tidyr)
library(tidyverse)
library(DESeq2)

getwd()
setwd("set your working directory")
abundance_files <- list.files(pattern = "\\.txt$")
abundance_files <- abundance_files[abundance_files != "metadata.txt"]
# Read and store each file into a list of data frames
RNA_input <- lapply(abundance_files, function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_id <- gsub("_abundance\\.txt|\\.txt", "", file)
  # Rename non-Locus columns
  colnames(df)[colnames(df) != "Locus"] <- paste0(colnames(df)[colnames(df) != "Locus"], "_", sample_id)
  return(df)
})
# Merge all data frames by row names
RNA_input <- Reduce(function(x, y) merge(x, y, by = "Locus", all = TRUE), RNA_input)
# Set row names and remove the now-extra row name column
rownames(RNA_input) <- RNA_input$Locus
RNA_input <- RNA_input %>%
  select(- "Locus")
colnames(RNA_input)

meta1 <- RNA_input %>%
  rownames_to_column(var = "Locus") %>%
  pivot_longer(!Locus, names_to = "SampleID", values_to = "reads") %>%
  select(SampleID) %>%
  unique() %>%
  mutate(SampleID2 = str_extract(SampleID, "^[^_]+")) %>%
  mutate(study = str_extract(SampleID, "(?<=_).*")) %>%
  select(SampleID2, study)

metadata <- read.delim("metadata.txt")
metadata$SampleID2 <- metadata$SRR_run

meta3 <- metadata %>% 
  inner_join(y = meta1, by = "SampleID2") %>%
  select(-study)

RNA_input_long <- RNA_input %>%
  rownames_to_column(var = "Locus") %>%
  pivot_longer(!Locus, names_to = "SampleID", values_to = "reads") %>%
  mutate(SampleID2 = str_extract(SampleID, "^[^_]+")) %>%
  filter(SampleID2 %in% meta3$SampleID2) %>%
  filter(reads >= 5)

RNA_input_wide <- RNA_input_long %>%
  select(Locus, reads, SampleID2) %>%
  pivot_wider(names_from = "SampleID2", values_from = "reads") %>%
  column_to_rownames(var = "Locus")
RNA_input_wide[is.na(RNA_input_wide)] <- 0
meta3 <- meta3 %>% column_to_rownames(var = "SRR_run")
all(colnames(RNA_input_wide) %in% rownames(meta3))
meta3 <- meta3[colnames(RNA_input_wide), ]
RNA_input_wide <- round(RNA_input_wide)
Deseq_RNA <- DESeqDataSetFromMatrix(
  countData = RNA_input_wide,
  colData = meta3,
  design = ~Genotype + tissue
)
Deseq_RNA
Deseq_RNA$Genotype <- relevel(Deseq_RNA$Genotype, ref = "Col_0")
Deseq_RNA <- DESeq(Deseq_RNA)
resultsNames(Deseq_RNA)
Genotype_ago1_27_vs_Col_0 <- as.data.frame(results(Deseq_RNA, name = "Genotype_ago1_27_vs_Col_0")) %>% rownames_to_column(var = "Locus")
Genotype_dcl1_5_vs_Col_0 <- as.data.frame(results(Deseq_RNA, name = "Genotype_dcl1_5_vs_Col_0")) %>% rownames_to_column(var = "Locus")
Genotype_dcl1_9_vs_Col_0 <- as.data.frame(results(Deseq_RNA, name = "Genotype_dcl1_9_vs_Col_0")) %>% rownames_to_column(var = "Locus")
Genotype_se_vs_Col_0 <- as.data.frame(results(Deseq_RNA, name = "Genotype_se_vs_Col_0")) %>% rownames_to_column(var = "Locus")
Genotype_se_1_vs_Col_0 <- as.data.frame(results(Deseq_RNA, name = "Genotype_se_1_vs_Col_0")) %>% rownames_to_column(var = "Locus")

rename_with_suffix <- function(df, suffix) {
  df %>%
    dplyr::rename_with(~ ifelse(.x != "Locus", paste0(.x, "_", suffix), .x))
}

merged_list_all <- list(
  ago1_27 = rename_with_suffix(Genotype_ago1_27_vs_Col_0, "ago1_27"),
  dcl1_5  = rename_with_suffix(Genotype_dcl1_5_vs_Col_0, "dcl1_5"),
  dcl1_9  = rename_with_suffix(Genotype_dcl1_9_vs_Col_0, "dcl1_9"),
  se      = rename_with_suffix(Genotype_se_vs_Col_0, "se"),
  se_1    = rename_with_suffix(Genotype_se_1_vs_Col_0, "se_1")
)
merged_list_all <- Reduce(function(x, y) merge(x, y, by = "Locus", all = TRUE), merged_list_all) %>% select(Locus, matches("log2|padj|base"))
#write.csv(merged_list_all, file = "merged_RNAseq_results_all.csv")

Genotype_ago1_27_vs_Col_0_sig <- Genotype_ago1_27_vs_Col_0 %>%
  filter(padj < 0.05)
Genotype_dcl1_5_vs_Col_0_sig <- Genotype_dcl1_5_vs_Col_0 %>%
  filter(padj < 0.05)
Genotype_dcl1_9_vs_Col_0_sig <- Genotype_dcl1_9_vs_Col_0 %>%
  filter(padj < 0.05)
Genotype_se_vs_Col_0_sig <- Genotype_se_vs_Col_0 %>%
  filter(padj < 0.05)
Genotype_se_1_vs_Col_0_sig <- Genotype_se_1_vs_Col_0 %>%
  filter(padj < 0.05)

rename_with_suffix <- function(df, suffix) {
  df %>%
    dplyr::rename_with(~ ifelse(.x != "Locus", paste0(.x, "_", suffix), .x))
}

sig_list <- list(
  ago1_27 = rename_with_suffix(Genotype_ago1_27_vs_Col_0_sig, "ago1_27"),
  dcl1_5  = rename_with_suffix(Genotype_dcl1_5_vs_Col_0_sig, "dcl1_5"),
  dcl1_9  = rename_with_suffix(Genotype_dcl1_9_vs_Col_0_sig, "dcl1_9"),
  se      = rename_with_suffix(Genotype_se_vs_Col_0_sig, "se"),
  se_1    = rename_with_suffix(Genotype_se_1_vs_Col_0_sig, "se_1")
)

merged_list <- Reduce(function(x, y) merge(x, y, by = "Locus", all = TRUE), sig_list) %>% 
  select(Locus, matches("log2|padj|base")) %>%
  filter(!grepl("ath-MIR844_chr2_9942051-9942435|ath-MIR400_chr1_11785836-11786137", Locus))
#write.csv(merged_list, file = "merged_RNAseq_results.csv")

known_targets <- read.delim("../../known_targets_all.txt")
novel_targets <- read.delim("../../novel_targets_all.txt")

sig_pri_miRNAs <- merged_list %>%
  filter(str_detect(Locus, "ath-MIR")) %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  mutate(miRNA = str_extract(Locus, "^[^_]+")) %>%
  mutate(miRNA = str_replace(miRNA, "^ath-MIR(\\d+)[a-z]?$", "ath-miR\\1")) %>%
  filter(!if_all(starts_with("log2FoldChange"), is.na)) %>%
  select(miRNA, matches ("log2FoldChange")) %>% 
  pivot_longer(!miRNA, names_to = "comparison", values_to = "log2FoldChange") %>%
  filter(!is.na(log2FoldChange)) %>%
  select(miRNA, comparison) %>% unique()

sig_pri_miRNAs_df <- merged_list %>%
  filter(str_detect(Locus, "ath-MIR")) %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  mutate(miRNA = str_extract(Locus, "^[^_]+")) %>%
  mutate(miRNA = str_replace(miRNA, "^ath-MIR(\\d+)[a-z]?$", "ath-miR\\1")) %>%
  filter(!if_all(starts_with("log2FoldChange"), is.na)) %>%
  select(Locus, miRNA, matches ("log2FoldChange"))
#write.csv(sig_pri_miRNAs_df, file = "sig_pri_miRNAs.csv")

#Figure 3A
merged_list %>%
  filter(str_detect(Locus, "ath-MIR")) %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  mutate(miRNA = str_extract(Locus, "^[^_]+")) %>%
  mutate(miRNA = str_replace(miRNA, "^ath-MIR(\\d+)[a-z]?$", "ath-miR\\1")) %>%
  #select(Locus, matches ("log2FoldChange")) %>% 
  pivot_longer(cols = matches("log2FoldChange"), names_to = "comparison", values_to = "log2FoldChange") %>%
  inner_join(sig_pri_miRNAs, by = c("miRNA", "comparison")) %>%
  filter(!is.na(log2FoldChange)) %>%
  ggplot(aes(x = comparison, y = log2FoldChange)) +
  geom_boxplot(width = 0.75, outlier.shape = NA, staplewidth = 0.05) +
  geom_jitter(aes(color = comparison), width = 0.15, alpha = 0.6, size = 3) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black") +
  #facet_wrap(~ comparison, ncol = 5) +
  scale_y_continuous(limits = c(-30, 15), breaks = seq(-30,15, 5)) +
  ylab("log2FC") +
  xlab("") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"), legend.position = "null")
#ggsave("priMIRNA_expression.png", width = 8, height = 6, bg = "white")

#Figure 3B
merged_list %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  pivot_longer(!Locus, names_to = "comparison", values_to = "log2FoldChange") %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(Locus = str_extract(Locus, "[^.]+")) %>% 
  inner_join(y = known_targets, by = "Locus", relationship = "many-to-many") %>%
  select(comparison, log2FoldChange) %>%
  unique() %>%
  #inner_join(sig_pri_miRNAs, by = c("miRNA", "comparison")) %>%
  ggplot(aes(x = comparison, y = log2FoldChange)) +
  geom_boxplot(width = 0.75, outlier.shape = NA, staplewidth = 0.05) +
  geom_jitter(aes(color = comparison), width = 0.15, alpha = 0.6, size = 3) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black") +
  #facet_wrap(~ comparison, ncol = 5) +
  scale_y_continuous(limits = c(-12, 12), breaks = seq(-12,12,3)) +
  ylab("log2FC") +
  xlab("") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"), legend.position = "null")
#ggsave("known_target_expression.png", width = 8, height = 6, bg = "white")

known_targets_metaanalysis_sig <- merged_list %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  pivot_longer(!Locus, names_to = "comparison", values_to = "log2FoldChange") %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(Locus = str_extract(Locus, "[^.]+")) %>% 
  inner_join(y = known_targets, by = "Locus", relationship = "many-to-many") %>%
#inner_join(sig_pri_miRNAs, by = c("miRNA", "comparison")) %>%
  select(-miRNA) %>% unique() %>%
  pivot_wider(names_from = comparison, values_from = log2FoldChange)
#write.csv(known_targets_metaanalysis_sig, file = "known_targets_metaanalysis_sig.csv")

novel_targets_metaanalysis_sig <- merged_list %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  pivot_longer(!Locus, names_to = "comparison", values_to = "log2FoldChange") %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(Locus = str_extract(Locus, "[^.]+")) %>% 
  inner_join(y = novel_targets, by = "Locus", relationship = "many-to-many") %>%
#inner_join(sig_pri_miRNAs, by = c("miRNA", "comparison")) %>%
  select(-miRNA) %>% unique() %>%
  pivot_wider(names_from = comparison, values_from = log2FoldChange)
#write.csv(novel_targets_metaanalysis_sig, file = "novel_targets_metaanalysis_sig.csv")

#Figure 3C
merged_list %>%
  select(Locus, matches ("log2FoldChange")) %>%
  select(-"log2FoldChange_ago1_27") %>%
  pivot_longer(!Locus, names_to = "comparison", values_to = "log2FoldChange") %>%
  filter(!is.na(log2FoldChange)) %>%
  mutate(Locus = str_extract(Locus, "[^.]+")) %>%
  inner_join(y = novel_targets, by = "Locus", relationship = "many-to-many") %>%
  select(comparison, log2FoldChange) %>%
  unique() %>%
  #inner_join(sig_pri_miRNAs, by = c("miRNA", "comparison")) %>% view()
  ggplot(aes(x = comparison, y = log2FoldChange)) +
  geom_boxplot(width = 0.75, outlier.shape = NA, staplewidth = 0.05) +
  geom_jitter(aes(color = comparison), width = 0.15, alpha = 0.6, size = 3) +
  stat_summary(fun = "mean",
               geom = "point",
               color = "black") +
  #facet_wrap(~ comparison, ncol = 5) +
  scale_y_continuous(limits = c(-14, 12), breaks = seq(-14,12,2)) +
  ylab("log2FC") +
  xlab("") +
  theme_minimal(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"), legend.position = "null")
#ggsave("novel_target_expression.png", width = 8, height = 6, bg = "white")


