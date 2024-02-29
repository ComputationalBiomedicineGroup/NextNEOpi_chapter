library(stringr)
library(readr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(DiversitySeq)
library(reshape2)

# Set path to nextNEOpi result folder, plots will be exported to this folder as well.
fastq_folder = "/path/to/your/nextNEOpi_result_folder/"

# Read in final result table containing the neoantigens
Sharma <- list.files(paste0(fastq_folder, "neoantigens/"), 
                     pattern = '*_MHC_Class_I_all_epitopes_ccf_ref_match.tsv', recursive = TRUE, full.names = TRUE)

Samples = as.data.frame(Sharma)
Samples = str_split_fixed(Samples$Sharma, "/", n = 9)[,7]

Input_mine = data.frame(Sample = Samples, Data = Sharma)

Sharma_extract <- function(input_1) {
  tempory <-read_tsv(input_1[2], col_names=TRUE)
  tempory$Sample = input_1[1]
  tempory = tempory[!duplicated(tempory[c("MT Epitope Seq")]), ]
  return(tempory)
}

All_combined_sharma = apply(Input_mine, 1, function(x) Sharma_extract(x))
All_combined_sharma = do.call(rbind, All_combined_sharma)

All_combined_sharma$Patient = str_split_fixed(All_combined_sharma$Sample, "-", n = 2)[,1]
All_combined_sharma = subset(All_combined_sharma, All_combined_sharma$`Best MT Percentile` <= 2)

All_combined_sharma_mean <- All_combined_sharma %>%
  group_by(Patient, `MT Epitope Seq`) %>%
  mutate(Shared_among_regions = n(),
         across(where(is.numeric), ~ mean(., na.rm = TRUE))) %>%
  ungroup()


All_combined_sharma_mean = All_combined_sharma_mean[!duplicated(All_combined_sharma_mean[c("Patient", "MT Epitope Seq")]), ]

All_combined_sharma_mean = subset(All_combined_sharma_mean, All_combined_sharma_mean$`Gene Expression` >= 1)
All_combined_sharma_mean = subset(All_combined_sharma_mean, All_combined_sharma_mean$`Tumor RNA VAF` >= 0.3)

# ------------------------------------------------------------------

# Violin

# ------------------------------------------------------------------

plot_violin <- function(input_1, input_2, input_3, input_4) {
  All_combined = subset(input_1, input_1$Shared_among_regions > input_3)
  All_combined = All_combined[!duplicated(All_combined[c("MT Epitope Seq")]), ]
  ggplot(All_combined, aes(y=get(input_2), x=Patient)) + 
    geom_violin(trim = FALSE) +
    geom_jitter(aes(color = as.factor(Shared_among_regions)), size = 2, alpha = 1) +
    ggtitle(input_4) +
    scale_color_manual(values = c("#99CCFF", "orange", "black", "red")) +
    theme(panel.spacing.y = unit(2, "mm"), strip.text.x = element_text(size = 16), panel.spacing.x = unit(2, "mm")) +
    labs(x = "", y = input_2, color = "Shared regions") +
    theme_minimal() +
    theme(axis.text=element_text(size=18), axis.text.y = element_text(angle = 90, hjust = 0.5),
          axis.title = element_text(size=20), plot.title = element_text(size=20), 
          legend.title = element_text(size=15), legend.text = element_text(size=15))
}

All_combined_sharma_mean$`Median Fold Change` = log10(as.numeric(All_combined_sharma_mean$`Median Fold Change`))
All_combined_sharma_mean$`Best MT IC50 Score` = log10(All_combined_sharma_mean$`Best MT IC50 Score`)
All_combined_sharma_mean$`Gene Expression` = log10(All_combined_sharma_mean$`Gene Expression`)
All_combined_sharma_mean$CCF.05 = round(All_combined_sharma_mean$CCF.05, 1)

Gene = plot_violin(All_combined_sharma_mean, "Gene Expression", 0, "")
Score = plot_violin(All_combined_sharma_mean, "Best MT IC50 Score", 0, "")
Perc = plot_violin(All_combined_sharma_mean, "Best MT Percentile", 0, "")
FC = plot_violin(All_combined_sharma_mean, "Median Fold Change", 0, "")
CCF = plot_violin(All_combined_sharma_mean, "CCF", 0, "")
CCF05 = plot_violin(All_combined_sharma_mean, "CCF.05", 0, "")
CCF95 = plot_violin(All_combined_sharma_mean, "CCF.95", 0, "")
pC = plot_violin(All_combined_sharma_mean, "pClonal", 0, "")
pS = plot_violin(All_combined_sharma_mean, "pSubclonal", 0, "")

wrap_plots(Gene, Score, Perc, FC, CCF, CCF05, CCF95, pC, pS, ncol =	3, nrow = 3, guides = "collect") + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))
ggsave(filename = paste0(fastq_folder, "violin.jpeg"), plot = last_plot(), width = 18, height = 12 )

# ------------------------------------------------------------------

# Barplot

# ------------------------------------------------------------------

All_TSAs = as.data.frame(table(All_combined_sharma$Sample))
All_TSAs$Var1 <- as.character(All_TSAs$Var1)
All_TSAs <- rbind(All_TSAs, c("P2-R1", 0))
All_TSAs$Freq <- as.numeric(All_TSAs$Freq)
All_TSAs$Patient = str_split_fixed(All_TSAs$Var1, "-", n = 2)[,1]

ggplot(All_TSAs, aes(y=Freq, x=Var1, fill = Patient)) + 
    geom_bar(position= "dodge", stat="identity", size = 0.3) +
    geom_text(size=12, aes(label = Freq), color = "black", vjust = -0.5) +
    ylab("Neoantigens") +
    xlab("Sample") +
    ylim(0,400) +
    scale_x_discrete(labels=c("R1","R2","R3", "R1","R2","R3", "R1","R2","R3", "R4")) +
    theme(axis.text=element_text(size=30), axis.title = element_text(size=25), plot.title = element_text(size=25), 
          legend.title = element_text(size=24), legend.text = element_text(size=24))
ggsave(filename = paste0(fastq_folder, "barplot.jpeg"), plot = last_plot(), width = 18, height = 12 )

# ------------------------------------------------------------------

# FUSIONS

# ------------------------------------------------------------------

Sharma_fusions <- list.files(paste0(fastq_folder, "neoantigens/"), pattern = '*_I_filtered_ref_match.tsv', recursive = T, full.names = T)

# extract the sample name from the file name
Samples = as.data.frame(Sharma_fusions)
Samples = str_split_fixed(Samples$Sharma_fusions, "/", n = 9)[,7]

# Input dataframe with file paths and sample ID
Input_mine_fusions = data.frame(Sample = Samples, Data = Sharma_fusions)

Sharma_fusions_extract <- function(input_1) {
  tempory <-read_tsv(input_1[2], col_names=TRUE)
  tempory$Sample = input_1[1]
  tempory = tempory[!duplicated(tempory[c("Fusion_Peptide")]), ]
  return(tempory)
}

All_combined_Sharma_fusions = apply(Input_mine_fusions, 1, function(x) Sharma_fusions_extract(x))
All_combined_Sharma_fusions = do.call(rbind, All_combined_Sharma_fusions)

All_combined_Sharma_fusions$Gene1_TPM = log(as.numeric(All_combined_Sharma_fusions$Gene1_TPM)+1)
All_combined_Sharma_fusions$Gene2_TPM = as.numeric(All_combined_Sharma_fusions$Gene2_TPM)
All_combined_Sharma_fusions[is.na(All_combined_Sharma_fusions)] <- 0
All_combined_Sharma_fusions$Gene2_TPM = log(All_combined_Sharma_fusions$Gene2_TPM+1)

All_combined_Sharma_fusions$Patient = str_split_fixed(All_combined_Sharma_fusions$Sample, "-", n = 2)[,1]

All_combined_Sharma_fusions <- All_combined_Sharma_fusions %>%
  group_by(Patient, Fusion_Peptide) %>%
  mutate(
    Shared_among_regions = n()
  ) %>%
  ungroup()

All_combined_Sharma_fusions$Fusion
ggplot(All_combined_Sharma_fusions, aes(y=Gene1_TPM, x=Gene2_TPM)) + 
  geom_point(size=4, aes(color = as.factor(Shared_among_regions))) +
  ggtitle("Fusion neoantigens") +
  #geom_text(size=4, aes(label = Fusion), color = "black") +
  theme(panel.spacing.y = unit(2, "mm"), strip.text.x = element_text(size = 16), panel.spacing.x = unit(2, "mm")) +
  facet_wrap(~Patient) +
  ylab("log(TPM+1) Gene1") +
  xlab("log(TPM+1) Gene2") +
  labs(color = "Shared regions") +
  theme(axis.text=element_text(size=22), axis.title = element_text(size=25), plot.title = element_text(size=25), 
        legend.title = element_text(size=18), legend.text = element_text(size=18))
ggsave(filename = paste0(fastq_folder, "fusions.jpeg"), plot = last_plot(), width = 18, height = 12 )

# ------------------------------------------------------------------

# TMB

# ------------------------------------------------------------------

Sharma_TMB <- list.files(paste0(fastq_folder, "analyses/"), pattern = '*_mutational_burden.txt', recursive = T, full.names = T)
Sharma_TMB_coding <- list.files(paste0(fastq_folder, "analyses/"), pattern = '*_mutational_burden_coding.txt', recursive = T, full.names = T)


# extract the sample name from the file name
Samples = as.data.frame(Sharma_TMB)
Samples = str_split_fixed(Samples$Sharma_TMB, "/", n = 9)[,7]

# Input dataframe with file paths and sample ID
Input_mine_TMB = data.frame(Sample = Samples, Data_1 = Sharma_TMB, Data_coding = Sharma_TMB_coding)

Sharma_TMB_extract <- function(input_1) {
  lines <- readLines(input_1[2])
  # Process the lines into key-value pairs
  key_value_pairs <- strsplit(lines, ":\\s*")
  data_list <- lapply(key_value_pairs, function(x) setNames(as.numeric(x[2]), x[1]))
  # Combne the data into a data frame
  data_df_TMB <- do.call(rbind, data_list)
  # If you want to transpose it so keys are columns and values are rows
  data_df_t_TMB <- as.data.frame(t(data_df_TMB))
  # Set the row names to NULL if you don't want row names
  colnames(data_df_t_TMB) <- sapply(key_value_pairs, `[`, 1)
  rownames(data_df_t_TMB) <- NULL
  data_df_t_TMB$Sample = input_1[1]
  
  lines <- readLines(input_1[3])
  # Process the lines into key-value pairs
  key_value_pairs <- strsplit(lines, ":\\s*")
  data_list <- lapply(key_value_pairs, function(x) setNames(as.numeric(x[2]), x[1]))
  # Combne the data into a data frame
  data_df_coding <- do.call(rbind, data_list)
  # If you want to transpose it so keys are columns and values are rows
  data_df_t_coding <- as.data.frame(t(data_df_coding))
  # Set the row names to NULL if you don't want row names
  colnames(data_df_t_coding) <- sapply(key_value_pairs, `[`, 1)
  rownames(data_df_t_coding) <- NULL
  data_df_t_coding$Sample = input_1[1]
  
  data_df_t = merge(data_df_t_TMB, data_df_t_coding, by = "Sample")
  
  return(data_df_t)
}

All_combined_Sharma_TMB = apply(Input_mine_TMB, 1, function(x) Sharma_TMB_extract(x))
All_combined_Sharma_TMB = do.call(rbind, All_combined_Sharma_TMB)

All_combined_Sharma_TMB$Patient = str_split_fixed(All_combined_Sharma_TMB$Sample, "-", n = 2)[,1]

# ------------------------------------------------------------------

# TCR

# ------------------------------------------------------------------

Sharma_TCR <- list.files(paste0(fastq_folder, "analyses/"), pattern = '*_tumor_RNA_mixcr.clones_TRAD.tsv', recursive = T, full.names = T)

Samples = as.data.frame(Sharma_TCR)
Samples = str_split_fixed(Samples$Sharma_TCR, "/", n = 9)[,7]

# Input dataframe with file paths and sample ID
Input_mine_TCR = data.frame(Sample = Samples, TCR = Sharma_TCR)

Sharma_TCR_extract <- function(input_1) {
  mixcr_output <- read_tsv(input_1[2])
  mixcr_output <- as.data.frame(mixcr_output)
  count_matrix <- mixcr_output[, c("cloneId", "readCount")]
  row.names(count_matrix) <- count_matrix$cloneId
  count_matrix <- count_matrix[, -1, drop = FALSE]
  Shannon = aindex(count_matrix,
         index = "Shannon")
  Richness = aindex(count_matrix,
         index = "Richness")
  Pielou = aindex(count_matrix,
         index = "Pielou")
  TCR = data.frame(Sample = input_1[1],
                  Shannon = Shannon$group1,
                  Richness = Richness$group1,
                  Eveness = Pielou$group1)
  return(TCR)
}

All_combined_Sharma_TCR = apply(Input_mine_TCR, 1, function(x) Sharma_TCR_extract(x))
All_combined_Sharma_TCR = do.call(rbind, All_combined_Sharma_TCR)
All_combined_Sharma_TCR$Patient = str_split_fixed(All_combined_Sharma_TCR$Sample, "-", n = 2)[,1]

# ------------------------------------------------------------------

# TMB + TCR visualization

# ------------------------------------------------------------------

All_combined_Sharma_TCR_TMB$`Mutational load (variants/Mbps).y`

All_combined_Sharma_TCR_TMB = merge(All_combined_Sharma_TCR, All_combined_Sharma_TMB, by = c("Sample", "Patient"))

plot_extra <- function(input_1, input_2) {
  max_value <- max(All_combined_Sharma_TCR_TMB[[input_1]], na.rm = TRUE)
  upper_limit <- max_value + (max_value * 0.08) # Increase by 25%
  
  ggplot(All_combined_Sharma_TCR_TMB, aes(y=get(input_1), x=Sample, fill = Patient)) + 
    geom_bar(position= "dodge", stat="identity", size = 0.3) +
    geom_text(size=6, aes(label = round(get(input_1), 2)), color = "black", vjust = -0.5) +
    ylab(input_2) +
    xlab("") +
    ylim(0, upper_limit) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    scale_x_discrete(labels=c("R1","R2","R3", "R1","R2","R3", "R1","R2","R3", "R4")) +
    theme_minimal() +
    theme(axis.text=element_text(size=22), axis.title = element_text(size=25), plot.title = element_text(size=25), 
          legend.title = element_text(size=18), legend.text = element_text(size=18))
}

Richness = plot_extra("Richness", "RNA T TCR Richness")
Shannon = plot_extra("Shannon", "RNA T TCR Shannon")
Eveness = plot_extra("Eveness", "RNA T TCR Eveness")
TMB = plot_extra("Mutational load (variants/Mbps).x", "TMB")
TMB_clonal = plot_extra("Mutational load clonal (variants/Mbps).x", "TMB clonal")
TMB_clonal_coding = plot_extra("Mutational load (variants/Mbps).y", "TMB clonal coding")

wrap_plots(TMB, TMB_clonal, TMB_clonal_coding, Richness, Shannon, Eveness, ncol =	3, nrow = 2, guides = "collect") + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 25))
ggsave(filename = paste0(fastq_folder, "patient_features.jpeg"), plot = last_plot(), width = 18, height = 12 )