library(dplyr)
library(readr)
library(stringr)

# Read in the metadata file SraRunTable.txt
SRA = read_csv("/path/to/SraRunTable.txt")
# Set the path to the folder where the downloaded FASTQ files are stored
fastq_folder = "/path/to/Sharma/FASTQ/folder/"

# Add the sampleType info required by nextNEOpi
SRA <- SRA %>%
  mutate(
	sampleType = ifelse(grepl("N_", `Sample Name`), "normal", "tumor") %>%
  	paste0("_", str_extract(`Sample Name`, "RNA|DNA"))
  )

# Remove normal_RNA entries if not removed earlier
SRA <- SRA %>%
  filter(sampleType != "normal_RNA")


# create a basic samplesheet
Sharma_samplesheet <- data.frame(sampleName = SRA$`Sample Name`,
                             	reads1 = paste(fastq_folder, SRA$Run, "_1.fastq.gz", sep = ""),
                             	reads2 = paste(fastq_folder, SRA$Run, "_2.fastq.gz", sep = ""),
                             	sampleType = SRA$sampleType,
                             	HLAfile = "",
                             	sex = SRA$sex)

Sharma_samplesheet_cancer = subset(Sharma_samplesheet, Sharma_samplesheet$sampleType != "normal_DNA")
Sharma_samplesheet_cancer$sampleName = sapply(strsplit(Sharma_samplesheet_cancer$sampleName, "_"), `[`, 1)

Sharma_samplesheet_normal = subset(Sharma_samplesheet, Sharma_samplesheet$sampleType == "normal_DNA")
Sharma_samplesheet_normal$sampleName = sapply(strsplit(Sharma_samplesheet_normal$sampleName, "-"), `[`, 1)

normal_freq = subset(Sharma_samplesheet, Sharma_samplesheet$sampleType != "normal_DNA")
normal_freq$sampleName = sapply(strsplit(normal_freq$sampleName, "-"), `[`, 1)
normal_freq = as.data.frame(table(normal_freq$sampleName))
normal_freq$Freq = normal_freq$Freq/2

Sharma_samplesheet_normal$numericSampleName <- as.numeric(sub("P", "", Sharma_samplesheet_normal$sampleName))
Sharma_samplesheet_normal <- Sharma_samplesheet_normal[order(Sharma_samplesheet_normal$numericSampleName), ]
Sharma_samplesheet_normal$numericSampleName <- NULL

Sharma_samplesheet_normal$Freq = normal_freq$Freq

expanded_df <- do.call(rbind, lapply(1:nrow(Sharma_samplesheet_normal), function(i) {
  replicate_df <- Sharma_samplesheet_normal[i, ]
  n_replicates <- Sharma_samplesheet_normal$Freq[i]
  replicate_df <- replicate_df[rep(1, n_replicates), ]
  replicate_df$sampleName <- paste0(Sharma_samplesheet_normal$sampleName[i], "-R", 1:n_replicates)
  return(replicate_df)
}))

rownames(expanded_df) <- NULL
expanded_df$Freq <- NULL

Sharma_samplesheet = rbind(Sharma_samplesheet_cancer, expanded_df)

# Write the samplesheet
write_csv(Sharma_samplesheet, "/path/to/out/Sharma_nextneopi_samplesheet.csv")