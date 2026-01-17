# ============================================================================ #
# Merge Different Epigenetic Clocks into One File
# ============================================================================ #
# File Paths
# ============================================================================ #
FILE_CASSIC_CLOCKS <- "data/ClassicClocks.csv"
FILE_GRIMAGE <- "data/DNAmGrimAgeOutput.csv"
FILE_GRIMAGEv2 <- "data/DNAmGrimAgev2Output.csv"
FILE_PC_CLOCKS <- "data/PCClocks.csv"
FILE_OTHER_CLOCKS <- "data/OtherClocks.csv"
FILE_OUTPUT <- "data/AllClocks_Merged.csv"
# ============================================================================ #
# Libraries
# ============================================================================ #
library(dplyr)
source("utils.R")
# ============================================================================ #
# Load Data
# ============================================================================ #
df_classic <- read.csv(FILE_CASSIC_CLOCKS)
df_grimage <- read.csv(FILE_GRIMAGE)
df_grimagev2 <- read.csv(FILE_GRIMAGEv2)
df_pcclocks <- read.csv(FILE_PC_CLOCKS)
df_other <- read.csv(FILE_OTHER_CLOCKS)
# ============================================================================ #
# For GrimAge and GrimAge v2 components, rename columns to avoid duplication
# ============================================================================ #

grimage_components <- c("DNAmGDF15", "DNAmB2M", "DNAmCystatinC", "DNAmTIMP1",
                        "DNAmADM", "DNAmPAI1", "DNAmLeptin", "DNAmPACKYRS")

df_grimage <- df_grimage %>%
  rename_with(~ paste0("GrimAge", .), all_of(intersect(names(.), grimage_components))) %>%
  rename(Sample = SampleID)

grimagev2_components <- c(grimage_components, "DNAmlogCRP", "DNAmlogA1C")
df_grimagev2 <- df_grimagev2 %>%
  rename_with(~ paste0("GrimAge2", .), all_of(intersect(names(.), grimagev2_components))) %>%
  rename(Sample = SampleID)
# ============================================================================ #
# Merge Clocks
# ============================================================================ #
df_merged <- df_classic %>%
    inner_join(df_grimage, by = "Sample") %>%
    inner_join(df_grimagev2, by = "Sample") %>%
    inner_join(df_pcclocks, by = "Sample") %>%
    inner_join(df_other, by = "Sample")

df_merged_clean <- collapse_duplicates(df_merged)
# ============================================================================ #
# Save Merged Results
# ============================================================================ #
write.table(df_merged_clean,
            file = FILE_OUTPUT,
            sep = ",",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
