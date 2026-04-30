# ============================================================================ #
# Merge Different Epigenetic Clocks into One File
# ============================================================================ #

# This script combines all epigenetic clock prediction outputs into a single
# sample-level table for downstream visualization and association analysis. It
# merges classical clocks, GrimAge v1, GrimAge v2, PC clocks, SystemsAge, and
# causal-age clocks by `Sample`.
#
# GrimAge v1 and v2 contain overlapping component names, so this script prefixes
# shared component columns with `GrimAge` or `GrimAge2` before merging. After all
# clock tables are joined, duplicate columns created by overlapping source files
# are reconciled with `collapse_duplicates()` from `utils.R`.
#
# Inputs:
#   - data/ClassicClocks.csv
#   - data/DNAmGrimAgeOutput.csv
#   - data/DNAmGrimAgev2Output.csv
#   - data/PCClocks.csv
#   - data/OtherClocks.csv
#
# Output:
#   - data/AllClocks_Merged.csv: one merged table containing all clock
#     predictions and clock components available for each shared sample.

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
