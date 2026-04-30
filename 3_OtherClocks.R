# ============================================================================ #
# SystemsAge and Causal Epigenetic Clock Prediction
# ============================================================================ #

# This script calculates additional DNA methylation age measures that are not
# produced by the classic-clock or PC-clock scripts. It uses the complete
# processed beta matrix from 2_prepare_complete_betas_for_clocks.qmd, transposes
# it to samples x CpGs, and applies CpG-imputation resources required by the
# SystemsAge and methylCIPHER clock implementations.
#
# Calculated outputs:
#   - SystemsAge outputs from `calcSystemsAge_cstm()`, including the renamed
#     chronological-age prediction column `SystemsChronAge_pred`
#   - AdaptAge from methylCIPHER::calcAdaptAge()
#   - DamAge from methylCIPHER::calcDamAge()
#   - CausAge from calcCausAge()
#
# Inputs:
#   - data/betas_complete_processed.RDS
#   - data/imputeMissingCpGs_SystemsClock.RDS
#   - data/SystemsAge_data.RData
#   - ~/public/SiyuanLuo/projects/LongAge/R/SystemsAge.R
#
# Output:
#   - data/OtherClocks.csv: SystemsAge, SystemsChronAge_pred, AdaptAge, DamAge,
#     and CausAge predictions for each sample.

# ============================================================================ #
# File Paths
# ============================================================================ #

# Input files
FILE_BETAS <- "data/betas_complete_processed.RDS"
FILE_METADATA <- "data/meta_3timepoints.csv"
FILE_IMPUTE_CPGS <- "data/imputeMissingCpGs_SystemsClock.RDS"
FILE_SystemsAge_RDATA <- "data/SystemsAge_data.RData"
# Source files
SYSTEMS_AGE_SOURCE <- "~/public/SiyuanLuo/projects/LongAge/R/SystemsAge.R"
# Ouput files
FILE_OUTPUT_AGES <- "data/OtherClocks.csv"
# ============================================================================ #
# Libraries
# ============================================================================ #

library(dplyr)
library(tidyr)
library(readxl)
library(methylCIPHER)

# ============================================================================ #
# Source External Scripts
# ============================================================================ #

source(SYSTEMS_AGE_SOURCE)

# ============================================================================ #
# Load Data
# ============================================================================ #
betas <- readRDS(FILE_BETAS)
datMeth <- t(betas)
# Load CpG imputation data
imputeMissingCpGs <- readRDS(file = FILE_IMPUTE_CPGS)

# ============================================================================ #
# Calculate Systems Age and Other Epigenetic Ages
# ============================================================================ #

# Systems Age
SystemsAge <- calcSystemsAge_cstm(datMeth, CpGImputation = imputeMissingCpGs, rdata_path=FILE_SystemsAge_RDATA)
SystemsAge <- SystemsAge %>% rename(SystemsChronAge_pred = Age_prediction)
SystemsAge$Sample <- rownames(SystemsAge)

df_age <- SystemsAge

# Calculate CausalAge ages
df_age$AdaptAge <- methylCIPHER::calcAdaptAge(datMeth, CpGImputation = imputeMissingCpGs)
df_age$DamAge <- methylCIPHER::calcDamAge(datMeth, CpGImputation = imputeMissingCpGs)
df_age$CausAge <- calcCausAge(datMeth, CpGImputation = imputeMissingCpGs)

# ============================================================================
# Save Results
# ============================================================================
write.table(df_age, 
            file = FILE_OUTPUT_AGES, 
            sep = ",", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)
