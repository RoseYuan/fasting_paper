# ============================================================================
# Title: Epigenetic PC clocks
# ============================================================================

# FILE PATHS - Update these at the front for easy modifications
# ============================================================================
clocks_dir <- "~/public/SiyuanLuo/projects/fasting_clinics/PC-Clocks/"
utils_file <- "utils.R"
betas_file <- "./data/betas_complete_processed.RDS"
metadata_file <- "data/meta_3timepoints.csv"
output_file <- "data/PCClocks.csv"

# ============================================================================
# Load Libraries
# ============================================================================
library(dplyr)
library(readxl)
library(tidyr)
# ============================================================================
# Load Custom Functions and PC Clocks Scripts
# ============================================================================
source(utils_file)
source(paste(clocks_dir, "run_calcPCClocks.R", sep = ""))
source(paste(clocks_dir, "run_calcPCClocks_Accel.R", sep = ""))

# ============================================================================
# Load Data
# ============================================================================
betas <- readRDS(betas_file)

# ============================================================================
# Prepare the Data
# ============================================================================
datMeth <- t(betas)

meta <- read.table(metadata_file, header = TRUE, sep = ",")
datPheno <- meta[,c("age","Gender","Sample")] %>% drop_na()
names(datPheno)[names(datPheno) == "age"] <- "Age"
datPheno$Female <- as.numeric(datPheno$Gender == "Female")
datMeth <- datMeth[datPheno$Sample, ]

# ============================================================================
# Calculate the Clocks
# ============================================================================
PCClock_DNAmAge <- calcPCClocks(path_to_PCClocks_directory = clocks_dir, 
                                 datMeth = datMeth, datPheno = datPheno)
PCClock_DNAmAge <- calcPCClocks_Accel(PCClock_DNAmAge)

# ============================================================================
# Save Results
# ============================================================================
write.table(PCClock_DNAmAge, 
            file = output_file, 
            sep = ",", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)
