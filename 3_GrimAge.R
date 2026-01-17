# ==============================================================================
# DNAm GrimAge and GrimAge v2 Calculation
# ==============================================================================
# 
# This script calculates DNA methylation-based GrimAge and GrimAge v2 
# epigenetic clocks along with their constituent protein surrogates
# 
# Input:  
#   - betas_complete_processed.RDS: Processed beta values from methylation array
#   - metadata/merged_3timepoints_meta.csv: Sample metadata with age and sex
#   - DNAmGrimAgeGitHub/input/: Reference files for GrimAge calculation
# 
# Output:
#   - DNAmGrimAgeOutput.csv: GrimAge and protein components
#   - myDNAmGrimAge2.csv: GrimAge v2 and protein components
# ==============================================================================

# Load required libraries
library(readr)
library(dplyr)

# Set options
options(stringsAsFactors = FALSE)

# ==============================================================================
# File Paths Configuration
# ==============================================================================

# Input files
PATH_BETA_CORRECTED <- "./data/betas_complete_processed.RDS"
PATH_METADATA <- "data/meta_3timepoints.csv"

# GrimAge reference files
PATH_GRIMAGE_V1_COEFS <- "../DNAmGrimAgeGitHub/input/ElasticNet_DNAmProtein_Vars_model4.csv"
PATH_GRIMAGE_V2_MODEL <- "../DNAmGrimAgeGitHub/input/DNAmGrimAge2_final.Rds"
PATH_GOLD_STANDARD <- "../DNAmGrimAgeGitHub/input/datMiniAnnotation3_Gold.csv"

# Intermediate and output files
PATH_INTERMEDIATE <- "./data/DNAmGrimAgeInput.csv"
PATH_OUTPUT_V1 <- "./data/DNAmGrimAgeOutput.csv"
PATH_OUTPUT_V2 <- "./data/DNAmGrimAgev2Output.csv"

# ==============================================================================
# Helper Functions
# ==============================================================================

#' Prepare methylation data with imputation for missing probes
#'
#' @param beta_matrix Matrix of beta values (probes x samples)
#' @param required_probes Character vector of required probe names
#' @param gold_standard_path Path to gold standard imputation values
#' @return Data frame with complete probe set (samples x probes)
prepare_methylation_data <- function(beta_matrix, required_probes, gold_standard_path) {
  
  # Identify common and missing probes
  common_probes <- intersect(rownames(beta_matrix), required_probes)
  message(sprintf("Found %d/%d required probes (%.1f%%)", 
                  length(common_probes), length(required_probes),
                  100 * length(common_probes) / length(required_probes)))
  
  # Extract available data
  mydata <- beta_matrix[common_probes, , drop = FALSE]
  
  # Impute missing probes using gold standard
  imputed_probes <- setdiff(required_probes, common_probes)
  if (length(imputed_probes) > 0) {
    message(sprintf("Imputing %d missing probes", length(imputed_probes)))
    
    datagold <- read_csv(gold_standard_path, show_col_types = FALSE)
    imputed_values <- datagold$gold[datagold$CpG %in% imputed_probes]
    
    imputed_mat <- matrix(imputed_values, 
                         nrow = length(imputed_probes), 
                         ncol = ncol(mydata))
    rownames(imputed_mat) <- imputed_probes
    colnames(imputed_mat) <- colnames(mydata)
    
    mydata <- rbind(mydata, imputed_mat)
  }
  
  # Transpose to samples x probes format and convert to data frame
  mydata <- as.data.frame(t(mydata))
  mydata$SampleID <- rownames(mydata)
  
  return(mydata)
}

#' Calculate DNAm protein surrogates
#'
#' @param dat_meth Data frame with methylation values
#' @param cpgs Data frame with CpG coefficients
#' @return Data frame with protein predictions added
calculate_protein_surrogates <- function(dat_meth, cpgs) {
  
  dat_meth$Intercept <- 1
  
  proteins <- unique(cpgs$Y.pred)
  message(sprintf("Calculating %d protein surrogates", length(proteins)))
  
  for (protein in proteins) {
    cpg_subset <- cpgs[cpgs$Y.pred == protein, ]
    probe_values <- as.matrix(dat_meth[, cpg_subset$var])
    dat_meth[[protein]] <- as.numeric(probe_values %*% cpg_subset$beta)
  }
  
  return(dat_meth)
}

#' Rename protein columns to standardized format
#'
#' @param data Data frame with protein columns
#' @param old_names Character vector of old column names
#' @param new_names Character vector of new column names
#' @return Data frame with renamed columns
rename_proteins <- function(data, old_names, new_names) {
  for (i in seq_along(old_names)) {
    col_idx <- which(names(data) == old_names[i])
    if (length(col_idx) > 0) {
      names(data)[col_idx] <- new_names[i]
    }
  }
  return(data)
}

# ==============================================================================
# GrimAge (Original Version)
# ==============================================================================

message("\n=== Calculating GrimAge (Original) ===\n")

# Load processed beta values
beta_corrected <- readRDS(PATH_BETA_CORRECTED)

# Load CpG coefficients for GrimAge
cpgs_v1 <- read.csv(PATH_GRIMAGE_V1_COEFS)
required_probes <- setdiff(unique(cpgs_v1$var), c("Intercept", "Age"))

# Prepare methylation data with imputation
meth_data <- prepare_methylation_data(
  beta_corrected, 
  required_probes,
  PATH_GOLD_STANDARD
)

# Merge with metadata
meta <- read.csv(PATH_METADATA)
meth_data <- merge(meth_data, 
                   meta[, c("Sample", "sex", "age")], 
                   by.x = "SampleID", 
                   by.y = "Sample")
meth_data <- meth_data %>% rename(Age = age, Female = sex)

# Save intermediate file
write_csv(meth_data, PATH_INTERMEDIATE)

# Calculate protein surrogates
dat_meth <- calculate_protein_surrogates(meth_data, cpgs_v1)

# Calculate raw GrimAge (COX variable)
proteins_v1 <- unique(cpgs_v1$Y.pred)
output_v1 <- dat_meth[, c('SampleID', 'Age', 'Female', proteins_v1)]

# Cox regression coefficients from original GrimAge paper
output_v1$COX <- 
  output_v1$DNAmGDF_15 * 0.000348777412272004 +
  output_v1$DNAmB2M * 4.59105969389204e-07 +
  output_v1$DNAmCystatin_C * 3.49816671441537e-06 +
  output_v1$DNAmTIMP_1 * 0.000143661105491888 +
  output_v1$DNAmadm * 0.00790270975255529 +
  output_v1$DNAmpai_1 * 2.55560382039825e-05 +
  output_v1$DNAmleptin * -7.32066983502079e-06 +
  output_v1$DNAmPACKYRS * 0.0303981613409142 +
  output_v1$Age * 0.0300823182194075 +
  output_v1$Female * -0.228468475622039

# Scale COX to age units (years)
mean_age <- 59.63951
sd_age <- 9.049608
mean_cox <- 13.20127
sd_cox <- 1.086805

cox_standardized <- (output_v1$COX - mean_cox) / sd_cox
output_v1$DNAmGrimAge <- (cox_standardized * sd_age) + mean_age
output_v1$COX <- NULL

# Calculate age acceleration (residuals)
output_v1$AgeAccelGrim <- residuals(
  lm(DNAmGrimAge ~ Age, data = output_v1, na.action = na.exclude)
)

# Standardize protein names
output_v1 <- rename_proteins(
  output_v1,
  old_names = c('DNAmadm', 'DNAmCystatin_C', 'DNAmGDF_15', 
                'DNAmleptin', 'DNAmpai_1', 'DNAmTIMP_1'),
  new_names = c('DNAmADM', 'DNAmCystatinC', 'DNAmGDF15', 
                'DNAmLeptin', 'DNAmPAI1', 'DNAmTIMP1')
)

# Save GrimAge results
write.csv(output_v1, PATH_OUTPUT_V1, row.names = FALSE, quote = FALSE)
message(sprintf("GrimAge results saved to %s", PATH_OUTPUT_V1))


# ==============================================================================
# GrimAge v2 (Updated Version)
# ==============================================================================

message("\n=== Calculating GrimAge v2 ===\n")

# Load GrimAge v2 model
grimage2_model <- readRDS(PATH_GRIMAGE_V2_MODEL)
cpgs_v2 <- grimage2_model[[1]]
glmnet_coefs <- grimage2_model[[2]]
gold_standard <- grimage2_model[[3]]

# Helper function to scale COX to age units
scale_cox_to_age <- function(data, cox_col, output_col, gold_standard) {
  cox_params <- gold_standard[gold_standard$var == 'COX', ]
  age_params <- gold_standard[gold_standard$var == 'Age', ]
  
  cox_standardized <- (data[[cox_col]] - cox_params$mean) / cox_params$sd
  data[[output_col]] <- (cox_standardized * age_params$sd) + age_params$mean
  
  return(data)
}

# Reload input data (already prepared)
dat_meth_v2 <- read_csv(PATH_INTERMEDIATE, show_col_types = FALSE)

# Calculate protein surrogates for GrimAge v2
dat_meth_v2 <- calculate_protein_surrogates(dat_meth_v2, cpgs_v2)

# Select relevant columns
proteins_v2 <- unique(cpgs_v2$Y.pred)
output_v2 <- dat_meth_v2[, c('SampleID', 'Age', 'Female', proteins_v2)]

# Calculate raw COX variable using glmnet coefficients
cox_predictors <- as.matrix(output_v2[, glmnet_coefs$var])
output_v2$COX <- as.numeric(cox_predictors %*% glmnet_coefs$beta)

# Scale COX to age units
output_v2 <- scale_cox_to_age(output_v2, 'COX', 'DNAmGrimAge2', gold_standard)

# Calculate age acceleration
output_v2$AgeAccelGrim2 <- residuals(
  lm(DNAmGrimAge2 ~ Age, data = output_v2, na.action = na.exclude)
)

# Remove temporary column
output_v2$COX <- NULL

# Standardize protein names
output_v2 <- rename_proteins(
  output_v2,
  old_names = c('DNAmadm', 'DNAmCystatin_C', 'DNAmGDF_15', 'DNAmleptin', 
                'DNAmpai_1', 'DNAmTIMP_1', 'DNAmlog.CRP', 'DNAmlog.A1C'),
  new_names = c('DNAmADM', 'DNAmCystatinC', 'DNAmGDF15', 'DNAmLeptin', 
                'DNAmPAI1', 'DNAmTIMP1', 'DNAmlogCRP', 'DNAmlogA1C')
)

# Save GrimAge v2 results
write.csv(output_v2, PATH_OUTPUT_V2, row.names = FALSE, quote = FALSE)
message(sprintf("GrimAge v2 results saved to %s", PATH_OUTPUT_V2))

message("\n=== GrimAge calculations completed ===\n")


