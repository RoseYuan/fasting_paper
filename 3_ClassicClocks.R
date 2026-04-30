# ============================================================================ #
# Classic Epigenetic Clock Prediction
# ============================================================================ #

# This script predicts classical DNA methylation clock values for the fasting
# cohort using batch-corrected, filtered beta matrices prepared in
# 2_qc_and_prepare_filtered_betas.qmd. It applies four external clock models
# (Horvath 353, Hannum, SkinBlood, and PhenoAge) through sesame::predictAge(),
# and computes DunedinPACE using DunedinPACE::PACEProjector().
#
# The probe-prioritized corrected beta matrix is used for the main prediction
# set because it preserves more CpGs for clock models. One sample
# (206467110003_R08C01) was excluded from that matrix during probe-prioritized
# filtering, so it is predicted separately from the sample-prioritized corrected
# matrix. For each clock, the script records both the predicted value and the
# number of model CpGs available in the beta matrix.
#
# Inputs:
#   - data/cleaned_beta_p_corrected.RDS: probe-prioritized corrected beta matrix
#   - data/cleaned_beta_s_corrected.RDS: sample-prioritized corrected beta matrix
#   - ../anno/Clock_*.rds: external classical clock model objects
#
# Output:
#   - data/ClassicClocks.csv: one row per sample with predicted clock values
#     and available-CpG counts for Horvath, Hannum, SkinBlood, PhenoAge, and
#     DunedinPACE.

# ============================================================================ #
# File Paths
# ============================================================================ #

# Input files
FILE_BETA_S_CORRECTED <- "./data/cleaned_beta_s_corrected.RDS"
FILE_BETA_P_CORRECTED <- "./data/cleaned_beta_p_corrected.RDS"
FILE_HORVATH353_MODEL <- "../anno/Clock_Horvath353.rds"
FILE_HANNUM_MODEL <- "../anno/Clock_Hannum.rds"
FILE_SKINBLOOD_MODEL <- "../anno/Clock_SkinBlood.rds"
FILE_PHENOAGE_MODEL <- "../anno/Clock_PhenoAge.rds"

# Output files
FILE_PREDICTED_AGE <- "data/ClassicClocks.csv"

# ============================================================================ #
# Libraries
# ============================================================================ #

library(tidyr)
library(dplyr)
library(sesame)
library(DunedinPACE)

# ============================================================================ #
# Load Data
# ============================================================================ #

# Load beta values
betass <- readRDS(FILE_BETA_S_CORRECTED)
betasp <- readRDS(FILE_BETA_P_CORRECTED)

# Load epigenetic clock models
model_Horvath353 <- readRDS(FILE_HORVATH353_MODEL)
model_Hannum <- readRDS(FILE_HANNUM_MODEL)
model_SkinBlood <- readRDS(FILE_SKINBLOOD_MODEL)
model_PhenoAge <- readRDS(FILE_PHENOAGE_MODEL)
model_ls <- list(
  Horvath = model_Horvath353,
  Hannum = model_Hannum,
  SkinBlood = model_SkinBlood,
  PhenoAge = model_PhenoAge
)

# ============================================================================ #
# Helper Functions
# ============================================================================ #

check_probes <- function(betas, model_list) {
  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    probes <- model$param$Probe_ID
    common_probes <- intersect(probes, rownames(betas))
    pct <- round(length(common_probes) / length(probes) * 100, 3)
    message(sprintf("%s: %d out of %d (%.3f%%) probes are preserved.",
                    model_name, length(common_probes), length(probes), pct))
  }
}

# ============================================================================ #
# Check Available Probes
# ============================================================================ #

check_probes(betasp, model_ls)

# Prepare data
betas <- betasp
los <- betass[, "206467110003_R08C01"]  # Left-out sample due to probe quality

# ============================================================================ #
# Predict Age Using Classical Epigenetic Clocks
# ============================================================================ #

# Minimum CpG requirements for each model
min_cpg_list <- list(
  Horvath = 300,
  Hannum = 50,
  SkinBlood = 300,
  PhenoAge = 450
)

# Predict age for all samples
df_age <- data.frame(Sample = character(), model = character(),
                     epi_age = numeric(), nCpG = integer())

for (model_name in names(model_ls)) {
  model <- model_ls[[model_name]]
  min_cpg <- min_cpg_list[[model_name]]
  
  # Predict for main samples
  for (pid in colnames(betas)) {
    profile <- betas[, pid]
    common_probes <- intersect(names(profile), model$param$Probe_ID)
    epi_age <- predictAge(profile, model, min_nonna = min_cpg)
    
    df_age <- rbind(df_age, data.frame(
      Sample = pid,
      model = model_name,
      epi_age = epi_age,
      nCpG = length(common_probes)
    ))
  }
  
  # Predict for left-out sample
  common_probes <- intersect(names(los), model$param$Probe_ID)
  epi_age <- predictAge(los, model, min_nonna = min_cpg)
  
  df_age <- rbind(df_age, data.frame(
    Sample = "206467110003_R08C01",
    model = model_name,
    epi_age = epi_age,
    nCpG = length(common_probes)
  ))
}

# ============================================================================ #
# DunedinPACE Clock
# ============================================================================ #

# Check probe availability
pbs_bg <- getRequiredProbes(backgroundList = TRUE)$DunedinPACE
n_common <- length(intersect(pbs_bg, rownames(betas)))
pct <- round(n_common / length(pbs_bg) * 100, 3)
message(sprintf("%d out of %d (%.3f%%) total probes are preserved.",
                n_common, length(pbs_bg), pct))

pbs <- getRequiredProbes(backgroundList = FALSE)$DunedinPACE
n_common <- length(intersect(pbs, rownames(betas)))
pct <- round(n_common / length(pbs) * 100, 3)
message(sprintf("%d out of %d (%.3f%%) modeling probes are preserved.",
                n_common, length(pbs), pct))

# Predict DunedinPACE for main samples
tmp <- PACEProjector(betas)
df_dunedin <- data.frame(
  Sample = names(tmp$DunedinPACE),
  model = "DunedinPACE",
  epi_age = tmp$DunedinPACE,
  nCpG = length(intersect(pbs, rownames(betas)))
)

# Predict DunedinPACE for left-out sample
tmp <- PACEProjector(betass, proportionOfProbesRequired = 0.7)
df_dunedin <- rbind(df_dunedin, data.frame(
  Sample = "206467110003_R08C01",
  model = "DunedinPACE",
  epi_age = tmp$DunedinPACE["206467110003_R08C01"],
  nCpG = length(intersect(pbs, names(los)))
))

# Combine all age predictions
df_age <- rbind(df_age, df_dunedin)

# ============================================================================ #
# Save Results
# ============================================================================ #

df_age <- df_age %>%
  pivot_wider(names_from = model, values_from = c(epi_age, nCpG))

write.table(df_age, file = FILE_PREDICTED_AGE, sep = ",",
            quote = FALSE, row.names = FALSE)
