# ==============================================================================
# Merge Fasting Study Metadata
# ==============================================================================

# This script builds the main metadata tables used throughout the fasting
# methylation and clinical analyses. It combines the IDAT sample sheet with
# daily clinical metadata, the three-timepoint clinical metadata export, and the
# AGE/glycation table. Timepoint labels are harmonized with project utilities,
# overlapping variables from the daily and three-timepoint sources are
# reconciled, and array/sample-sheet fields are retained for downstream
# methylation QC and batch modeling.
#
# Derived variables include:
#   - ID_Patient identifiers in the B__### format
#   - QC outlier flag for sample 4033964733
#   - BMI, baseline BMI, and baseline BMI category
#   - BMR estimates
#   - per-participant minimum weight/BMI/BMR2 and maximum loss/regain summaries
#
# Outputs:
#   - publication/data/meta_3timepoints.csv: merged D-1, D+10, and M+1 metadata
#     aligned to methylation sample-sheet fields and AGE/glycation variables.
#   - publication/data/meta_daily.csv: daily clinical metadata with harmonized
#     timepoints, BMI/BMR variables, and baseline BMI categories.
#
# The final block treats the daily table as a superset and reports disagreement
# or missingness for shared columns at the three methylation timepoints.
# ==============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  source("publication/utils.R")
})

# -------------------------------
# Paths (edit here only)
# -------------------------------
PATH_SAMPLESHEET <- "IDATs_sorted_Genesis/14.12.2022_sample_sheet_R&D_Buchinger_ALL.xlsx"
PATH_META1       <- "metadata/GENESIS_DATA_19.05.2022_blooddata.xlsx"      # all timepoints
PATH_META2       <- "metadata/GENESIS_DATA_FOR_GENKNOWME.xlsx"             # 3 timepoints (skip 12 rows)
PATH_AGE         <- "metadata/AGE-file.xlsx"

OUT_META_3TP     <- "publication/data/meta_3timepoints.csv"
OUT_META_DAILY   <- "publication/data/meta_daily.csv"


# ==============================================================================
# 1) Load sample sheet -> pheno
# ==============================================================================
samplesheet <- read_excel(PATH_SAMPLESHEET, sheet = 1) %>%
  mutate(Sample = Basename)

# Keep factors only where it makes sense; avoid converting everything to factor
pheno <- samplesheet %>%
  transmute(
    Sample,
    Sample_Name,
    Sample_Plate,
    ID_Patient = `ID Patient`,
    AMP_Plate,
    SentrixBarcode_A,
    SentrixPosition_A,
    Date_Scan,
    Gender,
    Timepoint
  ) %>%
  mutate(
    Sample = as.character(Sample),
    ID_Patient = as.character(ID_Patient),
    Timepoint = as.character(Timepoint),
    Gender = as.character(Gender),
    qc_outliers = Sample_Name == "4033964733"
  )

# ==============================================================================
# 2) Load metadata tables
# ==============================================================================
blmeta <- read_excel(PATH_META1, sheet = 1)

# Create ID_Patient-like field used later for merging (as in your original)
blmeta <- blmeta %>%
  mutate(
    ID_Patient = paste0("B__", sprintf("%03d", ID))
  )
blmeta2 <- read_excel(PATH_META2, sheet = 1, skip = 12)
age_gly    <- read_excel(PATH_AGE,  sheet = 1)

# Harmonise timepoints in all tables
pheno  <- pheno  %>% mutate(Timepoint = harmonise_timepoint(Timepoint, unknown = "error"))
blmeta <- blmeta %>% mutate(Timepoint = harmonise_timepoint(Timepoint, unknown = "error"))
blmeta2<- blmeta2%>% mutate(Timepoint = harmonise_timepoint(Timepoint, unknown = "error"))
age_gly<- age_gly %>% mutate(Timepoint = harmonise_timepoint(Timepoint, unknown = "error"))

# ==============================================================================
# 3) Add height from blmeta2 into blmeta + compute BMRs
# ==============================================================================
blmeta <- blmeta %>%
  left_join(
    blmeta2 %>% select(ID, Height) %>% distinct(),
    by = "ID"
  )

# Sex mapping from your original:
#   sex is integer-coded (0/1?), then +1 indexes c("male","female")
sex_chr <- c("male", "female")

tmp1 <- tibble(
  weight = as.numeric(blmeta$weight),
  age    = as.numeric(blmeta$age),
  sex    = sex_chr[as.integer(blmeta$sex) + 1L]
)
blmeta$BMR1 <- purrr::pmap_dbl(tmp1, BMR_Oxford)

tmp2 <- tibble(
  weight = as.numeric(blmeta$weight),
  age    = as.numeric(blmeta$age),
  sex    = sex_chr[as.integer(blmeta$sex) + 1L],
  height = as.numeric(blmeta$Height) * 0.01
)
blmeta$BMR2 <- purrr::pmap_dbl(tmp2, BMR_Oxford2)

# Quick NA diagnostics
na_rates <- sort(round(colMeans(is.na(as.data.frame(blmeta))), 3))
print(na_rates)

# find columns ~0.81 NA but with some non-NA rows
col_081 <- names(blmeta)[round(colMeans(is.na(as.data.frame(blmeta))), 2) == 0.81]
if (length(col_081) > 0) {
  message("Columns with ~0.81 NA rate: ", paste(col_081, collapse = ", "))
  print(
    blmeta %>%
      select(any_of(c("ID", "Timepoint", col_081))) %>%
      filter(if_any(all_of(col_081), ~ !is.na(.)))
  )
}

# ==============================================================================
# 4) Per-individual indices: max weight/BMI/BMR change (from daily blmeta)
# ==============================================================================
# ---- weight
weights_wide <- blmeta %>%
  select(ID, Timepoint, weight) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = weight,
    names_glue = "{.value}_{Timepoint}"
  )

weights_wide <- calc_min_and_changes(
  weights_wide,
  prefix = "weight",
  starting_col = "weight_D-2",
  followup_col = "weight_M+1"
)

# ---- BMI
blmeta <- blmeta %>%
  mutate(BMI = as.numeric(weight) / (as.numeric(Height) * 0.01)^2)

bmi_wide <- blmeta %>%
  select(ID, Timepoint, BMI) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = BMI,
    names_glue = "{.value}_{Timepoint}"
  )

bmi_wide <- calc_min_and_changes(
  bmi_wide,
  prefix = "BMI",
  starting_col = "BMI_D-2",
  followup_col = "BMI_M+1"
)

# ---- BMR2
bmr_wide <- blmeta %>%
  select(ID, Timepoint, BMR2) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = BMR2,
    names_glue = "{.value}_{Timepoint}"
  )

bmr_wide <- calc_min_and_changes(
  bmr_wide,
  prefix = "BMR2",
  starting_col = "BMR2_D-2",
  followup_col = "BMR2_M+1"
)

daily_meta <- blmeta
# ==============================================================================
# 5) Merge meta data: align to 3 timepoints (D-1/D+10/M+1)
# ==============================================================================
# Consistency check between blmeta and blmeta2 on common columns (optional but reproducible)
meta_com <- intersect(colnames(blmeta), colnames(blmeta2))

tmp1 <- blmeta %>%
  select(all_of(meta_com)) %>%
  filter(Timepoint %in% c("D-1", "D+10", "M+1")) 

tmp2 <- blmeta2 %>% select(all_of(meta_com))


merged_tmp <- merge(tmp1, tmp2, by = c("ID", "Timepoint"), all = TRUE)

# Reconcile all overlapping columns generically
common_vars <- setdiff(meta_com, c("ID", "Timepoint"))

for (v in common_vars) {
  merged_tmp[[v]] <- resolve_xy(
    merged_tmp[[paste0(v, ".x")]],
    merged_tmp[[paste0(v, ".y")]],
    name = v
  )
}
# drop the .x/.y columns
merged_tmp <- merged_tmp %>%
  select(ID, Timepoint, all_of(common_vars))

# Now build 3-timepoint blmeta and merge with blmeta2
blmeta2_only <- setdiff(colnames(blmeta2), colnames(blmeta))

merged_meta <- merged_tmp %>%
  left_join(
    blmeta2 %>%
      filter(Timepoint %in% c("D-1", "D+10", "M+1")) %>%
      select(ID, Timepoint, all_of(blmeta2_only)),
    by = c("ID", "Timepoint")
  )

# Now build 3-timepoint blmeta and merge with blmeta2
blmeta_only <- setdiff(colnames(blmeta), colnames(blmeta2))

merged_meta <- merged_meta %>%
  left_join(
    blmeta %>%
      filter(Timepoint %in% c("D-1", "D+10", "M+1")) %>%
      select(ID, Timepoint, all_of(blmeta_only)),
    by = c("ID", "Timepoint")
  )

# merge pheno
pheno <- pheno %>%
  mutate(
    ID = as.integer(sub("^B__", "", ID_Patient))
  )
meta <- merge(pheno, merged_meta, by = c("ID", "Timepoint"), all = TRUE)

# Merge AGE table (ID + Timepoint)
meta <- merge(meta, age_gly, by = c("ID", "Timepoint"), all = TRUE)

# Add per-individual indices
meta <- meta %>%
  left_join(weights_wide %>% select(ID, matches("^min_weight|max_weight_loss|max_weight_loss_perc|max_weight_regain|max_weight_regain_perc")),
            by = "ID") %>%
  left_join(bmi_wide %>% select(ID, matches("^min_BMI|max_BMI_loss|max_BMI_loss_perc|max_BMI_regain|max_BMI_regain_perc")),
            by = "ID") %>%
  left_join(bmr_wide %>% select(ID, matches("^min_BMR2|max_BMR2_loss|max_BMR2_loss_perc|max_BMR2_regain|max_BMR2_regain_perc")),
            by = "ID")

# Drop column if it exists (was in your QMD)
if ("others_describe" %in% names(meta)) {
  meta <- meta %>% select(-others_describe)
}

TP3 <- c("D-1", "D+10", "M+1")

daily_meta <- daily_meta %>%
  left_join(
    meta %>%
      filter(Timepoint %in% TP3) %>%
      select(ID, Timepoint, weight) %>%
      dplyr::rename(weight_meta = weight),
    by = c("ID", "Timepoint")
  ) %>%
  mutate(weight = coalesce(weight, weight_meta)) %>%
  select(-weight_meta)

# Finalise ID_Patient field
meta <- meta %>%
  mutate(ID_Patient = paste0("B__", sprintf("%03d", as.integer(ID)))) %>%
  select(-ID_Patient.x, -ID_Patient.y)

# Add BMI category at baseline (D-1) 
daily_meta <- daily_meta %>%
  mutate(BMI = as.numeric(BMI)) %>%
  group_by(ID) %>%
  mutate(BMI_baseline = BMI[Timepoint == "D-1"][1]) %>%
  mutate(BMI_baseline_category = BMI_category(BMI_baseline)) %>%
  ungroup()

meta <- meta %>%
  mutate(BMI = as.numeric(BMI)) %>%
  group_by(ID) %>%
  mutate(BMI_baseline = BMI[Timepoint == "D-1"][1]) %>%
  mutate(BMI_baseline_category = BMI_category(BMI_baseline)) %>%
  ungroup()
# ==============================================================================
# 6) Save
# ==============================================================================
# Keep daily_meta ID as character (avoid factor surprises)
daily_meta <- daily_meta %>% mutate(ID = as.character(ID))

write.table(meta, file = OUT_META_3TP, sep = ",", quote = FALSE, row.names = FALSE)
write.table(daily_meta, file = OUT_META_DAILY, sep = ",", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 7) treating daily_meta as the superset and verifying that all data in meta exists in daily_meta
# ==============================================================================
# Ensure same key types
meta <- meta %>% mutate(ID = as.integer(ID), Timepoint = as.character(Timepoint))
daily_meta <- daily_meta %>% mutate(ID = as.integer(ID), Timepoint = as.character(Timepoint))

# Join daily_meta onto meta (keep only the keys present in meta)
joined <- meta %>%
  select(ID, Timepoint, everything()) %>%
  left_join(
    daily_meta,
    by = c("ID", "Timepoint"),
    suffix = c(".meta", ".daily")
  )

shared_cols <- intersect(names(meta), names(daily_meta))
shared_cols <- setdiff(shared_cols, c("ID", "Timepoint"))

report <- bind_rows(lapply(shared_cols, \(c) check_shared_column(joined, c))) %>%
  arrange(desc(meta_non_na_daily_na), desc(n_disagree))

report
