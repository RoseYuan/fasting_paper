# -------------------------------
# Meta data processing
# -------------------------------

#' Categorize body mass index.
#'
#' Converts numeric BMI values into the project BMI categories used in metadata
#' and baseline plots.
#'
#' @param BMI Numeric vector of BMI values.
#' @return Character vector with BMI categories, or `NA` for missing BMI.
BMI_category <- function(BMI) {
  dplyr::case_when(
    is.na(BMI)            ~ NA_character_,
    BMI < 18.5            ~ "Underweight",
    BMI < 25              ~ "Healthy",
    BMI < 30              ~ "Overweight",
    BMI < 40              ~ "Obesity",
    TRUE                  ~ "Severe Obesity"
  )
}

#' Estimate basal metabolic rate with Oxford weight-only equations.
#'
#' Uses sex- and age-specific Oxford equations based on body weight. The output
#' is in MJ/day.
#'
#' @param weight Numeric body weight in kg.
#' @param age Numeric age in years.
#' @param sex Character scalar, either `"male"` or `"female"`.
#' @return Numeric BMR estimate in MJ/day, or `NA_real_` when inputs are missing
#'   or sex is not recognized.
BMR_Oxford <- function(weight, age, sex) {
  if (is.na(weight) || is.na(age) || is.na(sex)) return(NA_real_)
  if (sex == "male") {
    if (dplyr::between(age, 0, 3))      return(0.255  * weight - 0.141)
    if (dplyr::between(age, 3, 10))     return(0.0937 * weight + 2.15)
    if (dplyr::between(age, 10, 18))    return(0.0769 * weight + 2.43)
    if (dplyr::between(age, 18, 30))    return(0.0669 * weight + 2.28)
    if (dplyr::between(age, 30, 60))    return(0.0592 * weight + 2.48)
    if (dplyr::between(age, 60, 70))    return(0.0543 * weight + 2.37)
    if (age > 70)                       return(0.0573 * weight + 2.01)
  } else if (sex == "female") {
    if (dplyr::between(age, 0, 3))      return(0.246  * weight - 0.0965)
    if (dplyr::between(age, 3, 10))     return(0.0842 * weight + 2.12)
    if (dplyr::between(age, 10, 18))    return(0.0465 * weight + 3.18)
    if (dplyr::between(age, 18, 30))    return(0.0546 * weight + 2.33)
    if (dplyr::between(age, 30, 60))    return(0.0407 * weight + 2.90)
    if (dplyr::between(age, 60, 70))    return(0.0429 * weight + 2.39)
    if (age > 70)                       return(0.0417 * weight + 2.41)
  }
  NA_real_
}

#' Estimate basal metabolic rate with Oxford weight-height equations.
#'
#' Uses sex- and age-specific Oxford equations based on body weight and height.
#' The output is in MJ/day.
#'
#' @param weight Numeric body weight in kg.
#' @param height Numeric height in meters.
#' @param age Numeric age in years.
#' @param sex Character scalar, either `"male"` or `"female"`.
#' @return Numeric BMR estimate in MJ/day, or `NA_real_` when inputs are missing
#'   or sex is not recognized.
BMR_Oxford2 <- function(weight, height, age, sex) {
  # height in meters
  if (is.na(weight) || is.na(height) || is.na(age) || is.na(sex)) return(NA_real_)
  if (sex == "male") {
    if (dplyr::between(age, 0, 3))      return(0.118  * weight + 3.59 * height - 1.55)
    if (dplyr::between(age, 3, 10))     return(0.0632 * weight + 1.31 * height + 1.28)
    if (dplyr::between(age, 10, 18))    return(0.0651 * weight + 1.11 * height + 1.25)
    if (dplyr::between(age, 18, 30))    return(0.06   * weight + 1.31 * height + 0.473)
    if (dplyr::between(age, 30, 60))    return(0.0476 * weight + 2.26 * height - 0.574)
    if (age > 60)                       return(0.0478 * weight + 2.26 * height - 1.07)
  } else if (sex == "female") {
    if (dplyr::between(age, 0, 3))      return(0.127  * weight + 2.94 * height - 1.2)
    if (dplyr::between(age, 3, 10))     return(0.0666 * weight + 0.878 * height + 1.46)
    if (dplyr::between(age, 10, 18))    return(0.0393 * weight + 1.04 * height + 1.93)
    if (dplyr::between(age, 18, 30))    return(0.0433 * weight + 2.57 * height - 1.18)
    if (dplyr::between(age, 30, 60))    return(0.0342 * weight + 2.1  * height - 0.0486)
    if (age > 60)                       return(0.0356 * weight + 1.76 * height + 0.0448)
  }
  NA_real_
}

#' Calculate per-participant minimum value, loss, and regain summaries.
#'
#' Works on a wide table containing one row per participant and repeated
#' timepoint-specific columns such as `weight_D-2`, `weight_D+10`, and
#' `weight_M+1`. The function finds the minimum observed value across columns
#' matching `prefix`, calculates maximum loss from the starting column, and
#' calculates regain from the minimum to the follow-up column.
#'
#' @param wide_df Data frame containing an ID column and prefixed value columns.
#' @param prefix Column prefix to summarize, for example `"weight"` or `"BMI"`.
#' @param starting_col Baseline/start column used to calculate loss.
#' @param followup_col Follow-up column used to calculate regain.
#' @return `wide_df` with additional columns named `min_<prefix>`,
#'   `max_<prefix>_loss`, `max_<prefix>_loss_perc`,
#'   `max_<prefix>_regain`, and `max_<prefix>_regain_perc`.
calc_min_and_changes <- function(wide_df, prefix, starting_col, followup_col) {
  # wide_df: ID + columns like weight_Arrival, weight_FUP1...
  value_cols <- grep(paste0("^", prefix, "_"), names(wide_df), value = TRUE)
  
  if (!starting_col %in% names(wide_df)) {
    stop("Baseline column not found: ", starting_col)
  }
  if (!followup_col %in% names(wide_df)) {
    stop("Follow-up column not found: ", followup_col)
  }
  
  out <- wide_df %>%
    mutate(
      min_value = pmin(!!!rlang::syms(value_cols), na.rm = TRUE),
      max_loss = .data[[starting_col]] - min_value,
      max_loss_perc = (.data[[starting_col]] - min_value) / .data[[starting_col]],
      max_regain = .data[[followup_col]] - min_value,
      max_regain_perc = (.data[[followup_col]] - min_value) / min_value
    )
  
  # rename to match your original naming style
  out %>%
    rename_with(
      ~ paste0("min_", prefix),
      .cols = "min_value"
    ) %>%
    rename_with(
      ~ paste0("max_", prefix, "_loss"),
      .cols = "max_loss"
    ) %>%
    rename_with(
      ~ paste0("max_", prefix, "_loss_perc"),
      .cols = "max_loss_perc"
    ) %>%
    rename_with(
      ~ paste0("max_", prefix, "_regain"),
      .cols = "max_regain"
    ) %>%
    rename_with(
      ~ paste0("max_", prefix, "_regain_perc"),
      .cols = "max_regain_perc"
    )
}

#' Harmonize study timepoint labels.
#'
#' Converts heterogeneous timepoint labels from source metadata tables into the
#' project standard: `D-2`, `D-1`, `D0`, fasting days `D+1` to `D+12`,
#' reintroduction days `D+13` to `D+17`, and follow-up `M+1`.
#'
#' @param x Vector of raw timepoint labels.
#' @param unknown How to handle unmapped labels: `"keep"` preserves the original
#'   value with a warning, `"na"` replaces with `NA`, and `"error"` stops.
#' @return Character vector of harmonized timepoint labels.
harmonise_timepoint <- function(x, unknown = c("keep", "na", "error")) {
  unknown <- match.arg(unknown)
  
  x0 <- x
  s  <- as.character(x0)
  s  <- trimws(s)
  s  <- toupper(s)
  s  <- gsub("\\s+", "", s)
  
  out <- rep(NA_character_, length(s))
  
  
  # 2) Named special cases
  out[s %in% c("ARRIVAL", "Arrival")]    <- "D-2"
  out[s %in% c("BL","BASELINE","BEFORE","DIET","PRE","Diet",1)] <- "D-1" 
  out[s %in% c("TRANSITION","Transition")] <- "D0"
  out[s %in% c("END","AFTER",2)] <- "D+10" 
  out[s %in% c("FUP1",3)]   <- "M+1"
  
  # 3) Fasting days: F1..F12 -> D+1..D+12
  is_F <- grepl("^F0*[0-9]+$", s)
  if (any(is_F)) {
    k <- as.integer(sub("^F0*", "", s[is_F]))
    out[is_F] <- paste0("D+", k)
  }
  
  # 4) Refeed days: REF1..REF5 -> D+13..D+17
  is_REF <- grepl("^REF0*[0-9]+$", s)
  if (any(is_REF)) {
    r <- as.integer(sub("^REF0*", "", s[is_REF]))
    out[is_REF] <- paste0("D+", 12 + r)  # REF1->D+13
  }
  
  # Handle unknowns
  missing <- is.na(out) & !is.na(s) & s != ""
  if (any(missing)) {
    bad <- sort(unique(s[missing]))
    msg <- paste0("Unmapped Timepoint values: ", paste(bad, collapse = ", "))
    
    if (unknown == "error") stop(msg, call. = FALSE)
    if (unknown == "na") {
      warning(msg)
      out[missing] <- NA_character_
    } else { # keep
      warning(msg)
      out[missing] <- x0[missing]  # keep original for traceability
    }
  }
  
  out
}

#' Reconcile paired `.x` and `.y` columns after merging metadata tables.
#'
#' Compares two vectors representing the same variable from different sources.
#' Numeric-like values are compared within a tolerance; non-numeric values must
#' match exactly when both sources are non-missing. Missing values are filled
#' from the other source.
#'
#' @param x First source vector.
#' @param y Second source vector.
#' @param name Variable name used in conflict error messages.
#' @param tol Numeric tolerance for numeric-like comparisons.
#' @return A single reconciled vector, numeric when both inputs are numeric-like
#'   and character otherwise.
resolve_xy <- function(x, y, name = "var", tol = 1e-6) {
  # 1) Convert both to character, normalize common "NA" representations
  x_chr <- trimws(as.character(x))
  y_chr <- trimws(as.character(y))
  
  x_chr[x_chr %in% c("", "NA", "NaN", "N/A", "<NA>")] <- NA_character_
  y_chr[y_chr %in% c("", "NA", "NaN", "N/A", "<NA>")] <- NA_character_
  
  # 2) Try numeric comparison if both sides are parseable as numeric
  x_num <- suppressWarnings(as.numeric(x_chr))
  y_num <- suppressWarnings(as.numeric(y_chr))
  
  x_is_num <- !is.na(x_num) | is.na(x_chr)  # NA is allowed
  y_is_num <- !is.na(y_num) | is.na(y_chr)
  
  if (all(x_is_num) && all(y_is_num)) {
    conflict <- !is.na(x_num) & !is.na(y_num) & (abs(x_num - y_num) > tol)
    if (any(conflict)) {
      stop("Conflict in numeric field '", name, "'. مثال: ",
           paste(head(which(conflict), 3), collapse = ", "),
           call. = FALSE)
    }
    return(dplyr::coalesce(x_num, y_num))
  }
  
  # 3) Otherwise treat as string-like and require exact match (when both present)
  conflict <- !is.na(x_chr) & !is.na(y_chr) & (x_chr != y_chr)
  if (any(conflict)) {
    stop("Conflict in non-numeric field '", name, "'. Example values: '",
         x_chr[which(conflict)[1]], "' vs '", y_chr[which(conflict)[1]], "'",
         call. = FALSE)
  }
  
  dplyr::coalesce(x_chr, y_chr)
}

#' Summarize agreement for a shared metadata column.
#'
#' Used after joining the three-timepoint metadata table with the daily metadata
#' table. The input must contain paired columns named `<col>.meta` and
#' `<col>.daily`.
#'
#' @param df Joined metadata data frame containing `.meta` and `.daily` columns.
#' @param col Base column name to compare.
#' @param tol Numeric tolerance for numeric-like comparisons.
#' @return A one-row tibble with the column name, number of values present in
#'   `.meta` but missing in `.daily`, and number of disagreements.
check_shared_column <- function(df, col, tol = 1e-6) {
  a <- df[[paste0(col, ".meta")]]
  b <- df[[paste0(col, ".daily")]]
  
  # standardize "NA" strings to NA
  a_chr <- trimws(as.character(a)); a_chr[a_chr %in% c("", "NA", "<NA>", "NaN")] <- NA
  b_chr <- trimws(as.character(b)); b_chr[b_chr %in% c("", "NA", "<NA>", "NaN")] <- NA
  
  # try numeric comparison if both parse as numeric (allowing NA)
  a_num <- suppressWarnings(as.numeric(a_chr))
  b_num <- suppressWarnings(as.numeric(b_chr))
  a_numlike <- all(is.na(a_chr) | !is.na(a_num))
  b_numlike <- all(is.na(b_chr) | !is.na(b_num))
  
  meta_non_na_daily_na <- sum(!is.na(a_chr) & is.na(b_chr))
  
  if (a_numlike && b_numlike) {
    disagree <- !is.na(a_num) & !is.na(b_num) & abs(a_num - b_num) > tol
    n_disagree <- sum(disagree)
  } else {
    disagree <- !is.na(a_chr) & !is.na(b_chr) & (a_chr != b_chr)
    n_disagree <- sum(disagree)
  }
  
  tibble(
    column = col,
    meta_non_na_daily_na = meta_non_na_daily_na,
    n_disagree = n_disagree
  )
}

#' Collapse duplicated merge columns into one column per variable.
#'
#' Removes repeated `.x`/`.y` suffixes introduced by joins, checks duplicated
#' columns for row-level inconsistencies, and coalesces non-missing values into a
#' single output column.
#'
#' @param df Data frame containing a `Sample` column and possibly duplicated
#'   merge columns such as `Age.x`, `Age.y`, or `Age.x.x`.
#' @return Data frame with `Sample` plus one collapsed column per base variable.
collapse_duplicates <- function(df) {
  base_names <- names(df) %>%
    str_replace("\\.[xy](\\.[xy])*$", "") %>%  # remove .x, .y, .x.x, etc
    unique()

  out <- df %>% select(any_of("Sample"))

  for (nm in setdiff(base_names, "Sample")) {
    cols <- names(df)[str_replace(names(df), "\\.[xy](\\.[xy])*$", "") == nm]

    if (length(cols) == 1) {
      out[[nm]] <- df[[cols]]
    } else {
      # optional consistency check
      vals <- df[, cols, drop = FALSE]
      inconsistent <- apply(vals, 1, function(x) {
        ux <- unique(na.omit(x))
        length(ux) > 1
      })
      if (any(inconsistent)) {
        warning("Inconsistent values in duplicated column: ", nm)
      }
      out[[nm]] <- dplyr::coalesce(!!!vals)
    }
  }
  out
}

# -------------------------------
# Name mapping
# -------------------------------
var_to_full <- c(
  BASO = "Basophils (%)",
  BASO_abs = "Absolute basophil count (ABC)",
  Bilirubin = "Total bilirubin",
  EOSINOPHILE = "Eosinophils (%)",
  EOSINOPHILE_abs = "Absolute eosinophil count (AEC)",
  Erythrozytenverteilungsbreite = "Red cell distribution width (RDW)",
  GFR = "Estimated glomerular filtration rate (eGFR)",
  Homa = "HOMA-IR",
  `Hämoglobingehalt der Retis` = "Reticulocyte hemoglobin content (Ret-He)",
  Insulin = "Insulin",
  Ketosemia = "Blood ketone concentration (ketonemia)",
  LYMPHOZYTEN = "Lymphocytes (%)",
  Lymphozyten_absolut = "Absolute lymphocyte count (ALC)",
  MONOZYTEN = "Monocytes (%)",
  MONOZYTEN_abs = "Absolute monocyte count (AMC)",
  NEUTROPHILE = "Neutrophils (%)",
  NLR = "Neutrophil-to-lymphocyte ratio (NLR)",
  Neutrophile_absolut = "Absolute neutrophil count (ANC)",
  NonHDL.Cholesterin = "Non-HDL cholesterol (non-HDL-C)",
  ProMyeloMetamyelozyten = "Pro- and metamyelocytes (%)",
  ProMyeloMetamyelozyten_abs = "Absolute pro- and metamyelocyte count",
  Retikulozyten = "Reticulocytes (%)",
  `Retikulozyten ProduktionsIndex` = "Reticulocyte production index (RPI)",
  `Retikulozyten absolut` = "Absolute reticulocyte count",
  alk_phosphatase = "Alkaline phosphatase (ALP)",
  cholesterol = "Total cholesterol",
  creatinine = "Creatinine",
  crp_hs = "High-sensitivity C-reactive protein (hsCRP)",
  dbp = "Diastolic blood pressure (DBP)",
  erythrocytes = "Red blood cell count (RBC)",
  ggt = "Gamma-glutamyl transferase (GGT)",
  glucose = "Glucose",
  glykohemoglobin = "Glycated hemoglobin (HbA1c)",
  `got/AST` = "Aspartate aminotransferase (AST)",
  `gpt/ALT` = "Alanine aminotransferase (ALT)",
  hdl = "High-density lipoprotein cholesterol (HDL-C)",
  hematocrit = "Hematocrit (Hct)",
  hemoglobin = "Hemoglobin (Hb)",
  inr = "International normalized ratio (INR)",
  ldl = "Low-density lipoprotein cholesterol (LDL-C)",
  ldl_hdl_ratio = "LDL-to-HDL cholesterol ratio",
  leukocytes = "White blood cell count (WBC)",
  mch = "Mean corpuscular hemoglobin (MCH)",
  mchc = "Mean corpuscular hemoglobin concentration (MCHC)",
  mcv = "Mean corpuscular volume (MCV)",
  ptt = "Partial thromboplastin time (PTT)",
  pulse = "Heart rate",
  quick = "Prothrombin time (Quick)",
  sbp = "Systolic blood pressure (SBP)",
  thrombocytes = "Platelet count (PLT)",
  triglyceride = "Triglycerides",
  urea = "Urea",
  uric_acid = "Uric acid"
)

var_to_abbr <- c(
  hemoglobin = "Hb",
  hematocrit = "Hct",
  mcv = "MCV",
  mch = "MCH",
  mchc = "MCHC",
  Erythrozytenverteilungsbreite = "RDW",
  leukocytes = "WBC",
  erythrocytes = "RBC",
  thrombocytes = "PLT",
  crp_hs = "hsCRP",
  inr = "INR",
  ptt = "PTT",
  `got/AST` = "AST",
  `gpt/ALT` = "ALT",
  ggt = "GGT",
  alk_phosphatase = "ALP",
  hdl = "HDL-C",
  ldl = "LDL-C",
  NLR = "NLR",
  glykohemoglobin = "HbA1c",
  GFR = "eGFR"
)

var_to_short <- var_to_full
var_to_short[names(var_to_abbr)] <- var_to_abbr
# -------------------------------
# Visualization
# -------------------------------

###############
# PCA visualization
# https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
##############
suppressPackageStartupMessages({
library(PCAtools)
})

#' Create a scree plot with variance and cumulative-variance labels.
#'
#' Wraps `PCAtools::screeplot()` and adds text labels for per-component variance
#' explained and cumulative variance explained.
#'
#' @param pcaobj PCAtools PCA object, typically returned by `PCAtools::pca()`.
#' @param components Components to display. Defaults to all components returned
#'   by `PCAtools::getComponents(pcaobj)`.
#' @param ... Additional arguments passed to `PCAtools::screeplot()`.
#' @return A ggplot object.
myscreeplot <- function(pcaobj, components=getComponents(pcaobj), ...){
  plotobj <- data.frame(components, pcaobj$variance[components])
  colnames(plotobj) <- c('PC', 'Variance')
  plotobj$PC <- factor(plotobj$PC,
                       levels=plotobj$PC[seq_along(plotobj$PC)])
  cum_var_label <- paste0(round(cumsum(plotobj$Variance), 2), "%")
  cum_var_label[1] <- ""
  p <- screeplot(
    pcaobj,
    components= components,
    colBar = 'steelblue',
    colCumulativeSumLine = 'orange',
    colCumulativeSumPoints = 'orange',
    returnPlot = TRUE,
    ylim = c(0, 105),
    ...) +
    geom_text(data=plotobj, aes(y = Variance,
                                label = paste0(round(Variance, 2), "%")),
              vjust = -0.5, size = 5, color = "#000000") +
    geom_text(aes(y = cumsum(Variance),
                  label = cum_var_label),
              vjust = -0.5, size = 5, color = "#000000", nudge_x = 0.15)
  return(p)
}


#' Create a stable qualitative color palette for participant IDs.
#'
#' Builds a Polychrome palette and names each color by the corresponding ID so
#' participant colors remain stable across plots.
#'
#' @param ids Vector of participant IDs.
#' @return Named character vector of hex colors.
make_id_palette <- function(ids) {
  ids <- sort(unique(ids))
  seed <- c("#ff0000", "#00ff00", "#0000ff")
  pal <- Polychrome::createPalette(length(ids), seed)
  names(pal) <- as.character(ids)
  pal
}


palette("Polychrome 36")
theme_set(theme_minimal())

patient_colors <- setNames(palette("Polychrome 36")[1:33], 1:33)
labels_3timepoint <- c("D-1", "D+10", "M+1")
timepoint_numeric_values <- c(0, 2, 6)
