# -------------------------------
# Meta data processing
# -------------------------------

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

# Oxford equations (BMR in MJ/day)
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
# Visualization
# -------------------------------

###############
# PCA visualization
# https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
##############
suppressPackageStartupMessages({
library(PCAtools)
})

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


# helper: stable palette per ID
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
