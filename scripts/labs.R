#!/usr/bin/env Rscript
###############################################################################
# labs.R -- R port of 04_labs.sql.
#
# Output microbiology_Labs.csv: one row per
# (anon_id, pat_enc_csn_id_coded, order_proc_id_coded, Period_Day),
# Period_Day in {14, 30, 180} days. Per lab, 5 columns:
# Q25/Q75/median/first/last. Labs: wbc, neutrophils, lymphocytes, hgb, plt,
# na, hco3, bun, cr, lactate, procalcitonin. (crp is defined but never
# aggregated in the SQL, so it is omitted here too.)
#
# Run: Rscript labs.R  -- edit the CONFIGURATION block first.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})


# =============================================================================
# CONFIGURATION
# =============================================================================

RAW_DIR    <- "/PATH/TO/raw/edw/2024-12-31"
COHORT_CSV <- "/PATH/TO/microbiology_cultures_cohort.csv"
OUTPUT_CSV <- "microbiology_Labs.csv"

# ADAPT: deid. Set APPLY_DEID=TRUE to swap raw IDs for the project's deid
# integers before writing. Point MAPS_DIR at the folder holding
# anon_id_map.csv / enc_id_map.csv / order_proc_id_map.csv.
APPLY_DEID <- TRUE
MAPS_DIR   <- "/PATH/TO/maps"

# ADAPT: cohort column names (your file -> Stanford identifier).
COHORT_PATIENT_COL   <- "patientID"           # -> anon_id
COHORT_ENCOUNTER_COL <- "patientEncounterID"  # -> pat_enc_csn_id_coded
COHORT_ORDER_COL     <- "orderProcedureID"    # -> order_proc_id_coded
COHORT_TIME_COL      <- "cultureCollDTS"      # -> order_time_jittered_utc

# ADAPT: one CSV per lab (MGB EDW = one file per variable). Fix any filename
# that differs. crp omitted on purpose (see header).
LAB_FILES <- list(
  wbc           = file.path(RAW_DIR, "wbc.csv"),
  neutrophils   = file.path(RAW_DIR, "neutrophils.csv"),
  lymphocytes   = file.path(RAW_DIR, "lymphocytes.csv"),
  hgb           = file.path(RAW_DIR, "hgb.csv"),
  plt           = file.path(RAW_DIR, "plt.csv"),
  na            = file.path(RAW_DIR, "na.csv"),
  hco3          = file.path(RAW_DIR, "hco3.csv"),
  bun           = file.path(RAW_DIR, "bun.csv"),
  cr            = file.path(RAW_DIR, "cr.csv"),
  lactate       = file.path(RAW_DIR, "lactate.csv"),
  procalcitonin = file.path(RAW_DIR, "procalcitonin.csv")
)

# Lookback windows in days (matches SQL).
LOOKBACK_DAYS <- c(14L, 30L, 180L)


# =============================================================================
# SQL -> MGB SCHEMA MAPPING
#
# Stanford BigQuery                   MGB EDW CSV
# ---------------------------------   ----------------------------------------
# c.anon_id                           patientID (raw) / anon_id (deid)
# c.pat_enc_csn_id_coded              patientEncounterID (raw) / ..._coded (deid)
# c.order_proc_id_coded               orderProcedureID
# c.order_time_jittered_utc           cultureCollDTS
# lr.order_time_jittered_utc          orderDTS
# lr.base_name / lr.lab_name          implicit from FILE NAME (one lab per file)
# lr.reference_unit                   referenceRangeUnit
# lr.ord_value                        resultN (numeric) / resultTXT (e.g. '>10')
#
# ADAPT checklist -- verify against your MGB export:
# (A) UNIT_MAP is verbatim from the SQL (lowercased/trimmed). Extend an entry
#     if referenceRangeUnit has other spellings/casings (e.g. '10^3/uL').
# (B) If one CSV mixes lab types, set COMPONENT_NAME_FILTERS[[lab]].
# (C) resultN used first, resultTXT as fallback (USE_RESULTTXT_FALLBACK).
#
# Deliberate matches to the SQL (do NOT "fix"):
# * hgb g/dl -> mg/dl multiplies by 1000 (clinically odd, but matches SQL).
# * hco3 accepts meq/l or mmol/l; lactate accepts base_name 'lac'.
# * neutrophils/lymphocytes require unit '%' (absolute counts dropped).
# =============================================================================

# (A) Unit-string sets; compared lowercased/trimmed.
UNIT_MAP <- list(
  wbc = list(
    one_to_one  = c("thousand/ul","k/ul","10x3/ul","10*3/ul","x10e3/ul"),
    divide_1000 = c("/ul")
  ),
  neutrophils = list(pct_only = c("%")),
  lymphocytes = list(pct_only = c("%")),
  hgb = list(
    one_to_one    = c("mg/dl"),
    multiply_1000 = c("g/dl")
  ),
  plt = list(
    one_to_one  = c("x10e3/ul","10x3/ul","k/ul","10*3/ul","thousand/ul"),
    divide_1000 = c("ul")
  ),
  na            = list(one_to_one = c("mmol/l")),
  hco3          = list(one_to_one = c("meq/l","mmol/l")),
  bun           = list(one_to_one = c("mg/dl")),
  cr            = list(one_to_one = c("mg/dl")),
  lactate       = list(one_to_one = c("mmol/l","mmole/l")),
  procalcitonin = list(one_to_one = c("ng/ml"))
)

# (B) ADAPT: per-lab componentNM allow-set for files with mixed lab types.
# e.g. COMPONENT_NAME_FILTERS$wbc <- c("WBC", "WBC AUTO DIFFERENTIAL")
COMPONENT_NAME_FILTERS <- list(
)

# (C) TRUE = fall back to as.numeric(resultTXT) when resultN is NA.
USE_RESULTTXT_FALLBACK <- TRUE

# ADAPT: column names in every MGB lab CSV.
LAB_ID_COL        <- "patientID"
LAB_ENC_COL       <- "patientEncounterID"
LAB_TIME_COL      <- "orderDTS"
LAB_UNIT_COL      <- "referenceRangeUnit"
LAB_VALUE_N_COL   <- "resultN"
LAB_VALUE_TXT_COL <- "resultTXT"
LAB_COMPONENT_COL <- "componentNM"


# =============================================================================
# UNIT CONVERTERS -- return the value in the SQL's unified unit, else NA.
# =============================================================================

.norm <- function(x) tolower(trimws(as.character(x)))

convert_generic <- function(value, unit, rules) {
  u <- .norm(unit)
  v <- suppressWarnings(as.numeric(value))
  out <- rep(NA_real_, length(v))
  if (!is.null(rules$one_to_one)) {
    m <- u %chin% rules$one_to_one
    out[m] <- v[m]
  }
  if (!is.null(rules$multiply_1000)) {
    m <- u %chin% rules$multiply_1000
    out[m] <- v[m] * 1000
  }
  if (!is.null(rules$divide_1000)) {
    m <- u %chin% rules$divide_1000
    out[m] <- v[m] / 1000
  }
  if (!is.null(rules$pct_only)) {
    m <- u %chin% rules$pct_only
    out[m] <- v[m]
  }
  out
}

convert <- function(lab, value, unit) {
  convert_generic(value, unit, UNIT_MAP[[lab]])
}


# =============================================================================
# I/O
# =============================================================================

read_cohort <- function() {
  df <- fread(COHORT_CSV)
  required <- c(COHORT_PATIENT_COL, COHORT_ENCOUNTER_COL,
                COHORT_ORDER_COL,   COHORT_TIME_COL)
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop("Cohort missing columns: ", paste(missing, collapse=", "),
         "\nFix COHORT_*_COL constants at the top of the script.")
  }
  out <- df[, .(
    anon_id                 = as.character(get(COHORT_PATIENT_COL)),
    pat_enc_csn_id_coded    = as.character(get(COHORT_ENCOUNTER_COL)),
    order_proc_id_coded     = as.character(get(COHORT_ORDER_COL)),
    order_time_jittered_utc = as.POSIXct(get(COHORT_TIME_COL))
  )]
  out[]
}

read_lab <- function(lab, path) {
  if (!file.exists(path)) {
    message("  [skip] ", lab, ": file not found at ", path); return(NULL)
  }
  dt <- fread(path)
  required <- c(LAB_ID_COL, LAB_ENC_COL, LAB_TIME_COL, LAB_UNIT_COL, LAB_VALUE_N_COL)
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(path, " missing columns: ", paste(missing, collapse=", "))
  }

  # (B) Optional componentNM filter for files with mixed lab types.
  keep <- COMPONENT_NAME_FILTERS[[lab]]
  if (!is.null(keep) && LAB_COMPONENT_COL %in% names(dt)) {
    dt <- dt[get(LAB_COMPONENT_COL) %chin% keep]
  }

  # (C) resultN first, then optional resultTXT fallback.
  val <- suppressWarnings(as.numeric(dt[[LAB_VALUE_N_COL]]))
  if (USE_RESULTTXT_FALLBACK && LAB_VALUE_TXT_COL %in% names(dt)) {
    fallback <- suppressWarnings(as.numeric(dt[[LAB_VALUE_TXT_COL]]))
    val <- ifelse(is.na(val), fallback, val)
  }

  dt[, value := convert(lab, val, get(LAB_UNIT_COL))]
  dt <- dt[!is.na(value)]
  if (nrow(dt) == 0L) return(NULL)

  out <- dt[, .(
    anon_id              = as.character(get(LAB_ID_COL)),
    pat_enc_csn_id_coded = as.character(get(LAB_ENC_COL)),
    labtime              = as.POSIXct(get(LAB_TIME_COL)),
    value
  )]
  out <- out[!is.na(labtime)]
  out[]
}


# =============================================================================
# WINDOWED AGGREGATION -- per window: join lab stream to cohort on
# (anon_id, pat_enc_csn_id_coded), keep TIMESTAMP_DIFF(lab, cohort, DAY) in
# [-W, 0], then per (patient, enc, order) emit Q25/Q75/median/first/last.
# =============================================================================

agg_one_lab <- function(events, cohort, lookback_days, lab) {
  if (is.null(events) || nrow(events) == 0L) return(data.table())
  setkey(events, anon_id, pat_enc_csn_id_coded)
  setkey(cohort, anon_id, pat_enc_csn_id_coded)

  j <- cohort[events, on = c("anon_id","pat_enc_csn_id_coded"),
              allow.cartesian = TRUE, nomatch = NULL]
  if (nrow(j) == 0L) return(data.table())

  # TIMESTAMP_DIFF(..., DAY) truncates toward zero, so a lab 14.5d before the
  # culture yields -14 and is kept by BETWEEN -14 AND 0. trunc() matches.
  j[, dt_days := trunc(as.numeric(difftime(labtime, order_time_jittered_utc, units = "days")))]
  j <- j[dt_days >= -lookback_days & dt_days <= 0]
  if (nrow(j) == 0L) return(data.table())

  setorder(j, anon_id, pat_enc_csn_id_coded, order_proc_id_coded, labtime)

  # quantile(type = 7) approximates SQL APPROX_QUANTILES[OFFSET] to 2dp.
  agg <- j[, .(
    Q25    = round(quantile(value, 0.25, na.rm = TRUE, names = FALSE, type = 7), 2),
    Q75    = round(quantile(value, 0.75, na.rm = TRUE, names = FALSE, type = 7), 2),
    median = round(quantile(value, 0.50, na.rm = TRUE, names = FALSE, type = 7), 2),
    first  = round(value[1L],  2),
    last   = round(value[.N],  2)
  ), by = .(anon_id, pat_enc_csn_id_coded, order_proc_id_coded)]

  setnames(agg,
    old = c("Q25","Q75","median","first","last"),
    new = paste0(c("Q25_","Q75_","median_","first_","last_"), lab))
  agg[]
}

compute_window <- function(cohort, lookback_days, lab_streams) {
  message("  window: ", lookback_days, "d")
  parts <- lapply(names(lab_streams), function(lab) {
    message("    + ", lab)
    agg_one_lab(lab_streams[[lab]], cohort, lookback_days, lab)
  })
  parts <- Filter(function(d) nrow(d) > 0L, parts)
  if (length(parts) == 0L) return(data.table())
  out <- Reduce(function(a, b) merge(a, b, all = TRUE,
                                     by = c("anon_id","pat_enc_csn_id_coded",
                                            "order_proc_id_coded")),
                parts)
  out[, Period_Day := lookback_days]
  out[]
}


# =============================================================================
# DEIDENTIFICATION -- swap raw IDs for deid surrogates via the three ID maps;
# warn (don't drop) on unmapped rows. No date shift (output has no timestamps).
# =============================================================================

deidentify_output <- function(dt) {
  anon_map <- fread(file.path(MAPS_DIR, "anon_id_map.csv"),
                    colClasses = c(anon_id = "character"))
  enc_map  <- fread(file.path(MAPS_DIR, "enc_id_map.csv"),
                    colClasses = c(pat_enc_csn_id_coded = "character"))
  ord_map  <- fread(file.path(MAPS_DIR, "order_proc_id_map.csv"),
                    colClasses = c(order_proc_id_coded = "character"))

  anon_d <- setNames(anon_map$anon_id_deid,           anon_map$anon_id)
  enc_d  <- setNames(enc_map$pat_enc_csn_id_deid,     enc_map$pat_enc_csn_id_coded)
  ord_d  <- setNames(ord_map$order_proc_id_deid,      ord_map$order_proc_id_coded)

  dt <- copy(dt)
  dt[, anon_id              := unname(anon_d[as.character(anon_id)])]
  dt[, pat_enc_csn_id_coded := unname(enc_d [as.character(pat_enc_csn_id_coded)])]
  dt[, order_proc_id_coded  := unname(ord_d [as.character(order_proc_id_coded)])]

  for (col in c("anon_id","pat_enc_csn_id_coded","order_proc_id_coded")) {
    n_missing <- sum(is.na(dt[[col]]))
    if (n_missing > 0L) {
      message(sprintf("  WARN: %s rows have no %s mapping (left as NA)",
                      format(n_missing, big.mark=","), col))
    }
  }
  dt[]
}


# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  message("Reading cohort: ", COHORT_CSV)
  cohort <- read_cohort()
  message("Cohort rows: ", nrow(cohort))

  message("\nReading + unit-converting labs...")
  lab_streams <- list()
  for (lab in names(LAB_FILES)) {
    s <- read_lab(lab, LAB_FILES[[lab]])
    if (!is.null(s)) {
      lab_streams[[lab]] <- s
      message(sprintf("  [ok]   %-15s  %10s rows kept", lab, format(nrow(s), big.mark=",")))
    }
  }

  message("\nExtracting windowed features...")
  windows <- lapply(LOOKBACK_DAYS, function(w) compute_window(cohort, w, lab_streams))
  windows <- Filter(function(d) nrow(d) > 0L, windows)
  if (length(windows) == 0L) {
    message("No output rows."); return(invisible())
  }
  out <- unique(rbindlist(windows, use.names = TRUE, fill = TRUE))

  # Column order: keys first, then per-lab metrics grouped by lab.
  keys <- c("anon_id","pat_enc_csn_id_coded","order_proc_id_coded","Period_Day")
  ordered <- keys
  for (lab in names(LAB_FILES)) {
    for (stat in c("Q25","Q75","median","first","last")) {
      col <- paste0(stat, "_", lab)
      if (col %in% names(out)) ordered <- c(ordered, col)
    }
  }
  setcolorder(out, ordered)

  if (isTRUE(APPLY_DEID)) {
    message("\nDeidentifying with maps at ", MAPS_DIR, " ...")
    out <- deidentify_output(out)
  } else {
    message("\nWARNING: APPLY_DEID = FALSE -- output still contains raw ",
            "IDs. Do NOT share externally.")
  }

  fwrite(out, OUTPUT_CSV)
  message("\nWrote ", nrow(out), " rows -> ", OUTPUT_CSV)
}

main()
