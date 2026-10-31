#!/usr/bin/env Rscript
###############################################################################
# vitals.R -- R port of 05_vitals.sql.
#
# Output microbiology_Vitals.csv: one row per
# (anon_id, pat_enc_csn_id_coded, order_proc_id_coded), from the 48h before
# each culture order. Per vital, 5 columns: Q25/Q75/median/first/last. Vitals:
# heartrate, resprate, temp, sysbp, diasbp. (weight is filtered but never
# aggregated in the SQL, so it is omitted here too.)
#
# Run: Rscript vitals.R  -- edit the CONFIGURATION block first.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
})


# =============================================================================
# CONFIGURATION
# =============================================================================

RAW_DIR    <- "/PATH/TO/raw/edw/2024-12-31"
COHORT_CSV <- "/PATH/TO/microbiology_cultures_cohort.csv"
OUTPUT_CSV <- "microbiology_Vitals.csv"

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

# ADAPT: one CSV per vital. hr/resp/temp share a single-value schema; bp holds
# bpSystolic + bpDiastolic columns. Fix any filename that differs.
HRTEMP_FILES <- list(
  heartrate = file.path(RAW_DIR, "heartrate.csv"),
  resprate  = file.path(RAW_DIR, "resprate.csv"),
  temp      = file.path(RAW_DIR, "tempt.csv")
)
BP_FILE <- file.path(RAW_DIR, "bp.csv")

# Lookback window in hours (matches SQL).
LOOKBACK_HOURS <- 48L


# =============================================================================
# SQL -> MGB SCHEMA MAPPING
#
# Stanford BigQuery                   MGB EDW CSV
# ---------------------------------   ----------------------------------------
# c.anon_id                           patientID (raw) / anon_id (deid)
# c.pat_enc_csn_id_coded              patientEncounterID (raw) / ..._coded (deid)
# c.order_proc_id_coded               orderProcedureID
# c.order_time_jittered_utc           cultureCollDTS
# v.recorded_time_jittered_utc        recordedDTS
# v.row_disp_name                     implicit from FILE NAME (one vital / file)
# v.meas_value (numeric string)       resultN (hr/resp/temp)
# v.meas_value ("120/80" for BP)      bpSystolic + bpDiastolic (TWO columns)
#
# JOIN: the SQL joins vitals on v.anon_id = c.anon_id ONLY (not encounter), so
# a vital within 48h attaches to a culture regardless of encounter. We keep
# this; the output row's encounter is the cohort's, not the vital's.
#
# ADAPT checklist -- verify against your MGB export:
# (A) One file per vital, so no row_disp_name filter is needed. If a file mixes
#     types, add a FLOWSHEET_ID_FILTERS entry.
# (B) BP: SQL stores "SYS/DIA"; MGB stores numeric bpSystolic + bpDiastolic.
#     read_bp() handles it -- only touch it if the column names differ.
# (C) Temp unit: SQL's `>= 90` implies Fahrenheit. Confirm MGB temp is in F;
#     if Celsius, convert upstream or change TEMP_MIN (~32).
#
# Deliberate matches to the SQL (do NOT "fix"):
# * hr/resp/temp/bp pulled patient-only (no encounter constraint).
# * weight filtered but not aggregated => no weight columns emitted.
#
# MGB addition beyond the SQL: resultN preferred with as.numeric(measureTXT)
# fallback (USE via read_hrtemp; SQL uses the single meas_value column).
# =============================================================================

# (A) ADAPT: allowed flowsheet IDs if a file mixes measurement types.
# e.g. FLOWSHEET_ID_FILTERS$heartrate <- c("301070")
FLOWSHEET_ID_FILTERS <- list(
)

# (C) Value-range filters (verbatim from the SQL WHERE clause).
HEARTRATE_MIN <- 30
RESPRATE_MIN  <- 4
RESPRATE_MAX  <- 60
TEMP_MIN      <- 90    # Fahrenheit; see checklist item (C)
SYSBP_MIN     <- 40
DIASBP_MIN    <- 30

# ADAPT: column names in hr/resp/temp CSVs.
V_ID_COL        <- "patientID"
V_TIME_COL      <- "recordedDTS"
V_VALUE_N_COL   <- "resultN"
V_VALUE_TXT_COL <- "measureTXT"        # optional fallback
V_FMID_COL      <- "flowsheetMeasureID"

# ADAPT: column names in bp.csv.
BP_ID_COL   <- "patientID"
BP_TIME_COL <- "recordedDTS"
BP_SYS_COL  <- "bpSystolic"
BP_DIA_COL  <- "bpDiastolic"


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

.numeric_with_fallback <- function(primary, txt) {
  v <- suppressWarnings(as.numeric(primary))
  if (is.null(txt)) return(v)
  fb <- suppressWarnings(as.numeric(txt))
  ifelse(is.na(v), fb, v)
}

read_hrtemp <- function(label, path, lo = NULL, hi = NULL) {
  if (!file.exists(path)) {
    message("  [skip] ", label, ": file not found at ", path); return(NULL)
  }
  dt <- fread(path)
  required <- c(V_ID_COL, V_TIME_COL, V_VALUE_N_COL)
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(path, " missing columns: ", paste(missing, collapse=", "))
  }

  ids <- FLOWSHEET_ID_FILTERS[[label]]
  if (!is.null(ids) && V_FMID_COL %in% names(dt)) {
    dt <- dt[as.character(get(V_FMID_COL)) %chin% as.character(ids)]
  }

  txt <- if (V_VALUE_TXT_COL %in% names(dt)) dt[[V_VALUE_TXT_COL]] else NULL
  # SQL filters on the RAW cast, then rounds in projection. We match: filter on
  # raw, round only what we emit.
  val_raw <- .numeric_with_fallback(dt[[V_VALUE_N_COL]], txt)

  keep <- !is.na(val_raw)
  if (!is.null(lo)) keep <- keep & val_raw >= lo
  if (!is.null(hi)) keep <- keep & val_raw <= hi
  dt <- dt[keep]
  if (nrow(dt) == 0L) return(NULL)

  out <- dt[, .(
    anon_id   = as.character(get(V_ID_COL)),
    vitaltime = as.POSIXct(get(V_TIME_COL)),
    value     = round(val_raw[keep], 2)
  )]
  out <- out[!is.na(vitaltime)]
  out[]
}

read_bp <- function(path) {
  if (!file.exists(path)) {
    message("  [skip] bp: file not found at ", path); return(NULL)
  }
  dt <- fread(path)
  required <- c(BP_ID_COL, BP_TIME_COL, BP_SYS_COL, BP_DIA_COL)
  missing  <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(path, " missing columns: ", paste(missing, collapse=", "))
  }

  ids <- FLOWSHEET_ID_FILTERS[["bp"]]
  if (!is.null(ids) && V_FMID_COL %in% names(dt)) {
    dt <- dt[as.character(get(V_FMID_COL)) %chin% as.character(ids)]
  }

  # Filter on RAW cast (SQL semantics), round on emit.
  sys_raw <- suppressWarnings(as.numeric(dt[[BP_SYS_COL]]))
  dia_raw <- suppressWarnings(as.numeric(dt[[BP_DIA_COL]]))
  t       <- as.POSIXct(dt[[BP_TIME_COL]])
  pid     <- as.character(dt[[BP_ID_COL]])

  m_sys <- !is.na(sys_raw) & sys_raw >= SYSBP_MIN & !is.na(t)
  m_dia <- !is.na(dia_raw) & dia_raw >= DIASBP_MIN & !is.na(t)

  list(
    sysbp  = data.table(anon_id = pid[m_sys], vitaltime = t[m_sys],
                        value = round(sys_raw[m_sys], 2)),
    diasbp = data.table(anon_id = pid[m_dia], vitaltime = t[m_dia],
                        value = round(dia_raw[m_dia], 2))
  )
}


# =============================================================================
# 48h WINDOWED AGGREGATION -- per vital stream, join to cohort on anon_id alone
# (patient-only, per SQL), keep TIMESTAMP_DIFF(vital, cohort, HOUR) in [-48, 0],
# then aggregate per (anon_id, pat_enc_csn_id_coded, order_proc_id_coded).
# =============================================================================

agg_one_vital <- function(events, cohort, label) {
  if (is.null(events) || nrow(events) == 0L) return(data.table())
  setkey(events, anon_id)
  setkey(cohort, anon_id)

  j <- cohort[events, on = "anon_id", allow.cartesian = TRUE, nomatch = NULL]
  if (nrow(j) == 0L) return(data.table())

  # TIMESTAMP_DIFF(..., HOUR) truncates toward zero, so a vital 48.5h before the
  # culture yields -48 and is kept by BETWEEN -48 AND 0. trunc() matches.
  j[, dt_hr := trunc(as.numeric(difftime(vitaltime, order_time_jittered_utc, units = "hours")))]
  j <- j[dt_hr >= -LOOKBACK_HOURS & dt_hr <= 0]
  if (nrow(j) == 0L) return(data.table())

  setorder(j, anon_id, pat_enc_csn_id_coded, order_proc_id_coded, vitaltime)

  agg <- j[, .(
    Q25    = round(quantile(value, 0.25, na.rm = TRUE, names = FALSE, type = 7), 2),
    Q75    = round(quantile(value, 0.75, na.rm = TRUE, names = FALSE, type = 7), 2),
    median = round(quantile(value, 0.50, na.rm = TRUE, names = FALSE, type = 7), 2),
    first  = round(value[1L],  2),
    last   = round(value[.N],  2)
  ), by = .(anon_id, pat_enc_csn_id_coded, order_proc_id_coded)]

  setnames(agg,
    old = c("Q25","Q75","median","first","last"),
    new = paste0(c("Q25_","Q75_","median_","first_","last_"), label))
  agg[]
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

  message("\nReading hr / resp / temp ...")
  streams <- list()
  streams$heartrate <- read_hrtemp("heartrate", HRTEMP_FILES$heartrate,
                                   lo = HEARTRATE_MIN)
  streams$resprate  <- read_hrtemp("resprate",  HRTEMP_FILES$resprate,
                                   lo = RESPRATE_MIN, hi = RESPRATE_MAX)
  streams$temp      <- read_hrtemp("temp",      HRTEMP_FILES$temp,
                                   lo = TEMP_MIN)

  message("\nReading bp ...")
  bp <- read_bp(BP_FILE)
  if (!is.null(bp)) {
    streams$sysbp  <- bp$sysbp
    streams$diasbp <- bp$diasbp
  }

  message("\nExtracting 48h features ...")
  parts <- list()
  for (label in c("heartrate","resprate","temp","sysbp","diasbp")) {
    s <- streams[[label]]
    if (is.null(s) || nrow(s) == 0L) {
      message("  [skip] ", label, ": empty stream"); next
    }
    message("  + ", label)
    agg <- agg_one_vital(s, cohort, label)
    if (nrow(agg) > 0L) parts[[label]] <- agg
  }

  if (length(parts) == 0L) {
    message("No output rows."); return(invisible())
  }
  out <- Reduce(function(a, b) merge(a, b, all = TRUE,
                                     by = c("anon_id","pat_enc_csn_id_coded",
                                            "order_proc_id_coded")),
                parts)

  # Column order: keys first, then per-vital metrics grouped by vital.
  keys <- c("anon_id","pat_enc_csn_id_coded","order_proc_id_coded")
  ordered <- keys
  for (label in c("heartrate","resprate","temp","sysbp","diasbp")) {
    for (stat in c("Q25","Q75","median","first","last")) {
      col <- paste0(stat, "_", label)
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
