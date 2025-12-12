suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(stringr)
  library(glmnet)
  library(survival)
  library(Matrix)
  library(doParallel)
  library(foreach)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) stop("Usage: feature_selection_by_transition.R IN OUT INCLUDE EXCLUDE ALPHA NFOLDS MAXCAT SEED MAXFEATS")

IN       <- args[[1]]
OUT      <- args[[2]]
INCLUDE  <- args[[3]]     # supports "LIST::/abs/path/to/shard.list"
EXCLUDE  <- args[[4]]
ALPHA    <- as.numeric(args[[5]])
NFOLDS   <- as.integer(args[[6]])
MAXCAT   <- as.integer(args[[7]])
SEED     <- as.integer(args[[8]])
MAXFEATS <- as.integer(args[[9]])

REQUIRED <- c("person_id","start","stop","event")

logmsg <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(...))); flush.console() }
is_gcs <- function(path) grepl("^gs://", path)
ship <- function(local_path, subpath = basename(local_path), base_out = OUT) {
  if (is_gcs(base_out)) {
    dest <- file.path(base_out, subpath)
    cmd <- sprintf("gsutil -q cp %s %s", shQuote(local_path), shQuote(dest))
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  } else {
    dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
    file.copy(local_path, file.path(base_out, subpath), overwrite = TRUE)
  }
}
bt <- function(x) sprintf("`%s`", x)

get_threads <- function() {
  nt <- as.integer(Sys.getenv("N_THREADS", unset = "1"))
  if (is.na(nt) || nt < 1) nt <- 1L
  nt
}

# -------- helpers --------
make_numeric_safe <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  if (inherits(x, "Date")) return(as.numeric(x))
  if (inherits(x, "POSIXct")) return(as.numeric(x))
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(nx <- as.numeric(x))
  if (sum(!is.na(nx)) < length(nx) * 0.5) {
    y <- gsub("[^0-9eE+\\-.]", "", as.character(x))
    suppressWarnings(nx2 <- as.numeric(y))
    if (sum(!is.na(nx2)) >= sum(!is.na(nx))) nx <- nx2
  }
  nx
}

normalize_event01 <- function(x) {
  if (is.numeric(x)) return(as.integer(x > 0))
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) x <- as.character(x)
  xx <- tolower(trimws(as.character(x)))
  ones  <- c("1","true","t","yes","y","case","event","fail","failed","death","died")
  zeros <- c("0","false","f","no","n","control","censor","censored","alive")
  out <- rep(NA_integer_, length(xx))
  out[xx %in% ones]  <- 1L
  out[xx %in% zeros] <- 0L
  na <- is.na(out)
  if (any(na)) {
    suppressWarnings(num <- as.numeric(xx[na]))
    out[na] <- as.integer(!is.na(num) & num > 0)
  }
  out[is.na(out)] <- 0L
  out
}

coerce_numeric <- function(DT, cols, max_levels = 50L) {
  for (cn in cols) {
    v <- DT[[cn]]
    if (is.numeric(v)) next
    if (is.logical(v)) { DT[[cn]] <- as.numeric(v); next }
    if (is.factor(v)) {
      if (length(levels(v)) <= max_levels) DT[[cn]] <- as.numeric(v) - 1
      else {
        suppressWarnings(num <- as.numeric(as.character(v)))
        if (anyNA(num) && sum(!is.na(num)) < length(num) * 0.9) DT[[cn]] <- as.integer(as.integer(factor(v)) - 1)
        else DT[[cn]] <- num
      }
      next
    }
    if (is.character(v)) {
      suppressWarnings(num <- as.numeric(v))
      if (sum(!is.na(num)) >= length(num) * 0.9) DT[[cn]] <- num
      else {
        uq <- data.table::uniqueN(v)
        if (uq <= max_levels) DT[[cn]] <- as.integer(as.integer(factor(v)) - 1)
        else DT[[cn]] <- as.integer(factor(v)) - 1
      }
      next
    }
    suppressWarnings(num <- as.numeric(v))
    if (!all(is.na(num))) DT[[cn]] <- num
  }
  invisible(DT)
}

# to create an "empty but valid" result object when glmnet is skipped/failed
empty_result <- function(reason = "no_fit") {
  metrics   <- data.table(lambda = numeric(0), cve = numeric(0), cvsd = numeric(0), nzero = integer(0))
  coef_min  <- data.table(feature = character(0), coef = numeric(0))
  coef_1se  <- data.table(feature = character(0), coef = numeric(0))
  summary_txt <- sprintf("No glmnet fit performed. Reason: %s\n", reason)
  list(metrics = metrics, coef_min = coef_min, coef_1se = coef_1se, summary_txt = summary_txt)
}
# --------------------------------

# ---------- config ----------
set.seed(SEED)
ROWS_CAP <- as.integer(Sys.getenv("ROWS_CAP", unset="0"))          # 0 = no cap
FORCE_PER_TRANSITION <- as.integer(Sys.getenv("PER_TRANSITION", unset="1")) # 1=on if 'transition' exists
MIN_EVENTS <- as.integer(Sys.getenv("MIN_EVENTS", unset = "20"))   # minimum events required to attempt glmnet

logmsg("Input: %s", IN); logmsg("Out:   %s", OUT)
logmsg("include: %s | exclude: %s", INCLUDE, EXCLUDE)
logmsg("alpha=%.3f nfolds=%d max_cat=%d seed=%d max_feats=%d", ALPHA, NFOLDS, MAXCAT, SEED, MAXFEATS)

ds <- open_dataset(IN)
schema_cols <- names(ds$schema)
logmsg("RAW schema columns: %d", length(schema_cols))
logmsg("RAW head cols: %s", paste(utils::head(schema_cols, 50), collapse=", "))

missing_required <- setdiff(REQUIRED, schema_cols)
if (length(missing_required)) stop("Missing required columns: ", paste(missing_required, collapse=", "))
HAS_TRANSITION <- "transition" %in% schema_cols

# Determine candidate features (regex or list) without loading rows
determine_keep <- function(schema_cols) {
  keep <- character(0)
  if (startsWith(INCLUDE, "LIST::")) {
    list_path <- sub("^LIST::", "", INCLUDE)
    if (!file.exists(list_path)) stop("INCLUDE list file not found: ", list_path)
    keep <- readLines(list_path, warn = FALSE)
    keep <- intersect(keep, schema_cols)
  } else {
    keep <- schema_cols[stringr::str_detect(schema_cols, INCLUDE)]
    if (nzchar(EXCLUDE)) keep <- keep[!stringr::str_detect(schema_cols, EXCLUDE)]
  }
  keep <- setdiff(keep, REQUIRED)
  if (!length(keep)) stop("No candidate feature columns after filters.")
  if (!is.na(MAXFEATS) && MAXFEATS > 0L && length(keep) > MAXFEATS) {
    logmsg("Capping features from %d to MAXFEATS=%d", length(keep), MAXFEATS)
    keep <- keep[seq_len(MAXFEATS)]
  }
  keep
}

# ---------------- fitting ----------------
run_one_fit <- function(ds, keep, filter_transition = NULL, out_base = OUT) {
  cols_to_read <- unique(c(REQUIRED, if (HAS_TRANSITION) "transition", keep))
  s <- ds %>% dplyr::select(dplyr::all_of(cols_to_read))
  if (!is.null(filter_transition)) s <- s %>% dplyr::filter(transition == filter_transition)

  DT <- s %>% dplyr::collect() %>% as.data.table()
  n_in <- nrow(DT)
  if (!n_in) {
    logmsg("no rows in slice")
    return(NULL)
  }

  # Coerce required fields
  DT[, start := make_numeric_safe(start)]
  DT[, stop  := make_numeric_safe(stop)]
  DT[, event := normalize_event01(event)]

  # Drop non-finite outcome rows
  bad_nonfinite <- !(is.finite(DT$start) & is.finite(DT$stop) & is.finite(DT$event))
  n_bad_nonfinite <- sum(bad_nonfinite, na.rm=TRUE)
  if (n_bad_nonfinite) {
    logmsg("drop non-finite rows (start/stop/event): %d", n_bad_nonfinite)
    DT <- DT[!bad_nonfinite]
  }

  # Enforce strictly positive interval (stop > start)
  bad_interval <- !(DT$stop > DT$start)
  n_bad_interval <- sum(bad_interval, na.rm=TRUE)
  if (n_bad_interval) {
    logmsg("drop rows with non-positive interval (stop <= start): %d", n_bad_interval)
    DT <- DT[!bad_interval]
  }

  # Ensure event in {0,1}
  bad_event <- !(DT$event %in% c(0L, 1L))
  n_bad_event <- sum(bad_event, na.rm=TRUE)
  if (n_bad_event) {
    logmsg("drop rows with invalid event (not 0/1) after normalization: %d", n_bad_event)
    DT <- DT[!bad_event]
  }

  if (!nrow(DT)) {
    logmsg("no rows left after cleaning (was %d)", n_in)
    return(NULL)
  }

  # QC summary
  qc <- data.table(
    n_in = n_in,
    dropped_nonfinite = n_bad_nonfinite,
    dropped_nonpos_interval = n_bad_interval,
    dropped_bad_event = n_bad_event,
    n_out = nrow(DT)
  )
  qc_path <- file.path(tempdir(), "qc_cleaning.csv")
  fwrite(qc, qc_path)
  ship(qc_path, "qc_cleaning.csv", base_out = out_base)

  # Row cap
  if (ROWS_CAP > 0L && nrow(DT) > ROWS_CAP) {
    set.seed(SEED)
    DT <- DT[sample.int(nrow(DT), ROWS_CAP)]
    logmsg("Applied ROWS_CAP=%d, sampled rows: %d", ROWS_CAP, nrow(DT))
  }

  # Coerce feature columns numeric and clean them
  feat_cols <- setdiff(names(DT), c(REQUIRED, "transition"))
  DT <- coerce_numeric(DT, feat_cols, max_levels = MAXCAT)

  # Replace any non-finite feature values with 0
  for (cn in feat_cols) {
    v <- DT[[cn]]
    bad <- !is.finite(v)
    if (any(bad)) {
      data.table::set(DT, which(bad), cn, 0)
    }
  }

  # Drop zero-variance features *before* building design matrix
  if (length(feat_cols) > 0L) {
    nzv <- vapply(
      feat_cols,
      function(cn) {
        z <- DT[[cn]]
        stats::sd(z, na.rm = TRUE) > 0
      },
      logical(1L)
    )
    if (!all(nzv)) {
      logmsg("dropping %d zero-variance feature columns", sum(!nzv))
      feat_cols <- feat_cols[nzv]
    }
  }

  # Check if any features left
  if (!length(feat_cols)) {
    logmsg("no usable predictors after cleaning; skipping glmnet")
    return(empty_result("no_usable_predictors"))
  }

  # Design matrix and outcome
  form <- as.formula(paste("~", paste(bt(feat_cols), collapse="+"), "-1"))
  logmsg("design matrix: n=%d, p=%d (sparse)", nrow(DT), length(feat_cols))
  X <- Matrix::sparse.model.matrix(form, data = DT)
  y <- with(DT, Surv(start, stop, event))

  # Check number of events
  n_events <- sum(DT$event == 1L)
  if (n_events < MIN_EVENTS) {
    logmsg("Too few events (%d < MIN_EVENTS=%d); skipping glmnet", n_events, MIN_EVENTS)
    return(empty_result(sprintf("too_few_events_%d", n_events)))
  }

  # Parallel setup
  n_threads <- get_threads()
  if (n_threads > 1L) {
    cl <- parallel::makeCluster(n_threads)
    doParallel::registerDoParallel(cl)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    logmsg("doParallel enabled with %d workers", n_threads)
  } else {
    foreach::registerDoSEQ()
    logmsg("Running single-threaded")
  }

  # glmnet with tryCatch
  logmsg("glmnet(cv) starting…")
  t0 <- proc.time()
  fit <- tryCatch(
    {
      cv.glmnet(
        x = X, y = y, family = "cox",
        alpha = ALPHA, nfolds = NFOLDS, parallel = (n_threads > 1L),
        maxit = 100000, thresh = 1e-07, standardize = TRUE, type.measure = "deviance"
      )
    },
    error = function(e) {
      logmsg("cv.glmnet error: %s", conditionMessage(e))
      return(NULL)
    }
  )
  t1 <- proc.time()

  if (is.null(fit)) {
    logmsg("glmnet(cv) failed; returning empty result")
    return(empty_result("glmnet_error"))
  }

  logmsg("glmnet(cv) finished in %.1f sec", (t1 - t0)[["elapsed"]])

  lambda_seq <- fit$lambda
  metrics <- data.table(lambda = lambda_seq, cve = fit$cvm, cvsd = fit$cvsd, nzero = fit$nzero)

  coef_to_df <- function(m) {
    idx <- which(m != 0)
    if (!length(idx)) data.table(feature=character(0), coef=numeric(0))
    else data.table(feature = rownames(m)[idx], coef = as.numeric(m[idx]))
  }
  beta_min <- as.matrix(coef(fit, s = "lambda.min"))
  beta_1se <- as.matrix(coef(fit, s = "lambda.1se"))
  coef_min <- coef_to_df(beta_min)
  coef_1se <- coef_to_df(beta_1se)

  summary_txt <- paste0(
    "n_obs: ", nrow(DT), "\n",
    "p_feats: ", length(feat_cols), "\n",
    "n_events: ", n_events, "\n",
    "alpha: ", ALPHA, "\n",
    "nfolds: ", NFOLDS, "\n",
    "lambda.min: ", fit$lambda.min, "\n",
    "lambda.1se: ", fit$lambda.1se, "\n",
    "nonzero@min: ", nrow(coef_min), "\n",
    "nonzero@1se: ", nrow(coef_1se), "\n"
  )

  list(metrics=metrics, coef_min=coef_min, coef_1se=coef_1se, summary_txt=summary_txt)
}

write_outputs <- function(res, base_out) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  stg <- file.path(tempdir(), paste0("featsel_", stamp))
  dir.create(stg, showWarnings = FALSE, recursive = TRUE)
  f_metrics <- file.path(stg, "cv_metrics.csv")
  f_coefmin <- file.path(stg, "coef_lambda_min.csv")
  f_coef1se <- file.path(stg, "coef_lambda_1se.csv")
  f_sel1se  <- file.path(stg, "selected_features_1se.txt")
  f_sum     <- file.path(stg, "run_summary.txt")
  fwrite(res$metrics, f_metrics)
  fwrite(res$coef_min, f_coefmin)
  fwrite(res$coef_1se, f_coef1se)
  writeLines(res$coef_1se$feature, f_sel1se)
  writeLines(res$summary_txt, f_sum)
  ship(f_metrics, "cv_metrics.csv", base_out)
  ship(f_coefmin, "coef_lambda_min.csv", base_out)
  ship(f_coef1se, "coef_lambda_1se.csv", base_out)
  ship(f_sel1se,  "selected_features_1se.txt", base_out)
  ship(f_sum,     "run_summary.txt", base_out)
  logmsg("Wrote outputs under %s", base_out)
}

# ---- main ----
keep <- determine_keep(schema_cols)
logmsg("candidate feature columns: %d", length(keep))
logmsg("sample features: %s", paste(utils::head(keep, 30), collapse=", "))

if (HAS_TRANSITION && FORCE_PER_TRANSITION == 1L) {
  # collect unique transitions and drop bogus ones
  tr_vals <- (ds %>% dplyr::select(transition) %>% dplyr::collect() %>% unique())$transition
  bad_tr <- c(NA, "", "0", 0, "none", "NA", "na", "null", "NULL")
  tr_vals <- tr_vals[!(is.na(tr_vals) | tr_vals %in% bad_tr)]
  if (!length(tr_vals)) stop("transition column has no usable values.")

  logmsg("Per-transition mode. transitions: %s", paste(tr_vals, collapse=", "))

  for (tr in tr_vals) {
    sub_out <- file.path(OUT, sprintf("transition_%s", as.character(tr)))
    logmsg("=== transition=%s ===", as.character(tr))
    res <- run_one_fit(ds, keep, filter_transition = tr, out_base = sub_out)
    if (!is.null(res)) write_outputs(res, base_out = sub_out)
    rm(res); gc()
  }
} else {
  logmsg("Single-slice mode.")
  res <- run_one_fit(ds, keep, filter_transition = NULL, out_base = OUT)
  if (is.null(res)) stop("No data to fit.")
  write_outputs(res, base_out = OUT)
}
logmsg("DONE.")
