suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

source("experiments/appendix/runtime_benchmark/runtime_benchmark_helpers.R", encoding = "UTF-8")
configure_single_thread_math()
source_runtime_dependencies()

parse_flag_value <- function(arg, name) {
  sub(paste0("^--", name, "="), "", arg)
}

parse_runtime_args <- function(args) {
  cfg <- list(
    datasets = c("elec2", "aig"),
    methods = runtime_method_keys,
    windows = runtime_scale_windows,
    workers = default_runtime_workers,
    end_reps = runtime_end_to_end_reps,
    scale_reps = runtime_scale_reps,
    stage = "all",
    end_window = runtime_end_to_end_window,
    block_updates = runtime_scale_block_updates,
    end_timeout_sec = runtime_end_to_end_timeout_sec,
    scale_timeout_sec = runtime_scale_timeout_sec
  )

  for (arg in args) {
    if (grepl("^--datasets=", arg)) {
      value <- tolower(parse_flag_value(arg, "datasets"))
      cfg$datasets <- if (identical(value, "all")) {
        c("elec2", "aig")
      } else {
        trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      }
    } else if (grepl("^--methods=", arg)) {
      value <- tolower(parse_flag_value(arg, "methods"))
      cfg$methods <- if (identical(value, "all")) {
        runtime_method_keys
      } else {
        trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      }
    } else if (grepl("^--windows=", arg)) {
      cfg$windows <- as.integer(trimws(strsplit(parse_flag_value(arg, "windows"), ",", fixed = TRUE)[[1]]))
    } else if (grepl("^--workers=", arg)) {
      cfg$workers <- as.integer(parse_flag_value(arg, "workers"))
    } else if (grepl("^--end-reps=", arg)) {
      cfg$end_reps <- as.integer(parse_flag_value(arg, "end-reps"))
    } else if (grepl("^--scale-reps=", arg)) {
      cfg$scale_reps <- as.integer(parse_flag_value(arg, "scale-reps"))
    } else if (grepl("^--stage=", arg)) {
      cfg$stage <- tolower(parse_flag_value(arg, "stage"))
    } else if (grepl("^--end-window=", arg)) {
      cfg$end_window <- as.integer(parse_flag_value(arg, "end-window"))
    } else if (grepl("^--block-updates=", arg)) {
      cfg$block_updates <- as.integer(parse_flag_value(arg, "block-updates"))
    } else if (grepl("^--end-timeout-sec=", arg)) {
      cfg$end_timeout_sec <- as.numeric(parse_flag_value(arg, "end-timeout-sec"))
    } else if (grepl("^--scale-timeout-sec=", arg)) {
      cfg$scale_timeout_sec <- as.numeric(parse_flag_value(arg, "scale-timeout-sec"))
    } else if (nzchar(arg)) {
      stop("Unknown argument: ", arg)
    }
  }

  stopifnot(length(cfg$datasets) >= 1L)
  stopifnot(length(cfg$methods) >= 1L)
  stopifnot(length(cfg$windows) >= 1L)
  stopifnot(cfg$workers >= 1L)
  stopifnot(cfg$end_reps >= 1L)
  stopifnot(cfg$scale_reps >= 1L)
  stopifnot(cfg$end_window >= 1L)
  stopifnot(cfg$block_updates >= 1L)
  stopifnot(cfg$stage %in% c("all", "end_to_end", "update_only"))
  stopifnot(all(cfg$methods %in% runtime_method_keys))
  cfg
}

apply_runtime_overrides <- function(cfg) {
  runtime_end_to_end_window <<- as.integer(cfg$end_window)
  runtime_end_to_end_reps <<- as.integer(cfg$end_reps)
  runtime_scale_reps <<- as.integer(cfg$scale_reps)
  runtime_scale_block_updates <<- as.integer(cfg$block_updates)
  runtime_end_to_end_timeout_sec <<- as.numeric(cfg$end_timeout_sec)
  runtime_scale_timeout_sec <<- as.numeric(cfg$scale_timeout_sec)
  invisible(NULL)
}

setup_cluster <- function(workers, workdir, cfg) {
  workers <- max(1L, as.integer(workers))
  if (workers == 1L) {
    return(NULL)
  }

  cl <- parallel::makePSOCKcluster(workers)
  parallel::clusterCall(cl, function(workdir, cfg) {
    setwd(workdir)
    source("experiments/appendix/runtime_benchmark/runtime_benchmark_helpers.R", encoding = "UTF-8")
    configure_single_thread_math()
    source_runtime_dependencies()
    runtime_end_to_end_window <<- as.integer(cfg$end_window)
    runtime_end_to_end_reps <<- as.integer(cfg$end_reps)
    runtime_scale_reps <<- as.integer(cfg$scale_reps)
    runtime_scale_block_updates <<- as.integer(cfg$block_updates)
    runtime_end_to_end_timeout_sec <<- as.numeric(cfg$end_timeout_sec)
    runtime_scale_timeout_sec <<- as.numeric(cfg$scale_timeout_sec)
    invisible(NULL)
  }, workdir = workdir, cfg = cfg)
  cl
}

run_jobs <- function(jobs, worker_fun, workers, workdir, cfg) {
  if (length(jobs) == 0L) {
    return(list())
  }

  cl <- setup_cluster(workers = workers, workdir = workdir, cfg = cfg)
  if (is.null(cl)) {
    return(lapply(jobs, worker_fun))
  }

  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::parLapplyLB(cl, jobs, worker_fun)
}

make_end_jobs <- function(cfg) {
  grid <- expand.grid(
    dataset_slug = cfg$datasets,
    method_key = cfg$methods,
    repetition = seq_len(cfg$end_reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  lapply(seq_len(nrow(grid)), function(i) {
    list(
      dataset_slug = grid$dataset_slug[[i]],
      method_key = grid$method_key[[i]],
      repetition = grid$repetition[[i]]
    )
  })
}

make_update_jobs <- function(cfg) {
  grid <- expand.grid(
    dataset_slug = cfg$datasets,
    method_key = cfg$methods,
    window = as.integer(cfg$windows),
    repetition = seq_len(cfg$scale_reps),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  lapply(seq_len(nrow(grid)), function(i) {
    list(
      dataset_slug = grid$dataset_slug[[i]],
      method_key = grid$method_key[[i]],
      window = as.integer(grid$window[[i]]),
      repetition = grid$repetition[[i]]
    )
  })
}

end_job_worker <- function(job) {
  run_end_to_end_once(
    dataset_slug = job$dataset_slug,
    method_key = job$method_key,
    repetition = job$repetition
  )
}

update_job_worker <- function(job) {
  run_update_only_once(
    dataset_slug = job$dataset_slug,
    method_key = job$method_key,
    window = job$window,
    repetition = job$repetition
  )
}

cfg <- parse_runtime_args(commandArgs(trailingOnly = TRUE))
apply_runtime_overrides(cfg)

workdir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

merge_runtime_raw <- function(existing_df, new_df) {
  if (nrow(new_df) == 0L) {
    return(existing_df)
  }
  if (nrow(existing_df) == 0L) {
    return(new_df)
  }

  key_df <- new_df %>%
    select(any_of(c("benchmark", "dataset_slug", "method_key", "window", "repetition"))) %>%
    distinct()

  existing_keep <- dplyr::anti_join(
    existing_df,
    key_df,
    by = intersect(names(existing_df), names(key_df))
  )
  bind_rows(existing_keep, new_df)
}

end_to_end_raw <- read_runtime_csv_or_empty(end_to_end_raw_path)
update_only_raw <- read_runtime_csv_or_empty(update_only_raw_path)

if (cfg$stage %in% c("all", "end_to_end")) {
  message(
    sprintf(
      "Running end-to-end benchmark: datasets={%s}, methods=%d, reps=%d, w=%d, workers=%d",
      paste(cfg$datasets, collapse = ","),
      length(cfg$methods),
      cfg$end_reps,
      cfg$end_window,
      cfg$workers
    )
  )
  end_jobs <- make_end_jobs(cfg)
  end_res <- run_jobs(
    jobs = end_jobs,
    worker_fun = end_job_worker,
    workers = cfg$workers,
    workdir = workdir,
    cfg = cfg
  )
  end_to_end_raw <- merge_runtime_raw(end_to_end_raw, bind_rows(end_res))
}

if (cfg$stage %in% c("all", "update_only")) {
  message(
    sprintf(
      "Running update-only benchmark: datasets={%s}, methods=%d, windows={%s}, reps=%d, block_updates=%d, workers=%d",
      paste(cfg$datasets, collapse = ","),
      length(cfg$methods),
      paste(cfg$windows, collapse = ","),
      cfg$scale_reps,
      cfg$block_updates,
      cfg$workers
    )
  )
  update_jobs <- make_update_jobs(cfg)
  update_res <- run_jobs(
    jobs = update_jobs,
    worker_fun = update_job_worker,
    workers = cfg$workers,
    workdir = workdir,
    cfg = cfg
  )
  update_only_raw <- merge_runtime_raw(update_only_raw, bind_rows(update_res))
}

outputs <- save_runtime_outputs(
  end_to_end_raw = end_to_end_raw,
  update_only_raw = update_only_raw
)

cat("Saved end-to-end raw results to:", end_to_end_raw_path, "\n")
cat("Saved update-only raw results to:", update_only_raw_path, "\n")
cat("Saved end-to-end summary to:", end_to_end_summary_path, "\n")
cat("Saved update-only summary to:", update_only_summary_path, "\n")
cat("Saved main runtime table to:", main_table_path, "\n")
cat("Saved timeout summary to:", timeout_summary_path, "\n")
cat("Saved main LaTeX table to:", main_table_tex_path, "\n")
cat("Saved scalability LaTeX table to:", scale_table_tex_path, "\n")
cat("Saved protocol notes to:", protocol_note_path, "and", protocol_tex_path, "\n")
if (nrow(outputs$update_summary) > 0L) {
  cat("Saved scalability plot to:", scalability_plot_pdf_path, "and", scalability_plot_png_path, "\n")
}

print(outputs$main_table)
print(outputs$timeout_summary)
