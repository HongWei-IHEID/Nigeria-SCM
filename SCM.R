rm(list = ls())

# ===================== Packages =====================
library(tidyverse)
library(lubridate)
library(Synth)
library(future.apply)
library(openxlsx)

# ===================== Global settings =====================
treat_ccy  <- "NGN"
T0         <- as.Date("2021-02-01")
predictors <- c("cpi", "m2_growth", "reserves_gdp", "policy_rate", "kaopen")
nullfile   <- if (.Platform$OS.type == "windows") "NUL" else "/dev/null"

# Placebo pre-fit screen multiplier (not used in this paper, so set a very high number, most research will set it to 5 or 10)
pre_rmspe_mult <- 100

# ===================== Two specs in ONE loop =====================
specs <- list(
  level = list(
    data_path = "scm_panel_abokifx_full_sample.csv",
    y_var = "Y",
    seed = 88,
    out_root = "scm_results",
    suffix_prefix = "_abokifx_level_full",
    y_label = "Valuation Gap",
    is_demeaned = FALSE
  ),
  demeaning = list(
    data_path = "scm_panel_abokifx_demeaning_full_sample.csv",
    y_var = "Y_tilde",
    seed = 88,
    out_root = "scm_results",
    suffix_prefix = "_abokifx_demeaning_full",
    y_label = "Gap Deviation",
    is_demeaned = TRUE
  )
)

out_root <- unique(vapply(specs, \(x) x$out_root, character(1)))
if (length(out_root) != 1) stop("All specs must share the same out_root for a single-loop setup.")
out_root <- out_root[[1]]
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# ===================== Run =====================
for (tag in names(specs)) {

  cfg <- specs[[tag]]
  set.seed(cfg$seed)

  cat("\n\n============================================\n")
  cat("Running SCM spec:", tag, "\n")
  cat("Data:", cfg$data_path, "\n")
  cat("Outcome:", cfg$y_var, "\n")
  cat("============================================\n\n")

  suffix <- paste0(cfg$suffix_prefix, "_", tag)
  out_dir <- file.path(out_root, tag)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  plan(sequential)

  # --------------------- Data prep ---------------------
  panel <- read_csv(cfg$data_path, show_col_types = FALSE) %>%
    mutate(month = as.Date(month), time_id = as.integer(format(month, "%Y%m"))) %>%
    distinct(ccy, month, .keep_all = TRUE)

  # Required columns for port decomposition
  req_cols <- c("effective_rate", "crypto_shadow_rate")
  miss_cols <- setdiff(req_cols, names(panel))
  if (length(miss_cols) > 0) {
    stop("Missing columns in panel: ", paste(miss_cols, collapse = ", "),
         ". Please re-export with effective_rate and crypto_shadow_rate.")
  }

  month_seq <- seq.Date(as.Date("2019-01-01"), as.Date("2022-12-01"), by = "month")
  time_all  <- as.integer(format(month_seq, "%Y%m"))
  time_pre  <- time_all[month_seq < T0]

  pre_months <- month_seq[month_seq < T0]
  lag_ids <- as.integer(format(pre_months[unique(pmin(length(pre_months),
                                                      round(seq(1, length(pre_months), length.out = 8))))], "%Y%m"))

  # Special predictors: pre-mean + 8 evenly spaced lags (match the outcome variable)
  special_preds <- c(
    list(list(cfg$y_var, time_pre, "mean")),
    lapply(lag_ids, \(t) list(cfg$y_var, t, "mean"))
  )

  y_var <- cfg$y_var

  # Balanced panel screen (full window + no missing predictors/ports/outcome)
  valid_units <- panel %>%
    filter(time_id %in% time_all) %>%
    group_by(ccy) %>%
    summarise(
      ok = n_distinct(time_id) == length(time_all) &&
        !any(is.na(.data[[y_var]])) &&
        !any(is.na(effective_rate)) &&
        !any(is.na(crypto_shadow_rate)) &&
        !any((time_id %in% time_pre) & if_any(all_of(predictors), is.na)),
      .groups = "drop"
    ) %>%
    filter(ok) %>%
    pull(ccy)

  panel_scm <- panel %>%
    filter(ccy %in% valid_units, time_id %in% time_all) %>%
    arrange(ccy, month) %>%
    mutate(unit_id = as.integer(factor(ccy, levels = sort(unique(ccy)))))

  treat_id  <- panel_scm %>% filter(ccy == treat_ccy) %>% pull(unit_id) %>% unique()
  donor_ids <- panel_scm %>% filter(ccy != treat_ccy) %>% pull(unit_id) %>% unique()

  # Explicit donor exclusions, as an extension, not really used
  #exclude_ccy <- c("ARS", "CLP", "RUB")
  #exclude_ids <- panel_scm %>% filter(ccy %in% exclude_ccy) %>% pull(unit_id) %>% unique()
  #donor_ids   <- setdiff(donor_ids, exclude_ids)

  # ===================== SCM wrapper (ports standardized) =====================
  run_scm <- function(treat_uid, ctrl_uids, data_full, do_ports = FALSE) {

    dp <- dataprep(
      foo = as.data.frame(data_full),
      predictors = predictors,
      predictors.op = "mean",
      time.predictors.prior = time_pre,
      special.predictors = special_preds,
      dependent = cfg$y_var,
      unit.variable = "unit_id",
      time.variable = "time_id",
      treatment.identifier = treat_uid,
      controls.identifier = ctrl_uids,
      time.optimize.ssr = time_pre,
      time.plot = time_all,
      unit.names.variable = "ccy"
    )

    capture.output({ sc <- synth(dp) }, file = nullfile)

    y1  <- as.numeric(dp$Y1plot)
    y0  <- as.numeric(dp$Y0plot %*% sc$solution.w)
    gap <- y1 - y0
    is_pre <- dp$tag$time.plot < as.integer(format(T0, "%Y%m"))

    out <- list(
      w_mat = sc$solution.w, v_raw = sc$solution.v,
      X1 = dp$X1, X0 = dp$X0,
      y1 = y1, y0 = y0, gap = gap,
      pre_rmspe = sqrt(mean(gap[is_pre]^2)),
      post_rmspe = sqrt(mean(gap[!is_pre]^2))
    )

    if (isTRUE(do_ports)) {

      get_plot_series <- function(varname) {
        time_vec <- dp$tag$time.plot
        tu <- dp$tag$treatment.identifier
        cu <- dp$tag$controls.identifier

        dfv <- data_full %>%
          filter(time_id %in% time_vec, unit_id %in% c(tu, cu)) %>%
          select(time_id, unit_id, value = all_of(varname)) %>%
          distinct(time_id, unit_id, .keep_all = TRUE)

        wide <- dfv %>%
          pivot_wider(names_from = unit_id, values_from = value) %>%
          arrange(time_id)

        treat_vec <- as.numeric(wide[[as.character(tu)]])
        ctrl_mat  <- as.matrix(wide %>% select(all_of(as.character(cu))))

        list(treat = treat_vec, ctrl = ctrl_mat, time_id = wide$time_id)
      }

      E <- get_plot_series("effective_rate")
      S <- get_plot_series("crypto_shadow_rate")

      if (any(is.na(E$treat)) || any(is.na(S$treat)) || any(is.na(E$ctrl)) || any(is.na(S$ctrl))) {
        stop("Port decomposition failed due to NA in effective_rate/crypto_shadow_rate after screening.")
      }
      if (any(E$treat <= 0) || any(S$treat <= 0) || any((E$ctrl %*% sc$solution.w) <= 0) || any((S$ctrl %*% sc$solution.w) <= 0)) {
        stop("Port decomposition requires positive effective_rate and crypto_shadow_rate for logs.")
      }

      E1 <- E$treat
      E0 <- as.numeric(E$ctrl %*% sc$solution.w)
      S1 <- S$treat
      S0 <- as.numeric(S$ctrl %*% sc$solution.w)

      dE <- log(E1) - log(E0)
      dS <- log(S1) - log(S0)

      if (!isTRUE(cfg$is_demeaned)) {
        out$E1 <- E1; out$E0 <- E0
        out$S1 <- S1; out$S0 <- S0
        out$alpha_hat <- gap
        out$dE_port <- dE
        out$minus_dS_port <- -dS
        out$alpha_recon <- dE - dS
        out$alpha_diff_check <- max(abs(out$alpha_recon - gap), na.rm = TRUE)
      } else {
        alpha_level <- dE - dS
        mu_dE_pre  <- mean(dE[is_pre], na.rm = TRUE)
        mu_mds_pre <- mean((-dS)[is_pre], na.rm = TRUE)
        dE_tilde <- dE - mu_dE_pre
        minus_dS_tilde <- (-dS) - mu_mds_pre

        out$E1 <- E1; out$E0 <- E0
        out$S1 <- S1; out$S0 <- S0
        out$alpha_hat <- gap
        out$dE_port <- dE_tilde
        out$minus_dS_port <- minus_dS_tilde
        out$alpha_recon <- dE_tilde + minus_dS_tilde
        out$alpha_diff_check <- max(abs(out$alpha_recon - gap), na.rm = TRUE)
      }
    }

    out
  }

  # ===================== Main SCM + placebos =====================
  main_res <- run_scm(treat_id, donor_ids, panel_scm, do_ports = TRUE)
  main_res$ratio <- main_res$post_rmspe / main_res$pre_rmspe

  plan(multisession, workers = max(1, parallel::detectCores() - 1))

  placebo_list <- future_lapply(donor_ids, function(uid) {
    res <- run_scm(uid, setdiff(donor_ids, uid), panel_scm, do_ports = FALSE)
    ccy_name <- panel_scm %>% filter(unit_id == uid) %>% pull(ccy) %>% unique()
    tibble(
      ccy = ccy_name,
      pre_rmspe = res$pre_rmspe,
      post_rmspe = res$post_rmspe,
      ratio = res$post_rmspe / res$pre_rmspe,
      month = list(month_seq),
      gap = list(res$gap)
    )
  }, future.seed = TRUE) %>% bind_rows()

  plan(sequential)

  # Pre-fit screen for placebos
  main_pre <- main_res$pre_rmspe
  thr_pre  <- pre_rmspe_mult * max(main_pre, 1e-8)

  placebo_keep <- placebo_list %>%
    filter(!is.na(pre_rmspe), pre_rmspe <= thr_pre)

  dropped_n <- nrow(placebo_list) - nrow(placebo_keep)

  cat("\n[Placebo pre-fit screen]\n")
  cat("  Spec:", tag, "\n")
  cat("  Treated pre-RMSPE =", signif(main_pre, 6), "\n")
  cat("  Threshold (mult =", pre_rmspe_mult, "):", signif(thr_pre, 6), "\n")
  cat("  Kept placebos:", nrow(placebo_keep), "/", nrow(placebo_list),
      " (dropped ", dropped_n, ")\n", sep = "")

  k <- sum(placebo_keep$ratio >= main_res$ratio, na.rm = TRUE)
  J_eff <- sum(!is.na(placebo_keep$ratio))
  p_val <- (k + 1) / (J_eff + 1)

  # ===================== Tables =====================
  main_stat <- tibble(
    ccy = treat_ccy,
    pre_rmspe = main_res$pre_rmspe,
    post_rmspe = main_res$post_rmspe,
    ratio = main_res$ratio,
    p_value_right_tail = p_val
  )

  inference_tbl <- bind_rows(
    main_stat,
    placebo_list %>% select(ccy, pre_rmspe, post_rmspe, ratio) %>%
      mutate(p_value_right_tail = NA_real_)
  )

  id_map <- panel_scm %>% distinct(unit_id, ccy)

  w_tbl <- tibble(
    unit_id = suppressWarnings(as.integer(rownames(main_res$w_mat))),
    weight  = as.numeric(main_res$w_mat)
  ) %>%
    left_join(id_map, by = "unit_id") %>%
    select(ccy, weight) %>%
    arrange(desc(weight))

  bal_tbl <- tibble(
    predictor = rownames(main_res$X1),
    treated = as.numeric(main_res$X1),
    synthetic = as.numeric(main_res$X0 %*% main_res$w_mat),
    donor_mean = as.numeric(rowMeans(main_res$X0))
  )

  v_vec <- if (is.null(main_res$v_raw) || length(main_res$v_raw) != nrow(main_res$X1)) {
    rep(NA_real_, nrow(main_res$X1))
  } else {
    as.numeric(main_res$v_raw)
  }
  v_tbl <- tibble(predictor = rownames(main_res$X1), v = v_vec)

  fit_df <- tibble(
    month = month_seq,
    treated = main_res$y1,
    synthetic = main_res$y0,
    gap = main_res$gap
  )

  gap_placebo <- placebo_keep %>%
    select(ccy, month, gap) %>%
    unnest(c(month, gap))

  # ===================== Random donor samples (only serves as an extension, not actually reported in the paper) =====================
  B <- 1  #should be 1000
  J <- length(donor_ids)
  donor_draw_size <- max(5, ceiling(0.7 * J))

  plan(multisession, workers = max(1, parallel::detectCores() - 1))

  rand_list <- future_lapply(1:B, function(b) {
    ctrl_sub <- sample(donor_ids, size = donor_draw_size, replace = FALSE)

    res <- tryCatch(
      run_scm(treat_id, ctrl_sub, panel_scm, do_ports = FALSE),
      error = function(e) NULL
    )

    if (is.null(res)) {
      return(tibble(rep = b, month = month_seq, gap = NA_real_))
    }

    tibble(
      rep = b,
      month = month_seq,
      gap = res$gap
    )
  }, future.seed = TRUE) %>% bind_rows()

  plan(sequential)

  rand_list_ok <- rand_list %>%
    group_by(rep) %>%
    filter(!all(is.na(gap))) %>%
    ungroup()

  cat("\nRandom donor samples finished. Valid reps = ",
      n_distinct(rand_list_ok$rep), " / ", B, "\n", sep = "")

  # ===================== Port decomposition data (NGN) =====================
  if (!is.null(main_res$dE_port) && !is.null(main_res$minus_dS_port)) {

    decomp_ngn <- tibble(
      month = month_seq,
      alpha_hat = main_res$alpha_hat,
      dE_port = main_res$dE_port,
      minus_dS_port = main_res$minus_dS_port,
      alpha_recon = main_res$alpha_recon,
      eff_treated = main_res$E1,
      eff_synth = main_res$E0,
      crypto_treated = main_res$S1,
      crypto_synth = main_res$S0
    )

    cat("\n[Port decomposition check] max |alpha_recon - alpha_hat| = ",
        signif(main_res$alpha_diff_check, 6), "\n", sep = "")

    es_long <- decomp_ngn %>%
      transmute(
        month,
        `Effective (treated)` = log(eff_treated),
        `Effective (synthetic)` = log(eff_synth),
        `Crypto shadow (treated)` = log(crypto_treated),
        `Crypto shadow (synthetic)` = log(crypto_synth)
      ) %>%
      pivot_longer(-month, names_to = "series", values_to = "value")
  }

  # ===================== Figures =====================
  p1 <- ggplot(fit_df, aes(x = month)) +
    geom_line(aes(y = treated, color = "Treated"), linewidth = 0.8) +
    geom_line(aes(y = synthetic, color = "Synthetic"), linetype = "longdash", linewidth = 0.8) +
    geom_vline(xintercept = T0, linetype = "dotted") +
    scale_color_manual(values = c("Treated" = "grey40", "Synthetic" = "red")) +
    labs(x = NULL, y = cfg$y_label, color = NULL) +
    theme_bw() + theme(legend.position = "bottom")

  top_gap_country <- "RUB"
  gap_placebo_normal <- gap_placebo %>% filter(ccy != top_gap_country)
  gap_placebo_high   <- gap_placebo %>% filter(ccy == top_gap_country)
  last_point <- gap_placebo_high %>% filter(month == max(month))

  p2 <- ggplot() +
    geom_line(data = gap_placebo_normal,
              aes(x = month, y = gap, group = ccy, color = "Placebo Countries"),
              alpha = 0.4, linewidth = 0.5) +
    geom_line(data = gap_placebo_high,
              aes(x = month, y = gap, color = top_gap_country),
              linewidth = 1) +
    geom_line(data = fit_df,
              aes(x = month, y = gap, color = "Nigeria"),
              linewidth = 1.2) +
    geom_vline(xintercept = T0, linetype = "dotted") +
    geom_text(data = last_point, aes(x = month, y = gap, label = top_gap_country),
              color = "blue", hjust = -0.2, vjust = 0.5, fontface = "bold", size = 3) +
    scale_color_manual(values = setNames(c("grey70", "blue", "red"),
                                         c("Placebo Countries", top_gap_country, "Nigeria"))) +
    labs(x = NULL, y = "Prediction Error", color = NULL) +
    theme_bw() +
    theme(legend.position = "bottom") +
    scale_x_date(expand = expansion(mult = c(0.05, 0.1)))

  all_ratios <- bind_rows(
    main_stat %>% select(ccy, ratio),
    placebo_keep %>% select(ccy, ratio)
  )

  p3 <- ggplot(all_ratios, aes(x = reorder(ccy, ratio), y = ratio, fill = ccy == treat_ccy)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
    labs(x = NULL, y = "Ratio", fill = "Target") +
    theme_bw()

  p4 <- ggplot() +
    geom_line(
      data = rand_list_ok,
      aes(x = month, y = gap, group = rep),
      color = "grey70", alpha = 0.5, linewidth = 0.4
    ) +
    geom_line(
      data = fit_df,
      aes(x = month, y = gap),
      color = "red", linewidth = 1
    ) +
    geom_vline(xintercept = T0, linetype = "dotted") +
    labs(x = NULL, y = cfg$y_label) +
    theme_bw()

  if (exists("decomp_ngn")) {
    p5 <- ggplot(decomp_ngn, aes(x = month)) +
      geom_line(aes(y = alpha_hat, color = "SCM effect"), linewidth = 1.0) +
      geom_line(aes(y = dE_port, color = "Effective port"), linewidth = 0.8) +
      geom_line(aes(y = minus_dS_port, color = "Crypto port"), linewidth = 0.8) +
      geom_vline(xintercept = T0, linetype = "dotted") +
      labs(x = NULL, y = "Log points", color = NULL) +
      theme_bw() + theme(legend.position = "bottom")

    p6 <- ggplot(es_long, aes(x = month, y = value, linetype = series, color = series)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = T0, linetype = "dotted") +
      labs(x = NULL, y = "log(LC per USD)", color = NULL, linetype = NULL) +
      theme_bw() + theme(legend.position = "bottom")
  }

  print(p1); print(p2); print(p3); print(p4)
  if (exists("decomp_ngn")) { print(p5); print(p6) }

  ggsave(file.path(out_dir, paste0("scm_fit_path", suffix, ".png")), p1, width = 8.5, height = 4.8, dpi = 300)
  ggsave(file.path(out_dir, paste0("prediction_error_placebos", suffix, ".png")), p2, width = 8.5, height = 4.8, dpi = 300)
  ggsave(file.path(out_dir, paste0("rmspe_ratio_distribution", suffix, ".png")), p3, width = 8.5, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0("gap_random_donor_samples", suffix, ".png")), p4, width = 8.5, height = 4.8, dpi = 300)
  if (exists("decomp_ngn")) {
    ggsave(file.path(out_dir, paste0("port_decomposition", suffix, ".png")), p5, width = 8.5, height = 4.8, dpi = 300)
    ggsave(file.path(out_dir, paste0("effective_crypto_paths", suffix, ".png")), p6, width = 8.5, height = 4.8, dpi = 300)
  }

  # ===================== Export to Excel =====================
  wb <- createWorkbook()
  addWorksheet(wb, "weights");            writeData(wb, "weights", w_tbl)
  addWorksheet(wb, "predictor_balance");  writeData(wb, "predictor_balance", bal_tbl)
  addWorksheet(wb, "v_weights");          writeData(wb, "v_weights", v_tbl)
  addWorksheet(wb, "fit_series");         writeData(wb, "fit_series", fit_df)
  addWorksheet(wb, "inference_rmspe");    writeData(wb, "inference_rmspe", inference_tbl)

  if (exists("decomp_ngn")) {
    addWorksheet(wb, "port_decomp_ngn"); writeData(wb, "port_decomp_ngn", decomp_ngn)
    addWorksheet(wb, "port_paths_long"); writeData(wb, "port_paths_long", es_long)
  }

  saveWorkbook(wb, file.path(out_dir, paste0("scm_results", suffix, ".xlsx")), overwrite = TRUE)

  # ===================== Console output =====================
  cat("\n================ SCM Results ================\n")
  cat("\nSpec:", tag, "\n")
  cat("\n[Main model] Treated (NGN):\n")
  print(main_stat)
  cat("\n[Weights] Top donors:\n")
  print(w_tbl %>% head(10))
  cat("\n[Balance] Predictor balance:\n")
  print(bal_tbl)
  cat("\nResults saved to:", out_dir, "\n")
}
