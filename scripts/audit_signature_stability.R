#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(MutationalPatterns))

out_dir <- Sys.getenv("RCC_AUDIT_OUT", unset = "outputs")
load(file.path(out_dir, "pass_snv_signature_results.RData"))

cosine <- function(observed, reconstructed) {
  colSums(observed * reconstructed) /
    (sqrt(colSums(observed ^ 2)) * sqrt(colSums(reconstructed ^ 2)))
}

relative_contribution <- function(x) {
  sweep(x, 2, colSums(x), "/") * 100
}

extract_fit <- function(x) {
  if (!is.null(x$fit_res)) x$fit_res else x
}

summarize_fit <- function(name, fit_object, observed) {
  fit_result <- extract_fit(fit_object)
  rel <- relative_contribution(fit_result$contribution)
  get_signature <- function(signature) {
    if (signature %in% rownames(rel)) rel[signature, ] else rep(0, ncol(rel))
  }
  data.frame(
    sample = colnames(observed),
    fit = name,
    SBS2_percent = as.numeric(get_signature("SBS2")),
    SBS13_percent = as.numeric(get_signature("SBS13")),
    APOBEC_percent = as.numeric(get_signature("SBS2") + get_signature("SBS13")),
    reconstruction_cosine = as.numeric(cosine(observed, fit_result$reconstructed)),
    residual_fraction = as.numeric(colSums(abs(observed - fit_result$reconstructed)) /
      pmax(colSums(observed), 1)),
    active_signatures = as.numeric(colSums(fit_result$contribution > 0))
  )
}

set.seed(20260716)
cosmic <- get_known_signatures(muttype = "snv", source = "COSMIC_v3.2")
restricted_names <- intersect(
  c("SBS1", "SBS5", "SBS22", "SBS40", "SBS2", "SBS13", "SBS6", "SBS15"),
  colnames(cosmic)
)

restricted_regular <- fit_to_signatures(mut_mat, cosmic[, restricted_names, drop = FALSE])
full_regular <- fit_to_signatures(mut_mat, cosmic)
full_strict <- fit_to_signatures_strict(mut_mat, cosmic, max_delta = 0.004, method = "backwards")

fit_summary <- rbind(
  summarize_fit("restricted_8_regular", restricted_regular, mut_mat),
  summarize_fit("full_COSMIC_regular", full_regular, mut_mat),
  summarize_fit("full_COSMIC_strict", full_strict, mut_mat)
)
write.csv(fit_summary, file.path(out_dir, "signature_fit_stability_summary.csv"), row.names = FALSE)

# Bootstrap the full reference with strict refitting. Row names encode sample and
# replicate as <sample>_<replicate>; sample IDs in this cohort contain no '_'.
boot <- fit_to_signatures_bootstrapped(
  mut_mat,
  cosmic,
  n_boots = 100,
  max_delta = 0.004,
  method = "strict",
  verbose = FALSE
)
boot_samples <- sub("_[0-9]+$", "", rownames(boot))
boot_total <- pmax(rowSums(boot), 1)
boot_sbs2 <- if ("SBS2" %in% colnames(boot)) boot[, "SBS2"] / boot_total * 100 else rep(0, nrow(boot))
boot_sbs13 <- if ("SBS13" %in% colnames(boot)) boot[, "SBS13"] / boot_total * 100 else rep(0, nrow(boot))
boot_frame <- data.frame(
  sample = boot_samples,
  replicate = as.integer(sub("^.*_([0-9]+)$", "\\1", rownames(boot))),
  SBS2_percent = as.numeric(boot_sbs2),
  SBS13_percent = as.numeric(boot_sbs13),
  APOBEC_percent = as.numeric(boot_sbs2 + boot_sbs13)
)
write.csv(boot_frame, file.path(out_dir, "signature_full_strict_bootstrap_replicates.csv"), row.names = FALSE)

boot_summary <- do.call(rbind, lapply(split(boot_frame, boot_frame$sample), function(d) {
  data.frame(
    sample = d$sample[1],
    SBS2_median = median(d$SBS2_percent),
    SBS2_q025 = unname(quantile(d$SBS2_percent, 0.025)),
    SBS2_q975 = unname(quantile(d$SBS2_percent, 0.975)),
    SBS2_detection_fraction = mean(d$SBS2_percent > 0),
    APOBEC_median = median(d$APOBEC_percent),
    APOBEC_q025 = unname(quantile(d$APOBEC_percent, 0.025)),
    APOBEC_q975 = unname(quantile(d$APOBEC_percent, 0.975)),
    APOBEC_detection_fraction = mean(d$APOBEC_percent > 0)
  )
}))
write.csv(boot_summary, file.path(out_dir, "signature_full_strict_bootstrap_summary.csv"), row.names = FALSE)

save(
  restricted_regular,
  full_regular,
  full_strict,
  boot,
  file = file.path(out_dir, "signature_stability_results.RData")
)

cat("Completed signature stability audit for", ncol(mut_mat), "samples.\n")
cat("Median reconstruction cosine:\n")
print(tapply(fit_summary$reconstruction_cosine, fit_summary$fit, median))
cat("Strict-bootstrap SBS2 detection fractions for EIS01/LEN08/LEN15:\n")
print(boot_summary[boot_summary$sample %in% c("EIS01", "LEN08", "LEN15"), ])
