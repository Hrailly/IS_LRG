# 患者成组重抽样：同一患者及共享标识下的全部记录始终一起抽取。
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
options(stringsAsFactors = FALSE)
set.seed(20261010)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
module <- dirname(normalizePath(script_path, winslash = "/"))
root <- dirname(module)
input_dir <- file.path(root, "数据")
log_dir <- file.path(dirname(root), "logs")
dir.create(log_dir, showWarnings = FALSE)
objects <- qs::qread(file.path(input_dir, "01_supplement_objects.qs"))
x <- objects$biological_model_matrix
m <- objects$metadata
stopifnot(identical(colnames(x), m$Accession))
design <- cbind(Intercept = 1, Disease = as.integer(m$Condition == "IS"), Followup37587 = m$Followup37587)
patient <- factor(m$Patient_linked)
units <- levels(patient)
unit_condition <- vapply(units, function(id) {
  status <- unique(as.character(m$Condition[patient == id]))
  stopifnot(length(status) == 1)
  status
}, character(1))
repetitions <- 1000L
counts <- matrix(0L, nrow = repetitions, ncol = length(units), dimnames = list(NULL, units))
for (b in seq_len(repetitions)) {
  for (condition in c("Control", "IS")) {
    indices <- which(unit_condition == condition)
    drawn <- sample(indices, length(indices), replace = TRUE)
    counts[b, ] <- counts[b, ] + tabulate(drawn, nbins = length(units))
  }
}
stopifnot(all(rowSums(counts[, unit_condition == "IS", drop = FALSE]) == sum(unit_condition == "IS")),
          all(rowSums(counts[, unit_condition == "Control", drop = FALSE]) == sum(unit_condition == "Control")))
aggregate_units <- function(z) rowsum(z, group = patient, reorder = FALSE)[units, , drop = FALSE]

# 预先计算每名患者的充分统计量，避免反复展开相同样本矩阵。
unit_xtx <- lapply(seq_len(ncol(design)), function(j) {
  sapply(seq_len(ncol(design)), function(k) as.numeric(aggregate_units(matrix(design[, j] * design[, k], ncol = 1))))
})
xtx <- lapply(unit_xtx, function(z) counts %*% z)
inverse_rows <- matrix(NA_real_, nrow = repetitions, ncol = 3)
for (b in seq_len(repetitions)) {
  matrix_b <- rbind(xtx[[1]][b, ], xtx[[2]][b, ], xtx[[3]][b, ])
  stopifnot(qr(matrix_b)$rank == 3)
  inverse_rows[b, ] <- solve(matrix_b)[2, ]
}
bootstrap_effects <- matrix(0, nrow = repetitions, ncol = nrow(x))
message("计算1000次患者成组重抽样；表达预处理固定，疾病及时间系数逐次重新估计")
for (j in seq_len(ncol(design))) {
  unit_xty <- aggregate_units(t(x) * design[, j])
  bootstrap_effects <- bootstrap_effects + (counts * inverse_rows[, j]) %*% unit_xty
  message("  充分统计量 ", j, "/3 完成")
}
point <- as.numeric((solve(crossprod(design), crossprod(design, t(x))))[2, ])
intervals <- apply(bootstrap_effects, 2, quantile, probs = c(.025, .975), names = FALSE)
centered <- abs(sweep(bootstrap_effects, 2, point, "-"))
p_raw <- (1 + colSums(sweep(centered, 2, abs(point), ">="))) / (repetitions + 1)
results <- data.frame(Feature = rownames(x), Model = "Linked_patient_cluster_bootstrap", Contrast = "IS - Control; within-GSE37587 follow-up adjusted",
                      Effect = point, CI_low = intervals[1, ], CI_high = intervals[2, ], P_raw = p_raw, P_BH = p.adjust(p_raw, "BH"),
                      Patients_under_linkage_assumption = length(units), IS_units = sum(unit_condition == "IS"), Control_units = sum(unit_condition == "Control"),
                      Bootstrap_repetitions = repetitions, Matrix_refitted_in_bootstrap = FALSE)
write.csv(results, file.path(root, "结果表", "主要意见1", "13_patient_cluster_bootstrap.csv"), row.names = FALSE, fileEncoding = "UTF-8")
qs::qsave(list(seed = 20261010, patient_draw_counts = counts, patient_conditions = unit_condition, effects = bootstrap_effects,
               point = point, results = results), file.path(input_dir, "02_patient_bootstrap.qs"), preset = "high")
check <- qs::qread(file.path(input_dir, "02_patient_bootstrap.qs"))
stopifnot(identical(check$results, results), identical(check$patient_draw_counts, counts))
write.csv(data.frame(Check = c("Patient cluster preservation", "Case-control stratification", "QS reload", "Preprocessing uncertainty"),
                     Result = c("PASS: one multiplicity per linked patient for all records", "PASS: 59 linked case units and 44 control units per draw", "PASS", "Conditional on the already fitted expression matrix; ComBat is not refitted")),
           file.path(root, "结果表", "主要意见1", "13_patient_resampling_checks.csv"), row.names = FALSE, fileEncoding = "UTF-8")
capture.output(sessionInfo(), file = file.path(log_dir, "M17_major1_bootstrap_sessionInfo.txt"))
message("患者成组重抽样完成")
