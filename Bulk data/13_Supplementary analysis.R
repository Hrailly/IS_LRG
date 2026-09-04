# 审稿回复补充分析：队列内差异、跨队列荟萃分析与独立验证

set.seed(20260818)
options(stringsAsFactors = FALSE)

dir.create("output", showWarnings = FALSE)
dir.create("checkpoints", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

message("[1/8] 读取并清理输入数据")

read_expression <- function(matrix_file, sample_file) {
  expr <- read.csv(matrix_file, row.names = 1, check.names = FALSE)
  meta <- read.delim(sample_file, check.names = FALSE)
  rownames(expr) <- trimws(rownames(expr))
  expr <- expr[!grepl("///", rownames(expr), fixed = TRUE), , drop = FALSE]
  expr <- limma::avereps(as.matrix(expr))
  expr <- expr[, meta$Accession, drop = FALSE]
  list(expr = expr, meta = meta)
}

gse16561 <- read_expression("input/gse16561_matrix.csv", "input/gse16561_samples.txt")
gse22255 <- read_expression("input/gse22255_matrix.csv", "input/gse22255_samples.txt")
gse37587 <- read_expression("input/gse37587_matrix.csv", "input/gse37587_samples.txt")

lactylation_raw <- read.delim("input/lactylation_genes.txt", check.names = FALSE)$ID
lactylation_trimmed <- trimws(lactylation_raw)
symbol_corrections <- c(
  "HISTIH2BH" = "HIST1H2BH",
  "SPR14" = "SRP14",
  "PKM2" = "PKM"
)
lactylation_standardized <- lactylation_trimmed
replace_index <- lactylation_standardized %in% names(symbol_corrections)
lactylation_standardized[replace_index] <- unname(symbol_corrections[lactylation_standardized[replace_index]])
lactylation <- unique(lactylation_standardized)
original_24 <- read.delim("input/original_24_genes.txt", check.names = FALSE)[[1]]
original_six <- read.delim("input/original_six_genes.txt", check.names = FALSE)[[1]]

input_audit <- data.frame(
  cohort = c("GSE16561", "GSE22255", "GSE37587"),
  platform = c("GPL6883", "GPL570", "GPL6883"),
  cases = c(sum(gse16561$meta$condition == "IS"), sum(gse22255$meta$condition == "IS"), sum(gse37587$meta$condition == "IS")),
  controls = c(sum(gse16561$meta$condition == "Control"), sum(gse22255$meta$condition == "Control"), sum(gse37587$meta$condition == "Control")),
  genes_after_mapping = c(nrow(gse16561$expr), nrow(gse22255$expr), nrow(gse37587$expr)),
  minimum = c(min(gse16561$expr), min(gse22255$expr), min(gse37587$expr)),
  maximum = c(max(gse16561$expr), max(gse22255$expr), max(gse37587$expr))
)
write.csv(input_audit, "output/00_input_audit.csv", row.names = FALSE, fileEncoding = "UTF-8")

lactylation_audit <- data.frame(
  raw_rows = length(lactylation_raw),
  unique_after_trim = length(unique(lactylation_trimmed)),
  unique_after_symbol_standardization = length(lactylation),
  whitespace_affected_rows = sum(lactylation_raw != trimws(lactylation_raw)),
  duplicated_after_trim = sum(duplicated(lactylation_trimmed)),
  corrected_legacy_or_typographic_symbols = sum(replace_index)
)
write.csv(lactylation_audit, "output/00_lactylation_list_audit.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.table(data.frame(Gene = lactylation), "output/00_lactylation_genes_cleaned.tsv", sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

qs::qsave(list(gse16561 = gse16561, gse22255 = gse22255, gse37587 = gse37587, lactylation = lactylation), "checkpoints/01_clean_inputs.qs", preset = "high")

message("[2/8] 在两个病例-对照队列内分别估计差异")

run_limma <- function(cohort) {
  condition <- factor(cohort$meta$condition, levels = c("Control", "IS"))
  design <- stats::model.matrix(~ 0 + condition)
  colnames(design) <- c("Control", "IS")
  fit <- limma::lmFit(cohort$expr, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(IS - Control, levels = design))
  fit <- limma::eBayes(fit)
  out <- limma::topTable(fit, number = Inf, adjust.method = "BH", sort.by = "none")
  out$Gene <- rownames(out)
  out
}

de16561 <- run_limma(gse16561)
de22255 <- run_limma(gse22255)
write.csv(de16561, "output/01_GSE16561_limma.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(de22255, "output/02_GSE22255_limma.csv", row.names = FALSE, fileEncoding = "UTF-8")
qs::qsave(list(GSE16561 = de16561, GSE22255 = de22255), "checkpoints/02_cohort_limma.qs", preset = "high")

message("[3/8] 计算跨平台标准化效应与固定效应荟萃结果")

hedges_effect <- function(expr, meta) {
  case_values <- expr[, meta$Accession[meta$condition == "IS"], drop = FALSE]
  control_values <- expr[, meta$Accession[meta$condition == "Control"], drop = FALSE]
  n1 <- ncol(case_values)
  n0 <- ncol(control_values)
  m1 <- rowMeans(case_values)
  m0 <- rowMeans(control_values)
  s1 <- apply(case_values, 1, stats::sd)
  s0 <- apply(control_values, 1, stats::sd)
  pooled <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / (n1 + n0 - 2))
  d <- (m1 - m0) / pooled
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  g <- correction * d
  variance <- (n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0 - 2))
  data.frame(Gene = rownames(expr), g = g, variance = variance, stringsAsFactors = FALSE)
}

common_genes <- intersect(rownames(gse16561$expr), rownames(gse22255$expr))
effect16561 <- hedges_effect(gse16561$expr[common_genes, , drop = FALSE], gse16561$meta)
effect22255 <- hedges_effect(gse22255$expr[common_genes, , drop = FALSE], gse22255$meta)
names(effect16561)[2:3] <- c("g_GSE16561", "var_GSE16561")
names(effect22255)[2:3] <- c("g_GSE22255", "var_GSE22255")
meta_de <- merge(effect16561, effect22255, by = "Gene")
meta_de$w1 <- 1 / meta_de$var_GSE16561
meta_de$w2 <- 1 / meta_de$var_GSE22255
meta_de$g_meta <- (meta_de$w1 * meta_de$g_GSE16561 + meta_de$w2 * meta_de$g_GSE22255) / (meta_de$w1 + meta_de$w2)
meta_de$se_meta <- sqrt(1 / (meta_de$w1 + meta_de$w2))
meta_de$z_meta <- meta_de$g_meta / meta_de$se_meta
meta_de$p_meta <- 2 * stats::pnorm(-abs(meta_de$z_meta))
meta_de$fdr_meta <- stats::p.adjust(meta_de$p_meta, method = "BH")
meta_de$ci_low <- meta_de$g_meta - 1.96 * meta_de$se_meta
meta_de$ci_high <- meta_de$g_meta + 1.96 * meta_de$se_meta
meta_de$Q <- meta_de$w1 * (meta_de$g_GSE16561 - meta_de$g_meta)^2 + meta_de$w2 * (meta_de$g_GSE22255 - meta_de$g_meta)^2
meta_de$I2 <- pmax(0, (meta_de$Q - 1) / meta_de$Q) * 100
meta_de$direction_consistent <- sign(meta_de$g_GSE16561) == sign(meta_de$g_GSE22255)
meta_de$meta_significant <- meta_de$fdr_meta < 0.05 & meta_de$direction_consistent
meta_de <- meta_de[order(meta_de$fdr_meta, -abs(meta_de$g_meta)), ]

write.csv(meta_de, "output/03_meta_analysis_all_genes.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(meta_de[meta_de$Gene %in% lactylation, ], "output/04_lactylation_meta_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(meta_de[meta_de$Gene %in% original_24, ], "output/05_original_24_meta_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(meta_de[meta_de$Gene %in% original_six, ], "output/06_original_six_meta_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
qs::qsave(meta_de, "checkpoints/03_meta_analysis.qs", preset = "high")

message("[4/8] 生成联合分层表达分布与PCA图")

close_devices <- function() {
  while (grDevices::dev.cur() > 1) grDevices::dev.off()
}

expr_before <- cbind(gse16561$expr[common_genes, ], gse22255$expr[common_genes, ])
combined_meta <- rbind(
  data.frame(Accession = gse16561$meta$Accession, condition = gse16561$meta$condition, cohort = "GSE16561"),
  data.frame(Accession = gse22255$meta$Accession, condition = gse22255$meta$condition, cohort = "GSE22255")
)
combined_meta <- combined_meta[match(colnames(expr_before), combined_meta$Accession), ]
combat_mod <- stats::model.matrix(~ condition, data = combined_meta)
expr_after <- sva::ComBat(dat = expr_before, batch = combined_meta$cohort, mod = combat_mod, par.prior = TRUE, prior.plots = FALSE)

plot_expression_density <- function(expr, meta, file, title) {
  sampled_genes <- rownames(expr)[seq_len(min(1500, nrow(expr)))]
  long <- data.frame(
    expression = as.vector(expr[sampled_genes, , drop = FALSE]),
    group = rep(paste(meta$cohort, meta$condition, sep = " | "), each = length(sampled_genes))
  )
  close_devices()
  grDevices::pdf(file, width = 7.2, height = 5.4, onefile = FALSE)
  print(
    ggplot2::ggplot(long, ggplot2::aes(x = expression, color = group)) +
      ggplot2::geom_density(linewidth = 0.9, adjust = 1.1) +
      ggplot2::labs(title = title, x = "Expression", y = "Density", color = "Cohort | Group") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", color = "black", hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold", color = "black"),
        axis.text = ggplot2::element_text(color = "black"),
        legend.title = ggplot2::element_text(face = "bold", color = "black")
      )
  )
  grDevices::dev.off()
}

plot_pca <- function(expr, meta, file, title) {
  variances <- apply(expr, 1, stats::var)
  selected <- names(sort(variances, decreasing = TRUE))[seq_len(min(1000, length(variances)))]
  pca <- stats::prcomp(t(expr[selected, , drop = FALSE]), scale. = TRUE)
  variance_percent <- summary(pca)$importance[2, 1:2] * 100
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], cohort = meta$cohort, condition = meta$condition)
  close_devices()
  grDevices::pdf(file, width = 6.5, height = 5.4, onefile = FALSE)
  print(
    ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = condition, shape = cohort)) +
      ggplot2::geom_point(size = 2.8, alpha = 0.85) +
      ggplot2::scale_color_manual(values = c(Control = "#3C8DBC", IS = "#C83E4D")) +
      ggplot2::labs(
        title = title,
        x = sprintf("PC1 (%.1f%%)", variance_percent[1]),
        y = sprintf("PC2 (%.1f%%)", variance_percent[2]),
        color = "Group",
        shape = "Cohort"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", color = "black", hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold", color = "black"),
        axis.text = ggplot2::element_text(color = "black"),
        legend.title = ggplot2::element_text(face = "bold", color = "black")
      )
  )
  grDevices::dev.off()
}

plot_expression_density(expr_before, combined_meta, "output/01_Expression_Before.pdf", "Expression Before Harmonization")
plot_expression_density(expr_after, combined_meta, "output/02_Expression_After.pdf", "Expression After Harmonization")
plot_pca(expr_before, combined_meta, "output/03_PCA_Before.pdf", "PCA Before Harmonization")
plot_pca(expr_after, combined_meta, "output/04_PCA_After.pdf", "PCA After Harmonization")
qs::qsave(list(before = expr_before, after = expr_after, metadata = combined_meta), "checkpoints/04_visualization_matrices.qs", preset = "high")

message("[5/8] 在GSE16561发现并在GSE22255独立验证LASSO模型")

zscore_genes <- function(expr) {
  z <- t(scale(t(expr)))
  z[!is.finite(z)] <- 0
  z
}

model_genes <- Reduce(intersect, list(lactylation, rownames(gse16561$expr), rownames(gse22255$expr)))
train_expr <- zscore_genes(gse16561$expr[model_genes, , drop = FALSE])
valid_expr <- zscore_genes(gse22255$expr[model_genes, , drop = FALSE])
x_train <- t(train_expr)
x_valid <- t(valid_expr)
y_train <- ifelse(gse16561$meta$condition == "IS", 1, 0)
y_valid <- ifelse(gse22255$meta$condition == "IS", 1, 0)

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold_id <- integer(length(y))
  for (value in sort(unique(y))) {
    index <- which(y == value)
    fold_id[index] <- sample(rep(seq_len(k), length.out = length(index)))
  }
  fold_id
}

weights_train <- ifelse(y_train == 1, length(y_train) / (2 * sum(y_train == 1)), length(y_train) / (2 * sum(y_train == 0)))
fold_id <- make_stratified_folds(y_train, 10, 20260818)
lasso_cv <- glmnet::cv.glmnet(
  x = x_train,
  y = y_train,
  family = "binomial",
  alpha = 1,
  type.measure = "auc",
  nfolds = 10,
  foldid = fold_id,
  weights = weights_train,
  standardize = FALSE,
  keep = TRUE
)

lasso_coef <- as.matrix(stats::coef(lasso_cv, s = "lambda.min"))
selected_genes <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
selected_genes <- setdiff(selected_genes, "(Intercept)")
write.table(data.frame(Gene = selected_genes, Coefficient = lasso_coef[selected_genes, 1]), "output/07_discovery_lasso_model.tsv", sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

train_score <- as.numeric(stats::predict(lasso_cv, newx = x_train, s = "lambda.min", type = "response"))
valid_score <- as.numeric(stats::predict(lasso_cv, newx = x_valid, s = "lambda.min", type = "response"))

roc_metrics <- function(y, score, cohort, model_name, threshold = NULL) {
  roc_obj <- pROC::roc(response = factor(ifelse(y == 1, "IS", "Control"), levels = c("Control", "IS")), predictor = score, levels = c("Control", "IS"), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  if (is.null(threshold)) {
    coordinate <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
  } else {
    coordinate <- pROC::coords(roc_obj, x = threshold, input = "threshold", ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
  }
  data.frame(
    Cohort = cohort,
    Model = model_name,
    AUC = as.numeric(pROC::auc(roc_obj)),
    CI_low = ci[1],
    CI_high = ci[3],
    Threshold = as.numeric(coordinate$threshold),
    Sensitivity = as.numeric(coordinate$sensitivity),
    Specificity = as.numeric(coordinate$specificity),
    Brier = mean((score - y)^2),
    stringsAsFactors = FALSE
  )
}

lasso_train_metrics <- roc_metrics(y_train, train_score, "GSE16561", "LASSO panel")
lasso_valid_metrics <- roc_metrics(y_valid, valid_score, "GSE22255", "LASSO panel", threshold = lasso_train_metrics$Threshold)

single_gene_auc <- sapply(selected_genes, function(gene) {
  as.numeric(pROC::auc(pROC::roc(factor(ifelse(y_train == 1, "IS", "Control"), levels = c("Control", "IS")), x_train[, gene], levels = c("Control", "IS"), direction = "auto", quiet = TRUE)))
})
top_single_gene <- names(which.max(single_gene_auc))
single_direction <- if (mean(x_train[y_train == 1, top_single_gene]) >= mean(x_train[y_train == 0, top_single_gene])) 1 else -1
single_train_score <- single_direction * x_train[, top_single_gene]
single_valid_score <- single_direction * x_valid[, top_single_gene]
single_train_metrics <- roc_metrics(y_train, single_train_score, "GSE16561", paste0("Single gene: ", top_single_gene))
single_valid_metrics <- roc_metrics(y_valid, single_valid_score, "GSE22255", paste0("Single gene: ", top_single_gene), threshold = single_train_metrics$Threshold)

validation_metrics <- rbind(lasso_train_metrics, lasso_valid_metrics, single_train_metrics, single_valid_metrics)
write.csv(validation_metrics, "output/08_discovery_validation_metrics.csv", row.names = FALSE, fileEncoding = "UTF-8")

validation_logit <- stats::qlogis(pmin(pmax(valid_score, 1e-6), 1 - 1e-6))
calibration_model <- stats::glm(y_valid ~ validation_logit, family = stats::binomial())
calibration <- data.frame(
  Metric = c("Calibration intercept", "Calibration slope"),
  Estimate = c(stats::coef(calibration_model)[1], stats::coef(calibration_model)[2])
)
write.csv(calibration, "output/09_validation_calibration.csv", row.names = FALSE, fileEncoding = "UTF-8")

message("[6/8] 进行200次重复分层交叉验证并计算选择稳定性")

stability_count <- setNames(integer(length(model_genes)), model_genes)
for (iteration in seq_len(200)) {
  repeated_fold <- make_stratified_folds(y_train, 10, 20260818 + iteration)
  repeated_fit <- glmnet::cv.glmnet(
    x = x_train,
    y = y_train,
    family = "binomial",
    alpha = 1,
    type.measure = "auc",
    nfolds = 10,
    foldid = repeated_fold,
    weights = weights_train,
    standardize = FALSE
  )
  repeated_coef <- as.matrix(stats::coef(repeated_fit, s = "lambda.min"))
  repeated_selected <- setdiff(rownames(repeated_coef)[repeated_coef[, 1] != 0], "(Intercept)")
  stability_count[repeated_selected] <- stability_count[repeated_selected] + 1L
  if (iteration %% 20 == 0) message(sprintf("  稳定性评估：%d/200", iteration))
}

stability <- data.frame(Gene = names(stability_count), Selection_frequency = as.numeric(stability_count) / 200)
stability <- stability[order(stability$Selection_frequency, decreasing = TRUE), ]
write.csv(stability, "output/10_lasso_selection_stability.csv", row.names = FALSE, fileEncoding = "UTF-8")
qs::qsave(list(model = lasso_cv, selected_genes = selected_genes, stability = stability, metrics = validation_metrics, calibration = calibration), "checkpoints/05_lasso_validation.qs", preset = "high")

message("[7/8] 评估原六基因的队列内方向与AUC")

gene_auc_by_cohort <- function(cohort, cohort_name, genes) {
  y <- factor(cohort$meta$condition, levels = c("Control", "IS"))
  result <- lapply(genes, function(gene) {
    values <- as.numeric(cohort$expr[gene, ])
    roc_obj <- pROC::roc(y, values, levels = c("Control", "IS"), direction = "auto", quiet = TRUE)
    ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
    data.frame(
      Cohort = cohort_name,
      Gene = gene,
      Mean_difference_IS_minus_Control = mean(values[y == "IS"]) - mean(values[y == "Control"]),
      AUC = as.numeric(pROC::auc(roc_obj)),
      CI_low = ci[1],
      CI_high = ci[3]
    )
  })
  do.call(rbind, result)
}

six_auc <- rbind(
  gene_auc_by_cohort(gse16561, "GSE16561", original_six),
  gene_auc_by_cohort(gse22255, "GSE22255", original_six)
)
write.csv(six_auc, "output/11_original_six_cohort_auc.csv", row.names = FALSE, fileEncoding = "UTF-8")

message("[8/8] 将分型限定于单一发现队列并输出汇总")

subtype_genes <- intersect(original_six, rownames(gse16561$expr))
case_ids <- gse16561$meta$Accession[gse16561$meta$condition == "IS"]
subtype_matrix <- t(zscore_genes(gse16561$expr[subtype_genes, case_ids, drop = FALSE]))
subtype_fit <- cluster::pam(subtype_matrix, k = 2)
subtype_result <- data.frame(Accession = rownames(subtype_matrix), Subtype = paste0("Cluster ", subtype_fit$clustering))
write.csv(subtype_result, "output/12_single_cohort_subtypes.csv", row.names = FALSE, fileEncoding = "UTF-8")

original_clusters <- read.delim("input/original_clusters.txt", check.names = FALSE)
dataset_from_accession <- function(id) {
  ifelse(id %in% gse16561$meta$Accession, "GSE16561", ifelse(id %in% gse22255$meta$Accession, "GSE22255", "GSE37587"))
}
original_clusters$Cohort <- dataset_from_accession(original_clusters$ID)
original_subtype_table <- table(original_clusters$cluster, original_clusters$Cohort)
original_subtype_test <- stats::chisq.test(original_subtype_table)
write.csv(as.data.frame.matrix(original_subtype_table), "output/diagnostic_original_subtype_by_cohort.csv", fileEncoding = "UTF-8")
write.csv(data.frame(Statistic = unname(original_subtype_test$statistic), df = unname(original_subtype_test$parameter), P_value = original_subtype_test$p.value), "output/diagnostic_original_subtype_test.csv", row.names = FALSE, fileEncoding = "UTF-8")

summary_table <- data.frame(
  Item = c(
    "Primary case-control cohorts",
    "Primary case-control samples",
    "Common genes",
    "Meta-significant direction-consistent genes",
    "Cleaned lactylation genes",
    "Meta-significant lactylation genes",
    "Discovery LASSO selected genes",
    "Independent validation AUC",
    "Single-cohort subtype size",
    "Single-cohort silhouette"
  ),
  Value = c(
    "GSE16561 and GSE22255",
    ncol(gse16561$expr) + ncol(gse22255$expr),
    length(common_genes),
    sum(meta_de$meta_significant, na.rm = TRUE),
    length(lactylation),
    sum(meta_de$meta_significant & meta_de$Gene %in% lactylation, na.rm = TRUE),
    paste(selected_genes, collapse = "; "),
    sprintf("%.3f (95%% CI %.3f-%.3f)", lasso_valid_metrics$AUC, lasso_valid_metrics$CI_low, lasso_valid_metrics$CI_high),
    paste(table(subtype_result$Subtype), collapse = "; "),
    sprintf("%.3f", subtype_fit$silinfo$avg.width)
  )
)
write.csv(summary_table, "output/13_analysis_summary.csv", row.names = FALSE, fileEncoding = "UTF-8")
qs::qsave(list(summary = summary_table, subtype = subtype_result, subtype_fit = subtype_fit), "checkpoints/06_final_summary.qs", preset = "high")

message("分析完成。")
print(summary_table, row.names = FALSE)
