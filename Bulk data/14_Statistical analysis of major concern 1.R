# 主要问题1补充：保持原主分析，单独核查队列、校正模型、患者时间与原亚型。
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
options(stringsAsFactors = FALSE)
set.seed(20261010)
Sys.setenv(LANGUAGE = "en")
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
module <- dirname(normalizePath(script_path, winslash = "/"))
root <- dirname(module)
input_dir <- file.path(root, "数据")
log_dir <- file.path(dirname(root), "logs")
dir.create(log_dir, showWarnings = FALSE)
result_dir <- file.path(root, "结果表", "主要意见1")
figure_dir <- file.path(root, "补充图")
data_dir <- file.path(input_dir)
dir.create(result_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)
packages <- c("limma", "sva", "preprocessCore", "qs", "ggplot2", "lme4", "statmod")
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
save_table <- function(x, name) write.csv(x, file.path(result_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "NA")
read_tsv_matrix <- function(path) as.matrix(read.delim(path, row.names = 1, check.names = FALSE, fileEncoding = "UTF-8"))
fits <- list()
cohorts <- c("GSE16561", "GSE22255", "GSE37587")

message("[1/7] 读取原流程输入并建立样本及患者对应关系")
clinical <- read.delim(file.path(input_dir, "clinicaldata.txt"), check.names = FALSE)
metadata <- do.call(rbind, lapply(cohorts, function(cohort) {
  samples <- read.delim(file.path(input_dir, paste0(cohort, "样本.txt")))
  data.frame(Accession = samples$Accession, Cohort = cohort, Condition = samples$condition)
}))
stopifnot(!anyDuplicated(metadata$Accession), identical(sort(metadata$Accession), sort(clinical$Accession)))
info16561 <- read.delim(file.path(input_dir, "GSE16561_SampleInfo.txt"), check.names = FALSE)
info37587 <- read.delim(file.path(input_dir, "GSE37587_SampleInfo.txt"), check.names = FALSE)
metadata$Patient_within_series <- paste0(metadata$Cohort, ":", metadata$Accession)
metadata$Patient_linked <- metadata$Patient_within_series
metadata$Time <- "Not specified"
metadata$Followup37587 <- 0
index37587 <- match(info37587$Accession, metadata$Accession)
metadata$Patient_within_series[index37587] <- paste0("GSE37587:", info37587[["1:PATIENT NUMBER"]])
metadata$Patient_linked[index37587] <- metadata$Patient_within_series[index37587]
metadata$Time[index37587] <- info37587[["1:TIME"]]
metadata$Followup37587[index37587] <- as.integer(info37587[["1:TIME"]] == "Follow-Up")
stopifnot(length(unique(metadata$Patient_within_series[index37587])) == 34, all(table(metadata$Patient_within_series[index37587]) == 2))
alias_index <- match(sub("_Stroke$", "", info16561$Title), info37587$Title)
matched <- which(!is.na(alias_index))
aliases <- data.frame(
  Sample_GSE16561 = info16561$Accession[matched],
  Sample_GSE37587 = info37587$Accession[alias_index[matched]],
  Shared_title = info37587$Title[alias_index[matched]],
  Patient_GSE37587 = info37587[["1:PATIENT NUMBER"]][alias_index[matched]],
  Time_GSE37587 = info37587[["1:TIME"]][alias_index[matched]],
  Age_GSE16561 = info16561[["1:AGE"]][matched],
  Age_GSE37587 = info37587[["1:AGE"]][alias_index[matched]],
  Gender_GSE16561 = info16561[["1:GENDER"]][matched],
  Gender_GSE37587 = info37587[["1:GENDER"]][alias_index[matched]],
  Linkage_basis = "Identical specimen title, age and sex; conservative shared-participant assumption")
stopifnot(nrow(aliases) == 34, all(aliases$Age_GSE16561 == aliases$Age_GSE37587), all(aliases$Gender_GSE16561 == aliases$Gender_GSE37587))
metadata$Patient_linked[match(aliases$Sample_GSE16561, metadata$Accession)] <- paste0("GSE37587:", aliases$Patient_GSE37587)
save_table(metadata, "01_sample_patient_time.csv")
save_table(aliases, "02_cross_series_identifier_matches.csv")
input_matrices <- lapply(cohorts, function(cohort) read_tsv_matrix(file.path(input_dir, paste0(cohort, ".txt"))))
names(input_matrices) <- cohorts
common_genes <- Reduce(intersect, lapply(input_matrices, rownames))
preprocessed <- lapply(input_matrices, function(x) {
  x <- limma::avereps(x)
  quantiles <- as.numeric(quantile(x, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
  take_log <- quantiles[5] > 100 || (quantiles[6] - quantiles[1] > 50 && quantiles[2] > 0)
  if (take_log) { x[x < 0] <- 0; x <- log2(x + 1) }
  limma::normalizeBetweenArrays(x)
})
cohort_matrices <- lapply(preprocessed, function(x) x[common_genes, , drop = FALSE])
combined <- do.call(cbind, cohort_matrices)
meta_combined <- metadata[match(colnames(combined), metadata$Accession), ]
rownames(meta_combined) <- meta_combined$Accession
meta_combined$Condition <- factor(meta_combined$Condition, levels = c("Control", "IS"))
meta_combined$Cohort <- factor(meta_combined$Cohort, levels = cohorts)
sample_audit <- do.call(rbind, lapply(cohorts, function(cohort) {
  m <- meta_combined[meta_combined$Cohort == cohort, ]
  data.frame(Cohort = cohort, Samples = nrow(m), IS = sum(m$Condition == "IS"), Control = sum(m$Condition == "Control"),
             Within_series_units = length(unique(m$Patient_within_series)), Shared_identifier_matches = sum(m$Accession %in% aliases$Sample_GSE16561 | m$Accession %in% aliases$Sample_GSE37587))
}))
save_table(sample_audit, "01_cohort_design.csv")

# 所有效应表保留全基因范围；显示基因仅对应上一轮回复已声称的一致性。
claim_genes <- c("HIST2H2BE", "EMG1", "RBM14", "SRP14")
save_table(data.frame(Setting = c("Disease contrast", "Subtype contrast", "Gene universe", "Displayed concordance genes", "Display source", "Multiplicity family", "Cross-series linkage"),
  Value = c("IS - Control", "B - A (original explicit contrast)", paste(length(common_genes), "original common genes"), paste(claim_genes, collapse = "; "),
            "Previous response paragraph 10; not a new candidate selection", "BH across every tested row within each model", "Matching specimen title, age and sex are linked only as a conservative sensitivity assumption")), "00_analysis_settings.csv")

fit_effects <- function(x, design, contrast, label, block = NULL, feature_type = "Gene") {
  stopifnot(ncol(x) == nrow(design), qr(design)$rank == ncol(design))
  correlation <- NA_real_
  if (is.null(block)) {
    fit <- limma::lmFit(x, design)
  } else {
    stopifnot(length(block) == ncol(x))
    correlation <- limma::duplicateCorrelation(x, design, block = block)$consensus.correlation
    fit <- limma::lmFit(x, design, block = block, correlation = correlation)
  }
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(contrasts = contrast, levels = design))
  fit <- limma::eBayes(fit)
  out <- limma::topTable(fit, number = Inf, sort.by = "none", adjust.method = "BH", confint = .95)
  table <- data.frame(Feature = rownames(out), Model = label, Contrast = contrast, Feature_type = feature_type,
                      Effect = out$logFC, CI_low = out$CI.L, CI_high = out$CI.R, P_raw = out$P.Value, P_BH = out$adj.P.Val,
                      SE_moderated = as.numeric(fit$stdev.unscaled[, 1] * sqrt(fit$s2.post)), DF_moderated = fit$df.total,
                      Samples = ncol(x), Independent_units = if (is.null(block)) ncol(x) else length(unique(block)),
                      Within_patient_correlation = correlation)
  fits[[label]] <<- fit
  list(table = table, fit = fit, correlation = correlation)
}
disease_design <- function(m, time = FALSE) {
  disease <- as.integer(m$Condition == "IS")
  if (time) cbind(Intercept = 1, Disease = disease, Followup37587 = m$Followup37587) else cbind(Intercept = 1, Disease = disease)
}
quantile_matrix <- function(x) {
  out <- preprocessCore::normalize.quantiles(as.matrix(x))
  dimnames(out) <- dimnames(x)
  out
}

message("[2/7] 重现原数值并估计两个队列内的疾病效应")
reconstructed_combat <- sva::ComBat(combined, batch = meta_combined$Cohort, par.prior = TRUE)
original_environment <- new.env()
load(file.path(input_dir, "merge_all.RDATA"), envir = original_environment)
preserved_combat <- as.matrix(original_environment$outTab)[common_genes, colnames(combined), drop = FALSE]
preserved_main <- as.matrix(read.csv(file.path(input_dir, "merge矩阵new.csv"), row.names = 1, check.names = FALSE))
preserved_main <- preserved_main[common_genes, clinical$Accession, drop = FALSE]
reconstructed_main <- quantile_matrix(reconstructed_combat[, clinical$Accession, drop = FALSE])
original_meta <- meta_combined[clinical$Accession, ]
original_fit <- fit_effects(preserved_main, disease_design(original_meta), "Disease", "Original_published")
published_de <- read.csv(file.path(input_dir, "差异分析结果.csv"), row.names = 1, check.names = FALSE)
original_effect <- original_fit$table$Effect[match(rownames(published_de), original_fit$table$Feature)]
original_p <- original_fit$table$P_raw[match(rownames(published_de), original_fit$table$Feature)]
reproduction <- data.frame(Check = c("Text-input ComBat replay", "Text-input final matrix replay", "Published effects", "Published raw P", "Saved ComBat to final matrix"),
                           Maximum_absolute_difference = c(max(abs(reconstructed_combat - preserved_combat)), max(abs(reconstructed_main - preserved_main)),
                                                           max(abs(original_effect - published_de$logFC)), max(abs(original_p - published_de$P.Value)),
                                                           max(abs(quantile_matrix(preserved_combat[, clinical$Accession, drop = FALSE]) - preserved_main))),
                           Check_role = c("Input replay trace", "Input replay trace", "Preserved result verification", "Preserved result verification", "Preserved result verification"))
reproduction$PASS <- reproduction$Maximum_absolute_difference < 1e-7
save_table(reproduction, "00_original_numeric_reproduction.csv")
stopifnot(all(reproduction$PASS[reproduction$Check_role == "Preserved result verification"]))
replay_fit <- fit_effects(reconstructed_main, disease_design(original_meta), "Disease", "Text_input_replay")
save_table(replay_fit$table, "00_text_input_replay_effects.csv")
cohort_effects <- do.call(rbind, lapply(cohorts[1:2], function(cohort) {
  x <- cohort_matrices[[cohort]]
  m <- meta_combined[colnames(x), ]
  table <- fit_effects(x, disease_design(m), "Disease", cohort)$table
  table$Cohort <- cohort
  table$N_IS <- sum(m$Condition == "IS")
  table$N_Control <- sum(m$Condition == "Control")
  case <- x[, m$Condition == "IS", drop = FALSE]
  control <- x[, m$Condition == "Control", drop = FALSE]
  n1 <- ncol(case); n0 <- ncol(control)
  pooled_sd <- sqrt(((n1 - 1) * apply(case, 1, sd)^2 + (n0 - 1) * apply(control, 1, sd)^2) / (n1 + n0 - 2))
  correction <- 1 - 3 / (4 * (n1 + n0) - 9)
  table$Hedges_g <- correction * (rowMeans(case) - rowMeans(control)) / pooled_sd
  table$Hedges_SE <- sqrt((n1 + n0) / (n1 * n0) + table$Hedges_g^2 / (2 * (n1 + n0 - 2)))
  table$Hedges_CI_low <- table$Hedges_g - 1.96 * table$Hedges_SE
  table$Hedges_CI_high <- table$Hedges_g + 1.96 * table$Hedges_SE
  table
}))
save_table(cohort_effects, "03_cohort_specific_effects.csv")
claim_table <- cohort_effects[cohort_effects$Feature %in% claim_genes, ]
claim_table$Display_basis <- "Direction claim in previous response paragraph 10"
save_table(claim_table, "04_previous_claim_concordance.csv")

message("[3/7] 排除GSE37587并计算疾病模型ComBat补充对照")
keep <- meta_combined$Cohort != "GSE37587"
exclusion_raw <- combined[, keep, drop = FALSE]
exclusion_meta <- droplevels(meta_combined[keep, ])
exclusion_combat <- sva::ComBat(exclusion_raw, batch = exclusion_meta$Cohort, par.prior = TRUE)
exclusion_matrix <- quantile_matrix(exclusion_combat)
exclusion_fit <- fit_effects(exclusion_matrix, disease_design(exclusion_meta), "Disease", "Exclude_GSE37587")
save_table(exclusion_fit$table, "05_exclusion_effects.csv")
combat_mod <- model.matrix(~ Condition + Followup37587, data = meta_combined)
preserved_disease_combat <- sva::ComBat(combined, batch = meta_combined$Cohort, mod = combat_mod, par.prior = TRUE)
disease_matrix <- quantile_matrix(preserved_disease_combat)
old_matrix_ordered <- preserved_main[, colnames(combined), drop = FALSE]
# 仅比较校正方式的均值效应，不给忽略重复采样的新推断P值。
mean_difference <- function(x, m) rowMeans(x[, m$Condition == "IS", drop = FALSE]) - rowMeans(x[, m$Condition == "Control", drop = FALSE])
combat_comparison <- data.frame(Feature = common_genes, Original_effect = mean_difference(old_matrix_ordered, meta_combined),
                                Biological_model_effect = mean_difference(disease_matrix, meta_combined))
combat_comparison$Same_direction <- sign(combat_comparison$Original_effect) == sign(combat_comparison$Biological_model_effect)
save_table(combat_comparison, "05_combat_model_effect_comparison.csv")

message("[4/7] 患者相关性及时间调整；分别计算队列内患者与共享标识敏感性")
time_design <- disease_design(meta_combined, time = TRUE)
own_patient_fit <- fit_effects(disease_matrix, time_design, "Disease", "Biological_model_within_series_patients", block = meta_combined$Patient_within_series)
linked_patient_fit <- fit_effects(disease_matrix, time_design, "Disease", "Biological_model_linked_patients", block = meta_combined$Patient_linked)
original_linked_fit <- fit_effects(old_matrix_ordered, time_design, "Disease", "Original_matrix_linked_patients", block = meta_combined$Patient_linked)
patient_effects <- rbind(own_patient_fit$table, linked_patient_fit$table, original_linked_fit$table)
save_table(patient_effects, "06_patient_time_effects.csv")
models <- rbind(original_fit$table, replay_fit$table, exclusion_fit$table, patient_effects)
save_table(models, "07_disease_effects_all_models.csv")
baseline <- original_fit$table
baseline <- baseline[match(common_genes, baseline$Feature), ]
comparisons <- do.call(rbind, lapply(unique(models$Model)[-1], function(label) {
  x <- models[models$Model == label, ]
  x <- x[match(baseline$Feature, x$Feature), ]
  selected <- baseline$P_BH < .05
  data.frame(Model = label, Genes = nrow(x), Original_significant = sum(selected), Model_significant = sum(x$P_BH < .05),
             Spearman_effect = cor(baseline$Effect, x$Effect, method = "spearman"), Pearson_effect = cor(baseline$Effect, x$Effect),
             All_direction_agreement = mean(sign(baseline$Effect) == sign(x$Effect)),
             Original_significant_direction_agreement = mean(sign(baseline$Effect[selected]) == sign(x$Effect[selected])),
             Original_significant_retained_BH = sum(selected & x$P_BH < .05),
             Original_significant_retained_same_direction_BH = sum(selected & x$P_BH < .05 & sign(baseline$Effect) == sign(x$Effect)),
             Samples = x$Samples[1], Units = x$Independent_units[1], Patient_correlation = x$Within_patient_correlation[1])
}))
save_table(comparisons, "08_disease_sensitivity_summary.csv")

message("[5/7] 原亚型按队列和时间分层；患者分组的亚型构成检验")
clusters <- read.delim(file.path(input_dir, "Cluster.txt"), check.names = FALSE)
sub_meta <- meta_combined[clusters$ID, ]
stopifnot(all(sub_meta$Condition == "IS"), nrow(sub_meta) == 127)
sub_meta$Subtype <- factor(clusters$cluster, levels = c("A", "B"))
sub_meta$SubtypeB <- as.integer(sub_meta$Subtype == "B")
save_table(sub_meta, "09_original_subtype_sample_metadata.csv")
distribution <- as.data.frame(table(Cohort = sub_meta$Cohort, Subtype = sub_meta$Subtype))
names(distribution)[3] <- "Samples"
distribution$Percent <- 100 * distribution$Samples / ave(distribution$Samples, distribution$Cohort, FUN = sum)
save_table(distribution, "09_subtype_by_cohort.csv")
longitudinal <- sub_meta[sub_meta$Cohort == "GSE37587", ]
time_distribution <- as.data.frame(table(Time = longitudinal$Time, Subtype = longitudinal$Subtype))
names(time_distribution)[3] <- "Samples"
time_distribution$Percent <- 100 * time_distribution$Samples / ave(time_distribution$Samples, time_distribution$Time, FUN = sum)
save_table(time_distribution, "09_subtype_by_time.csv")
long_pairs <- merge(longitudinal[longitudinal$Time == "Baseline", c("Patient_within_series", "Accession", "Subtype")],
                    longitudinal[longitudinal$Time == "Follow-Up", c("Patient_within_series", "Accession", "Subtype")], by = "Patient_within_series", suffixes = c("_Baseline", "_Followup"))
stopifnot(nrow(long_pairs) == 34)
save_table(long_pairs, "09_patient_subtype_transitions.csv")
transitions <- table(long_pairs$Subtype_Baseline, long_pairs$Subtype_Followup)
discordant <- transitions[1, 2] + transitions[2, 1]
paired_p <- if (discordant == 0) 1 else binom.test(transitions[1, 2], discordant, p = .5)$p.value
save_table(data.frame(Test = "Exact paired McNemar (binomial discordant pairs)", Patients = 34, A_to_B = transitions[1, 2], B_to_A = transitions[2, 1], P_raw = paired_p), "09_subtype_time_test.csv")
subtype_glmm <- lme4::glmer(SubtypeB ~ Cohort + Followup37587 + (1 | Patient_linked), data = sub_meta, family = binomial(),
                            control = lme4::glmerControl(optimizer = "bobyqa"))
without_cohort <- lme4::glmer(SubtypeB ~ Followup37587 + (1 | Patient_linked), data = sub_meta, family = binomial(),
                            control = lme4::glmerControl(optimizer = "bobyqa"))
cohort_lrt <- anova(without_cohort, subtype_glmm, test = "Chisq")
cohort_p <- cohort_lrt$`Pr(>Chisq)`[2]
coef_glmm <- coef(summary(subtype_glmm))
save_table(data.frame(Term = rownames(coef_glmm), Log_odds_B = coef_glmm[, 1], SE = coef_glmm[, 2],
                      CI_low = coef_glmm[, 1] - 1.96 * coef_glmm[, 2], CI_high = coef_glmm[, 1] + 1.96 * coef_glmm[, 2], P_raw = coef_glmm[, 4]), "09_subtype_cohort_model.csv")
save_table(data.frame(Test = "Patient-grouped logistic mixed-model cohort LRT", Statistic = cohort_lrt$Chisq[2], DF = cohort_lrt$Df[2], P_raw = cohort_p,
                      Singular_fit = lme4::isSingular(subtype_glmm), Linked_units = length(unique(sub_meta$Patient_linked))), "09_subtype_cohort_test.csv")
aliases$Subtype_GSE16561 <- as.character(sub_meta[aliases$Sample_GSE16561, "Subtype"])
aliases$Subtype_GSE37587 <- as.character(sub_meta[aliases$Sample_GSE37587, "Subtype"])
aliases$Same_subtype <- aliases$Subtype_GSE16561 == aliases$Subtype_GSE37587
save_table(aliases, "02_cross_series_identifier_matches.csv")

message("[6/7] 原亚型表达和已有GSVA分数的队列内效应")
subtype_design <- function(m, cohort = FALSE, time = FALSE) {
  d <- cbind(Intercept = 1, SubtypeB = m$SubtypeB)
  if (cohort) d <- cbind(d, Cohort22255 = as.integer(m$Cohort == "GSE22255"), Cohort37587 = as.integer(m$Cohort == "GSE37587"))
  if (time) d <- cbind(d, Followup37587 = m$Followup37587)
  d
}
subtype_stratified <- function(x, feature_type, prefix) {
  out <- lapply(cohorts, function(cohort) {
    m <- sub_meta[sub_meta$Cohort == cohort, ]
    paired <- cohort == "GSE37587"
    r <- fit_effects(x[, m$Accession, drop = FALSE], subtype_design(m, time = paired), "SubtypeB", paste0(prefix, "_", cohort),
                     block = if (paired) m$Patient_within_series else NULL, feature_type = feature_type)$table
    r$Cohort <- cohort
    r$N_A <- sum(m$Subtype == "A"); r$N_B <- sum(m$Subtype == "B")
    r
  })
  r <- fit_effects(x[, sub_meta$Accession, drop = FALSE], subtype_design(sub_meta, cohort = TRUE, time = TRUE), "SubtypeB", paste0(prefix, "_All_cohorts_adjusted"),
                   block = sub_meta$Patient_linked, feature_type = feature_type)$table
  r$Cohort <- "All cohorts (adjusted)"; r$N_A <- sum(sub_meta$Subtype == "A"); r$N_B <- sum(sub_meta$Subtype == "B")
  do.call(rbind, c(out, list(r)))
}
# 原分型使用ComBat后、跨队列二次分位数归一化前的表达矩阵。
sub_expression <- preserved_combat[, sub_meta$Accession, drop = FALSE]
sub_gene_effects <- subtype_stratified(sub_expression, "Gene", "Subtype_expression")
save_table(sub_gene_effects, "10_subtype_gene_effects.csv")
pathway_sets <- c(Hallmark = "hallmark", KEGG = "KEGG", Reactome = "Reactome")
pathway_effects <- list()
pathway_reproduction <- list()
pathway_original <- list()
for (set in names(pathway_sets)) {
  message("  原有GSVA分数：", set)
  folder <- input_dir
  scores <- read_tsv_matrix(file.path(folder, paste0(set, "_gsvaOut.txt")))
  scores <- scores[, sub_meta$Accession, drop = FALSE]
  original <- fit_effects(scores, subtype_design(sub_meta), "SubtypeB", paste0(set, "_Original"), feature_type = set)$table
  old <- read.delim(file.path(folder, paste0(set, "_B-A.all.txt")), row.names = 1, check.names = FALSE)
  delta <- max(abs(original$Effect[match(rownames(old), original$Feature)] - old$logFC))
  pathway_reproduction[[set]] <- data.frame(Set = set, Pathways = nrow(scores), Maximum_effect_difference = delta, PASS = delta < 1e-7)
  stopifnot(delta < 1e-7)
  original$Cohort <- "Original pooled"; original$N_A <- sum(sub_meta$Subtype == "A"); original$N_B <- sum(sub_meta$Subtype == "B")
  pathway_original[[set]] <- original
  pathway_effects[[set]] <- subtype_stratified(scores, set, paste0("Subtype_", set))
}
save_table(do.call(rbind, pathway_effects), "11_subtype_pathway_effects.csv")
save_table(do.call(rbind, pathway_original), "11_original_pathway_effects_reproduced.csv")
save_table(do.call(rbind, pathway_reproduction), "00_original_pathway_reproduction.csv")
path_summary <- do.call(rbind, lapply(names(pathway_effects), function(set) {
  original <- pathway_original[[set]]
  original <- original[order(original$Feature), ]
  do.call(rbind, lapply(unique(pathway_effects[[set]]$Model), function(label) {
    x <- pathway_effects[[set]][pathway_effects[[set]]$Model == label, ]
    x <- x[match(original$Feature, x$Feature), ]
    selected <- original$P_BH < .05 & abs(original$Effect) > .1
    data.frame(Set = set, Model = label, Original_selected = sum(selected),
               Original_selected_direction_agreement = mean(sign(original$Effect[selected]) == sign(x$Effect[selected])),
               Original_selected_retained_BH = sum(selected & x$P_BH < .05 & sign(original$Effect) == sign(x$Effect)))
  }))
}))
save_table(path_summary, "12_subtype_pathway_summary.csv")

message("[7/7] 保存对象与运行信息")
qs::qsave(list(metadata = meta_combined, genes = common_genes, original_matrix = preserved_main, original_combat_matrix = preserved_combat,
               exclusion_matrix = exclusion_matrix, biological_model_matrix = disease_matrix, combat_mod = combat_mod, subtype_metadata = sub_meta,
               subtype_glmm = subtype_glmm, fits = fits, comparisons = comparisons), file.path(data_dir, "01_supplement_objects.qs"), preset = "high")
reloaded <- qs::qread(file.path(data_dir, "01_supplement_objects.qs"))
stopifnot(identical(reloaded$genes, common_genes), identical(reloaded$metadata$Accession, meta_combined$Accession),
          identical(reloaded$comparisons, comparisons))
capture.output(sessionInfo(), file = file.path(log_dir, "M17_major1_sessionInfo.txt"))
save_table(data.frame(Package = packages, Version = vapply(packages, function(p) as.character(packageVersion(p)), character(1))), "00_package_versions.csv")
save_table(data.frame(Check = c("QS reload", "Gene universe", "GSE37587 paired samples", "Original labels preserved"),
                      Value = c("PASS", as.character(length(common_genes)), "34 patients; 68 samples", "127 original subtype labels")), "19_analysis_completion.csv")
message("统计分析完成；接续患者重抽样、图表生成与交付核查")
