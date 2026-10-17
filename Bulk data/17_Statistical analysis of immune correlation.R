# 次要意见2：计算23种免疫细胞的完整相关性并逐基因BH校正。
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
options(stringsAsFactors = FALSE)
Sys.setenv(LANGUAGE = "en")
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
module <- dirname(normalizePath(script_path, winslash = "/"))
root <- dirname(module)
input_dir <- file.path(root, "数据")
log_dir <- file.path(dirname(root), "logs")
dir.create(log_dir, showWarnings = FALSE)
source_dir <- input_dir
result_dir <- file.path(root, "结果表", "次要意见2")
save_table <- function(x, name) write.csv(x, file.path(result_dir, name), row.names = FALSE, fileEncoding = "UTF-8", na = "NA")
stopifnot(requireNamespace("qs", quietly = TRUE))

message("[1/3] 核查23种免疫细胞、原绘图基因与样本对应")
scores <- as.matrix(read.delim(file.path(source_dir, "merge_ssGSEA.result.txt"), row.names = 1, check.names = FALSE))
expression <- as.matrix(read.csv(file.path(source_dir, "merge矩阵new.csv"), row.names = 1, check.names = FALSE))
clinical <- read.delim(file.path(source_dir, "clinicaldata.txt"), check.names = FALSE)
genes <- read.delim(file.path(source_dir, "机器学习交集.txt"), check.names = FALSE)[[1]]
gmt <- strsplit(readLines(file.path(input_dir, "immune.gmt"), encoding = "UTF-8"), "\t")
gmt_ids <- vapply(gmt, `[`, character(1), 1)
stopifnot(nrow(scores) == 23, length(genes) == 6, !anyDuplicated(genes), !anyDuplicated(rownames(scores)),
          identical(gmt_ids, rownames(scores)), all(genes %in% rownames(expression)),
          setequal(clinical$Accession, colnames(scores)), setequal(colnames(scores), colnames(expression)), ncol(scores) == 171)
samples <- colnames(scores)
expression <- expression[genes, samples, drop = FALSE]
clinical <- clinical[match(samples, clinical$Accession), ]
stopifnot(all(is.finite(scores)), all(is.finite(expression)), !anyDuplicated(samples),
          all(apply(scores, 1, sd) > 0), all(apply(expression, 1, sd) > 0))

# 按评分ID规范显示名称，保持与原评分的逐行对应。
cell_labels <- gsub("\\.+", " ", sub("na$", "", rownames(scores)))
legacy_labels <- gsub("\\.", " ", gsub("na", "", rownames(scores), fixed = TRUE))
legacy_first <- read.csv(file.path(source_dir, paste0(genes[1], "相关性分析结果.csv")), check.names = FALSE)
audit <- data.frame(Score_row = seq_len(nrow(scores)), Score_ID = rownames(scores), Immune_cell = cell_labels,
                     Legacy_cell_label = legacy_labels, Legacy_included = legacy_labels %in% legacy_first$cell,
                     GMT_gene_entries = vapply(gmt, function(x) length(x) - 2L, integer(1)))
stopifnot(sum(audit$Legacy_included) == 22, audit$Immune_cell[!audit$Legacy_included] == "Activated B cell")
save_table(audit, "02_immune_cell_type_audit.csv")
save_table(data.frame(Accession = samples, Condition = clinical$condition), "01_correlation_sample_scope.csv")
save_table(data.frame(Setting = c("Correlation", "Exact", "Alternative", "Samples", "Immune cell types", "BH family", "BH family size", "Total tests", "Figure scope"),
                      Value = c("Spearman", "FALSE (original cor.test setting)", "two.sided", "171 existing records: 127 IS and 44 Control", "23 existing ssGSEA score rows",
                                "Separate family for each displayed gene", "23 tests per gene; no raw-P prefilter", "6 x 23 = 138", "Only gene-immune correlation figure; original A/B and all original graphics retained")),
           "00_analysis_settings.csv")

message("[2/3] 计算138项完整相关性并逐基因进行23项BH校正")
complete <- do.call(rbind, lapply(genes, function(gene) {
  table <- do.call(rbind, lapply(seq_len(nrow(scores)), function(i) {
    test <- stats::cor.test(as.numeric(scores[i, ]), as.numeric(expression[gene, ]), method = "spearman", exact = FALSE, alternative = "two.sided")
    data.frame(Gene = gene, Score_ID = rownames(scores)[i], Immune_cell = cell_labels[i], Spearman_rho = unname(test$estimate),
               N = ncol(scores), P_raw = test$p.value, Legacy_included = audit$Legacy_included[i])
  }))
  table$P_BH <- stats::p.adjust(table$P_raw, method = "BH")
  table$BH_family <- paste0(gene, ": all 23 immune score types")
  table$BH_family_size <- 23L
  table$Significant_raw <- table$P_raw < .05
  table$Significant_BH <- table$P_BH < .05
  table
}))
stopifnot(nrow(complete) == 138, !anyDuplicated(complete[c("Gene", "Score_ID")]), all(is.finite(complete$P_BH)))
save_table(complete, "03_complete_immune_correlations.csv")
save_table(complete[complete$Immune_cell == "Activated B cell", ], "06_activated_B_correlations.csv")

message("[3/3] 逐项核对原132项数值、统计家族与存储重载")
reproduction <- do.call(rbind, lapply(genes, function(gene) {
  previous <- read.csv(file.path(source_dir, paste0(gene, "相关性分析结果.csv")), check.names = FALSE)
  current <- complete[complete$Gene == gene & complete$Legacy_included, ]
  current$Legacy_cell_label <- legacy_labels[match(current$Score_ID, rownames(scores))]
  current <- current[match(previous$cell, current$Legacy_cell_label), ]
  stopifnot(nrow(previous) == 22, !anyNA(current$Score_ID))
  delta_rho <- max(abs(previous$cor - current$Spearman_rho))
  delta_p <- max(abs(previous$p.value - current$P_raw))
  data.frame(Gene = gene, Original_tests = nrow(previous), Max_rho_difference = delta_rho, Max_raw_P_difference = delta_p,
             PASS = delta_rho < 1e-12 & delta_p < 1e-12)
}))
stopifnot(all(reproduction$PASS))
save_table(reproduction, "05_original_correlation_reproduction.csv")
family_summary <- do.call(rbind, lapply(genes, function(gene) {
  current <- complete[complete$Gene == gene, ]
  previous <- read.csv(file.path(source_dir, paste0(gene, "相关性分析结果.csv")), check.names = FALSE)
  data.frame(Gene = gene, Original_tests = 22L, Complete_tests = 23L,
             Original_raw_P_below_0_05 = sum(previous$p.value < .05), Complete_raw_P_below_0_05 = sum(current$P_raw < .05),
             Original_22_test_BH_below_0_05 = sum(stats::p.adjust(previous$p.value, "BH") < .05),
             Complete_23_test_BH_below_0_05 = sum(current$P_BH < .05),
             Complete_BH_significant_in_original_22 = sum(current$P_BH < .05 & current$Legacy_included))
}))
save_table(family_summary, "04_BH_family_summary.csv")
qs::qsave(list(scores = scores, expression = expression, samples = clinical, genes = genes, audit = audit,
               complete_correlations = complete, reproduction = reproduction, family_summary = family_summary,
               method = "Spearman; exact=FALSE; BH separately over every gene's 23 tests"),
           file.path(input_dir, "01_immune_correlations.qs"), preset = "high")
reload <- qs::qread(file.path(input_dir, "01_immune_correlations.qs"))
stopifnot(identical(reload$complete_correlations, complete), identical(reload$scores, scores), identical(reload$expression, expression))
save_table(data.frame(Check = c("Complete 138 tests", "23 cell types", "Original 132 results", "QS reload"), Result = "PASS"), "07_analysis_checks.csv")
capture.output(sessionInfo(), file = file.path(log_dir, "M18_immune_sessionInfo.txt"))
print(family_summary)
message("完整相关性及BH校正完成；接续绘图与交付核查")
