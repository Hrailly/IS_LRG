# 从已计算的统计表生成单页图，保证点估计、置信区间与P值使用同一模型。
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
options(stringsAsFactors = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
module <- dirname(normalizePath(script_path, winslash = "/"))
root <- dirname(module)
input_dir <- file.path(root, "数据")
log_dir <- file.path(dirname(root), "logs")
dir.create(log_dir, showWarnings = FALSE)
result_dir <- file.path(root, "结果表", "主要意见1")
figure_dir <- file.path(root, "补充图")
read_result <- function(name) read.csv(file.path(result_dir, name), check.names = FALSE)
cohorts <- c("GSE16561", "GSE22255", "GSE37587")
claim_genes <- c("HIST2H2BE", "EMG1", "RBM14", "SRP14")
claim_table <- read_result("04_previous_claim_concordance.csv")
all_models <- read_result("07_disease_effects_all_models.csv")
baseline <- all_models[all_models$Model == "Original_published", ]
exclusion_fit <- list(table = read_result("05_exclusion_effects.csv"))
combat_comparison <- read_result("05_combat_model_effect_comparison.csv")
linked_patient_fit <- list(table = all_models[all_models$Model == "Biological_model_linked_patients", ])
distribution <- read_result("09_subtype_by_cohort.csv")
time_distribution <- read_result("09_subtype_by_time.csv")
cohort_p <- read_result("09_subtype_cohort_test.csv")$P_raw
paired_p <- read_result("09_subtype_time_test.csv")$P_raw
theme_science <- ggplot2::theme_classic(base_size = 11, base_family = "Arial") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold", color = "black", size = 14),
                 axis.title = ggplot2::element_text(face = "bold", color = "black"), axis.text = ggplot2::element_text(color = "black"),
                 legend.title = ggplot2::element_text(face = "bold"), legend.position = "top", plot.margin = ggplot2::margin(9, 12, 9, 9))
save_pdf <- function(plot, name, width, height) {
  while (grDevices::dev.cur() > 1) grDevices::dev.off()
  grDevices::cairo_pdf(file.path(figure_dir, name), width = width, height = height, family = "Arial")
  print(plot)
  grDevices::dev.off()
}
format_p <- function(p) if (p < .001) format(p, scientific = TRUE, digits = 2) else sprintf("%.3f", p)
forest <- claim_table
forest$Feature <- factor(forest$Feature, levels = rev(claim_genes))
forest$Cohort <- factor(forest$Cohort, levels = cohorts[1:2])
forest$P_label <- paste0("P = ", vapply(forest$P_raw, format_p, character(1)))
effect_span <- max(forest$CI_high) - min(forest$CI_low)
label_position <- max(forest$CI_high) + .10 * effect_span
position <- ggplot2::position_dodge(width = .55)
p <- ggplot2::ggplot(forest, ggplot2::aes(x = Effect, y = Feature, color = Cohort, group = Cohort)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  ggplot2::geom_errorbar(ggplot2::aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = .12, position = position) +
  ggplot2::geom_point(size = 2.7, position = position) +
  ggplot2::geom_text(ggplot2::aes(x = label_position, label = P_label), position = position, hjust = 0, size = 3.0, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = c("#315A7D", "#B66B46")) +
  ggplot2::scale_x_continuous(limits = c(min(forest$CI_low) - .08 * effect_span, max(forest$CI_high) + .45 * effect_span)) +
  ggplot2::labs(title = "Disease Effects by Cohort", x = "Expression difference (IS - Control)", y = NULL, color = "Cohort") + theme_science
save_pdf(p, "01_Cohort_Disease_Effects.pdf", 8.8, 4.6)

scatter_effect <- function(x, y, title, xlab, ylab, file) {
  d <- data.frame(Original = x, Supplement = y, Published = baseline$P_BH < .05)
  rho <- cor(x, y, method = "spearman")
  p <- ggplot2::ggplot(d, ggplot2::aes(Original, Supplement)) +
    ggplot2::geom_point(ggplot2::aes(color = Published), size = .65, alpha = .38) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = c("FALSE" = "#AEB7C0", "TRUE" = "#315A7D"), labels = c("Other genes", "Original BH < 0.05")) +
    ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -.03, vjust = 1.3, label = sprintf("Spearman rho = %.3f\nGenes = %s", rho, format(nrow(d), big.mark = ",")), size = 3.5) +
    ggplot2::labs(title = title, x = xlab, y = ylab, color = NULL) + theme_science
  save_pdf(p, file, 6.8, 6.3)
}
ordered_exclusion <- exclusion_fit$table[match(baseline$Feature, exclusion_fit$table$Feature), ]
scatter_effect(baseline$Effect, ordered_exclusion$Effect, "GSE37587 Exclusion Sensitivity", "Original effect (IS - Control)", "Exclusion effect (IS - Control)", "02_GSE37587_Exclusion.pdf")
combat_comparison <- combat_comparison[match(baseline$Feature, combat_comparison$Feature), ]
scatter_effect(combat_comparison$Original_effect, combat_comparison$Biological_model_effect, "ComBat Biological Model Sensitivity", "Original mean difference (IS - Control)", "Biological model mean difference (IS - Control)", "03_ComBat_Model_Comparison.pdf")
linked <- linked_patient_fit$table[match(baseline$Feature, linked_patient_fit$table$Feature), ]
scatter_effect(baseline$Effect, linked$Effect, "Patient and Time Sensitivity", "Original effect (IS - Control)", "Patient/time-adjusted effect (IS - Control)", "04_Patient_Time_Sensitivity.pdf")

distribution$Label <- paste0(distribution$Samples, "\n", sprintf("%.1f%%", distribution$Percent))
p <- ggplot2::ggplot(distribution, ggplot2::aes(Cohort, Percent, fill = Subtype)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(.7), width = .6) +
  ggplot2::geom_text(ggplot2::aes(label = Label), position = ggplot2::position_dodge(.7), vjust = -.25, size = 3.4) +
  ggplot2::scale_fill_manual(values = c(A = "#315A7D", B = "#B66B46")) +
  ggplot2::scale_y_continuous(limits = c(0, 85), expand = ggplot2::expansion(mult = c(0, .01))) +
  ggplot2::annotate("text", x = 2, y = 81, label = paste0("Patient-grouped cohort test: P = ", format_p(cohort_p)), size = 3.4) +
  ggplot2::labs(title = "Subtype Composition by Cohort", x = NULL, y = "Samples within cohort (%)", fill = "Subtype") +
  theme_science + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_pdf(p, "05_Subtype_by_Cohort.pdf", 6.4, 5.5)
time_distribution$Time <- factor(time_distribution$Time, levels = c("Baseline", "Follow-Up"))
time_distribution$Label <- paste0(time_distribution$Samples, "\n", sprintf("%.1f%%", time_distribution$Percent))
p <- ggplot2::ggplot(time_distribution, ggplot2::aes(Time, Percent, fill = Subtype)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(.7), width = .6) +
  ggplot2::geom_text(ggplot2::aes(label = Label), position = ggplot2::position_dodge(.7), vjust = -.25, size = 3.5) +
  ggplot2::scale_fill_manual(values = c(A = "#315A7D", B = "#B66B46")) +
  ggplot2::scale_y_continuous(limits = c(0, 90), expand = ggplot2::expansion(mult = c(0, .01))) +
  ggplot2::annotate("text", x = 1.5, y = 85, label = paste0("Paired exact test: P = ", format_p(paired_p)), size = 3.4) +
  ggplot2::labs(title = "GSE37587 Subtypes over Time", x = NULL, y = "Patients at each time (%)", fill = "Subtype") +
  theme_science + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_pdf(p, "06_Subtype_by_Time.pdf", 5.5, 5.3)


message("六张单页PDF生成完成")
