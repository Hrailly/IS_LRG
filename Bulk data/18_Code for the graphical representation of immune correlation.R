# 仅绘制新的基因—免疫相关性图；所有原图及原A/B面板保留。
Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8")
options(stringsAsFactors = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
module <- dirname(normalizePath(script_path, winslash = "/"))
root <- dirname(module)
input_dir <- file.path(root, "数据")
log_dir <- file.path(dirname(root), "logs")
dir.create(log_dir, showWarnings = FALSE)
result_dir <- file.path(root, "结果表", "次要意见2")
figure_dir <- file.path(root, "补充图")
objects <- qs::qread(file.path(input_dir, "01_immune_correlations.qs"))
plot_data <- objects$complete_correlations

# 沿用原免疫组图的基因显示顺序；23种细胞按原评分矩阵顺序显示。
gene_order <- c("CREBBP", "EMG1", "HIST2H2BE", "RBM14", "S100A4", "SRP14")
stopifnot(setequal(gene_order, objects$genes))
plot_data$Gene <- factor(plot_data$Gene, levels = gene_order)
plot_data$Immune_cell <- factor(plot_data$Immune_cell, levels = rev(objects$audit$Immune_cell))
plot_data$BH_label <- ifelse(plot_data$P_BH < .001, "***", ifelse(plot_data$P_BH < .01, "**", ifelse(plot_data$P_BH < .05, "*", "ns")))
plot_data$Text_color <- ifelse(abs(plot_data$Spearman_rho) > .45, "white", "black")
write.csv(plot_data, file.path(result_dir, "08_figure4_plot_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")

plot <- ggplot2::ggplot(plot_data, ggplot2::aes(Gene, Immune_cell)) +
  ggplot2::geom_tile(ggplot2::aes(fill = Spearman_rho), color = "white", linewidth = .55) +
  ggplot2::geom_text(ggplot2::aes(label = BH_label, color = Text_color), size = 3.4, fontface = "bold", show.legend = FALSE) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradient2(low = "#315E7D", mid = "#F5F4EF", high = "#A45E43", midpoint = 0,
                               limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1), name = "Spearman rho") +
  ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(title = "Gene and Immune Correlations", subtitle = "BH adjusted P: * < 0.05; ** < 0.01; *** < 0.001; ns >= 0.05", x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 11, base_family = "Arial") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5, face = "bold", color = "black", size = 15),
                 plot.subtitle = ggplot2::element_text(hjust = .5, color = "black", size = 9.2),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black", size = 11),
                 axis.text.y = ggplot2::element_text(color = "black", size = 10.5),
                 panel.grid = ggplot2::element_blank(), legend.position = "bottom", legend.title = ggplot2::element_text(face = "bold", color = "black"),
                 legend.text = ggplot2::element_text(color = "black"), plot.margin = ggplot2::margin(10, 12, 7, 8))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
grDevices::cairo_pdf(file.path(figure_dir, "07_Immune_Correlations_BH.pdf"), width = 8.5, height = 8.6, family = "Arial")
print(plot)
grDevices::dev.off()
qs::qsave(list(plot_data = plot_data, plot = plot, gene_order = gene_order), file.path(input_dir, "02_figure4_objects.qs"), preset = "high")
reload <- qs::qread(file.path(input_dir, "02_figure4_objects.qs"))
stopifnot(identical(reload$plot_data, plot_data))
capture.output(sessionInfo(), file = file.path(log_dir, "M18_immune_figures_sessionInfo.txt"))
message("独立Figure 4相关性PDF已生成；原A/B及原图未改")
