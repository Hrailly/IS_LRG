library(Seurat)
library(limma)
library(Matrix)

outdir <- file.path(getwd(), "pseudobulk_DE_output")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

voom_de <- function(m, design) {
  lib <- colSums(m)
  m <- m[, lib > 0, drop = FALSE]
  dd <- design[lib > 0, , drop = FALSE]
  cpm <- t(t(m) / colSums(m)) * 1e6
  m <- m[rowSums(cpm > 1) >= 3, , drop = FALSE]
  y <- voom(m, dd, lib.size = colSums(m))
  fit <- eBayes(lmFit(y, dd))
  tt <- topTable(fit, coef = "groupStroke", number = Inf, sort.by = "P")
  tt$n_samples <- nrow(dd)
  tt$df_residual <- fit$df.residual
  tt
}

run_block <- function(obj, column, tag) {
  meta <- obj@meta.data
  cnt <- LayerData(obj, assay = "RNA", layer = "counts")

  meta$group <- factor(sub("\\(.*", "", as.character(meta$treatment)), levels = c("Sham", "Stroke"))
  meta$time <- factor(ifelse(grepl("D14", as.character(meta$treatment)), "14d", "2d"),
                      levels = c("2d", "14d"))
  meta$age <- factor(ifelse(as.integer(substr(as.character(meta$sample_id), 8, 10)) >= 827,
                            "aged", "young"), levels = c("young", "aged"))

  lv <- sort(unique(as.character(meta$sample_id)))
  sm <- meta[match(lv, meta$sample_id), ]
  d_primary <- model.matrix(~ group + time + age, data = sm)
  d_sens <- model.matrix(~ group + time, data = sm)
  ind <- sparseMatrix(i = seq_len(nrow(meta)), j = match(meta$sample_id, lv), x = 1,
                      dims = c(nrow(meta), length(lv)))

  ids <- split(seq_len(nrow(meta)), as.character(meta[[column]]))
  ids <- ids[order(names(ids), decreasing = TRUE)]
  if (all(grepl("^[0-9]+$", names(ids)))) ids <- ids[order(as.numeric(names(ids)))]

  summ <- list()
  full <- list()
  for (nm in names(ids)) {
    m <- as.matrix(cnt[, ids[[nm]], drop = FALSE] %*% ind[ids[[nm]], , drop = FALSE])
    tp <- voom_de(m, d_primary)
    ts <- voom_de(m, d_sens)
    pass <- tp$adj.P.Val < 0.05
    keep <- if (sum(pass) >= 100) pass else seq_len(min(100, nrow(tp)))
    full[[nm]] <- data.frame(
      feature = nm,
      gene = rownames(tp)[keep],
      tp[keep, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "n_samples")],
      Passes_FDR_0.05 = ifelse(pass[keep], "yes", "no"),
      row.names = NULL, check.names = FALSE
    )
    summ[[nm]] <- data.frame(
      Population = nm,
      Libraries = tp$n_samples[1],
      Genes_tested = nrow(tp),
      FDR_primary = sum(pass),
      FDR_sensitivity = sum(ts$adj.P.Val < 0.05),
      Top_gene = rownames(tp)[1],
      Top_gene_log2FC = tp$logFC[1],
      Top_gene_FDR = tp$adj.P.Val[1],
      row.names = NULL, check.names = FALSE
    )
  }

  s <- do.call(rbind, summ)
  f <- do.call(rbind, full)
  write.csv(s, file.path(outdir, paste0(tag, "_summary.csv")), row.names = FALSE)
  write.csv(f, file.path(outdir, paste0(tag, "_full_list.csv")), row.names = FALSE)
  print(s, digits = 4)
  invisible(list(summary = s, full = f))
}

res_R5 <- run_block(readRDS(file.choose()), "cellType_1", "Table_R5_celltype")
invisible(gc())

res_R6 <- run_block(readRDS(file.choose()), "Mics_subcluster", "Table_R6_subcluster")
invisible(gc())

cat("\nwritten to ", outdir, "\n", sep = "")
