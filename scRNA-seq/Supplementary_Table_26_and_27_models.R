ct <- read.csv(file.choose(), stringsAsFactors = FALSE)
sc <- read.csv(file.choose(), stringsAsFactors = FALSE)

names(ct)[names(ct) == "cellType_1"] <- "feature"
names(sc)[names(sc) == "Mics_subcluster"] <- "feature"

prep <- function(d) {
  d$age <- factor(ifelse(as.integer(substr(d$sample_id, 8, 10)) >= 827, "aged", "young"),
                  levels = c("young", "aged"))
  d$group <- factor(d$group, levels = c("Sham", "Stroke"))
  d$time <- factor(d$time, levels = c("2d", "14d"))
  d$p <- d$cells / d$sample_total_cells
  d$y <- qlogis((d$cells + 0.5) / (d$sample_total_cells + 1))
  d
}

models <- list(primary = y ~ group + time + age,
               sensitivity = y ~ group + time)

fit_all <- function(d) {
  do.call(rbind, lapply(split(d, d$feature), function(s) {
    do.call(rbind, lapply(names(models), function(mn) {
      m <- lm(models[[mn]], data = s)
      cf <- summary(m)$coefficients
      k <- grep("groupStroke", rownames(cf))
      data.frame(feature = s$feature[1], model = mn, n = nrow(s), df = m$df.residual,
                 sham_mean = mean(s$p[s$group == "Sham"]),
                 stroke_mean = mean(s$p[s$group == "Stroke"]),
                 estimate = cf[k, 1], se = cf[k, 2], t = cf[k, 3], p = cf[k, 4],
                 row.names = NULL)
    }))
  }))
}

res_ct <- fit_all(prep(ct))
res_sc <- fit_all(prep(sc))

res_ct$fdr <- ave(res_ct$p, res_ct$model, FUN = function(x) p.adjust(x, "BH"))
res_sc$fdr <- ave(res_sc$p, res_sc$model, FUN = function(x) p.adjust(x, "BH"))

print(res_ct, digits = 4)
print(res_sc, digits = 4)
