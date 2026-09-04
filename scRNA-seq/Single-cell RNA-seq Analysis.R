# Single-cell transcriptome analysis of brain and blood in the acute and subacute phases of experimental stroke ----------------------------------------------
# GSE225948_Single-cell quality control and annotation --------------------
setwd("D:/实验性脑卒中急性期和亚急性期脑和血液单细胞转录组分析/Rscript")
rm(list=ls());gc()
setwd("D:/实验性脑卒中急性期和亚急性期脑和血液单细胞转录组分析/Rscript")
dir.create("1.Seurat", showWarnings = F); setwd("1.Seurat")

library(Seurat)
library(Matrix)
library(data.table)
library(dplyr)
library(R.utils)

if(FALSE){base_dir <- "D:/实验性脑卒中急性期和亚急性期脑和血液单细胞转录组分析/Rscript/0.rawdata/GSE225948"
raw_dir <- file.path(base_dir, "GSE225948_RAW")


# 读取单个样本
read_sample <- function(counts_gz) {
  sample_id <- sub("_counts\\.csv\\.gz$", "", basename(counts_gz))
  
  # 读取 counts
  csv_counts <- sub("\\.gz$", "", counts_gz); gunzip(counts_gz, destname = csv_counts, remove = F)
  df <- fread(csv_counts, header = T, data.table = F); file.remove(csv_counts)
  genes <- df[[1]]; counts_mat <- as.matrix(df[, -1, drop = F])
  orig_bc <- colnames(df)[-1]; colnames(counts_mat) <- orig_bc; rownames(counts_mat) <- genes
  
  # 读取 metadata
  meta_gz <- sub("_counts\\.csv\\.gz$", "_metadata.csv.gz", counts_gz)
  if (file.exists(meta_gz)) {
    csv_meta <- sub("\\.gz$", "", meta_gz); gunzip(meta_gz, destname = csv_meta, remove = F)
    mdf <- fread(csv_meta, header = T, data.table = F); file.remove(csv_meta)
    rn <- mdf[[1]]; mdf <- as.data.frame(mdf[, -1, drop = F]); rownames(mdf) <- rn
    if (!all(colnames(counts_mat) == rownames(mdf))) mdf <- mdf[colnames(counts_mat), , drop = F]
    meta <- mdf
  } else {
    meta <- NULL
  }
  
  list(sample_id = sample_id, counts = as(counts_mat, "dgCMatrix"), meta = meta)
}

# 读取所有样本
count_files <- sort(list.files(raw_dir, pattern = "_counts\\.csv\\.gz$", full.names = T))
samples <- lapply(count_files, read_sample)

# 合并基因和 counts
all_genes <- unique(unlist(lapply(samples, function(x) rownames(x$counts))))
aligned <- lapply(samples, function(s) {
  m <- s$counts
  M <- Matrix(0, nrow = length(all_genes), ncol = ncol(m), sparse = T)
  idx <- match(rownames(m), all_genes); M[idx, ] <- m
  colnames(M) <- paste0(s$sample_id, ".", colnames(m))
  M
})
merged_counts <- do.call(cbind, aligned); rownames(merged_counts) <- all_genes

# 合并 metadata
metas_prefix <- lapply(samples, function(s) {
  if (is.null(s$meta)) return(NULL)
  meta <- s$meta
  rownames(meta) <- paste0(s$sample_id, ".", rownames(meta))
  meta$sample_id <- s$sample_id
  meta
})
metas_prefix <- metas_prefix[!sapply(metas_prefix, is.null)]
cell_meta <- do.call(rbind, metas_prefix)

# 验证并排序
if (!all(colnames(merged_counts) %in% rownames(cell_meta))) {
  stop("Metadata 缺少细胞信息")
}
cell_meta <- cell_meta[colnames(merged_counts), , drop = F]

# 创建 Seurat 对象
scRNA <- CreateSeuratObject(counts = merged_counts, meta.data = cell_meta, 
                            project = "GSE225948", min.cells = 3, min.features = 200)

# 保存结果
saveRDS(scRNA, file.path(base_dir, "scRNA_GSE225948_loaded.rds"))
write.csv(data.frame(cell = colnames(scRNA), scRNA@meta.data), file.path(base_dir, "CellMetadata_GSE225948.csv"))

cat("Done: cells =", ncol(scRNA), "| genes =", nrow(scRNA), 
    "| samples =", length(unique(scRNA$sample_id)), "\n")}

#scRNA

# QC-cluster ----------------------------------------------------------------------
if(FALSE){source("./refdata/my_clean_functions.R")
table(scRNA@meta.data$cellType_1)
scRNA <- readRDS("../Rscript/0.rawdata/GSE225948/scRNA_GSE225948_loaded.rds")
dim(scRNA)
#[1]  21946 102282
refdata_path <- "../Rscript/Rawdata/refdata/Homo/"
load(paste0(refdata_path, "BlueprintEncode.Rdata"))
load(paste0(refdata_path, "DatabaseImmuneCell.se.Rdata"))
load(paste0(refdata_path, "HumanPrimaryCellAtlas.Rdata"))
load(paste0(refdata_path, "MonacoImmune.Rdata"))
load(paste0(refdata_path, "NovershternHematopoietic.Rdata"))
table(data.merge.harmony$disease_state.y,data.merge.harmony$sample_id)

scRNA <- scRNAAutoAnno(
  Path = "./",
  SeuratObject = scRNA,
  DoubletFinder = T,
  Multi = TRUE,
  ref = list(
    BlueprintEncode,
    DatabaseImmuneCell.se,
    NovershternHematopoietic,
    HumanPrimaryCellAtlas,
    MonacoImmune
  ),
  labels = list(
    BlueprintEncode$label.main,
    DatabaseImmuneCell.se$label.main,
    NovershternHematopoietic$label.main,
    HumanPrimaryCellAtlas$label.main,
    MonacoImmune$label.main
  ),
  Idents = "RNA_snn_res.0.2")}

scRNA <- readRDS("./0.rawdata/GSE225948/CellAnnoted.cellType_1.rds")
dim(scRNA)
#[1] 21946 42251
Idents(scRNA) <- "cellType_1"
scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05,]
write.csv(scRNA.markers, "./3.CellAnnotate/Findall.markers.cellType_1.csv")
RePlot(scRNA = scRNA,Ident = "RNA_snn_res.0.2",Multi = TRUE)

# CellAnnote ------------------------------------------------
# 参数设置  -----------------------------------------------------------------

input_rds      <- "./0.rawdata/GSE225948/data.merge.harmony.2000.rds"  # 输入RDS路径
output_rds     <- "./0.rawdata/GSE225948/CellAnnoted.cellType_1.rds"   # 输出RDS路径
output_dir     <- "03.CellAnnote"                                       # 输出目录
cluster_col    <- "RNA_snn_res.0.2"                                     # 聚类列名
celltype_col   <- "cellType_1"                                          # 注释列名
treatment_col  <- "treatment"                                           # 分组列名

# 标记基因（按细胞类型分组，注释说明对应cluster编号）
features <- c(
  "Ccl5","Trbc2","Nkg7",      # T cells         - cluster 4
  "Igkc","Ms4a1","Cd79a",     # B cells         - cluster 8
  "S100a8","S100a9","Cxcr2",  # Neutrophils     - cluster 5
  "Lyve1","Mrc1","Cd163",     # Macrophages     - cluster 6
  "Chil3","Lyz2","Lgals3",    # Monocytes       - cluster 1
  "Hexb","Mef2c","Gpr34",     # Microglial cells- cluster 0,7
  "H2-Ab1","H2-Aa","Cd74",    # APC             - cluster 3
  "Cldn5","Esam","Slc2a1",    # Endothelial     - cluster 2
  "Plp1","Mbp","Mag"          # Oligodendrocytes- cluster 9
)

# cluster编号 → 细胞类型映射（左：cluster编号字符串，右：细胞类型名称）
cluster_map <- c(
  "0" = "Microglial cells",
  "1" = "Monocytes",
  "2" = "Endothelial cells",
  "3" = "Antigen-presenting cells",
  "4" = "T cells",
  "5" = "Neutrophils",
  "6" = "Macrophages",
  "7" = "Microglial cells",
  "8" = "B cells",
  "9" = "Oligodendrocytes"
)

# cellType因子顺序
celltype_levels <- c(
  "T cells", "B cells", "Neutrophils", "Macrophages", "Monocytes",
  "Microglial cells", "Antigen-presenting cells", "Endothelial cells", "Oligodendrocytes"
)

# 配色
my_pal2 <- c(
  "#D4477D","#D24B27","#4DBBD5","#6387C5","#6E4B9E",
  "#C10020","#1E78B4","#FCBF6E","#83AD00","#9ebcda",
  "#74a9cf","#fbdf72","#FF8E00","#F37B7D","#CF4A31","#F37B7D","#FF8E00"
)

# 读取数据
data.merge.harmony <- readRDS(input_rds)
dim(data.merge.harmony)
#single_cell
# 0     1     2     3     4     5     6     7     8     9
# Hepatocytes        19310     0     0     0     0     0     0     0     0     0
# Neurons               0  9102     0  3379  2202     0     0     0     0     0
# Platelets             0     0  4847     0     0     0     0     0     0     0
# unknow                0     0     0     0     0  1099     0     0     0     0
# Endothelial_cells     0     0     0     0     0     0   896     0     0     0
# Dendritic cells       0     0     0     0     0     0     0   806     0     0
# Melanocytes           0     0     0     0     0     0     0     0   386     0
# Endothelial cells     0     0     0     0     0     0     0     0     0   224

# 1. 各基因 FeaturePlot
for (i in features) {
  p <- FeaturePlot(data.merge.harmony, cols = c("lightgrey","red"), features = i, reduction = "umap")
  ggsave_fun(filename = paste0(output_dir, "/01.Featureplot/", i, ".featureplot"), plot = p, width = 7, height = 7)
}


# 2. 聚类 DotPlot（注释前）
p <- DotPlot(data.merge.harmony, features = features, cols = "RdYlBu", group.by = cluster_col) +
  scale_size_continuous(range = c(0, 10)) +
  theme(
    panel.border     = element_rect(colour = "black"),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position  = "top",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    text             = element_text(size = 10)
  )
ggsave_fun(filename = paste0(output_dir, "/Cluster"), plot = p, width = 10, height = 10)


# 3. Cluster → 细胞类型注释
cluster_label <- as.character(data.merge.harmony@meta.data[[cluster_col]])
for (k in names(cluster_map)) {
  cluster_label <- gsub(paste0("^", k, "$"), cluster_map[k], cluster_label)
}
data.merge.harmony <- AddMetaData(data.merge.harmony, cluster_label, col.name = celltype_col)
data.merge.harmony[[celltype_col]] <- factor(data.merge.harmony[[celltype_col]], levels = celltype_levels)


# 4. UMAP 可视化
p <- DimPlot(data.merge.harmony, group.by = celltype_col, raster = FALSE, label = TRUE, reduction = "umap") + NoLegend()
ggsave_fun(filename = paste0(output_dir, "/03.umap.plot"), plot = p, width = 7, height = 7)


# 5. 注释后 DotPlot
p <- DotPlot(data.merge.harmony, features = features, cols = "RdYlBu", group.by = celltype_col) +
  scale_size_continuous(range = c(0, 10)) +
  theme(
    panel.border      = element_rect(colour = "black"),
    axis.text.x       = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position   = "top",
    legend.key.height = unit(0.3, "cm"),
    legend.key.width  = unit(0.8, "cm"),
    text              = element_text(size = 10)
  )
ggsave_fun(filename = paste0(output_dir, "/04.dotplot"), plot = p, width = 14, height = 7)


# 6. 细胞比例图
library(scRNAtoolVis)
data.merge.harmony[[treatment_col]] <- gsub("\\(D\\d+\\)", "", data.merge.harmony[[treatment_col]])
print(table(data.merge.harmony[[treatment_col]]))

cellRatio <- cellRatioPlot(
  object        = data.merge.harmony,
  sample.name   = treatment_col,
  celltype.name = celltype_col,
  flow.curve    = 0.5,
  fill.col      = my_pal2
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave_fun(filename = paste0(output_dir, "/10.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 15, height = 15)


# 7. 保存 & 统计
saveRDS(data.merge.harmony, output_rds)
print(table(data.merge.harmony[[celltype_col]]))
# Expression of keygene ---------------------------------------------------
source("./refdata/my_clean_functions.R")

scRNA <- readRDS("../Rscript/0.rawdata/GSE225948/CellAnnoted.cellType_1.rds")
key_genes <- c("Crebbp","Emg1","Hist2h2be","Rbm14","S100a4","Srp14")

# 运行7
scRNA <- scRNAAutoAnno(
  Path = ".",
  SeuratObject = scRNA,
  KeyGene = key_genes
)
