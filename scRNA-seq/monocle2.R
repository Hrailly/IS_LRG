# Microglial cells --------------------------------------------------------
setwd("D:/实验性脑卒中急性期和亚急性期脑和血液单细胞转录组分析/Rscript")
# Monocle2 ----------------------------------------------------------------
# 参数设置 --------------------------------------------------------------------
source("./refdata/my_clean_functions.R")
input_rds       <- "./0.rawdata/GSE225948/CellAnnoted.cellType_1.rds"  # 输入RDS文件路径
output_dir      <- "./2.Monocle2/"                                      # 输出目录
cell_type_col   <- "cellType_1"                                         # 细胞类型列名
cell_type_value <- "Microglial cells"                                   # 目标细胞类型
sub_col_name    <- "Mics_subcluster"                                    # 亚群列名
output_rds      <- "scRNA_Mics_reclustered.rds"                         # 输出RDS文件名
output_dimplot  <- "Mics_reclustering.pdf"                              # UMAP图文件名
output_heatmap  <- "Mics_subcluster_markers.pdf"                        # 热图文件名
output_markers  <- "Mics_markers.csv"                                   # 标记基因文件名
cluster_res     <- 0.2                                                  # 选用的聚类分辨率
nfeatures_hvg   <- 3000                                                 # 高变基因数（整合用）
npcs            <- 50                                                   # PCA维数
res_range       <- seq(0.2, 1.2, by = 0.1)                             # 分辨率范围
nfeatures_var   <- 2000                                                 # FindVariableFeatures基因数
min_pct         <- 0.25                                                 # FindAllMarkers最小表达比例
logfc_threshold <- 0.25                                                 # FindAllMarkers logFC阈值
top_n_markers   <- 5                                                    # 热图展示每群Top N基因数


# 1. 提取和聚类
scRNA <- readRDS(input_rds)
scRNA_sub <- subset(scRNA, !!sym(cell_type_col) == cell_type_value)

scRNA_sub <- Harmony.integration.reduceDimension(
  seurat.object = scRNA_sub,
  assay          = "RNA",
  set.resolutions = res_range,
  nfeatures      = nfeatures_hvg,
  npcs           = npcs
)

res_col <- paste0("RNA_snn_res.", cluster_res)
scRNA_sub[[sub_col_name]] <- scRNA_sub[[res_col]]
Idents(scRNA_sub) <- sub_col_name

scRNA_sub <- FindVariableFeatures(scRNA_sub, nfeatures = nfeatures_var)
saveRDS(scRNA_sub, file.path(output_dir, output_rds))


# 2. 聚类分布
cat(sub_col_name, "亚群分布：\n")
print(table(scRNA_sub[[sub_col_name, drop = TRUE]]))


# 3. 可视化
p <- DimPlot(scRNA_sub, group.by = sub_col_name, label = TRUE) +
  ggtitle(paste(sub_col_name, "subclusters")) +
  theme(legend.position = "right")
ggsave(file.path(output_dir, output_dimplot), p, width = 8, height = 6)

#scRNA_sub <- ScaleData(scRNA_sub, features = rownames(scRNA_sub))
# 4. 标记基因
markers <- FindAllMarkers(
  scRNA_sub,
  only.pos       = TRUE,
  min.pct        = min_pct,
  logfc.threshold = logfc_threshold
)
write.csv(markers, file.path(output_dir, output_markers), row.names = FALSE)

top_genes <- markers %>%
  group_by(cluster) %>%
  top_n(top_n_markers, avg_log2FC)

p <- DoHeatmap(scRNA_sub, features = top_genes$gene) + NoLegend()
ggsave(file.path(output_dir, output_heatmap), p, width = 10, height = 8)

# # 5. 起始亚群选择（基于小胶质细胞静息标志物）
# resting_markers <- c("P2ry12", "Tmem119", "Cx3cr1", "Gpr34")
# avg_resting <- AverageExpression(scRNA_Mics, 
#                                  features = resting_markers,
#                                  group.by = "Mics_subcluster")
# resting_index <- colMeans(avg_resting$RNA)
# root_cluster <- names(resting_index)[which.max(resting_index)]
# cat("起始亚群（静息小胶质细胞标志物高表达）：", root_cluster, "\n")
# print(resting_index)


setwd("~/zbq/实验性脑卒中急性期和亚急性期脑和血液单细胞转录组分析/Rscript")
dir.create("2.Monocle2", showWarnings = F); setwd("2.Monocle2")

source("../refdata/my_clean_functions.R")

#scRNA <- readRDS("./0.rawdata/GSE225948/CellAnnoted.cellType_1.rds")

scRNA_Mics <- readRDS(file.path(output_dir, output_rds))

OutPath <- "./"
group_by <- "Mics_subcluster"   #group
ShowGene <- c("Crebbp","Emg1","Hist2h2be","Rbm14","S100a4","Srp14")


cds <-
  monocle2_module(
    object = scRNA_Mics,
    OutPath = OutPath,
    group_by = group_by,
    ShowGene = ShowGene,
    OnlyPlot = F
  )
monocle2_branch_module(
  cds = cds,
  OutPath = OutPath,
  branch_point = 1,#看需修改
  species = "mouse",
  OnlyPlot = F
)

# #查看细胞分化命运对应的states
# # 先重建 branch_cds（内部数据）
# branch_cds <- buildBranchCellDataSet(
#   cds,
#   branch_point = 1,
#   progenitor_method = 'duplicate'#还可共享原始细胞，参数为overlap
# )
# head(branch_cds)
# 
# # branch_cds 里直接有 Branch 列
# branch_pdata <- pData(branch_cds)
# 
# branch_pdata %>%
#   as.data.frame() %>%
#   #dplyr::filter(!grepl("duplicate", rownames(.))) %>%
#   dplyr::select(State, Branch) %>%
#   dplyr::group_by(State, Branch) %>%
#   dplyr::summarise(n = dplyr::n())
# 
# levels(branch_pdata$Branch)
# # [1] "Y_18" "Y_34"
# #      ↑         ↑
# #   第1位       第2位
# # Cell fate 1  Cell fate 2
# # state2       state3
