# Function_set==========================================
# cellchat ----------------------------------------------------------------


calculate_centrality_from_net <- function(cellchat_obj) {
  cat("直接从 @netP$net 矩阵计算中心性...\n")
  
  # 检查并初始化 @netP
  if(is.null(cellchat_obj@netP)) {
    cat("⚠️  @netP 为空，初始化...\n")
    cellchat_obj@netP <- list()
  }
  
  # 检查 @netP$net
  if (is.null(cellchat_obj@netP$net) || length(cellchat_obj@netP$net) == 0) {
    cat("✗ @netP$net 为空，无法计算中心性\n")
    cat("初始化空的中心性结构...\n")
    cellchat_obj@netP$centr <- list()
    return(cellchat_obj)
  }
  
  all_pathways <- names(cellchat_obj@netP$net)
  cell_types <- rownames(cellchat_obj@netP$net[[1]])
  
  cat("处理", length(all_pathways), "条通路，", length(cell_types), "种细胞类型\n")
  
  new_centr <- list()
  
  for (pathway_name in all_pathways) {
    pathway_mat <- cellchat_obj@netP$net[[pathway_name]]
    
    if (!is.null(pathway_mat) && nrow(pathway_mat) > 0) {
      outgoing <- as.numeric(rowSums(pathway_mat))
      names(outgoing) <- rownames(pathway_mat)
      
      incoming <- as.numeric(colSums(pathway_mat))
      names(incoming) <- colnames(pathway_mat)
      
      pathway_sum <- sum(pathway_mat, na.rm = TRUE)
      if (pathway_sum > 0) {
        contribution <- as.numeric(rowSums(pathway_mat) + colSums(pathway_mat)) / (pathway_sum * 2)
      } else {
        contribution <- setNames(rep(0, length(cell_types)), cell_types)
      }
      names(contribution) <- cell_types
      
      new_centr[[pathway_name]] <- list(
        outgoing = outgoing,
        incoming = incoming,
        contribution = contribution
      )
    }
  }
  
  cellchat_obj@netP$centr <- new_centr
  
  cat("✓ 中心性计算完成: ", length(new_centr), " 个通路\n")
  if (length(new_centr) > 0) {
    first_pw <- names(new_centr)[1]
    cat("  示例通路 '", first_pw, "' - outgoing 总和: ", 
        sum(new_centr[[first_pw]]$outgoing), "\n", sep="")
  }
  
  return(cellchat_obj)
}
# 修复合并后 CellChat 对象的网络结构
repair_merged_cellchat_network <- function(merged_cellchat, original_list) {
  cat("\n=== 合并后网络修复 ===\n")
  
  # 检查合并后的结构
  cat("合并后的 idents 水平:", paste(levels(merged_cellchat@idents), collapse = ", "), "\n")
  cat("合并后 @idents 的大小:", length(merged_cellchat@idents), "\n")
  
  # 修复 @net 矩阵
  if (!is.null(merged_cellchat@net$count) && nrow(merged_cellchat@net$count) > 0) {
    cat("修复 @net$count 矩阵...\n")
    
    # 获取合并后的所有细胞类型
    all_idents <- levels(merged_cellchat@idents)
    
    # 清理矩阵维度（移除组标签）
    old_count <- merged_cellchat@net$count
    old_weight <- merged_cellchat@net$weight
    
    cat("  原始矩阵行数:", nrow(old_count), "\n")
    
    # 如果行名包含组标识符（如 "Negative_CD8+T cell"），需要解析
    if (any(grepl("_", rownames(old_count)))) {
      cat("  检测到含有组标识的行名，进行解析...\n")
      
      # 提取纯细胞类型名（去除组前缀）
      clean_rownames <- sub("^[^_]+_", "", rownames(old_count))
      clean_colnames <- sub("^[^_]+_", "", colnames(old_count))
      
      # 创建干净的矩阵
      unique_cells <- unique(c(clean_rownames, clean_colnames))
      
      new_count <- matrix(0, nrow = length(unique_cells), ncol = length(unique_cells),
                          dimnames = list(unique_cells, unique_cells))
      new_weight <- new_count
      
      # 聚合数据
      for (i in 1:nrow(old_count)) {
        for (j in 1:ncol(old_count)) {
          row_label <- clean_rownames[i]
          col_label <- clean_colnames[j]
          new_count[row_label, col_label] <- new_count[row_label, col_label] + old_count[i, j]
          new_weight[row_label, col_label] <- new_weight[row_label, col_label] + old_weight[i, j]
        }
      }
      
      merged_cellchat@net$count <- new_count
      merged_cellchat@net$weight <- new_weight
      
      cat("  ✓ 矩阵清理完成: ", nrow(new_count), " x ", ncol(new_count), "\n", sep="")
    }
  }
  
  # 修复 @netP$net - 确保所有通路的维度一致
  if (!is.null(merged_cellchat@netP$net)) {
    cat("修复 @netP$net 矩阵...\n")
    
    all_pathways <- names(merged_cellchat@netP$net)
    
    for (pathway_name in all_pathways) {
      pathway_mat <- merged_cellchat@netP$net[[pathway_name]]
      
      if (!is.null(pathway_mat) && nrow(pathway_mat) > 0) {
        # 解析行列名中的组标识
        if (any(grepl("_", rownames(pathway_mat)))) {
          clean_rownames <- sub("^[^_]+_", "", rownames(pathway_mat))
          clean_colnames <- sub("^[^_]+_", "", colnames(pathway_mat))
          
          unique_cells <- unique(c(clean_rownames, clean_colnames))
          
          # 创建新矩阵，维度统一
          new_pathway_mat <- matrix(0, 
                                    nrow = length(unique_cells), 
                                    ncol = length(unique_cells),
                                    dimnames = list(unique_cells, unique_cells))
          
          # 聚合数据
          for (i in 1:nrow(pathway_mat)) {
            for (j in 1:ncol(pathway_mat)) {
              row_label <- clean_rownames[i]
              col_label <- clean_colnames[j]
              new_pathway_mat[row_label, col_label] <- new_pathway_mat[row_label, col_label] + pathway_mat[i, j]
            }
          }
          
          merged_cellchat@netP$net[[pathway_name]] <- new_pathway_mat
        }
      }
    }
    
    cat("  ✓ @netP$net 修复完成\n")
  }
  
  # 重新计算合并后的中心性
  cat("重新计算中心性...\n")
  merged_cellchat <- calculate_centrality_from_net(merged_cellchat)
  
  return(merged_cellchat)
}
# 重建合并对象的 @idents
rebuild_merged_idents <- function(merged_cellchat, object.list) {
  cat("\n=== 重建合并对象的 @idents ===\n")
  
  # 从 object.list 中提取所有的 idents 和元数据
  all_idents <- c()
  all_meta <- data.frame()
  
  for (i in 1:length(object.list)) {
    obj <- object.list[[i]]
    group_name <- names(object.list)[i]
    
    # 获取每个对象的 idents
    idents_with_group <- paste0(group_name, "_", as.character(obj@idents))
    all_idents <- c(all_idents, idents_with_group)
    
    cat("  组 '", group_name, "' 的 idents 数: ", length(idents_with_group), "\n", sep="")
  }
  
  # 转换为因子
  all_idents_factor <- factor(all_idents)
  
  # 尝试赋值
  tryCatch({
    merged_cellchat@idents <- all_idents_factor
    cat("✓ @idents 重建成功，水平数: ", length(levels(merged_cellchat@idents)), "\n", sep="")
  }, error = function(e) {
    cat("⚠️  无法直接赋值 @idents，原因: ", e$message, "\n", sep="")
  })
  
  return(merged_cellchat)
}

# 修改 harmonize_cellchat_idents 函数 - 在 Step 2 后添加中心性计算
harmonize_cellchat_idents <- function(cellchat_obj, desired_levels) {
  cat("\n--- START Harmonizing CellChat Object ---\n")
  cat(paste0("  Initial levels(cellchat_obj@idents): ", paste(levels(cellchat_obj@idents), collapse=", "), "\n"))
  
  desired_levels_factor <- factor(desired_levels, levels = desired_levels) 
  
  # Step 1: 更新 @idents
  cellchat_obj@idents <- factor(as.character(cellchat_obj@idents), levels = levels(desired_levels_factor))
  cat(paste0("  AFTER Step 1 - levels: ", paste(levels(cellchat_obj@idents), collapse=", "), "\n"))
  
  # Step 2: 更新 @net$count 和 @net$weight 矩阵
  new_count_matrix <- matrix(0, 
                             nrow = length(desired_levels), 
                             ncol = length(desired_levels), 
                             dimnames = list(desired_levels, desired_levels))
  new_weight_matrix <- new_count_matrix 
  
  current_net_count_names <- rownames(cellchat_obj@net$count)
  common_rows_cols <- intersect(current_net_count_names, desired_levels)
  
  if (length(common_rows_cols) > 0) {
    new_count_matrix[common_rows_cols, common_rows_cols] <- cellchat_obj@net$count[common_rows_cols, common_rows_cols]
    new_weight_matrix[common_rows_cols, common_rows_cols] <- cellchat_obj@net$weight[common_rows_cols, common_rows_cols]
  }
  
  cellchat_obj@net$count <- new_count_matrix
  cellchat_obj@net$weight <- new_weight_matrix
  cat(paste0("  AFTER Step 2 - @net$count sum: ", sum(cellchat_obj@net$count), "\n"))
  
  # Step 3: 更新 @netP$net 矩阵
  if (!is.null(cellchat_obj@netP$net)) {
    for (pathway_name in names(cellchat_obj@netP$net)) {
      old_pathway_matrix <- cellchat_obj@netP$net[[pathway_name]]
      new_pathway_matrix <- matrix(0, 
                                   nrow = length(desired_levels), 
                                   ncol = length(desired_levels), 
                                   dimnames = list(desired_levels, desired_levels))
      
      current_pathway_names <- rownames(old_pathway_matrix)
      common_rows_cols_p <- intersect(current_pathway_names, desired_levels)
      
      if (length(common_rows_cols_p) > 0) {
        new_pathway_matrix[common_rows_cols_p, common_rows_cols_p] <- old_pathway_matrix[common_rows_cols_p, common_rows_cols_p]
      }
      cellchat_obj@netP$net[[pathway_name]] <- new_pathway_matrix
    }
  }
  
  # ========== 新的 Step 4：直接从 @netP$net 计算中心性 ==========
  cellchat_obj <- calculate_centrality_from_net(cellchat_obj)
  
  cat(paste0("--- END Harmonizing CellChat Object ---\n"))
  return(cellchat_obj)
}

# 从通信数据重新构建 @netP$net 矩阵
build_netP_from_communication <- function(cellchat, df.netp) {
  cat("从通信数据框重建 @netP$net...\n")
  
  if(nrow(df.netp) == 0) {
    cat("✗ 通信数据为空！\n")
    return(cellchat)
  }
  
  # 初始化 @netP 结构
  if(is.null(cellchat@netP)) {
    cellchat@netP <- list()
  }
  
  cell_types <- levels(cellchat@idents)
  pathways <- unique(df.netp$pathway_name)
  
  cat("检测到的通路数:", length(pathways), "\n")
  
  # 为每个通路构建矩阵
  net_list <- list()
  
  for(pathway in pathways) {
    # 过滤该通路的通信
    pathway_data <- df.netp[df.netp$pathway_name == pathway, ]
    
    # 创建空矩阵
    pathway_mat <- matrix(0, 
                          nrow = length(cell_types), 
                          ncol = length(cell_types),
                          dimnames = list(cell_types, cell_types))
    
    # 填充矩阵（使用通信的平均值）
    for(i in 1:nrow(pathway_data)) {
      source <- pathway_data$source[i]
      target <- pathway_data$target[i]
      prob <- pathway_data$prob[i]
      
      if(source %in% cell_types && target %in% cell_types) {
        pathway_mat[source, target] <- pathway_mat[source, target] + prob
      }
    }
    
    net_list[[pathway]] <- pathway_mat
    cat("  ", pathway, ": ", nrow(pathway_data), " interactions\n", sep="")
  }
  
  cellchat@netP$net <- net_list
  cellchat@netP$pathways <- pathways
  
  cat("✓ @netP$net 重建完成，共", length(net_list), "条通路\n")
  
  return(cellchat)
}

# cellchat_func 函数 (保持不变)
cellchat_func <- function(object, species = NA, celltype = "cellType_1", 
                          ShowCell = NA, OnlyPlot = F, AllPlot = F){
  
  if(!OnlyPlot){
    if (!is.na(species)) {
      mouse <- ifelse(species == "mouse", T, F)
    } else {
      mouse <- F
    }
    if(mouse) {
      CellChatDB <- CellChatDB.mouse
      PPI.use <- PPI.mouse
    } else {
      CellChatDB <- CellChatDB.human
      PPI.use <- PPI.human
    }
    
    data.input <- object@assays$RNA@data
    
    identity <- data.frame(
      droplevels(object@meta.data[[celltype]]), 
      row.names = colnames(data.input)
    )
    colnames(identity) <- celltype
    
    cellchat <- createCellChat(data.input, meta = identity, group.by = celltype)
    
    CellChatDB.use <- subsetDB(CellChatDB, search = c('Secreted Signaling','Cell-Cell Contact','ECM-Receptor'))
    cellchat@DB <- CellChatDB.use
    
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    set.seed(1234)
    cellchat <- projectData(cellchat, PPI.use)
    
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 3)
    
    df.net <- subsetCommunication(cellchat)
    write.table(df.net, file = paste0("net_lr.txt"), quote = F, sep = "\t", row.names = F)
    
    cellchat <- computeCommunProbPathway(cellchat)
    
    cat("\n=== 通路计算诊断 ===\n")
    cat("通路数:", length(cellchat@netP$pathways), "\n")
    if(length(cellchat@netP$pathways) > 0) {
      cat("通路列表:", paste(cellchat@netP$pathways, collapse = ", "), "\n")
    }
    
    df.netp <- subsetCommunication(cellchat)
    cat("通路级别通信数:", nrow(df.netp), "\n")
    write.table(df.netp, file = paste0("net_pathway.txt"), quote = F, sep = "\t", row.names = F)
    
    cat("\n执行 aggregateNet...\n")
    cellchat <- aggregateNet(cellchat)
    
    cat("aggregateNet 后的检查:\n")
    cat("  @net$count 存在:", !is.null(cellchat@net$count), "\n")
    cat("  @net$weight 存在:", !is.null(cellchat@net$weight), "\n")
    cat("  @netP$net 存在:", !is.null(cellchat@netP$net), "\n")
    
    if(!is.null(cellchat@netP$net)) {
      cat("  @netP$net 中的键:", paste(names(cellchat@netP$net), collapse = ", "), "\n")
    }
    
    if(is.null(cellchat@netP$net) || length(cellchat@netP$net) == 0) {
      cat("\n⚠️  WARNING: @netP$net 为空，尝试从通路数据重新构建...\n")
      cellchat <- build_netP_from_communication(cellchat, df.netp)
      cat("重建后 @netP$net 的键:", paste(names(cellchat@netP$net), collapse = ", "), "\n")
    }
    
    cat("\n计算中心性...\n")
    cellchat <- calculate_centrality_from_net(cellchat)
    
    pathways.shows <- cellchat@netP$pathways
    head(cellchat@LR$LRsig)
    
    saveRDS(cellchat, file = paste0('cellchat.rds'))
  }
  
  # ========== 读取和绘图部分 ==========
  cellchat <- readRDS('cellchat.rds')
  groupSize <- as.numeric(table(cellchat@idents))
  pathways.shows <- cellchat@netP$pathways
  
  # 汇总【fig01】
  tryCatch({
    fig1 <- c("./fig01_Net_number_circle_barplot")
    dir.create(fig1, recursive = T, showWarnings = FALSE) 
    pdf(file = paste0(fig1, "/fig01_Net_number_strength.pdf"), width = 12, height = 7)
    par(mfrow = c(1,2))
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
    dev.off()
    if(file.exists(paste0(fig1, "/fig01_Net_number_strength.pdf"))) {
      file.copy(paste0(fig1, "/fig01_Net_number_strength.pdf"), "1.Net_number_strength.pdf", overwrite = TRUE)
    }
    
    df.net <- read.table(paste0("net_lr.txt"), sep = "\t", check.names = F, header = T)
    data <- as.data.frame(table(c(
      df.net$source,
      df.net$target
    )))
    colnames(data) <- c("Cell_Type", "all_sum")
    data <- data[order(data$all_sum, decreasing = T), ]
    data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)
    
    sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A") 
    
    if(length(levels(data$Cell_Type)) > length(sample_color)) {
      warning("Not enough sample_color defined for all cell types. Using default ggplot colors.")
      pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
        geom_bar(stat = "identity", width = 0.8) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x = "Cell Type", y = "Counts") +
        theme(axis.text.y = element_text(size = 12, colour = "black")) +
        theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
    } else {
      sample_color2 <- ifelse(data$all_sum > summary(data$all_sum)[2], 
                              ifelse(data$all_sum > summary(data$all_sum)[3], 
                                     ifelse(data$all_sum > summary(data$all_sum)[5], 
                                            sample_color[1], sample_color[2]), 
                                     sample_color[3]), 
                              sample_color[4])
      sample_color1 <- ifelse(data$all_sum > summary(data$all_sum)[2], 
                              ifelse(data$all_sum > summary(data$all_sum)[3], 
                                     ifelse(data$all_sum > summary(data$all_sum)[5], 
                                            alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), 
                                     alpha(sample_color[3], 0.9)), 
                              alpha(sample_color[4], 0.9))
      pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_manual(values = sample_color1) +
        scale_colour_manual(values = sample_color2) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x = "Cell Type", y = "Counts") +
        theme(axis.text.y = element_text(size = 12, colour = "black")) +
        theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
    }
    
    ggsave(paste0(fig1, "/fig01_Interaction Count.pdf"), plot = pB2, width = 10, height = 8)
    if(file.exists(paste0(fig1, "/fig01_Interaction Count.pdf"))) {
      file.copy(paste0(fig1, "/fig01_Interaction Count.pdf"), "2.Interaction Count.pdf", overwrite = TRUE)
    }
    
    if(!is.na(ShowCell)){
      ident1_idx <- which(levels(cellchat@idents) %in% ShowCell)
      ident2_idx <- which(!levels(cellchat@idents) %in% ShowCell)
      
      if(length(ident1_idx) == 0) {
        warning(paste("ShowCell '", ShowCell, "' not found in current cellchat object's cell types for single bubble plot. Skipping.", sep=""))
      } else {
        gg1 <- netVisual_bubble(cellchat, sources.use = ident1_idx, targets.use = ident2_idx, angle.x = 45, 
                                remove.isolate = F, font.size = 12, font.size.title = 15, return.data = T,
                                title.name = paste0('Signaling from ', ShowCell))
        gg2 <- netVisual_bubble(cellchat, sources.use = ident2_idx, targets.use = ident1_idx, angle.x = 45, 
                                remove.isolate = F, font.size = 12, font.size.title = 15, return.data = T,
                                title.name = paste0('Signaling to ', ShowCell))
        
        bubble_plot_width <- 5 + length(levels(cellchat@idents)) * 1
        bubble_plot_height <- 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05
        
        ggsave(paste0(fig1, "/fig01_Bubble.pdf"), 
               plot = gg1$gg.obj + gg2$gg.obj, 
               width = bubble_plot_width, 
               height = bubble_plot_height)
        if(file.exists(paste0(fig1, "/fig01_Bubble.pdf"))) {
          file.copy(paste0(fig1, "/fig01_Bubble.pdf"), "3.Bubble.pdf", overwrite = TRUE)
        }
      }
    }
    message("✓ Fig01 生成成功")
  }, error = function(e) {
    warning("Fig01 生成失败: ", e$message)
  })
  
  if(AllPlot){
    # 细胞圈图【fig02】
    tryCatch({
      fig2 <- c("./fig02_NetVisual_circle")
      dir.create(fig2, recursive = T, showWarnings = FALSE)
      
      mat <- cellchat@net$weight
      
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        
        pdf(file = paste0(fig2, '/fig02_netVisual_weight_', rownames(mat)[i], '.pdf'), 
            width = 10, height = 9)
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, 
                         edge.weight.max = max(mat), title.name = rownames(mat)[i])
        dev.off()
      }
      message("✓ Fig02 生成成功")
    }, error = function(e) {
      warning("Fig02 生成失败: ", e$message)
    })
    
    # pathway圈图【fig03】
    tryCatch({
      fig3 <- c("./fig03_NetVisual_aggregate")
      dir.create(fig3, recursive = T, showWarnings = FALSE)
      
      for (pathways.show in pathways.shows) {
        tryCatch({
          pdf(file = paste0(fig3, '/fig03_netVisual_aggregate_', pathways.show, '.pdf'), 
              width = 10, height = 9)
          netVisual_aggregate(cellchat,
                              signaling = pathways.show,
                              layout = "circle",
                              pt.title = 50,
                              vertex.label.cex = 1)        
          dev.off()
        }, error = function(e) {
          cat("警告: 通路 ", pathways.show, " 圈图生成失败: ", e$message, "\n", sep="")
        })
      }
      message("✓ Fig03 生成成功")
    }, error = function(e) {
      warning("Fig03 生成失败: ", e$message)
    })
    
    # pathway热图【fig04】
    tryCatch({
      fig4 <- c("./fig04_NetVisual_heatmap")
      dir.create(fig4, recursive = T, showWarnings = FALSE)
      for (pathways.show in pathways.shows) {
        tryCatch({
          pdf(file = paste0(fig4, '/fig04_netVisual_heatmap_', pathways.show, '.pdf'), 
              width = 10, height = 9)
          print(netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", 
                                  font.size = 14, font.size.title = 20))
          dev.off()
        }, error = function(e) {
          cat("警告: 通路 ", pathways.show, " 热图生成失败: ", e$message, "\n", sep="")
        })
      }
      message("✓ Fig04 生成成功")
    }, error = function(e) {
      warning("Fig04 生成失败: ", e$message)
    })
    
    # pathway贡献图【fig05】
    tryCatch({
      fig5 <- c("./fig05_NetAnalysis_contribution")
      dir.create(fig5, recursive = T, showWarnings = FALSE)
      for (pathways.show in pathways.shows) {
        tryCatch({
          res <- extractEnrichedLR(cellchat, signaling = pathways.show, 
                                   geneLR.return = TRUE, enriched.only = T)
          p <- netAnalysis_contribution(cellchat, signaling = pathways.show,
                                        font.size = 15, font.size.title = 15,
                                        title = paste0("Contribution of each L-R pair in ", pathways.show)) +
            theme(text = element_text(face = 'bold'))
          ggsave(filename = paste0(fig5, '/fig05_netAnalysis_contribution_', pathways.show, ".pdf"), 
                 plot = p, width = 10, height = 8)
        }, error = function(e) {
          cat("警告: 通路 ", pathways.show, " 贡献图生成失败: ", e$message, "\n", sep="")
        })
      }
      message("✓ Fig05 生成成功")
    }, error = function(e) {
      warning("Fig05 生成失败: ", e$message)
    })
    
    # 气泡图【fig06】
    tryCatch({
      fig6 <- c("./fig06_NetVisual_bubble")
      dir.create(fig6, recursive = T, showWarnings = FALSE)
      
      # from 
      for (sources.use in levels(cellchat@idents)) {
        tryCatch({
          targets.use <- levels(cellchat@idents)
          p <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, 
                                angle.x = 45, remove.isolate = F,
                                font.size = 10, font.size.title = 15, return.data = T,
                                title.name = paste0('Signaling from ', sources.use))
          ggsave(filename = paste0(fig6, '/fig06_netVisual_bubble_from_', sources.use, ".pdf"), 
                 plot = p$gg.obj, width = 8, height = 8 + nrow(p$communication) * 0.05)
        }, error = function(e) {
          cat("警告: ", sources.use, " from 气泡图生成失败\n", sep="")
        })
      }
      
      # to
      for (targets.use in levels(cellchat@idents)) {
        tryCatch({
          sources.use <- levels(cellchat@idents)
          p <- netVisual_bubble(cellchat, sources.use = sources.use, targets.use = targets.use, 
                                angle.x = 45, remove.isolate = F,
                                font.size = 10, font.size.title = 15, return.data = T,
                                title.name = paste0('Signaling to ', targets.use))
          ggsave(filename = paste0(fig6, '/fig06_netVisual_bubble_to_', targets.use, ".pdf"), 
                 plot = p$gg.obj, width = 8, height = 8 + nrow(p$communication) * 0.05)
        }, error = function(e) {
          cat("警告: ", targets.use, " to 气泡图生成失败\n", sep="")
        })
      }
      message("✓ Fig06 生成成功")
    }, error = function(e) {
      warning("Fig06 生成失败: ", e$message)
    })
    
    # 互作热图和散点图【fig07】
    tryCatch({
      fig7 <- c("./fig07_NetAnalysis_signalingRole")
      dir.create(fig7, recursive = T, showWarnings = FALSE)
      
      # 散点图
      tryCatch({
        gg1 <- netAnalysis_signalingRole_scatter(cellchat)
        ggsave(filename = paste0(fig7, '/fig07_netAnalysis_signalingRole_scatter.pdf'), 
               plot = gg1, width = 8, height = 8)
      }, error = function(e) {
        cat("警告: 散点图生成失败: ", e$message, "\n", sep="")
      })
      
      # 热图
      tryCatch({
        ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                                 height = 7 + length(cellchat@netP$pathways) * 0.125)
        ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", 
                                                 height = 7 + length(cellchat@netP$pathways) * 0.125)
        pdf(paste0(fig7, '/fig07_netAnalysis_signalingRole_heatmap.pdf'), 
            width = 15, height = 7 + length(cellchat@netP$pathways) * 0.125)
        print(ht1 + ht2)
        dev.off()
      }, error = function(e) {
        cat("警告: 热图生成失败: ", e$message, "\n", sep="")
      })
      
      message("✓ Fig07 生成成功")
    }, error = function(e) {
      warning("Fig07 生成失败: ", e$message)
    })
    
    # 受配体小提琴图【fig08】
    tryCatch({
      fig8 <- c("./fig08_PlotGeneExpression")
      dir.create(fig8, recursive = T, showWarnings = FALSE)
      for (pathways.show in pathways.shows) {
        tryCatch({
          res <- extractEnrichedLR(cellchat, signaling = pathways.show, 
                                   geneLR.return = TRUE, enriched.only = T)
          p <- plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = TRUE, type = "violin")
          ggsave(filename = paste0(fig8, '/fig08_PlotGeneExpression_', pathways.show, ".pdf"), 
                 plot = p, width = 10, height = 7 + length(res$geneLR) * 0.225)
        }, error = function(e) {
          cat("警告: 通路 ", pathways.show, " 小提琴图生成失败: ", e$message, "\n", sep="")
        })
      }
      message("✓ Fig08 生成成功")
    }, error = function(e) {
      warning("Fig08 生成失败: ", e$message)
    })
  }
  
  return(cellchat)
}

# cellchat_module 函数 (修改了 Signaling Role Scatter Plot 的生成方式)
cellchat_module <- function(object, OutPath, group_by = "all", celltype = "cellType_1", 
                            idents = NULL, ShowCell = NA, species = NA, OnlyPlot = F, AllPlot = F){
  
  if (class(object) != 'Seurat') {
    stop('object should be a SeuratObject')
  }
  
  # ⭐ 第一步：保存原始工作目录
  original_wd <- getwd()
  cat("原始工作目录:", original_wd, "\n")
  
  # ⭐ 第二步：将 OutPath 转换为绝对路径
  OutPath <- normalizePath(OutPath, winslash = "/", mustWork = FALSE)
  cat("输出路径（绝对）:", OutPath, "\n")
  
  # ⭐ 第三步：添加错误处理和最终清理
  tryCatch({
    
    if(group_by == "all" & is.null(idents)){
      
      out_path_single <- file.path(OutPath, "CellChat")
      if(!dir.exists(out_path_single)) {
        dir.create(out_path_single, recursive = T)
      }
      setwd(out_path_single)
      cellchat <- cellchat_func(object = object, species = species, celltype = celltype, 
                                ShowCell = ShowCell, OnlyPlot = OnlyPlot, AllPlot = AllPlot)
      
    } else if(group_by != "all" & length(idents) == 2){
      
      message("正在独立分析每个组的所有细胞类型...")
      
      # ⭐ 使用绝对路径
      cellchat_base_dir <- file.path(OutPath)  # 确保是绝对路径
      
      # 验证和创建目录
      if(!dir.exists(cellchat_base_dir)) {
        dir.create(cellchat_base_dir, recursive = T, showWarnings = FALSE)
      }
      
      cat("CellChat 基础目录:", cellchat_base_dir, "\n")
      
      # ⭐ 只在这里 setwd() 一次，到基础目录
      setwd(cellchat_base_dir)
      
      all_cell_types_union <- levels(droplevels(object@meta.data[[celltype]]))
      
      # 独立运行每个组的 CellChat 分析
      for (ident_name in idents) { 
        # ⭐ 使用绝对路径，不依赖当前工作目录
        current_group_output_dir <- file.path(cellchat_base_dir, ident_name)
        
        if(!dir.exists(current_group_output_dir)) {
          dir.create(current_group_output_dir, recursive = T, showWarnings = FALSE)
        }
        
        cat("为组 '", ident_name, "' 创建目录: ", current_group_output_dir, "\n", sep="")
        
        # ⭐ setwd() 到该组的目录
        setwd(current_group_output_dir)
        
        sub.object <- object[, object@meta.data[[group_by]] %in% ident_name]
        
        message(paste("正在为组:", ident_name, "运行 CellChat，包含", ncol(sub.object), "个细胞。"))
        cellchat_func(object = sub.object, species = species, celltype = celltype, 
                      ShowCell = ShowCell, OnlyPlot = OnlyPlot, AllPlot = AllPlot)
      }
      
      # ⭐ 返回到基础目录，而不是相对路径切换
      setwd(cellchat_base_dir)
      
      # 读取各个组的 CellChat 对象（使用绝对路径）
      cellchat.group1 <- readRDS(file.path(cellchat_base_dir, idents[1], "cellchat.rds"))
      cellchat.group2 <- readRDS(file.path(cellchat_base_dir, idents[2], "cellchat.rds"))
      
      # 后续分析...
      cellchat.group1 <- calculate_centrality_from_net(cellchat.group1)
      cellchat.group2 <- calculate_centrality_from_net(cellchat.group2)
      
      message("正在协调 CellChat 对象...")
      cellchat.group1_harmonized <- harmonize_cellchat_idents(cellchat.group1, all_cell_types_union)
      cellchat.group2_harmonized <- harmonize_cellchat_idents(cellchat.group2, all_cell_types_union)
      
      object.list <- setNames(list(cellchat.group1_harmonized, cellchat.group2_harmonized), idents[1:2])
      message("正在合并 CellChat 对象...")
      cellchat <- mergeCellChat(object.list, add.names = names(object.list))
      save(cellchat, file = file.path(cellchat_base_dir, paste0("cellchat_merged_", idents[1], "_", idents[2], ".RData")))
      
      message("正在修复合并后的网络结构...")
      cellchat <- repair_merged_cellchat_network(cellchat, object.list)
      
      # 创建比较结果输出目录
      dir <- file.path(cellchat_base_dir, "comparison")  # ⭐ 使用绝对路径
      dir.create(dir, recursive = T, showWarnings = FALSE) 
      
      # ⭐ 后续所有文件操作都使用绝对路径
      message("正在生成比较图表...")
      
      message("重建合并对象的 @idents...")
      cellchat <- rebuild_merged_idents(cellchat, object.list)
      
      # 1. 总体相互作用比较
      tryCatch({
        gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2))
        gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1, 2), measure = "weight")
        ggsave(file.path(dir, "01.CompareInteractions.pdf"), plot = gg1 + gg2, width = 8, height = 6)
        message("✓ 生成了 01.CompareInteractions.pdf")
      }, error = function(e) {
        warning("生成 CompareInteractions 时出错: ", e$message)
      })
      
      # 2. 网络可视化 (差异)
      tryCatch({
        cat("\n生成网络差异可视化...\n")
        
        # 计算 count 的最大值
        count.max <- sapply(object.list, function(x) {
          max(x@net$count)
        })
        
        # 计算 weight 的最大值
        weight.max <- sapply(object.list, function(x) {
          max(x@net$weight)
        })
        
        # 验证
        cat("Count 最大值:", count.max, "\n")
        cat("Weight 最大值:", weight.max, "\n")        
        # ========== Count 版本 ==========
        cat("  生成交互数量对比...\n")
        pdf(file.path(dir, "02.NetVisual_count.pdf"), width = 16, height = 8)
        par(mfrow = c(1, 2), xpd = TRUE)
        
        for (i in 1:length(object.list)) {
          current_groupSize <- as.numeric(table(object.list[[i]]@idents))
          netVisual_circle(object.list[[i]]@net$count, 
                           vertex.weight = current_groupSize,
                           weight.scale = T, 
                           label.edge = F,
                           edge.weight.max = count.max[i],    # ✅ 改成 count.max[i]！
                           edge.width.max = 12,
                           title.name = paste0("Number of Interactions - ", names(object.list)[i]))
        }
        
        par(mfrow = c(1, 1))
        dev.off()
        message("✓ 生成了 02.NetVisual_count.pdf")
        
        # ========== Weight 版本 ==========
        cat("  生成交互强度对比...\n")
        pdf(file.path(dir, "03.NetVisual_weight.pdf"), width = 16, height = 8)
        par(mfrow = c(1, 2), xpd = TRUE)
        
        for (i in 1:length(object.list)) {
          current_groupSize <- as.numeric(table(object.list[[i]]@idents))
          netVisual_circle(object.list[[i]]@net$weight, 
                           vertex.weight = current_groupSize,
                           weight.scale = T, 
                           label.edge = F,
                           edge.weight.max = weight.max[i],    # ✅ 改成 weight.max[i]！
                           edge.width.max = 12,
                           title.name = paste0("Interaction Weights - ", names(object.list)[i]))
        }
        
        par(mfrow = c(1, 1))
        dev.off()
        message("✓ 生成了 03.NetVisual_weight.pdf")
        
      }, error = function(e) {
        cat("❌ 生成失败:", e$message, "\n")
      })
      
      # 3. 热图比较
      tryCatch({
        pdf(file.path(dir, "04.NetVisualHeatmap.pdf"), width = 12, height = 6)
        gg1 <- netVisual_heatmap(cellchat)
        gg2 <- netVisual_heatmap(cellchat, measure = "weight")
        print(gg1 + gg2)
        dev.off()
        message("✓ 生成了 04.NetVisualHeatmap.pdf")
      }, error = function(e) {
        warning("生成 NetVisualHeatmap 失败，跳过此步骤")
        cat("原因:", e$message, "\n")
      })
      
      
      
      # 4. 信号角色散点图 (改进版 - 带完整错误检查)
      # 改进的散点图生成函数
      message("开始生成信号角色散点图...")
      
      # 获取所有细胞类型
      all_cell_types_global <- unique(unlist(lapply(object.list, function(x) {
        rownames(x@net$count)
      })))
      
      gg <- list()
      
      for (i in 1:length(object.list)) {
        tryCatch({
          current_obj <- object.list[[i]]
          group_name <- names(object.list)[i]
          
          cat(paste0("\n--- 为组 '", group_name, "' 生成散点图 ---\n"))
          
          # 从 @net$weight 或 @net$count 直接计算
          if (!is.null(current_obj@net$weight) && sum(current_obj@net$weight) > 0) {
            net_matrix <- current_obj@net$weight
          } else if (!is.null(current_obj@net$count) && sum(current_obj@net$count) > 0) {
            net_matrix <- current_obj@net$count
          } else {
            warning(paste("组", group_name, "没有网络数据"))
            gg[[i]] <- NULL
            next
          }
          
          cat(paste0("  网络矩阵大小: ", nrow(net_matrix), " x ", ncol(net_matrix), "\n"))
          cat(paste0("  网络矩阵总和: ", sum(net_matrix, na.rm = TRUE), "\n"))
          
          # 计算 outgoing (行求和) 和 incoming (列求和)
          outgoing_all <- rowSums(net_matrix, na.rm = TRUE)
          incoming_all <- colSums(net_matrix, na.rm = TRUE)
          
          # 确保所有细胞类型都在矩阵中
          outgoing_final <- setNames(rep(0, length(all_cell_types_global)), all_cell_types_global)
          incoming_final <- setNames(rep(0, length(all_cell_types_global)), all_cell_types_global)
          
          # 填充有数据的细胞类型
          common_cells_out <- intersect(names(outgoing_all), all_cell_types_global)
          common_cells_in <- intersect(names(incoming_all), all_cell_types_global)
          
          outgoing_final[common_cells_out] <- outgoing_all[common_cells_out]
          incoming_final[common_cells_in] <- incoming_all[common_cells_in]
          
          cat(paste0("  Outgoing cells: ", paste(names(outgoing_final[outgoing_final > 0]), collapse = ", "), "\n"))
          cat(paste0("  Outgoing sum: ", sum(outgoing_final), "\n"))
          cat(paste0("  Incoming sum: ", sum(incoming_final), "\n"))
          
          # 构造数据框
          df_plot <- data.frame(
            x = outgoing_final,
            y = incoming_final,
            labels = names(outgoing_final),
            stringsAsFactors = FALSE,
            row.names = NULL
          )
          
          df_plot$activity.total <- df_plot$x + df_plot$y
          
          # 移除所有活动都为0的细胞类型
          df_plot <- df_plot[df_plot$activity.total > 0, ]
          
          cat(paste0("  有活动的细胞类型数: ", nrow(df_plot), "\n"))
          
          if (nrow(df_plot) == 0) {
            warning(paste("没有发现信号活动，跳过组:", group_name))
            gg[[i]] <- NULL
            next
          }
          
          # 绘图
          p <- ggplot(df_plot, aes(x = x, y = y, label = labels, color = activity.total)) +
            geom_point(aes(size = activity.total), alpha = 0.8) +
            ggrepel::geom_text_repel(
              size = 3.5, 
              max.overlaps = Inf,
              box.padding = unit(0.5, "lines"),
              point.padding = unit(0.3, "lines"),
              segment.color = "grey50",
              segment.alpha = 0.5
            ) +
            scale_color_gradientn(
              colors = RColorBrewer::brewer.pal(n = 10, name = "RdYlBu"),
              name = "Total Activity"
            ) +
            scale_size_area(max_size = 15, name = "Total Activity") +
            labs(
              x = "Outgoing Signaling",
              y = "Incoming Signaling",
              title = paste0("Signaling Role - ", group_name)
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 10),
              legend.position = "right",
              panel.border = element_rect(fill = NA, colour = "black", size = 0.5)
            )
          
          gg[[i]] <- p
          cat(paste0("  ✓ 散点图生成成功\n"))
          
        }, error = function(e) {
          cat(paste0("  ✗ 错误: ", e$message, "\n"))
          gg[[i]] <<- NULL
        })
      }
      
      # 过滤掉 NULL 元素
      gg <- gg[!sapply(gg, is.null)]
      
      cat(paste0("\n成功生成的散点图数: ", length(gg), "\n"))
      
      if (length(gg) > 0) {
        # 统一坐标轴范围
        all_x_vals <- c()
        all_y_vals <- c()
        
        for (plot in gg) {
          tryCatch({
            plot_data <- ggplot_build(plot)$data[[1]]
            all_x_vals <- c(all_x_vals, plot_data$x)
            all_y_vals <- c(all_y_vals, plot_data$y)
          }, error = function(e) {
            cat(paste0("警告: 无法提取坐标: ", e$message, "\n"))
          })
        }
        
        # 清理数据
        all_x_vals <- all_x_vals[!is.na(all_x_vals) & is.finite(all_x_vals)]
        all_y_vals <- all_y_vals[!is.na(all_y_vals) & is.finite(all_y_vals)]
        
        if (length(all_x_vals) > 0 && length(all_y_vals) > 0) {
          max_x <- max(all_x_vals, na.rm = TRUE) * 1.15
          max_y <- max(all_y_vals, na.rm = TRUE) * 1.15
          
          cat(paste0("坐标轴范围 - X: [0, ", round(max_x, 2), "], Y: [0, ", round(max_y, 2), "]\n"))
          
          gg <- lapply(gg, function(p) {
            p + xlim(0, max_x) + ylim(0, max_y)
          })
        }
        
        # 确定布局
        n_plots <- length(gg)
        ncol_layout <- min(2, n_plots)
        nrow_layout <- ceiling(n_plots / ncol_layout)
        
        cat(paste0("布局: ", nrow_layout, " 行 x ", ncol_layout, " 列\n"))
        
        combined_plot <- patchwork::wrap_plots(
          plots = gg,
          ncol = ncol_layout,
          nrow = nrow_layout
        )
        
        pdf_width <- 6 * ncol_layout
        pdf_height <- 6 * nrow_layout
        
        ggsave(
          file.path(dir, "05.SignalingRole.pdf"),
          plot = combined_plot,
          width = pdf_width,
          height = pdf_height,
          dpi = 300
        )
        
        message("✓ 信号角色散点图已保存到: ", file.path(dir, "05.SignalingRole.pdf"))
        cat(paste0("  图片尺寸: ", pdf_width, " x ", pdf_height, " inches\n"))
        
      } else {
        warning("没有成功生成任何散点图!")
      }
      # 6. 气泡图比较（改进版）
      if(!is.na(ShowCell)) {
        tryCatch({
          # 重建后的 @idents 可能仍有问题，直接从 object.list 获取细胞类型
          all_celltypes <- unique(unlist(lapply(object.list, function(x) levels(x@idents))))
          ident1_idx_merged <- which(all_celltypes %in% ShowCell)
          ident2_idx_merged <- which(!all_celltypes %in% ShowCell)
          
          if(length(ident1_idx_merged) == 0) {
            warning(paste("ShowCell '", ShowCell, "'未在细胞类型中找到。", sep=""))
          } else {
            gg1 <- netVisual_bubble(cellchat, sources.use = ident1_idx_merged, targets.use = ident2_idx_merged, 
                                    comparison = c(1, 2), angle.x = 45, remove.isolate = F,
                                    font.size = 12, font.size.title = 15, return.data = T,
                                    title.name = paste0('Signaling from ', ShowCell))
            gg2 <- netVisual_bubble(cellchat, sources.use = ident2_idx_merged, targets.use = ident1_idx_merged, 
                                    comparison = c(1, 2), angle.x = 45, remove.isolate = F,
                                    font.size = 12, font.size.title = 15, return.data = T,
                                    title.name = paste0('Signaling to ', ShowCell))
            
            bubble_plot_width <- 5 + length(all_celltypes) * 1 
            bubble_plot_height <- 8 + max(nrow(gg1$communication), nrow(gg2$communication)) * 0.05
            
            ggsave(file.path(dir, "6.Bubble.pdf"), plot = gg1$gg.obj + gg2$gg.obj, 
                   width = bubble_plot_width, height = bubble_plot_height)
            message("✓ 生成了 6.Bubble.pdf")
          }
        }, error = function(e) {
          warning("生成 Bubble 失败: ", e$message)
        })
      }  # ← ⭐ 关闭 if(!is.na(ShowCell))
      
    }  # ← ⭐ 关闭 else if(group_by != "all" & length(idents) == 2)
    
  }, error = function(e) {  # ← ⭐ 最外层 tryCatch 的 error
    cat("❌ 错误发生:", e$message, "\n")
    cat("正在恢复工作目录...\n")
    setwd(original_wd)
    stop(e)
  }, finally = {
    setwd(original_wd)
    cat("✓ 已恢复原始工作目录:", original_wd, "\n")
  })  # ← ⭐ 关闭最外层 tryCatch
  
  return(cellchat)  # ← ⭐ 只有一个 )
}  # ← ⭐ 关闭函数


# AUCcell -----------------------------------------------------------------


# --- Function: advanced_filter_gene_sets ---
advanced_filter_gene_sets <- function(gene_set, 
                                      scRNA_data,      # Seurat对象的RNA表达数据
                                      top_percentage = 0.3,  # 选择Top 30%的基因
                                      min_expression_threshold = 0.1) {
  # 1. 确保基因在数据集中存在
  available_genes <- intersect(gene_set, rownames(scRNA_data))
  
  if(length(available_genes) == 0) {
    return(NULL)
  }
  
  # 2. 计算每个基因在数据集中的平均表达
  gene_mean_expr <- rowMeans(scRNA_data[available_genes, ])
  
  # 3. 按平均表达排序
  sorted_genes <- names(sort(gene_mean_expr, decreasing = TRUE))
  
  # 4. 选择Top N%的基因
  top_n <- max(round(length(sorted_genes) * top_percentage), 1)
  top_genes <- sorted_genes[1:top_n]
  
  # 5. 过滤掉表达非常低的基因
  high_expr_genes <- top_genes[gene_mean_expr[top_genes] > min_expression_threshold]
  
  return(high_expr_genes)
}


# --- Function: extract_auc_threshold ---
extract_auc_threshold <- function(assignment) {
  # 检查assignment是否为NULL或不包含预期结构
  if (is.null(assignment) || 
      length(assignment) == 0 || 
      !is.list(assignment)) {
    return(NULL)
  }
  
  # 提取具体的基因集名称（因为有嵌套结构）
  gene_set_name <- names(assignment)
  
  # 如果gene_set_name为空，返回NULL
  if (length(gene_set_name) == 0) {
    return(NULL)
  }
  
  # 获取具体的assignment
  specific_assignment <- assignment[[gene_set_name]]
  
  # 提取阈值
  auc_thr <- specific_assignment$aucThr$selected
  
  return(auc_thr)
}


# ggsave ------------------------------------------------------------------


# --- Function: ggsave_fun ---
ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}



# Harmony.integration.reduceDimension -------------------------------------


# --- Function: Harmony.integration.reduceDimension ---
Harmony.integration.reduceDimension <- function(seurat.object, set.resolutions, assay = "RNA", nfeatures = 3000, PC = 50, npcs = 100){
  # require(Seurat)
  require(harmony)
  require(clustree)
  # require(dplyr)
  #OK nfeatures
  if (!dir.exists("02.Cluster")) {
    dir.create("02.Cluster")
  } else {
    print("Dir already exists!")
  }
  
  DefaultAssay(seurat.object) <- assay
  
  seurat.object <- NormalizeData(object = seurat.object)
  seurat.object <- FindVariableFeatures(object = seurat.object, selection.method = "vst", nfeatures = nfeatures)
  top10 <- head(VariableFeatures(seurat.object), 10)
  plot1 <- VariableFeaturePlot(seurat.object)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  ggsave_fun(filename = paste0("02.Cluster/01.VariableFeaturePlot"),plot = plot1 + plot2,width = 14)
  seurat.object <- ScaleData(object = seurat.object)
  seurat.object <- RunPCA(seurat.object, verbose = T, npcs = npcs,features = VariableFeatures(object = seurat.object))
  p <- ElbowPlot(object = seurat.object, ndims = npcs)
  ggsave_fun(filename = paste0("02.Cluster/02.ElbowPlot"),plot = p)
  #If dims.use is not specified, all PCs will be used by default
  #assay.use defaults to RNA. If the SCTransform standardization method is used, you need to specify assay.use="SCT" in RunHarmony
  #Compare the differences between the two modes and find that the default is better and more complete
  seurat.object <- RunHarmony(object = seurat.object, group.by.vars = "orig.ident", assay.use=assay, verbose = FALSE)
  seurat.object <- RunUMAP(seurat.object, reduction = "harmony", dims = 1:PC, verbose = T)
  seurat.object <- RunTSNE(seurat.object, reduction = "harmony", dims = 1:PC, verbose = T)
  seurat.object <- FindNeighbors(seurat.object, dims = 1:PC, reduction = "harmony", verbose = T) #使用harmony替代PCA
  seurat.object <- FindClusters(seurat.object, resolution = set.resolutions, verbose = T)
  p <- DimPlot(object = seurat.object, reduction = 'pca',label = F, group.by = "orig.ident")
  ggsave_fun(filename = paste0("02.Cluster/03.PcaPlot"),plot = p)
  p <- DimPlot(object = seurat.object, reduction = 'harmony',label = F, group.by = "orig.ident")
  ggsave_fun(filename = paste0("02.Cluster/04.HarmonyPlot"),plot = p)
  pdf("02.Cluster/05.Data.merge.harmony.pdf")
  p <- clustree(seurat.object)
  print(p)
  p <- DimPlot(object = seurat.object, reduction = 'umap',label = TRUE, group.by = "orig.ident")
  print(p)
  merge.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = seurat.object, reduction = 'umap',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
    print(p)
  })
  p <- DimPlot(object = seurat.object, reduction = 'tsne',label = TRUE, group.by = "orig.ident")
  print(p)
  merge.res <- sapply(set.resolutions, function(x){
    p <- DimPlot(object = seurat.object, reduction = 'tsne',label = TRUE, group.by = paste0(assay, "_snn_res.", x)) + NoLegend()
    print(p)
  })
  dev.off()
  return(seurat.object)
}




# human2mouse -------------------------------------------------------------


# --- Function: human2mouse ---
human2mouse <- function(x) {
  x1 <- substr(x, 1, 1)
  x2 <- substr(x, 2, nchar(x))
  return(paste0(x1, tolower(x2)))
}



# KeyGenePipeline --------------------------------------------------------------------


# --- Function: gseaplot3 ---
gseaplot3 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                       rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                       ES_geom = "line"){
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- enrichplot:::gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8),  
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
                                     plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                     axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") {
    title <- unlist(strsplit(title, split = "_"))
    # title2 <- paste(str_to_title(title[2:length(title)]), collapse = " ")
    # title <- paste(title[1],title2,collapse = " ") %>% stringr::str_wrap(., width = 40)
    title <- paste(title,collapse = " ") %>% stringr::str_wrap(., width = 40)
    p.res <- p.res + ggtitle(title)
  }
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust",'enrichmentScore','NES')]
    rownames(pd) <- ''
    pd <- pd[, -1]
    colnames(pd) <- c( "pvalue", "p.adjust",'ES','NES')
    pd <- round(pd, 4)
    #tp <- enrichplot:::tableGrob2(pd, p.res)
    tp <- enrichplot:::tableGrob2(pd)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), 
                                                                         xmax = quantile(p.res$data$x, 0.95), 
                                                                         ymin = quantile(p.res$data$runningScore,0.75), 
                                                                         ymax = quantile(p.res$data$runningScore,0.9))
  }
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}


# --- Function: gsea_dotplot ---
gsea_dotplot <- function(object,x = "geneRatio",color = "p.adjust",showCategory = 10,
                         size = NULL,split = T,title = "",
                         col_low = "#D15047", col_high =  "#3793DE", scales = 'free_x',
                         orderBy = "x",decreasing = TRUE,y_wrap = 40) {
  colorBy <- match.arg(color, c("pvalue", "p.adjust", "qvalue"))
  if (x == "geneRatio" || x == "GeneRatio") {
    x <- "GeneRatio"
    if (is.null(size)) 
      size <- "Count"
  }else if (x == "count" || x == "Count") {
    x <- "Count"
    if (is.null(size)) 
      size <- "GeneRatio"
  }else if (is(x, "formula")) {
    x <- as.character(x)[2]
    if (is.null(size)) 
      size <- "Count"
  }else {
    if (is.null(size)) 
      size <- "Count"
  }
  
  df <- fortify(object, showCategory = showCategory, split = ".sign")
  if (orderBy != "x" && !orderBy %in% colnames(df)) {
    message("wrong orderBy parameter; set to default `orderBy = \"x\"`")
    orderBy <- "x"
  }
  if (orderBy == "x") {
    df <- dplyr::mutate(df, x = eval(parse(text = x)))
  }
  
  label_func <- function(str,label_format=y_wrap){str <- gsub("_", " ", str);stringr::str_wrap(str, label_format)}
  
  idx <- order(df[[orderBy]], decreasing = T)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))
  result <<- df
  p <- ggplot(df, aes_string(x = x, y = "Description", size = size, fill = colorBy)) + 
    geom_point(shape=21) + 
    scale_fill_continuous(low = col_low, high =  col_high, name = color) + 
    scale_y_discrete(labels = label_func) + ylab(NULL) + 
    ggtitle(title) + DOSE::theme_dose(12) + scale_size(range = c(3,9))+
    guides(fill = guide_colorbar(order = 1,reverse = TRUE),
           size = guide_legend(order = 2,override.aes = list(fill='grey40')))+
    theme(text = element_text(face = 'bold',color = 'black'),
          axis.line.x.top  = element_line(size=0.8),
          axis.ticks = element_blank(),
          strip.background = element_rect(color = 'black',size = 1.3),
          panel.border = element_rect(color = 'black',size = 1.3),
          panel.grid.major = element_line(linetype = 5),
          panel.grid.minor = element_blank())
  
  if(split){
    p <- p + facet_wrap(~.sign,scales=scales) + theme(strip.text=element_text(size = 15,face = 'bold'))
  }
  return(p)
}


# --- Function: getScatterplot ---
getScatterplot <- function(object, gene1, gene2, cor.method = "pearson", jitter.num = 0.15, pos = TRUE) {
  if (!gene1 %in% rownames(object)) {
    print("gene1 was not found")
    if (!gene2 %in% rownames(object)) {
      print("gene2 was not found")
    }
  } else {
    exp.mat <- GetAssayData(object = object, assay = "RNA") %>%
      .[c(gene1, gene2), ] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()
    if (pos) {
      if (nrow(exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]) > (nrow(exp.mat) * 0.01)) {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]
      } else {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 | exp.mat[, 2] > 0), ]
      }
    }
    colnames(exp.mat) <- c("Var1", "Var2")
    plots <- ggplot(data = exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) +
      geom_smooth(method = "lm", se = T, color = "red", size = 1) +
      stat_cor(method = cor.method) +
      labs(x = gene1, y = gene2) +
      geom_jitter(width = jitter.num, height = jitter.num, color = "black", size = 1, alpha = 1) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.text = element_text(colour = "black", size = 10),
        axis.text = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(size = 1, linetype = "solid", colour = "black"),
        panel.background = element_rect(fill = "white")
      )
    return(plots)
  }
}


# --- Function: gather_graph_node ---
gather_graph_node <- function(df, index = NULL, value = tail(colnames(df), 1), root = NULL) {
  require(dplyr)
  if (length(index) < 2) {
    stop("please specify at least two index column(s)")
  } else {
    list <- lapply(seq_along(index), function(i) {
      dots <- index[1:i]
      df %>%
        group_by(.dots = dots) %>%
        summarise(
          node.size = sum(.data[[value]]),
          node.level = index[[i]],
          node.count = n()
        ) %>%
        mutate(
          node.short_name = as.character(.data[[dots[[length(dots)]]]]),
          node.branch = as.character(.data[[dots[[1]]]])
        ) %>%
        tidyr::unite(node.name, dots, sep = "/")
    })
    data <- do.call("rbind", list) %>% as_tibble()
    data$node.level <- factor(data$node.level, levels = index)
    
    if (is.null(root)) {
      return(data)
    } else {
      root_data <- data.frame(
        node.name = root,
        node.size = sum(df[[value]]),
        node.level = root,
        node.count = 1,
        node.short_name = root,
        node.branch = root,
        stringsAsFactors = F
      )
      data <- rbind(root_data, data)
      data$node.level <- factor(data$node.level, levels = c(root, index))
      return(data)
    }
  }
}


# --- Function: gather_graph_edge ---
gather_graph_edge <- function(df, index = NULL, root = NULL) {
  require(dplyr)
  if (length(index) < 2) {
    stop("please specify at least two index column(s)")
  } else if (length(index) == 2) {
    data <- df %>%
      mutate(from = .data[[index[[1]]]]) %>%
      tidyr::unite(to, index, sep = "/") %>%
      dplyr::select(from, to) %>%
      mutate_at(c("from", "to"), as.character)
  } else {
    list <- lapply(seq(2, length(index)), function(i) {
      dots <- index[1:i]
      df %>%
        tidyr::unite(from, dots[-length(dots)], sep = "/", remove = F) %>%
        tidyr::unite(to, dots, sep = "/") %>%
        dplyr::select(from, to) %>%
        mutate_at(c("from", "to"), as.character)
    })
    data <- do.call("rbind", list)
  }
  data <- as_tibble(data)
  if (is.null(root)) {
    return(data)
  } else {
    root_data <- df %>%
      group_by(.dots = index[[1]]) %>%
      summarise(count = n()) %>%
      mutate(from = root, to = as.character(.data[[index[[1]]]])) %>%
      dplyr::select(from, to)
    rbind(root_data, data)
  }
}


# --- Function: CoreAlg ---
CoreAlg <- function(X, y) {
  # try different values of nu
  svn_itor <- 3
  
  res <- function(i) {
    if (i == 1) {
      nus <- 0.25
    }
    if (i == 2) {
      nus <- 0.5
    }
    if (i == 3) {
      nus <- 0.75
    }
    model <- svm(X, y, type = "nu-regression", kernel = "linear", nu = nus, scale = F)
    model
  }
  
  if (Sys.info()["sysname"] == "Windows") {
    out <- mclapply(1:svn_itor, res, mc.cores = 1)
  } else {
    out <- mclapply(1:svn_itor, res, mc.cores = svn_itor)
  }
  
  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)
  
  # do cibersort
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  # pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  # get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  w <- (q / sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}


# --- Function: doPerm ---
doPerm <- function(perm, X, Y) {
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while (itor <= perm) {
    # print(itor)
    
    # random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
    
    # standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    # run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    # store correlation
    if (itor == 1) {
      dist <- mix_r
    } else {
      dist <- rbind(dist, mix_r)
    }
    
    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}


# --- Function: CIBERSORT ---
CIBERSORT <- function(dir, sig_matrix, mixture_file, perm = 0, QN = TRUE) {
  library(e1071)
  library(parallel)
  library(preprocessCore)
  
  # read in data
  X <- read.table(sig_matrix, header = T, sep = "\t", row.names = 1, check.names = F)
  Y <- mixture_file
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  # order
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  
  P <- perm # number of permutations
  
  # anti-log if max < 50 in mixture file
  if (max(Y) < 50) {
    Y <- 2^Y
  }
  
  # quantile normalization of mixture file
  if (QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  if (substr(Sys.Date(), 6, 7) > 5) {
    next
  }
  
  # intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, ]
  
  # standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  # empirical null distribution of correlation coefficients
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y)$dist)
  }
  
  # print(nulldist)
  
  header <- c("Mixture", colnames(X), "P-value", "Correlation", "RMSE")
  # print(header)
  
  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  # iterate through mixtures
  while (itor <= mixtures) {
    y <- Y[, itor]
    
    # standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    # run SVR core algorithm
    result <- CoreAlg(X, y)
    
    if (substr(Sys.Date(), 1, 4) > 3021) {
      next
    }
    
    # get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    # calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }
    
    # print output
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }
    
    itor <- itor + 1
  }
  
  # save results
  write.table(rbind(header, output), file = paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
  
  # return matrix object containing all results
  obj <- rbind(header, output)
  obj <- obj[, -1]
  obj <- obj[-1, ]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  obj
}


# --- Function: CorMatrix ---
CorMatrix <- function(cor, p) {
  ut <- upper.tri(cor)
  data.frame(
    row = rownames(cor)[row(cor)[ut]],
    column = rownames(cor)[col(cor)[ut]],
    cor = (cor)[ut],
    p = p[ut]
  )
}


# --- Function: capitalize_first ---
capitalize_first <- function(string) {
  words <- strsplit(string, " ")[[1]]
  capitalized_words <- toupper(substring(words, 1, 1))
  rest_of_words <- substring(words, 2)
  capitalized_string <- paste0(capitalized_words, rest_of_words, collapse = " ")
  return(capitalized_string)
}


# --- Function: KeyGenePipeline ---
# ==================== 核心包 ====================
library(limma)          # 差异表达分析
library(clusterProfiler) # GSEA/富集分析
library(org.Hs.eg.db)  # 人类基因数据库
#ibrary(org.Mm.eg.db)  # 小鼠基因数据库
library(GSVA)          # GSVA分析
library(GSEABase)      # GSEA数据格式
# ==================== 数据处理 ====================
library(dplyr)         # 数据框操作
library(tidyr)         # 数据整形
library(reshape2)      # reshape功能
# ==================== 相关性分析 ====================
library(Hmisc)         # rcorr函数
# ==================== 绘图 ====================
library(ggplot2)       # 基础绘图
library(ggrepel)       # 文字排斥
library(cowplot)       # 拼图
library(gridExtra)     # 网格排列
library(corrplot)      # 相关性热图
library(aplot)         # insert_top功能
#library(ggprism)       # 美化主题
library(ggsignif)      # 显著性标记
library(ggpubr)        # 比较检验
# ==================== 网络/可视化 ====================
library(igraph)        # 网络图
library(tidygraph)     # tbl_graph
library(ggraph)        # ggraph绘图
# ==================== GOplot ====================
library(GOplot)        # circle_dat, gather_graph_node等
# ==================== 主函数 ====================
library(stringr)       # 字符串处理
library(ComplexHeatmap) # 复杂热图
KeyGenePipeline <- function(Path = "", mat, species = NA, KeyGene = NA,all.pathway = F) {
  cat("1.GSEA\n")
  cat("2.GSVA\n")
  cat("3.免疫浸润\n")
  cat("4.免疫调控因子\n")
  cat("5.差异表达\n")
  cat("6.疾病相关性\n")
  cat("7.Motif转录因子\n")
  cat("8.非编码 RNA 网络\n")
  cat("\033[0;34m请选择要分析的内容，按序号填写，英文逗号分隔")
  AnalysisIndex <- scan(sep = ",", quiet = TRUE)
  
  if (!is.na(species)) {
    mouse <- ifelse(species == "mouse", T, F)
  } else {
    mouse <- F
  }
  
  mat <- as.matrix(mat)
  n <- 0
  
  for (step in AnalysisIndex) {
    n <- n + 1
    # GSEA
    if (step == 1) {
      message("*************************GSEA*************************")
      
      dir <- paste0(Path, "/.GSEA")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      cat("1.KEGG\n")
      cat("2.GO\n")
      cat("\033[0;34m请选择要分析的内容，按序号填写")
      isKEGGIndex <- scan(what = "numeric",quiet = TRUE)
      if(isKEGGIndex == 1){
        isKEGG <- TRUE
      }else{
        isKEGG <- FALSE
      }
      
      for (GeneIndex in 1:length(KeyGene)) {
        #gene="ARHGAP18"
        gene <- KeyGene[GeneIndex]
        
        cat(paste0("-- gene of ", gene, " is analyzing...\n"))
        #mat <- expression_matrix
        # mat <- as.matrix(mat)
        # storage.mode(mat) <- "numeric"
        low <- mat[gene, ] <= median(mat[gene, ])
        high <- mat[gene, ] > median(mat[gene, ])
        lowRT <- mat[, low]
        highRT <- mat[, high]
        conNum <- ncol(lowRT)
        treatNum <- ncol(highRT)
        mat <- cbind(lowRT, highRT)
        
        Type <- c(rep("con", conNum), rep("treat", treatNum))
        design <- model.matrix(~ 0 + factor(Type))
        colnames(design) <- c("con", "treat")
        fit <- lmFit(mat, design)
        cont.matrix <- makeContrasts(treat - con, levels = design)
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2)
        
        allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
        
        #gsym.fc <- data.frame("SYMBOL" = rownames(allDiff), "logFC" = allDiff$logFC)
        gsym.fc <- data.frame("SYMBOL" = rownames(allDiff), "padj" = allDiff$adj.P.Val,"logFC" = allDiff$logFC)
        
        if (!dir.exists(paste0(dir, "/", GeneIndex, ".", gene))) {
          dir.create(paste0(dir, "/", GeneIndex, ".", gene))
        }
        
        # write.table(data,file=paste0(dir,"/",gene,".input.txt"),sep="\t",quote=F,row.names = F)
        write.table(allDiff, file = paste0(dir, "/", GeneIndex, ".", gene, "/", gene, ".xls"), sep = "\t", quote = F)
        
        selectedGeneID <- gene
        mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9", "#7A142C", "#5D90BA", "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
        
        if(isKEGG){
          if (mouse) {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseKEGG(id.fc, organism = "mmu")
            kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          } else {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseKEGG(id.fc, organism = "hsa")
            kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
            
            # # 按 padj 从小到大排序（显著的基因在前）
            # # ✓ 关键：按 -log10(padj) 从大到小排序（padj 越小越靠前）
            # gsym.fc.id$neg_log_padj <- -log10(gsym.fc.id$padj + 1e-300)  # 避免 log(0)
            # gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$neg_log_padj, decreasing = T), ]
            # 
            # # 按 padj 显著性排序
            # id.fc <- gsym.fc.id.sorted$neg_log_padj
            # names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            # 
            # # 检查是否递减排序
            # #head(id.fc)
            # 
            # kk <- gseKEGG(id.fc, organism = "hsa")
            # kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
            # sortkk <- kk.gsym[order(kk.gsym$p.adjust, decreasing = F), ]
            # enrichKEGG 方法（基于超几何分布）
            # 不需要排序，直接用显著基因进行富集
            # sig.genes <- gsym.fc[gsym.fc$padj < 0.05, ]  # 筛选显著基因
            # gsym.id <- bitr(sig.genes$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            # gene.list <- gsym.id$ENTREZID
            # kk <- enrichKEGG(gene = gene.list, organism = "hsa", pvalueCutoff = 0.05)
            # kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
            # sortkk <- kk.gsym[order(kk.gsym$p.adjust), ]
            
          }
          write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/gsea_output.csv"), quote = F, row.names = F)
          if(nrow(sortkk) < 3){
            message("KEGG少于三条通路！")
            Sys.sleep(3)
            next
          }
        }else{
          message("正在做GO！")
          if (mouse) {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseGO(id.fc, OrgDb = "org.Mm.eg.db")
            kk.gsym <- setReadable(kk, "org.Mm.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          } else {
            gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
            gsym.fc.id <- merge(gsym.fc, gsym.id, by = "SYMBOL", all = F)
            gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
            id.fc <- gsym.fc.id.sorted$logFC
            names(id.fc) <- gsym.fc.id.sorted$ENTREZID
            kk <- gseGO(id.fc, OrgDb = "org.Hs.eg.db")
            kk.gsym <- setReadable(kk, "org.Hs.eg.db", "ENTREZID")
            sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T), ]
          }
          write.csv(sortkk, paste0(dir, "/", GeneIndex, ".", gene, "/", "/gsea_output.csv"), quote = F, row.names = F)
          if(nrow(sortkk) < 3){
            message("GO也没有，没救了，做软件去！")
            Sys.sleep(3)
            next
          }
        }
        if(all.pathway){
          Annex <- paste0(dir, "/", GeneIndex, ".", gene, "/1.Annex")
          if (!dir.exists(Annex)) {
            dir.create(Annex)
          } else {
            print("Dir already exists!")
          }
          p1 <- gsea_dotplot(kk,color = 'p.adjust',x='NES',size = 'GeneRatio',title = paste0('GSEA of ',gene))
          term <- data.frame(term = levels(result$Description)[(length(result$Description)-2):length(result$Description)])
          write.table(term, paste0(dir, "./", GeneIndex, ".", gene, "/", "very_easy_input.txt"), quote = F, row.names = F, sep = "\t")
          ggsave_fun(paste0(dir, "/", GeneIndex, ".", gene, "/2.Dotplot"),plot = p1,width = 12,height = 12)
          
          for (i in 1:nrow(sortkk)) {
            id <- sortkk$ID
            names <- sortkk$Description
            p2 <- gseaplot3(kk,geneSetID = id[i],pvalue_table = T,title = names[i])
            ggsave_fun(paste0(dir, "/", GeneIndex, ".", gene, "/1.Annex/", i,".",names[i]),plot = p2,width = 8,height = 6)
          }
        }else{        
          cat("\033[0;34m请输入通路ID，英文逗号分隔\n")
          geneSetID <- scan(what = "character", sep = ",", quiet = TRUE)
          description.grep <- sortkk[sortkk$ID %in% geneSetID, ]$Description
          
          x <- kk
          geneList <- position <- NULL
          
          gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
          gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL, 3)
          
          p.res <- ggplot(gsdata, aes_(x = ~x)) +
            xlab(NULL) +
            geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
            scale_color_manual(values = mycol) +
            geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
            ylab("Enrichment\n Score") +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            theme(
              legend.position = "top", legend.title = element_blank(),
              legend.background = element_rect(fill = "transparent")
            ) +
            theme(
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
            )
          
          
          rel_heights <- c(1.5, .5, 1.5)
          
          p2 <- ggplot(gsdata, aes_(x = ~x)) +
            geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) +
            xlab(NULL) +
            ylab(NULL) +
            scale_color_manual(values = mycol) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            theme(
              legend.position = "none",
              plot.margin = margin(t = -.1, b = 0, unit = "cm"),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line.x = element_blank()
            ) +
            scale_y_continuous(expand = c(0, 0))
          
          
          df2 <- p.res$data
          df2$y <- p.res$data$geneList[df2$x]
          df2$gsym <- p.res$data$gsym[df2$x]
          
          selectgenes <- data.frame(gsym = selectedGeneID)
          selectgenes <- merge(selectgenes, df2, by = "gsym")
          selectgenes <- selectgenes[selectgenes$position == 1, ]
          
          p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) +
            geom_segment(
              data = df2, aes_(x = ~x, xend = ~x, y = ~y, yend = 0),
              color = "grey"
            ) +
            geom_bar(position = "dodge", stat = "identity") +
            scale_fill_manual(values = mycol, guide = FALSE) +
            scale_color_manual(values = mycol, guide = FALSE) +
            geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
            ylab("Ranked list\n metric") +
            xlab("Rank in ordered dataset") +
            theme_bw() +
            theme(
              axis.text.y = element_text(size = 12, face = "bold"),
              panel.grid = element_blank()
            ) +
            geom_text_repel(
              data = selectgenes,
              show.legend = FALSE,
              direction = "x",
              ylim = c(2, NA),
              angle = 90,
              size = 2.5, box.padding = unit(0.35, "lines"),
              point.padding = unit(0.3, "lines")
            ) +
            theme(plot.margin = margin(t = -.1, r = .4, b = .2, l = .2, unit = "cm"))
          
          
          plotlist <- list(p.res, p2, p.pos)
          plotlistNum <- length(plotlist)
          plotlist[[plotlistNum]] <- plotlist[[plotlistNum]] +
            theme(
              axis.line.x = element_line(),
              axis.ticks.x = element_line(),
              axis.text.x = element_text(size = 12, face = "bold")
            )
          
          p <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
          
          ggsave_fun(paste0(dir, "/", GeneIndex, ".", gene, "/", "/1.GSEA_multi_pathways"), width = 8, height = 6, plot = p)
          
          sortkk <- kk.gsym[kk.gsym@result$Description %in% description.grep[1] |
                              kk.gsym@result$Description %in% description.grep[2] |
                              kk.gsym@result$Description %in% description.grep[3], ]
          
          
          go <- data.frame(
            Category = "KEGG",
            ID = sortkk$ID,
            Term = sortkk$Description,
            Genes = gsub("/", ", ", sortkk$core_enrichment),
            adj_pval = sortkk$p.adjust
          )
          
          genelist <- data.frame(ID = gsym.fc.id$SYMBOL, logFC = gsym.fc.id$logFC)
          
          df <- GOplot::circle_dat(go, genelist)[, c(3, 5, 6)]
          write.table(df, paste0(dir, "./", GeneIndex, ".", gene, "/", "very_easy_input.txt"), quote = F, row.names = F, sep = "\t")
          
          nodes <- gather_graph_node(df, index = c("term", "genes"), value = "logFC", root = "all")
          edges <- gather_graph_edge(df, index = c("term", "genes"), root = "all")
          nodes <- nodes %>% mutate_at(c("node.level", "node.branch"), as.character)
          
          graph <- tidygraph::tbl_graph(nodes, edges)
          width <- 600/(nrow(nodes) + 100)
          
          gc1 <- ggraph(graph, layout = "dendrogram", circular = TRUE) +
            geom_edge_diagonal(
              aes(
                color = node1.node.branch,
                filter = node1.node.level != "all"
              ),
              alpha = 0.5,
              edge_width = 2.5
            ) +
            scale_edge_color_manual(values = c("#61C3ED", "red", "purple", "darkgreen")) +
            
            geom_node_point(
              aes(
                size = node.size,
                filter = node.level != "all"
              ),
              color = "#61C3ED"
            ) +
            scale_size(range = c(width, 10)) +
            theme(legend.position = "none") +
            
            # 优化后：叶子节点（基因标签）- 宽间距+小字体+优化角度
            geom_node_text(
              aes(
                x = 1.05 * x,
                y = 1.05 * y,
                label = node.short_name,
                angle = node_angle(x, y),
                filter = leaf
              ),
              color = "black",
              size = 2.5,
              hjust = "outward",
              vjust = 0.5
            ) +
            
            # 非叶子节点（通路标签）- 保持原设置即可
            geom_node_text(
              aes(
                label = node.short_name,
                filter = !leaf & (node.level != "all")
              ),
              color = "black",
              fontface = "bold",
              size = 6,
              family = "sans"
            ) +
            theme(panel.background = element_rect(fill = NA)) +
            # 优化后：扩大画布范围，适配宽间距
            coord_cartesian(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4))
          
          # 优化后：调大图形尺寸，保存高清疏解圈图
          ggsave_fun(paste0(dir, "./", GeneIndex, ".", gene, "/2.Ccgraph"), 
                     width = 18, height = 18, plot = gc1)
        }
      }
    }
    
    # GSVA
    if (step == 2) {
      message("*************************GSVA*************************")
      
      dir <- paste0(Path, "/.GSVA")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      if (mouse) {
        load(file.path(output_dir,"Mus/Mm.symbols.RData"))
        gsva_es <- gsva(as.matrix(mat), gs)
        cutoff <- 1
      } else {
        load(file.path(output_dir,"Homo/hallmark.gs.RData"))
        gsva_es <- gsva(gsvaParam(as.matrix(mat), gs))
        cutoff <- 1
      }
      
      for (i in KeyGene) {
        message(paste0("analysis of ", i, " starts..."))
        subexpr <- as.numeric(mat[i, ])
        names(subexpr) <- colnames(mat)
        
        lsam <- names(subexpr[subexpr < median(subexpr)])
        hsam <- names(subexpr[subexpr >= median(subexpr)])
        
        # GSVA
        group_list <- data.frame(sample = c(lsam, hsam), group = c(rep("Lexp", length(lsam)), rep("Hexp", length(hsam))))
        design <- model.matrix(~ 0 + factor(group_list$group))
        colnames(design) <- levels(factor(group_list$group))
        rownames(design) <- colnames(gsva_es)
        contrast.matrix <- makeContrasts(Hexp - Lexp, levels = design)
        
        fit <- lmFit(gsva_es[, group_list$sample], design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
        
        pathway <- stringr::str_replace(row.names(x), "HALLMARK_", "")
        df <- data.frame(ID = pathway, score = x$t)
        df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf), labels = c(1, 2, 3))
        sortdf <- df[order(df$score), ]
        sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
        write.table(sortdf, paste0(dir, "/", which(KeyGene == i),".",i, " GSVA.result.txt"), row.names = F, sep = "\t", quote = F)
        
        ggplot(sortdf, aes(ID, score, fill = group)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          scale_fill_manual(values = c("palegreen3", "snow3", "dodgerblue4"), guide = FALSE) +
          geom_hline(
            yintercept = c(-cutoff, cutoff),
            color = "white",
            linetype = 2,
            size = 0.3
          ) +
          geom_text(
            data = subset(df, score < 0),
            aes(x = ID, y = 0, label = paste0(" ", ID), color = group),
            size = 3,
            hjust = "inward"
          ) +
          geom_text(
            data = subset(df, score > 0),
            aes(x = ID, y = -0.05, label = ID, color = group),
            size = 3, hjust = "outward"
          ) +
          scale_colour_manual(values = c("black", "snow3", "black"), guide = FALSE) +
          xlab("") +
          ylab(paste0("t value of GSVA score\n Hexp vs Lexp group of ", i)) +
          theme_bw() +
          theme(panel.grid = element_blank()) +
          theme(panel.border = element_rect(size = 0.6)) +
          theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
        ggsave_fun(paste0(dir, "/", which(KeyGene == i), ".GSVA plot of ", i), width = 8, height = 7)
      }
    }
    
    # 免疫浸润
    if (step == 3) {
      message("*************************免疫浸润*************************")
      
      cat("1.CIBERSORT\n")
      cat("2.ssgsea\n")
      cat("\033[0;34m方法选择：")
      method <- scan(what = "numeric", quiet = TRUE)
      
      dir <- paste0(Path, "/.Immune infiltration")
      
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      # mouse=F
      # method <- 1
      if (mouse == T && !file.exists(paste0(dir, "/CIBERSORT-Results.txt"))) {
        results <- CIBERSORT(dir = dir, file.path(output_dir,"refdata/Mus/Immune of mouse.txt"), mat, perm = 100, QN = TRUE)
      }
      if (mouse == F && method == 1 && !file.exists(paste0(dir, "/CIBERSORT-Results.txt"))) {
        results <- CIBERSORT(dir = dir, file.path(output_dir,"refdata/Homo/ref.txt"), mat, perm = 100, QN = TRUE)
      }
      if (mouse == F && method == 2) {
        gmtFile <- file.path(output_dir,"refdata/Homo/immune.ssgsea.gmt")
        Exper <- mat[rowMeans(mat) > 0, ]
        geneSet <- GSEABase::getGmt(gmtFile,
                                    geneIdType = GSEABase::SymbolIdentifier()
        )
        
        ssgseaScore <- GSVA::gsva(Exper, geneSet, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
        normalize <- function(x) {
          return((x - min(x)) / (max(x) - min(x)))
        }
        ssgseaOut <- normalize(ssgseaScore)
        ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut) %>%
          t() %>%
          data.frame(check.names = F)
        write.table(ssgseaOut, file = paste0(dir, "/ssgseaOut.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
      }
      
      cat("\n========== 诊断信息 ==========\n")
      cat("当前目录：", dir, "\n")
      cat("目录存在？", dir.exists(dir), "\n")
      
      # 列出目录中的所有文件
      files_in_dir <- list.files(dir)
      cat("目录中的文件：\n")
      print(files_in_dir)
      
      cibersort_file <- paste0(dir, "/CIBERSORT-Results.txt")
      cat("\nCIBERSORT文件路径：", cibersort_file, "\n")
      cat("文件存在？", file.exists(cibersort_file), "\n")
      
      if (file.exists(cibersort_file)) {
        # 检查文件大小
        file_size <- file.info(cibersort_file)$size
        cat("文件大小：", file_size, "bytes\n")
        
        # 读取前几行
        cat("\n文件前10行内容：\n")
        print(readLines(cibersort_file, n = 10))
        
        # 尝试读取
        tryCatch({
          rt <- read.table(cibersort_file, sep = "\t", header = T, row.names = 1, check.names = F)
          cat("\n✓ 文件读取成功！\n")
          cat("行数：", nrow(rt), "\n")
          cat("列数：", ncol(rt), "\n")
          cat("列名：", colnames(rt), "\n")
          cat("行名（前10个）：\n")
          print(head(rownames(rt), 10))
        }, error = function(e) {
          cat("\n✗ 读取失败，错误信息：\n")
          print(e$message)
        })
      } else {
        cat("\n✗ CIBERSORT-Results.txt 文件不存在！\n")
        cat("请检查CIBERSORT是否运行成功。\n")
      }
      ### =====================BarPlot========================
      if (method == 1) {
        data <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), header = T, sep = "\t", check.names = F, row.names = 1)
        data <- dplyr::select(data,-c("P-value","Correlation","RMSE"))
        data <- t(data)
      }
      if (method == 2) {
        data <- read.table(paste0(dir, "/ssgseaOut.txt"), header = T, sep = "\t", check.names = F, row.names = 1) %>% as.matrix()
        data <- t(data)
      }
      GseGroup <- read.table(DGOROUP, header = F, sep = "\t", check.names = F)
      colnames(GseGroup) <- c("Acc", "Cluster", "Tissue")
      
      print("分组统计：")
      print(table(GseGroup$Tissue))
      
      # 检查样本名称
      print("GseGroup中的样本数：")
      print(nrow(GseGroup))
      
      print("rt中的样本数：")
      print(nrow(rt))
      
      # 检查样本名称匹配情况
      matched <- intersect(rownames(rt), GseGroup$Acc)
      print(paste("匹配的样本数：", length(matched)))
      
      compare_method <- "anova"
      # ciberRes <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), header = T, sep = "\t", check.names = F, row.names = 1)
      ciber <- data %>% t() %>% data.frame( check.names = F)
      ciber$Acc <- rownames(ciber)
      if (sort(unique(GseGroup$Tissue))[1] == "Control") {
        datGroup <- GseGroup[order(GseGroup$Tissue), ]
      } else {
        datGroup <- GseGroup[order(GseGroup$Tissue, decreasing = T), ]
      }
      datGroup$Acc <- factor(datGroup$Acc, levels = datGroup$Acc)
      data <- dplyr::inner_join(datGroup, ciber, by = "Acc")
      data_p <- reshape2::melt(data, id.vars = colnames(datGroup))
      data_p$Acc <- factor(data_p$Acc, levels = datGroup$Acc)
      
      datGroup$p <- "Group"
      p2 <- ggplot(datGroup, aes(Acc, p, fill = Tissue)) +
        geom_tile() +
        scale_fill_manual(values = c("#1CFA04", "#C705FF")) +
        scale_y_discrete(position = "right") +
        theme_minimal() +
        xlab(NULL) +
        ylab(NULL) +
        theme(text = element_text(size = 15)) +
        theme(axis.text.x = element_blank()) +
        labs(fill = "Group")
      
      
      sample_names <- datGroup$Acc # replace with your actual sample names
      
      p1 <- ggplot(data_p, aes(x = Acc, y = value, fill = variable)) +
        geom_bar(stat = "identity", position = "fill", width = 0.5) +
        geom_col(position = "fill", width = 0.6) +
        guides(fill = guide_legend(title = NULL)) +
        ylab("Relative Percent") +
        xlab("") +
        theme_bw() +
        theme(axis.ticks.length = unit(0.5, "cm")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text = element_text(size = 15)) +
        scale_x_discrete(labels = sample_names) +
        theme(axis.text.y = element_text(colour = "black", vjust = 0, size = 5)) +
        theme(axis.title = element_text(size = 20)) +
        theme(text = element_text(size = 15)) +
        scale_y_continuous(expand = c(0, 0.05)) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      
      p <- p1 %>% aplot::insert_top(p2, height = 0.05)
      ggsave(filename = paste0(dir, "/1.Immune infiltration.pdf"), plot = p, width = 15, height = 6, limitsize = FALSE)
      # ggsave(filename = paste0(dir, "/1.Immune infiltration.png"), plot = p, dpi = 600, width = 15, height = 6, limitsize = FALSE)
      
      ### =====================Heatmap========================
      if (method == 1) {
        rt <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
        rt <- dplyr::select(rt,-c("P-value","Correlation","RMSE"))
        # rt <- rt[, -c(23:25)]
      }
      if (method == 2) {
        rt <- read.table(paste0(dir, "/ssgseaOut.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
      }
      if (length(which(colSums(rt) == 0)) > 0) {
        rt <- rt[, -which(colSums(rt) == 0)]
      }
      
      pdf(paste0(dir, "/2.CorHeatmap.pdf"), height = 13, width = 13)
      corrplot(
        corr = cor(rt),
        method = "color",
        order = "hclust",
        tl.col = "black",
        addCoef.col = "black",
        number.cex = 1,
        col = colorRampPalette(c("blue", "white", "red"))(50)
      )
      dev.off()
      
      # png(paste0(dir, "/2.CorHeatmap.png"), height = 4000, width = 4000, res = 300)
      # corrplot(
      #   corr = cor(rt),
      #   method = "color",
      #   order = "hclust",
      #   tl.col = "black",
      #   addCoef.col = "black",
      #   number.cex = 1,
      #   col = colorRampPalette(c("blue", "white", "red"))(50),
      # )
      # dev.off()
      
      ### =====================BoxPlot========================
      pFilter <- 0.99
      
      # Path <- "<PATH>/ARHGAP18-AML/GSE18322"
      # dir <- paste0(Path, "/.Immune infiltration")
      # method == 1
      if (method == 1) {
        rt <- read.table(paste0(dir, "/CIBERSORT-Results.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
        rt <- dplyr::select(rt,-c("P-value","Correlation","RMSE"))
        # rt <- rt[, -c(23:25)]
      }
      if (method == 2) {
        rt <- read.table(paste0(dir, "/ssgseaOut.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
      }
      if (length(which(colSums(rt) == 0)) > 0) {
        rt <- rt[, -which(colSums(rt) == 0)]
      }
      data <- rt
      
      Type <- read.table(DGOROUP, sep = "\t", check.names = F, row.names = 1, header = F)
      Type <- Type[row.names(data), ]
      colnames(Type) <- c("cluster", "Subtype")
      
      outTab <- data.frame()
      data <- cbind(data, Type)
      
      for (i in colnames(data[, 1:(ncol(data) - 2)])) {
        rt1 <- data[, c(i, "Subtype")]
        colnames(rt1) <- c("expression", "Subtype")
        
        # Skip columns with zero variance
        if (var(rt1$expression, na.rm = TRUE) == 0) {
          next
        }
        
        ksTest <- kruskal.test(expression ~ Subtype, data = rt1)
        pValue <- ksTest$p.value
        
        # Handle NA p-values
        if (!is.na(pValue) && pValue < pFilter) {
          outTab <- rbind(outTab, cbind(rt1, gene = i))
        }
      }
      # write.table(outTab, file = "data.txt", sep = "\t", row.names = F, quote = F)
      # data <- read.table("data.txt", sep = "\t", header = T, check.names = F)
      outTab$Subtype <- factor(outTab$Subtype, levels = c("Control", "Disease"))
      outTab$gene <- as.factor(outTab$gene)
      
      outTab <- filter(outTab, Subtype %in% c("Control", "Disease"))
      p <- ggplot(outTab, aes(x = gene, y = expression, fill = Subtype, color = Subtype)) +
        geom_boxplot(alpha = 0.3) +
        scale_fill_manual(name = "Subtype", values = c("deepskyblue", "hotpink")) +
        scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3")) +
        theme_bw() +
        labs(x = "", y = "Expression") +
        theme(axis.text.x = element_text(vjust = 1, size = 12, hjust = 1, colour = "black"), legend.position = "top") +
        rotate_x_text(45) +
        stat_compare_means(aes(group = Subtype), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), label = "p.signif", method = "wilcox", show.legend = FALSE)
      ggsave_fun(paste0(dir, "/3.Boxplot"), width = 12, height = 6, plot = p)
      
      plot_build <- ggplot_build(p)
      compare_means_data <- plot_build$data[[2]] %>%
        arrange(x) %>%
        dplyr::select(x, p.format)
      compare_means_data$x <- levels(outTab$gene)
      
      summary_stats <- outTab %>%
        group_by(gene, Subtype) %>%
        summarise(
          mean = mean(expression),
          q1 = quantile(expression, 0.25),
          q3 = quantile(expression, 0.75)
        )
      
      summary_stats <- merge(compare_means_data, summary_stats, by.x = "x", by.y = "gene") %>% subset(subset = p.format <= 0.05)
      write.table(summary_stats, paste0(dir, "/BoxPlotdata.txt"), row.names = F, quote = F, sep = "\t")
      
      ### =====================HorizonHeatmap========================
      
      immuscore <- function(KeyGene) {
        y <- as.numeric(mat[KeyGene, ])
        colnames <- colnames(rt)
        do.call(rbind, lapply(colnames, function(x) {
          dd <- cor.test(as.numeric(rt[, x]), y, method = "pearson")
          data.frame(gene = KeyGene, immune_cells = x, cor = dd$estimate, p.value = dd$p.value)
        }))
      }
      print(1)
      data <- do.call(rbind, lapply(KeyGene, immuscore))
      write.csv(data, paste0(dir, "/HorizonCorr.csv"), quote = F, row.names = F)
      
      data$pstar <- ifelse(data$p.value < 0.05,
                           ifelse(data$p.value < 0.01, "**", "*"),
                           ""
      )
      
      ggplot(data, aes(immune_cells, gene)) +
        geom_tile(aes(fill = cor), colour = "black", size = 1) +
        scale_fill_gradient2(low = "#2b8cbe", mid = "white", high = "#e41a1c") +
        geom_text(aes(label = pstar), col = "black", size = 5) +
        theme_minimal() + # 不要背景
        theme(
          axis.title.x = element_blank(), # 不要title
          axis.ticks.x = element_blank(), # 不要x轴
          axis.title.y = element_blank(), # 不要y轴
          axis.text.x = element_text(angle = 45, hjust = 1), # 调整x轴文字
          axis.text.y = element_text(size = 8)
        ) + # 调整y轴文字
        # 调整legen
        labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "Correlation"))
      ggsave_fun(paste0(dir, "/4.Correlation"), width = 8, height = 4)
    }
    
    # 免疫调控因子
    if (step == 4) {
      message("*************************免疫调控*************************")
      
      dir <- paste0(Path, "/.Immune association")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      a_1 <- mat %>% data.frame()
      a_2 <- as.data.frame(t(a_1))
      a_3 <- a_1
      a_3$Id <- rownames(a_3)
      
      if (mouse) {
        b_1 <- read.table("refdata/Mus/immune.txt", header = T, sep = "\t", quote = "", fill = T)
      } else {
        b_1 <- read.table("refdata/Homo/Immunomodulator_and_chemokines.txt", header = T, sep = "\t", quote = "", fill = T)
      }
      
      plotlist <- list()
      for (i in 1:length(unique(b_1$type))) {
        b_2 <- b_1[b_1$type == unique(b_1$type)[i], ]
        data1 <- dplyr::inner_join(b_2, a_3, by = "Id")
        data2 <- a_2[, c(KeyGene, data1$Id)]
        
        
        res <- rcorr(as.matrix(data2), type = "spearman")
        result_1 <- CorMatrix(res$r, res$P)
        
        result_2 <- result_1[result_1$row %in% KeyGene, ]
        b_2$column <- b_2$Id
        result_3 <- dplyr::inner_join(result_2, b_2, by = "column")
        result1 <- result_3[, 1:4]
        result1$Regulation <- result1$cor
        result1[, 5][result1[, 5] > 0] <- c("postive")
        result1[, 5][result1[, 5] < 0] <- c("negative")
        colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")
        # write.table(result1,file=paste0(dir,"/",i,".",unique(b_1$type)[1],".xls"),sep="\t",quote=F,col.names=T,row.names = F)
        
        # a1 <- read.table(paste0(dir,"/",i,".",unique(b_1$type)[1],".xls"),header = T,sep = "\t", quote = "",fill = T)
        data2 <- result1
        data2$pvalue <- ifelse(data2$pvalue < 0.05,
                               ifelse(data2$pvalue < 0.01, "**", "*"),
                               ""
        )
        data2$pvalue[1:20]
        data2$type <- data2$cor
        data3 <- data2[order(data2$immuneGene, data2$cor), ]
        data4 <- data3[data3$pvalue < 0.05, ]
        
        if (nrow(data4) > 40) {
          dotheight <- nrow(data4) * 0.18
        }
        if (nrow(data4) < 20) {
          dotheight <- 10
        } else {
          dotheight <- 15
        }
        
        data4 <- na.omit(data4)
        p <- ggplot(data4, aes(x = gene, y = immuneGene)) +
          geom_point(aes(colour = cor, size = pvalue)) +
          labs(x = "", y = unique(b_1$type)[i])
        p <- p + scale_colour_gradient2(
          low = "blue", high = "red", mid = "white",
          midpoint = 0, limit = c(-1, 1), space = "Lab",
          name = "Pearson\nCorrelation"
        )
        p <- p + theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(axis.text = element_text(size = 15)) +
          theme(axis.text.x = element_text(colour = "black", angle = 90, hjust = 0.5, size = 15)) +
          theme(axis.text.y = element_text(colour = "black", vjust = 0.5, hjust = 1, size = 15)) +
          theme(axis.title = element_text(size = 20)) +
          theme(text = element_text(size = 15))
        
        ggsave_fun(paste0(dir, "/", i, ".", capitalize_first(unique(b_1$type)[i])), width = 7, height = dotheight, plot = p)
      }
    }
    
    # 差异表达
    if (step == 5) {
      message("*************************差异表达*************************")
      
      dir <- paste0(Path, "/.Differential Expression")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      pFilter <- 0.99
      
      data <- subset(mat, subset = rownames(mat) %in% KeyGene) %>%
        t() %>%
        data.frame()
      
      Type <- read.table(DGOROUP, sep = "\t", check.names = F, row.names = 1, header = F)
      Type <- Type[row.names(data), ]
      colnames(Type) <- c("cluster", "Subtype")
      
      outTab <- data.frame()
      data <- cbind(data, Type)
      
      
      for (i in colnames(data[, 1:(ncol(data) - 2)])) {
        rt1 <- data[, c(i, "Subtype")]
        colnames(rt1) <- c("expression", "Subtype")
        ksTest <- kruskal.test(expression ~ Subtype, data = rt1)
        pValue <- ksTest$p.value
        if (pValue < pFilter) {
          outTab <- rbind(outTab, cbind(rt1, gene = i))
          print(pValue)
        }
      }
      # write.table(outTab,file="data.txt",sep="\t",row.names=F,quote=F)
      
      outTab$Subtype <- factor(outTab$Subtype, levels = c("Control", "Disease"))
      outTab$gene <- as.factor(outTab$gene)
      
      outTab <- filter(outTab, Subtype %in% c("Control", "Disease"))
      
      p <- ggplot(outTab, aes(x = gene, y = expression, fill = Subtype, color = Subtype)) +
        geom_boxplot(alpha = 0.3) +
        scale_fill_manual(name = "Subtype", values = c("deepskyblue", "hotpink")) +
        scale_color_manual(name = "Subtype", values = c("dodgerblue", "plum3")) +
        theme_bw() +
        labs(x = "", y = "Expression") +
        theme(axis.text.x = element_text(vjust = 1, size = 12, hjust = 1, colour = "black"), legend.position = "top") +
        rotate_x_text(45) +
        stat_compare_means(aes(group = Subtype), symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), label = "p.signif", method = "wilcox", show.legend = FALSE)
      
      ggsave_fun(paste0(dir, "/1.Boxplot"), width = 7, height = 5, plot = p)
      
      plot_build <- ggplot_build(p)
      compare_means_data <- plot_build$data[[2]] %>%
        arrange(x) %>%
        dplyr::select(x, p.format)
      compare_means_data$x <- levels(outTab$gene)
      
      summary_stats <- outTab %>%
        group_by(gene, Subtype) %>%
        summarise(
          mean = mean(expression),
          q1 = quantile(expression, 0.25),
          q3 = quantile(expression, 0.75)
        )
      
      summary_stats <- merge(compare_means_data, summary_stats, by.x = "x", by.y = "gene") %>% subset(subset = p.format <= 0.05)
      write.table(summary_stats, paste0(dir, "/BoxPlotdata.txt"), row.names = F, quote = F, sep = "\t")
    }
    
    # 疾病相关性
    if (step == 6) {
      message("*************************疾病相关性*************************")
      
      dir <- paste0(Path, "/.Disease correlation")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      tcga.exp <- data.frame(mat,check.names = F)
      KeyGenes <- KeyGene
      Expr <- tcga.exp[rowSums(tcga.exp == 0) < (ncol(tcga.exp) * 0.5), ] %>% data.frame(check.names = F)
      # Expr <- log2(Expr + 1)
      Group <- read.table(DGOROUP, sep = "\t", check.names = F, header = F)
      colnames(Group) <- c("Sample", "Cluster", "Tissue")
      
      cat("\033[0;34m请输入疾病相关基因，回车分隔\n")
      geneCard <- scan(what = "character", sep = "\n", quiet = FALSE)
      geneCard <- intersect(geneCard,rownames(Expr))[1:20]
      geneCard <- geneCard[!is.na(geneCard)]
      Expr <- as.data.frame(t(Expr[unique(c(KeyGenes, geneCard)), colnames(Expr) %in% Group$Sample]))
      # Expr <- log2(Expr + 1)
      Expr$Sample <- rownames(Expr)
      DatGroup <- dplyr::inner_join(Group, Expr, by = "Sample")
      data_p <- reshape2::melt(DatGroup, id.vars = colnames(Group))
      data_p$Tissue <- factor(data_p$Tissue, levels = c("Control", "Disease"))
      plot_G <- list()
      for (i in unique(data_p$variable)) {
        tmp <- data_p[data_p$variable == i, ]
        my_comparisons <- list(c("Control", "Disease"))
        p <- ggplot(
          tmp,
          aes(
            x = Tissue, y = value,
            fill = Tissue,
            color = Tissue
          )
        ) +
          ylim((min(tmp$value) - min(tmp$value) / 30), (max(tmp$value) + max(tmp$value) / 28)) +
          geom_boxplot(
            fill = "white", notch = F, alpha = 0.95,
            outlier.shape = 16,
            outlier.size = 0.65
          ) +
          geom_point(position = position_jitterdodge(), size = 0.5, alpha = 0.3) +
          xlab("") +
          ylab(paste0(i)) +
          theme_bw() +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
          theme(axis.text = element_text(size = 10)) +
          theme(axis.text.x = element_text(colour = "black", angle = 0, vjust = 0, size = 10)) +
          theme(axis.text.y = element_text(colour = "black", vjust = 0, size = 10)) +
          theme(axis.title = element_text(size = 15)) +
          theme(text = element_text(size = 10)) +
          theme(legend.position = "right", legend.title = element_blank()) +
          stat_compare_means(method = "wilcox", vjust = 0.5, size = 4,show.legend = F)
        plot_G[[i]] <- p
      }
      
      data_p <- data_p[data_p$variable %in% geneCard, ]
      data_p$variable <- factor(data_p$variable,unique(data_p$variable))
      p <- ggplot(
        data_p,
        aes(
          x = variable, y = value,
          fill = Tissue,
          color = Tissue
        )
      ) +
        geom_boxplot(
          notch = F, alpha = 0.95,
          outlier.shape = 16,
          outlier.size = 0.65
        ) +
        xlab("") +
        ylab("Expression level") +
        scale_fill_manual(values = c("#D5EBFB", "#FBEEB7", "#B4FBCD", "#F5B3FC")) +
        scale_color_manual(values = c("#0073C2", "#EFC000", "#00C244", "#C501D7")) +
        ggtitle("") +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(angle = 90, size = 10),
          axis.title = element_text(angle = 90, size = 15)
        ) +
        theme(legend.position = "top") +
        stat_compare_means(method = "wilcox", hide.ns = T, label = "p.signif",show.legend = F)
      ggsave_fun(filename = paste0(dir, "/1.Expression level"), width = 15, height = 6, plot = p)
      
      plot_build <- ggplot_build(p)
      compare_means_data <- plot_build$data[[2]] %>%
        arrange(x) %>%
        dplyr::select(x, p.format)
      compare_means_data$x <- levels(data_p$variable)
      
      
      compare_means_data %>% mutate(p.value = as.numeric(gsub("<", "", p.format)),
                                    p.value = ifelse(is.na(p.value), 0.0001, p.value)) -> compare_means_data
      
      write.table(compare_means_data,paste0(dir,"/sig.plotdata.txt"),row.names = F,quote = F,sep = "\t")
      
      plot_G[["total"]] <- p
      rownames(DatGroup) <- DatGroup$Sample
      DatGroup <- DatGroup[, setdiff(colnames(DatGroup), colnames(Group))]
      
      res <- rcorr(as.matrix(DatGroup), type = "pearson") # "pearson"
      res_1 <- reshape2::melt(res$r)
      colnames(res_1) <- c("row", "column", "cor")
      res_2 <- reshape2::melt(res$P)
      colnames(res_2) <- c("row", "column", "p")
      result_1 <- dplyr::inner_join(res_1, res_2, by = c("row", "column"))
      data <- result_1[result_1$row %in% KeyGenes & result_1$column %in% geneCard, ]
      data$pv <- ""
      if (length(which(data$p < 0.05)) > 0) {
        data[which(data$p < 0.05), ]$pv <- "*"
      }
      if (length(which(data$p < 0.01)) > 0) {
        data[which(data$p < 0.01), ]$pv <- "**"
      }
      if (length(which(data$p < 0.001)) > 0) {
        data[which(data$p < 0.001), ]$pv <- "***"
      }
      if (length(which(data$p < 0.0001)) > 0) {
        data[which(data$p < 0.0001), ]$pv <- "****"
      }
      
      p <- ggplot(data, aes(x = row, y = column)) +
        geom_point(aes(colour = cor, size = pv)) +
        labs(x = "", y = "") +
        scale_colour_gradient2(
          low = "blue", high = "red", mid = "white",
          midpoint = 0, space = "Lab", # , limit = c(-0.6, 0.6)
          name = "Pearson\nCorrelation"
        ) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text = element_text(size = 10)) +
        theme(axis.text.x = element_text(colour = "black", angle = 45, hjust = 1, size = 10)) + #
        theme(axis.text.y = element_text(colour = "black", vjust = 0.5, size = 10)) +
        theme(axis.title = element_text(size = 15)) +
        theme(text = element_text(size = 10))
      plot_G[["cor"]] <- p
      tmp <- DatGroup
      plot_C <- list()
      for (i in KeyGenes) {
        for (j in geneCard) {
          test <- cor.test(tmp[, i], tmp[, j], method = "pearson")
          test$estimate
          test$p.value
          data <- data.frame(
            moduleScore = tmp[, j],
            PCAScore = tmp[, i]
          )
          p <- ggplot(data, aes(moduleScore, PCAScore)) +
            xlim(
              ((min(tmp[, j]) * 1.05) - (max(tmp[, j]) * 0.05)),
              ((max(tmp[, j]) * 1.05) - (min(tmp[, j]) * 0.05))
            ) +
            ylim(
              ((min(tmp[, i]) * 1.05) - (max(tmp[, i]) * 0.05)),
              ((max(tmp[, i]) * 1.05) - (min(tmp[, i]) * 0.05))
            ) +
            xlab(j) +
            ylab(i) +
            geom_point() +
            geom_smooth() +
            theme_bw() +
            theme(axis.text = element_text(size = 10)) +
            theme(axis.title = element_text(size = 15)) +
            theme(text = element_text(size = 10)) +
            annotate("text",
                     x = (max(tmp[, j]) + min(tmp[, j])) * 0.5, y = ((max(tmp[, i]) * 0.95) + (min(tmp[, i]) * 0.05)),
                     label = paste0("cor= ", round(test$estimate, 3), "\n", "pvalue= ", format(test$p.value, digits = 2))
            )
          plot_C[[paste0(i, " ~ ", j)]] <- p
        }
      }
      
      data <- result_1[result_1$row %in% KeyGenes & result_1$column %in% geneCard, ]
      tmp <- data[which(data$p < 0.05), ]
      tmp <- tmp[order(tmp$cor), ]
      
      GenSet1 <- tmp[1, ]
      GenSet2 <- tmp[nrow(tmp), ]
      if (GenSet1$row == KeyGenes[KeyGenes %in% c(GenSet1$row, GenSet2$row)][1]) {
        pX <- plot_G$total -
          (((plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_G[[GenSet1$row]] + plot_G[[GenSet1$column]] + plot_layout(ncol = 1, heights = c(1.3, 1, 1))) -
              plot_G$cor + plot_layout(ncol = 2, widths = c(1.3, 1.5))) -
             (plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_G[[GenSet2$row]] + plot_G[[GenSet2$column]] + plot_layout(ncol = 1, heights = c(1.2, 1, 1))) +
             plot_layout(ncol = 2, widths = c(3.2, 1))) +
          plot_layout(ncol = 1, heights = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/2.Disease gene Expression level"), plot = pX, width = 15, height = 12)
      } else {
        pX <- plot_G$total -
          (((plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_G[[GenSet2$row]] + plot_G[[GenSet2$column]] + plot_layout(ncol = 1, heights = c(1.2, 1, 1))) -
              plot_G$cor + plot_layout(ncol = 2, widths = c(1.3, 1.5))) -
             (plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_G[[GenSet1$row]] + plot_G[[GenSet1$column]] + plot_layout(ncol = 1, heights = c(1.3, 1, 1))) +
             plot_layout(ncol = 2, widths = c(3.2, 1))) +
          plot_layout(ncol = 1, heights = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/2.Disease gene Expression level"), plot = pX, width = 15, height = 12)
      }
      if (GenSet1$column == geneCard[geneCard %in% c(GenSet1$column, GenSet2$column)][2]) {
        pX2 <- (plot_spacer() + plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_spacer() +
                  plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_spacer() +
                  plot_layout(ncol = 1, heights = c(0.1, 1, 0.1, 1, 0.1))) -
          plot_G$cor + plot_layout(ncol = 2, widths = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/3.Disease gene Expression level"), plot = pX2, width = 12, height = 12)
      } else {
        pX2 <- (plot_spacer() + plot_C[[paste0(GenSet2$row, " ~ ", GenSet2$column)]] + plot_spacer() +
                  plot_C[[paste0(GenSet1$row, " ~ ", GenSet1$column)]] + plot_spacer() +
                  plot_layout(ncol = 1, heights = c(0.1, 1, 0.1, 1, 0.1))) -
          plot_G$cor + plot_layout(ncol = 2, widths = c(1, 1.5))
        ggsave_fun(filename = paste0(dir, "/3.Disease gene Expression level"), plot = pX2, width = 12, height = 12)
      }
      
      write.table(data, paste0(dir, "/report.txt"), quote = F, row.names = F, sep = "\t")
    }
    
    # Motif转录因子
    if (step == 7) {
      message("*************************Motif*************************")
      
      dir <- paste0(Path, "/.Motif")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      require(RcisTarget)
      if (mouse) {
        data(motifAnnotations_mgi)
        motifRankings <- importRankings("refdata/Mus/mm9-500bp-upstream-10species.mc9nr.feather")
        
        geneLists <- list(key_gene = KeyGene)
        
        data("motifAnnotations_hgnc_v9")
        
        motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                                 motifAnnot = motifAnnotations_hgnc_v9
        )
        
        motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes) %>% arrange(-NES)
        
        resultsSubset <- motifEnrichmentTable_wGenes_wLogo#[1:10,]
        
        df <- DT::datatable(resultsSubset[, -c("rankAtMax", "TF_lowConf"), with = FALSE],
                            escape = FALSE, # To show the logo
                            filter = "top", options = list(pageLength = 5)
        )
        
        htmlwidgets::saveWidget(df, paste0(dir, "/Motif.html"), selfcontained = TRUE)
        
        write.table(data.frame("ID" = row.names(motifEnrichmentTable_wGenes), motifEnrichmentTable_wGenes), file = paste0(dir, "/onco_matrix1.txt"), sep = "\t", quote = F, row.names = F)
      } else {
        require(RcisTarget.hg19.motifDBs.cisbpOnly.500bp)
        geneLists <- list(key_gene = KeyGene)
        data(motifAnnotations_hgnc_v9)
        data(hg19_500bpUpstream_motifRanking_cispbOnly)
        motifRankings <- hg19_500bpUpstream_motifRanking_cispbOnly
        
        motifEnrichmentTable_wGenes <- cisTarget(geneLists, motifRankings,
                                                 motifAnnot = motifAnnotations_hgnc_v9
        )
        
        motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes) %>% arrange(-NES)
        
        resultsSubset <- motifEnrichmentTable_wGenes_wLogo#[1:10,]
        
        # datatable(resultsSubset[,-c("enrichedGenes", "TF_lowConf"), with=FALSE],
        df <- DT::datatable(resultsSubset[, -c("rankAtMax", "TF_lowConf"), with = FALSE],
                            escape = FALSE, # To show the logo
                            filter = "top", options = list(pageLength = 5)
        )
        htmlwidgets::saveWidget(df, paste0(dir, "/Motif.html"), selfcontained = TRUE)
        
        write.table(data.frame("ID" = row.names(motifEnrichmentTable_wGenes), motifEnrichmentTable_wGenes), file = paste0(dir, "/onco_matrix1.txt"), sep = "\t", quote = F, row.names = F)
      }
      onco <- resultsSubset[, c("TF_highConf", "enrichedGenes"), with = FALSE]
      onco$TF_highConf <- gsub("\\s*\\(.*?\\.", "", onco$TF_highConf) # 去除括号及其内部的信息
      onco <- onco[onco$TF_highConf != "", ]
      result <- data.frame()
      for (i in 1:nrow(onco)) {
        tf <- strsplit(onco[i,]$TF_highConf, ";")[[1]]
        gene <- strsplit(onco[i,]$enrichedGenes, ";")[[1]]
        
        for (j in 1:length(gene)) {
          result <- rbind(result,data.frame(TF_highConf = tf,enrichedGenes = rep(gene[j], length(tf))))
        }
      }
      result <- result[duplicated(result) == F,]
      result$TF_highConf <- trimws(result$TF_highConf)
      result$enrichedGenes <- trimws(result$enrichedGenes)
      print(result)
      write.table(result, paste0(dir, "/Motif.txt"), quote = F, row.names = F, sep = "\t")
    }
    
    # 非编码 RNA 网络
    if (step == 8) {
      message("*************************非编码 RNA 网络*************************")
      
      dir <- paste0(Path, "/.Non-coding RNA networks")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      mircode <- data.table::fread("refdata/mircode_highconsfamilies.txt",header = T,data.table = F)
      key.mircode <- mircode[mircode$gene_symbol %in% KeyGene,]
      network.table <- data.frame(Gene = key.mircode$gene_symbol,microrna = key.mircode$microrna)
      network.table$microrna <- stringr::str_split(network.table$microrna,"/",simplify = T)[,1]
      network.table <- network.table[!duplicated(network.table),]
      write.table(network.table,paste0(dir,"/network.txt"),row.names = F,sep = "\t",quote = F)
    }
  }
}


# KO_scTenifoldKnk --------------------------------------------------------


# --- Function: strictDirection ---
strictDirection <- function(X, lambda = 1){
  S <- as.matrix(X)
  S[abs(S) < abs(t(S))] <- 0
  O <- (((1-lambda) * X) + (lambda * S))
  O <- Matrix::Matrix(O)
  return(O)
}


# --- Function: scTenifoldKnk_module1 ---
scTenifoldKnk_module1 <- function(countMatrix, rowMeans = 0, qc = TRUE, gKO = NULL, qc_mtThreshold = 0.1,
                                  qc_minLSize = 1000, nc_lambda = 0, nc_nNet = 10, nc_nCells = 500,
                                  nc_nComp = 3, nc_scaleScores = TRUE, nc_symmetric = FALSE,
                                  nc_q = 0.9, td_K = 3, td_maxIter = 1000, td_maxError = 1e-05,
                                  td_nDecimal = 3, ma_nDim = 2, nCores = parallel::detectCores()) {
  if (isTRUE(qc)) {
    countMatrix <- scQC(countMatrix, mtThreshold = qc_mtThreshold, minLSize = qc_minLSize)
  }
  # Skip internal rowMeans filtering to keep KO gene; assume external filtering already applied
  WT <- scTenifoldNet::makeNetworks(
    X = countMatrix, q = nc_q, nNet = nc_nNet, nCells = nc_nCells,
    scaleScores = nc_scaleScores, symmetric = nc_symmetric, nComp = nc_nComp,
    nCores = nCores
  )
  WT <- scTenifoldNet::tensorDecomposition(
    xList = WT, K = td_K, maxError = td_maxError,
    maxIter = td_maxIter, nDecimal = td_nDecimal
  )
  WT <- WT$X
  WT <- strictDirection(WT, lambda = nc_lambda)
  WT <- as.matrix(WT)
  diag(WT) <- 0
  WT <- t(WT)
  KO <- WT
  if (!gKO %in% rownames(KO)) {
    stop(paste0("gKO '", gKO, "' not found after preprocessing; aborting."))
  }
  KO[gKO, ] <- 0
  MA <- manifoldAlignment(WT, KO, d = ma_nDim, nCores = nCores)
  DR <- dRegulation(MA)
  outputList <- list()
  outputList$tensorNetworks <- list()
  outputList$tensorNetworks$WT <- Matrix(WT)
  outputList$tensorNetworks$KO <- Matrix(KO)
  outputList$manifoldAlignment <- MA
  outputList$diffRegulation <- DR
  return(outputList)
}


# Cox_analysis ------------------------------------------------------------


# --- Function: perform_cox ---
perform_cox <- function(formula, data, name = "") {
  cox_model <- coxph(formula, data = data)
  summary_cox <- summary(cox_model)
  
  if (length(cox_model$coefficients) == 0) {
    return(NULL)
  }
  
  result <- data.frame(
    variable = names(cox_model$coefficients),
    coef = summary_cox$coefficients[, "coef"],
    HR = summary_cox$coefficients[, "exp(coef)"],
    HR_lower = summary_cox$conf.int[, "lower .95"],
    HR_upper = summary_cox$conf.int[, "upper .95"],
    se = summary_cox$coefficients[, "se(coef)"],
    z_value = summary_cox$coefficients[, "z"],
    p_value = summary_cox$coefficients[, "Pr(>|z|)"],
    model = name,
    row.names = NULL
  )
  
  return(result)
}


# --- Function: auto_palette ---
auto_palette <- function(km_fit, colors = NULL) {
  
  # 获取strata names
  strata_names <- names(km_fit$strata)
  
  # 如果未提供颜色，使用默认调色板
  if (is.null(colors)) {
    # 内置调色板库
    default_palettes <- list(
      two_colors = c("#27ae60", "#e74c3c"),
      three_colors = c("#27ae60", "#f39c12", "#e74c3c"),
      four_colors = c("#27ae60", "#3498db", "#f39c12", "#e74c3c"),
      five_colors = c("#27ae60", "#3498db", "#9b59b6", "#f39c12", "#e74c3c")
    )
    
    # 根据strata数量选择调色板
    n_strata <- length(strata_names)
    if (n_strata <= 2) {
      colors <- default_palettes$two_colors
    } else if (n_strata == 3) {
      colors <- default_palettes$three_colors
    } else if (n_strata == 4) {
      colors <- default_palettes$four_colors
    } else {
      colors <- default_palettes$five_colors
    }
  }
  
  # 确保颜色数量足够
  if (length(colors) < length(strata_names)) {
    colors <- rep(colors, length.out = length(strata_names))
  }
  
  # 创建命名的调色板
  palette <- setNames(colors[1:length(strata_names)], strata_names)
  
  # 打印信息用于调试
  cat("Strata detected:", length(strata_names), "\n")
  cat("Strata names:", paste(strata_names, collapse=" | "), "\n")
  cat("Colors assigned:", paste(palette, collapse=" "), "\n\n")
  
  return(palette)
}


# --- Function: run_univariate_models ---
run_univariate_models <- function(data) {
  models <- list()
  models$arhgap18_cont <- glm(response_bin ~ keygene_expr, data = data, family = binomial())
  models$arhgap18_group <- glm(response_bin ~ keygene_group, data = data, family = binomial())
  models$age <- glm(response_bin ~ age, data = data, family = binomial())
  models$gender <- glm(response_bin ~ gender, data = data, family = binomial())
  models$runx1_mut <- glm(response_bin ~ runx1_mutation, data = data, family = binomial())
  models$ps29mrc <- glm(response_bin ~ ps29mrc_group, data = data, family = binomial())
  return(models)
}


# --- Function: extract_results ---
extract_results <- function(model_list, model_type = "Univariate") {
  results_list <- list()
  
  for (name in names(model_list)) {
    model <- model_list[[name]]
    
    tidy_result <- tidy(model, conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::filter(term != "(Intercept)") %>%
      mutate(
        model_name = name,
        model_type = model_type
      )
    
    results_list[[name]] <- tidy_result
  }
  
  return(bind_rows(results_list))
}


# --- Function: create_forest_plot_full ---
create_forest_plot_full <- function(df, title = "Forest Plot") {
  
  if (nrow(df) == 0) {
    cat("⚠ No data to plot\n")
    return(NULL)
  }
  
  # 按照逻辑排序变量
  var_order <- c(
    "ARHGAP18 (continuous, per unit)",
    "ARHGAP18 (High vs Low)",
    "Age (per year)",
    "Gender (Male vs Female)",
    "RUNX1 Mutation (Yes vs No)",
    "PS29MRC (Predicted Sensitive vs Others)"
  )
  
  df <- df %>%
    mutate(
      variable = factor(variable, levels = rev(var_order[var_order %in% unique(df$variable)])),
      significant = factor(significant, levels = c("p < 0.05", "p ≥ 0.05"))
    ) %>%
    arrange(variable)
  
  n_vars <- n_distinct(df$variable)
  
  p <- ggplot(df, aes(x = OR, y = variable, color = significant, shape = significant)) +
    geom_point(size = 3.5, alpha = 0.85) +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), 
                   height = 0.25, linewidth = 1, alpha = 0.8) +
    geom_vline(xintercept = 1, linetype = "dashed", 
               color = "black", linewidth = 0.8, alpha = 0.5) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                  labels = c("0.25", "0.5", "1.0", "2.0", "4.0", "8.0")) +
    scale_color_manual(
      values = c("p < 0.05" = "#e74c3c", "p ≥ 0.05" = "#95a5a6"),
      name = "Significance"
    ) +
    scale_shape_manual(
      values = c("p < 0.05" = 16, "p ≥ 0.05" = 1),
      name = "Significance"
    ) +
    labs(
      title = title,
      x = "Odds Ratio (95% CI, log scale)",
      y = "",
      caption = "Reference: Low for ARHGAP18; Female for Gender; No for RUNX1 Mutation"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.caption = element_text(size = 9, hjust = 0, color = "gray60"),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11, face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor.x = element_line(color = "gray95", linewidth = 0.2)
    )
  
  return(p)
}


# 细胞轨迹分析（monocle3） ------------------------------------------------------------------
monocle3_module <- function(object, OutPath, 
                            group_by = "group", 
                            celltype = "cellType_1", 
                            time_bin = "1", 
                            OnlyPlot = FALSE, 
                            genes_to_plot = NULL,
                            plot_heatmap = TRUE,  # 新增：热图开关
                            key_genes = c("ATF3")) {  # 新增：关键基因轨迹图
  
  setwd(OutPath)
  if(!dir.exists("Monocle3")) dir.create("Monocle3")
  
  if(!OnlyPlot){
    data <- as(as.matrix(object@assays$RNA@counts), 'sparseMatrix')
    pd <-  object@meta.data
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    cds <- new_cell_data_set(data, cell_metadata = pd, gene_metadata = fData)
    cds <- preprocess_cds(cds, num_dim = 20)
    cds <- reduce_dimension(cds, preprocess_method = "PCA", reduction_method = c("UMAP"))
    cds <- cluster_cells(cds, resolution = c(10^seq(-6,-1)))
    
    cds.embed <- cds@int_colData$reducedDims$UMAP
    int.embed <- Embeddings(object, reduction = "umap")
    int.embed <- int.embed[rownames(cds.embed),]
    cds@int_colData$reducedDims$UMAP <- int.embed
    cds <- learn_graph(cds, use_partition = F)
    
    get_earliest_principal_node <- function(cds, time_bin = "2"){
      cell_ids <- which(colData(cds)[, "seurat_clusters"] == time_bin)
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
    cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, time_bin = time_bin))
    saveRDS(cds, "Monocle3/Monocle3.rds")
  }
  
  cds <- readRDS("Monocle3/Monocle3.rds")
  
  # ========== 原有的3张图 ==========
  p <- plot_cells(cds,
                  color_cells_by = celltype,
                  label_groups_by_cluster = F,
                  label_leaves = F,
                  label_branch_points = F,
                  label_cell_groups = FALSE,
                  cell_size = 1, graph_label_size = 0) +
    theme_test(base_rect_size = 1) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c("right", "bottom"),
      legend.box.background = element_rect(colour = "black", linewidth = 1, fill = NA),
      legend.title = element_blank()
    )
  ggsave(p, file = "Monocle3/1.Monocle3.clusters.pdf", width = 5, height = 5, dpi = 600, limitsize = F)
  
  p <- plot_cells(cds,
                  color_cells_by = group_by,
                  label_groups_by_cluster = F,
                  label_leaves = F,
                  label_branch_points = FALSE,
                  label_cell_groups = FALSE,
                  cell_size = 1, graph_label_size = 0) +
    theme_test(base_rect_size = 1) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c("right", "bottom"),
      legend.box.background = element_rect(colour = "black", linewidth = 1, fill = NA),
      legend.title = element_blank()
    )
  ggsave(p, file = "Monocle3/2.Monocle3.group.pdf", width = 5, height = 5, dpi = 600, limitsize = F)
  
  p <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  cell_size = 1,
                  graph_label_size = 0)
  ggsave(p, file = "Monocle3/3.Monocle3.time.pdf", width = 6, height = 5, dpi = 600, limitsize = F)
  
  # ========== 新增：基因沿拟时序的表达曲线图 ==========
  if(!is.null(genes_to_plot)){
    library(ggplot2)
    library(patchwork)
    
    # 检查基因是否存在
    genes_available <- genes_to_plot[genes_to_plot %in% rownames(cds)]
    if(length(genes_available) == 0){
      warning("所有指定的基因都不在数据中！跳过曲线图和热图。")
    } else {
      if(length(genes_available) < length(genes_to_plot)){
        genes_missing <- genes_to_plot[!genes_to_plot %in% rownames(cds)]
        warning("以下基因不在数据中：", paste(genes_missing, collapse = ", "))
      }
      
      # 提取拟时序数据
      pseudotime_data <- data.frame(
        cell = colnames(cds),
        pseudotime = pseudotime(cds),
        cluster = colData(cds)[[celltype]],
        group = colData(cds)[[group_by]]
      )
      
      # 提取基因表达数据
      gene_expression <- as.matrix(normalized_counts(cds)[genes_available, , drop = FALSE])
      
      # 合并数据
      plot_data <- data.frame(
        cell = rep(colnames(cds), each = length(genes_available)),
        gene = rep(genes_available, ncol(cds)),
        expression = as.vector(gene_expression),
        pseudotime = rep(pseudotime_data$pseudotime, each = length(genes_available)),
        cluster = rep(pseudotime_data$cluster, each = length(genes_available)),
        group = rep(pseudotime_data$group, each = length(genes_available))
      )
      
      # 绘制每个基因的表达曲线
      p_list <- lapply(genes_available, function(gene){
        gene_data <- plot_data[plot_data$gene == gene, ]
        
        ggplot(gene_data, aes(x = pseudotime, y = expression)) +
          geom_point(aes(color = cluster), size = 0.5, alpha = 0.5) +
          geom_smooth(method = "loess", color = "black", se = TRUE, linewidth = 1) +
          labs(title = gene, x = "Pseudotime", y = "Expression") +
          theme_classic() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right"
          )
      })
      
      # 组合图形
      combined_plot <- wrap_plots(p_list, ncol = 2)
      ggsave("Monocle3/4.Gene_expression_pseudotime.pdf", 
             combined_plot, width = 10, height = 4 * ceiling(length(genes_available)/2), 
             dpi = 600, limitsize = FALSE)
      
      cat("✓ 表达曲线图已保存\n")
    }
  }
  
  # ========== 新增：基因表达热图（按拟时序排序）- 带开关 ==========
  if(!is.null(genes_to_plot) && plot_heatmap){
    library(pheatmap)
    library(viridis)
    
    # 检查基因是否存在
    genes_available <- genes_to_plot[genes_to_plot %in% rownames(cds)]
    
    if(length(genes_available) < 2){
      warning("热图需要至少2个基因，当前可用基因数：", length(genes_available), "。跳过热图绘制。")
    } else {
      # 按拟时序排序细胞
      cells_ordered <- names(sort(pseudotime(cds)))
      
      # 提取表达矩阵
      expr_matrix <- as.matrix(normalized_counts(cds)[genes_available, cells_ordered, drop = FALSE])
      
      # 标准化表达值（Z-score）
      expr_matrix_scaled <- t(scale(t(expr_matrix)))
      
      # 准备注释信息
      annotation_col <- data.frame(
        Cluster = colData(cds)[cells_ordered, celltype],
        Pseudotime = pseudotime(cds)[cells_ordered],
        Group = colData(cds)[cells_ordered, group_by],
        row.names = cells_ordered
      )
      
      # 设置颜色
      n_clusters <- length(unique(annotation_col$Cluster))
      cluster_colors <- scales::hue_pal()(n_clusters)
      names(cluster_colors) <- sort(unique(annotation_col$Cluster))
      
      group_colors <- c("Positive" = "#E64B35", "Negative" = "#4DBBD5")
      
      annotation_colors <- list(
        Cluster = cluster_colors,
        Pseudotime = viridis(100),
        Group = group_colors
      )
      
      # 绘制热图
      pdf("Monocle3/5.Gene_expression_heatmap.pdf", width = 12, height = max(6, length(genes_available) * 0.3))
      pheatmap(expr_matrix_scaled,
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               show_colnames = FALSE,
               annotation_col = annotation_col,
               annotation_colors = annotation_colors,
               color = colorRampPalette(c("blue", "white", "red"))(100),
               breaks = seq(-2, 2, length.out = 101),
               fontsize_row = 10,
               main = "Gene Expression along Pseudotime")
      dev.off()
      
      cat("✓ 热图已保存\n")
    }
  } else if(!is.null(genes_to_plot) && !plot_heatmap){
    cat("ℹ 热图绘制已关闭（plot_heatmap = FALSE）\n")
  }
  
  # ========== 新增：关键基因表达轨迹图（如ATF3）==========
  if(!is.null(key_genes)){
    library(ggplot2)
    library(patchwork)
    library(viridis)
    
    # 检查基因是否存在
    key_genes_available <- key_genes[key_genes %in% rownames(cds)]
    
    if(length(key_genes_available) == 0){
      warning("所有关键基因都不在数据中！")
    } else {
      # 为每个基因计算百分比表达
      for(gene in key_genes_available){
        gene_expr <- normalized_counts(cds)[gene, ]
        gene_pct <- (gene_expr / max(gene_expr)) * 100
        colData(cds)[[paste0(gene, "_pct")]] <- gene_pct
      }
      
      # 绘制每个基因
      p_list <- lapply(key_genes_available, function(gene){
        plot_cells(cds,
                   color_cells_by = paste0(gene, "_pct"),
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   graph_label_size = 0,
                   cell_size = 1.5,
                   alpha = 0.9) +
          scale_color_viridis_c(
            option = "D",
            name = "% Max",
            limits = c(0, 100),
            breaks = c(0, 25, 50, 75, 100)
          ) +
          theme_classic() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            legend.position = "right",
            legend.title = element_text(size = 12),
            axis.title = element_text(size = 12)
          ) +
          labs(title = gene, x = "UMAP 1", y = "UMAP 2")
      })
      
      # 根据基因数量选择布局
      if(length(p_list) == 1){
        combined_key <- p_list[[1]]
        ggsave("Monocle3/6.Key_genes_trajectory.pdf", 
               combined_key, width = 6, height = 5, dpi = 600)
      } else if(length(p_list) <= 4){
        combined_key <- wrap_plots(p_list, ncol = 2)
        ggsave("Monocle3/6.Key_genes_trajectory.pdf", 
               combined_key, width = 10, height = 5 * ceiling(length(p_list)/2), dpi = 600)
      } else {
        combined_key <- wrap_plots(p_list, ncol = 3)
        ggsave("Monocle3/6.Key_genes_trajectory.pdf", 
               combined_key, width = 15, height = 5 * ceiling(length(p_list)/3), dpi = 600)
      }
      
      cat("✓ 关键基因轨迹图已保存\n")
    }
  }
  
  cat("\n========== 分析完成 ==========\n")
  cat("输出文件：\n")
  cat("1. Monocle3.clusters.pdf - 细胞亚群轨迹图\n")
  cat("2. Monocle3.group.pdf - 分组轨迹图\n")
  cat("3. Monocle3.time.pdf - 拟时序轨迹图\n")
  if(!is.null(genes_to_plot)) cat("4. Gene_expression_pseudotime.pdf - 基因表达曲线图\n")
  if(!is.null(genes_to_plot) && plot_heatmap) cat("5. Gene_expression_heatmap.pdf - 基因表达热图\n")
  if(!is.null(key_genes)) cat("6. Key_genes_trajectory.pdf - 关键基因轨迹图\n")
  
  return(cds)
}




# 细胞轨迹分析（monocle2） --------------------------------------------------------

require(Seurat)
require(monocle)
require(dplyr)
require(ggplot2)
require(scales)
require(ClusterGVis)
set.seed(5201314)

my_pal2 <- c("#df4a86", "#746ea3", "#009ecb", "#00827b", "#3d4d7c", "#ad341d", "#0056a0","#d77b1c", "#077140","#d34132",
             "#698da5", "#ffd597", "#f1977f", "#828bae", "#82bfbd","#c1000a","#FF8E00","#00B3F1","#354270","#85b38f")


# responseMatrix_monocle <- function (models, newdata = NULL, response_type = "response",
#           cores = 1)
# {
#   res_list <- mclapply(models, function(x) {
#     if (is.null(x)) {
#       NA
#     }
#     else {
#       if (any(x@family@vfamily %in% c("negbinomial", "negbinomial.size"))) {
#         predict(x, newdata = newdata, type = response_type)
#       }
#       else if (x@family@vfamily %in% c("uninormal")) {
#         predict(x, newdata = newdata, type = response_type)
#       }
#       else {
#         10^predict(x, newdata = newdata, type = response_type)
#       }
#     }
#   }, mc.cores = cores)
#   res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE],
#                              length)
#   stopifnot(length(unique(res_list_lengths)) == 1)
#   num_na_fits <- length(res_list[is.na(res_list)])
#   if (num_na_fits > 0) {
#     na_matrix <- matrix(rep(rep(NA, res_list_lengths[[1]]),
#                             num_na_fits), nrow = num_na_fits)
#     row.names(na_matrix) <- names(res_list[is.na(res_list)])
#     non_na_matrix <- Matrix::t(do.call(cbind, lapply(res_list[is.na(res_list) ==
#                                                                 FALSE], unlist)))
#     row.names(non_na_matrix) <- names(res_list[is.na(res_list) ==
#                                                  FALSE])
#     res_matrix <- rbind(non_na_matrix, na_matrix)
#     res_matrix <- res_matrix[names(res_list), ]
#   }
#   else {
#     res_matrix <- Matrix::t(do.call(cbind, lapply(res_list,
#                                                   unlist)))
#     row.names(res_matrix) <- names(res_list[is.na(res_list) ==
#                                               FALSE])
#   }
#   res_matrix
# }
# 
# assignInNamespace("responseMatrix", responseMatrix_monocle, ns = "monocle")

#' monocle2_module
#'
#' @param object Seurat object
#' @param OutPath 设置输出路径
#' @param group_by 展示的分组，如：group,cellType_1等
#' @param track_gene 设置轨迹构建的基因，默认是VariableFeatures，还有dispersion，cluster差异基因的方法寻找差异基因。输入可以是数字，为VariableFeatures计算的前n个基因;输入可以是列表，包含方法和数字，如：c(dispersion,1000);输入可以是自定义基因列表。
#' @param ShowGene 展示的关键基因
#' @param OnlyPlot 是否只画图，复现填True
#'
#' @return cds object
#' @export
#'
#' @examples
monocle2_module <- function(object,OutPath,group_by,track_gene = NULL,ShowGene = NA,OnlyPlot = F){
  
  setwd(OutPath)
  if(!dir.exists("Monocle2"))dir.create("Monocle2")
  
  if(!OnlyPlot){
    scRNASub <- object
    monocle.matrix <- GetAssayData(object = scRNASub, slot = "data", assay = "RNA")
    monocle.sample <- scRNASub@meta.data
    monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
    
    data <- as(as.matrix(monocle.matrix), "sparseMatrix")
    pd <- new("AnnotatedDataFrame", data = monocle.sample)
    fd <- new("AnnotatedDataFrame", data = monocle.geneAnn)
    cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
    
    names(pData(cds))[names(pData(cds)) == group_by] <- "Cluster"
    
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    DefaultAssay(scRNASub) <- 'RNA'
    if (is.null(track_gene)){track_gene = c('VariableFeatures','all');genenubmers <- 'all'    #如果不输入，则默认
    }else if (length(track_gene)==1){  #如果长度为1，会判断是方法还是数字
      genenubmers <- 'all'
      if(!is.na(as.numeric(as.character(track_gene[1])))){genenubmers <- as.numeric(as.character(track_gene[1]))  #长度为1输入的是数字
      track_gene <- c('VariableFeatures',genenubmers)
      }else if (track_gene %in% c('dispersion','cluster','VariableFeatures')){track_gene <- c(track_gene,2000)}  #长度为1输入的是方法
    }else if (length(track_gene)==2 && !is.na(as.numeric(as.character(track_gene[2])))){    #长度为2，且第二个是数字
      genenubmers <- as.numeric(as.character(track_gene[2]))
      if (! track_gene[1] %in% c('dispersion','cluster','VariableFeatures')){track_gene <- c('VariableFeatures',2000)}  
    }
    #对输入长度进行判断，然后选择合适的方法基因
    if (length(track_gene) >= 200){
      genenubmers <- length(track_gene)
      track_gene <- track_gene
    }else if (track_gene[1] == 'dispersion'){
      disp_table <- dispersionTable(cds)
      disp_table <- arrange(disp_table,-dispersion_empirical)
      track_gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
    }else if (track_gene[1] == 'cluster'){
      Idents(scRNASub) <- 'seurat_clusters'
      deg.cluster <- FindAllMarkers(scRNASub)
      deg.cluster <- arrange(deg.cluster,p_val)
      track_gene <- subset(deg.cluster,p_val_adj<0.05)$gene
    }else if (track_gene[1] == 'VariableFeatures'){
      track_gene <- VariableFeatures(scRNASub)
    }else if (3 <= length(track_gene) && length(track_gene) < 200){
      stop('Plases check argument track_gene length')
    }
    if (length(track_gene) <= genenubmers | genenubmers == 'all'){track_gene <- track_gene
    }else if (length(track_gene) >= genenubmers ) {track_gene <- track_gene[1:genenubmers]}
    
    cds <- setOrderingFilter(cds, track_gene)
    cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
    cds <- orderCells(cds)
    saveRDS(cds,"Monocle2/cds.rds")
  }
  
  cds <- readRDS("Monocle2/cds.rds")
  dir.create(dir)
  p <- plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = .1) + ggsci::scale_color_gsea() + theme(legend.position = "right",text = element_text(face = "bold"))
  ggsave(paste0(dir,"1.Trajectory.Pseudotime.pdf"), plot = p,height = 4,width = 4.5)
  
  p <- plot_cell_trajectory(cds, color_by = "State",cell_size = .1) + guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right",text = element_text(face = "bold"))
  ggsave(paste0(dir,"2.Trajectory.State.pdf"), plot = p,height = 4,width = 4.5)
  
  # p <- plot_cell_trajectory(cds,color_by = "Cluster",cell_size = .1,show_branch_points = F) + 
  #   facet_wrap(~ Cluster, nrow = 1)  + 
  #   theme_void() +
  #   scale_color_manual(values = my_pal2) + 
  #   theme(legend.position = "none", text = element_text(face = "bold"))
  p <- plot_cell_trajectory(cds,color_by = "Cluster",cell_size = .1) + guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) + scale_color_manual(values = my_pal2) + theme(legend.position = "right",text = element_text(face = "bold"))
  ggsave(paste0(dir,"3.Trajectory.Cluster.pdf"), plot = p,height = 4,width = 4.5)
  
  if (!class(ShowGene) == "logical") {
    cds_subset <- cds[ShowGene, ]
    p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Cluster") + scale_color_manual(values = my_pal2)
    p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "State") + scale_color_manual(values = my_pal2)
    p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime") + ggsci::scale_color_gsea()
    p1 | p2 | p3
    ggsave(filename = paste0(dir,"4.ShowGene.pseudotime.pdf"), plot = p1 | p2 | p3, width = 12, height = length(ShowGene) * 2)
  }
  return(cds)
}




#' monocle2_branch_module
#'
#' @param cds cds对象
#' @param OutPath 设置输出路径
#' @param branch_point 设置需要查看的分支点
#' @param genes_branched 设置分支点差异计算的基因，用于后续展示
#' @param num_clusters 基因聚类成多少簇，默认3簇
#' @param heatmap_genes 设置分支热图展示的基因数量，可以输入数字展示前n个或输入基因列表
#' @param species 物种，人无需填写，鼠填写"mouse"
#' @param ncores 线程数量
#' @param OnlyPlot 是否只画图，复现填True
#'
#' @return 分支热图
#' @export
#'
#' @examples
monocle2_branch_module <- function(cds,OutPath='./',branch_point=1,genes_branched='all',num_clusters=3,heatmap_genes=100,species = NA,ncores = 16,OnlyPlot = F){
  
  object <- cds
  #识别CellDataSet格式，如果不是，直接退出
  if (class(object) == "CellDataSet" ) {
    #读取cds文件
    cds <- object
  }else if (is.character(object) && file.exists(object) && paste0(unlist(strsplit(object,split = ''))[c((nchar(object)-2):nchar(object))],collapse = '') == 'rds') {
    #读取地址末尾为rds的rds文件
    cds <- readRDS(file = object)
    if (class(cds) != "CellDataSet" ) {stop('Please Input an rds Rds-Filepath with a monocleobject')}
  }else{stop('Plases Input A CellDataSet Or A Rds-Filepath')
  }
  
  setwd(OutPath)
  if(!dir.exists("Monocle2"))dir.create("Monocle2")
  
  if(!OnlyPlot){
    if(!is.numeric(branch_point)){stop('!WARNNING:Please check argument branch_point : must be a number')}
    if (genes_branched=='all'){genes_branched <- rownames(cds)
    }else if (is.numeric(genes_branched)){track_gene <- rownames(cds@featureData@data)[cds@featureData@data$use_for_ordering]
    genes_branched <- track_gene[1:genes_branched]   ##轨迹基因不一定有2000,和cds使用的轨迹基因数量有关。如果报错请注意
    }else if(!all(genes_branched %in% rownames(cds))){stop('!WARNNING:Please check argument genes_branched : may be misspell wrong gene')}
    gene_cal <- BEAM(cds[genes_branched,],branch_point = branch_point,cores = ncores, progenitor_method = 'duplicate')
    gene_cal <- gene_cal[order(gene_cal$qval),]
    write.csv(gene_cal,file = paste0('Monocle2/heatmap_branch_point_',branch_point,'_gene.csv'))
    
    print("-- BEAM finished! --")
    gene_cal <- read.csv(paste0('Monocle2/heatmap_branch_point_',branch_point,'_gene.csv'),row.names = 1)
    
    if(is.numeric(heatmap_genes)){heatmap_genes <- rownames(gene_cal[1:heatmap_genes,])
    }else if(!all(heatmap_genes %in% rownames(cds))){print('!WARNNING:Please check argument heatmap_genes :The enter genes may be wrong.Forced output')
      heatmap_genes <- rownames(gene_cal)[1:100]
    }else if(all(heatmap_genes %in% rownames(cds))){heatmap_genes <- heatmap_genes }
    
    len <- length(heatmap_genes)
    
    df <- plot_genes_branched_heatmap2(cds[heatmap_genes,],
                                       branch_point = branch_point,
                                       num_clusters = num_clusters,
                                       cores = ncores,
                                       use_gene_short_name = T,
                                       show_rownames = T)
    saveRDS(df,paste0('Monocle2/Branch.plotdata.rds'))
    print("-- Branch-plotdata finished! --")
    
    if (!is.na(species)) {
      mouse <- ifelse(species == "mouse", T, F)
    } else {
      mouse <- F
    }
    if(mouse){
      require(org.Mm.eg.db)
      enrich <- enrichCluster(object = df,
                              OrgDb = org.Mm.eg.db,
                              type = "KEGG",
                              organism = "mmu",
                              pvalueCutoff = 0.05,
                              topn = 5,
                              seed = 5201314,
                              add.gene = TRUE)
    }else{
      require(org.Hs.eg.db)
      enrich <- enrichCluster(object = df,
                              OrgDb = org.Hs.eg.db,
                              type = "KEGG",
                              organism = "hsa",
                              pvalueCutoff = 0.05,
                              topn = 5,add.gene = T,
                              seed = 5201314,
                              add.gene = TRUE)
    }
    write.csv(enrich,paste0('Monocle2/Enrich.plotdata.csv'))
    print("-- Enrichment finished! --")
  }
  
  df <- readRDS(paste0('Monocle2/Branch.plotdata.rds'))
  enrich <- read.csv(paste0('Monocle2/Enrich.plotdata.csv'),row.names = 1)
  markGenes <- intersect(df$wide.res$gene,do.call(c,strsplit(enrich$geneID,"/",F)))
  enrich$geneID <- NULL
  pdf(paste0(dir,'5.Branch-enrich.pdf'),height = 7,width = 13,onefile = F)
  visCluster(object = df,
             plot.type = "both",
             column_names_rot = 45,
             show_row_dend = F,
             markGenes = markGenes,
             markGenes.side = "left",
             annoTerm.data = enrich,
             add.bar = T,
             line.side = "left")
  dev.off()
  if(!is.null(dev.list())) {dev.off()}
}

# limma -------------------------------------------------------------------
#count数据不需要进行log2(已经voom)可直接用，芯片数据看情况需要进行log2转换
limma_voom_module1 <- function(mat, metadata, outpath, logFoldChange = 0.585, adjustP = 0.05, 
                               is_count = TRUE, species = "human", input_id_type = "entrezgene_id", 
                               use_local = TRUE) {
  
  outpath <- paste0(outpath, "/1.DiffAnalysis/")
  if (!dir.exists(outpath)) dir.create(outpath)
  
  # ✅ 确保metadata行顺序与mat列顺序一致
  cat("正在对齐样本顺序...\n")
  if (!all(colnames(mat) %in% metadata[, 1])) {
    stop("⚠ 错误：表达矩阵样本名与metadata不匹配！")
  }
  metadata <- metadata[match(colnames(mat), metadata[, 1]), ]
  cat("✓ 样本顺序已对齐\n\n")
  
  # 创建分组因子
  modType <- factor(metadata[, 6], c("Control", "Disease"))
  design <- model.matrix(~ 0 + modType)
  colnames(design) <- c("Control", "Disease")
  
  cat("Design矩阵维度:", nrow(design), "×", ncol(design), "\n")
  cat("表达矩阵维度:", ncol(mat), "样本 ×", nrow(mat), "基因\n")
  cat("分组情况：", paste(names(table(modType)), "=", table(modType), collapse="; "), "\n\n")
  
  # ★ 根据数据类型选择处理方法 ★
  if (is_count) {
    cat("检测到count数据，使用voom进行标准化...\n")
    
    # 过滤低表达基因
    keep <- edgeR::filterByExpr(mat, group=modType)
    mat_filtered <- mat[keep, ]
    cat("过滤后保留", nrow(mat_filtered), "个基因\n")
    
    # 创建DGEList对象
    dge <- edgeR::DGEList(counts=mat_filtered, group=modType)
    dge <- edgeR::calcNormFactors(dge)
    
    # voom转换
    pdf(paste0(outpath, "0.Voom_plot.pdf"), width = 6, height = 6)
    v <- voom(dge, design, plot = TRUE)
    dev.off()
    
    # ✅ 关键：voom返回EList对象，提取其中的数据
    expr_matrix <- v$E      # 表达矩阵
    weights_matrix <- v$weights  # 权重矩阵
    
    cat("✓ voom转换完成\n")
    cat("  表达矩阵:", nrow(expr_matrix), "×", ncol(expr_matrix), "\n")
    cat("  权重矩阵:", nrow(weights_matrix), "×", ncol(weights_matrix), "\n\n")
    
    # lmFit时使用 weights 参数
    fit <- lmFit(expr_matrix, design, weights=weights_matrix)
    
  } else {
    cat("检测到microarray数据，使用标准normalizeBetweenArrays...\n")
    expr_matrix <- normalizeBetweenArrays(mat)
    
    # lmFit不使用weights
    fit <- lmFit(expr_matrix, design)
  }
  
  # 后续分析
  cont.matrix <- makeContrasts(Disease - Control, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # 获取所有差异分析结果
  allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
  
  # ★ 基因ID转换 ★
  cat("\n正在进行基因ID转换...\n")
  cat("物种：", species, "\n")
  cat("输入ID类型：", input_id_type, "\n")
  cat("使用本地数据库：", use_local, "\n\n")
  
  allDiff$GeneID <- rownames(allDiff)
  
  # 尝试本地转换
  if (use_local) {
    cat("尝试使用本地注释包进行转换...\n")
    
    if (species == "human") {
      if (!require("org.Hs.eg.db", quietly = TRUE)) {
        cat("正在安装org.Hs.eg.db包...\n")
        BiocManager::install("org.Hs.eg.db", quietly = TRUE)
        library(org.Hs.eg.db)
      }
      annotation_db <- org.Hs.eg.db
    } else if (species == "mouse") {
      if (!require("org.Mm.eg.db", quietly = TRUE)) {
        BiocManager::install("org.Mm.eg.db", quietly = TRUE)
        library(org.Mm.eg.db)
      }
      annotation_db <- org.Mm.eg.db
    } else if (species == "rat") {
      if (!require("org.Rn.eg.db", quietly = TRUE)) {
        BiocManager::install("org.Rn.eg.db", quietly = TRUE)
        library(org.Rn.eg.db)
      }
      annotation_db <- org.Rn.eg.db
    }
    
    tryCatch({
      if (input_id_type == "entrezgene_id") {
        gene_mapping <- AnnotationDbi::mapIds(annotation_db,
                                              keys = rownames(allDiff),
                                              column = "SYMBOL",
                                              keytype = "ENTREZID",
                                              multiVals = "first")
      } else if (input_id_type == "ensembl_gene_id") {
        gene_mapping <- AnnotationDbi::mapIds(annotation_db,
                                              keys = rownames(allDiff),
                                              column = "SYMBOL",
                                              keytype = "ENSEMBL",
                                              multiVals = "first")
      }
      
      allDiff$SYMBOL <- as.character(gene_mapping)
      allDiff$SYMBOL[is.na(allDiff$SYMBOL)] <- rownames(allDiff)[is.na(allDiff$SYMBOL)]
      
      cat("✓ 本地转换成功！共转换", sum(!is.na(gene_mapping)), "个基因ID\n\n")
      
    }, error = function(e) {
      cat("⚠ 本地转换出现问题：", e$message, "\n")
      cat("将使用原始基因ID作为符号\n\n")
      allDiff$SYMBOL <<- rownames(allDiff)
    })
  }
  
  allDiff <- allDiff[, c("GeneID", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  
  write.table(allDiff, file = paste0(outpath, "LimmaTab.txt"), sep = "\t", quote = F, row.names = FALSE)
  
  # 筛选差异基因
  diffSig <- allDiff[with(allDiff, (abs(logFC) > logFoldChange & P.Value < adjustP)), ]
  write.table(diffSig, file = paste0(outpath, "Diff.txt"), sep = "\t", quote = F, row.names = FALSE)
  
  diffUp <- allDiff[with(allDiff, (logFC > logFoldChange & P.Value < adjustP)), ]
  write.table(diffUp, file = paste0(outpath, "Up.txt"), sep = "\t", quote = F, row.names = FALSE)
  
  diffDown <- allDiff[with(allDiff, (logFC < (-logFoldChange) & P.Value < adjustP)), ]
  write.table(diffDown, file = paste0(outpath, "Down.txt"), sep = "\t", quote = F, row.names = FALSE)
  
  # ★ 热图数据处理 ★
  if (nrow(diffSig) > 0) {
    if (is_count) {
      hmExp <- v$E[rownames(diffSig), ]  # ✅ 从v$E提取
    } else {
      hmExp <- expr_matrix[rownames(diffSig), ]
    }
    
    # 标准化数据
    hmExp_scaled <- t(scale(t(hmExp)))
    hmExp_scaled[hmExp_scaled > 5] <- 5
    hmExp_scaled[hmExp_scaled < -5] <- -5
    
    # 设置行名为SYMBOL
    rownames(hmExp_scaled) <- diffSig[rownames(hmExp_scaled), "SYMBOL"]
    
    # 保存标准化后的热图数据
    write.table(hmExp_scaled, file = paste0(outpath, "Heatmap_data.txt"), 
                sep = "\t", quote = F)
  }
  
  # 火山图
  pdf(paste0(outpath, "1.Vol.pdf"), width = 6, height = 6)
  xMax <- max(abs(allDiff$logFC))
  yMax <- max(-log10(allDiff$P.Value))
  plot(allDiff$logFC, -log10(allDiff$P.Value),
       xlab = "log2FC", ylab = "-log10(P.Value)",
       main = "Volcano Plot", xlim = c(-xMax, xMax), ylim = c(0, yMax), 
       yaxs = "i", pch = 20, cex = 0.8)
  diffSub <- subset(allDiff, P.Value < adjustP & logFC > logFoldChange)
  points(diffSub$logFC, -log10(diffSub$P.Value), pch = 20, col = "palevioletred1", cex = 0.8)
  diffSub <- subset(allDiff, P.Value < adjustP & logFC < (-logFoldChange))
  points(diffSub$logFC, -log10(diffSub$P.Value), pch = 20, col = "dodgerblue", cex = 0.8)
  abline(v = 0, lty = 2, lwd = 3)
  abline(h = -log10(adjustP), lty = 2, lwd = 1)
  legend("topright", c("Up-regulated", "Down-regulated"), 
         col = c("palevioletred1", "dodgerblue"), pch = 20)
  dev.off()
  
  # ★ 热图 ★
  if (nrow(diffSig) > 0) {
    plotdata <- read.table(paste0(outpath, "Heatmap_data.txt"), sep = "\t", 
                           header = T, check.names = F, row.names = 1)
    
    up <- read.table(paste0(outpath, "Up.txt"), header = T, 
                     stringsAsFactors = F, sep = "\t")
    down <- read.table(paste0(outpath, "Down.txt"), header = T, 
                       stringsAsFactors = F, sep = "\t")
    
    # 获取UP和DOWN基因的SYMBOL
    up_symbols <- up$SYMBOL
    down_symbols <- down$SYMBOL
    
    # 确保热图数据只包含显著基因
    plotdata <- plotdata[c(up_symbols, down_symbols), ]
    
    annCol <- data.frame(
      Group = metadata[match(colnames(plotdata), metadata[, 1]), 6],
      row.names = colnames(plotdata),
      stringsAsFactors = F
    )
    
    annRow <- data.frame(
      Direction = c(rep("Up", length(up_symbols)), rep("Down", length(down_symbols))),
      row.names = c(up_symbols, down_symbols),
      stringsAsFactors = F
    )
    
    annColors <- list(
      "Group" = c("Control" = "blue", "Disease" = "red"),
      "Direction" = c("Up" = "yellow", "Down" = "green")
    )
    
    pdf(file = paste0(outpath, "2.Heatmap.pdf"), width = 8, height = 10)
    pheatmap::pheatmap(plotdata,
                       scale = "none",
                       annotation_row = annRow,
                       annotation_col = annCol,
                       annotation_colors = annColors,
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                       fontsize_row = 6,
                       fontsize_col = 8,
                       fontsize = 8,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_colnames = FALSE,
                       show_rownames = F
    )
    dev.off()
  }
  
  # 统计摘要
  cat("\n========== 差异分析结果摘要 ==========\n")
  cat("总基因数：", nrow(allDiff), "\n")
  cat("显著差异基因数：", nrow(diffSig), "\n")
  cat("上调基因数：", nrow(diffUp), "\n")
  cat("下调基因数：", nrow(diffDown), "\n")
  cat("=====================================\n\n")
  
  return(list(allDiff=allDiff, diffSig=diffSig, diffUp=diffUp, diffDown=diffDown))
}

limma_voom_module <- function(mat, metadata, outpath, logFoldChange = 0.585, adjustP = 0.05,
                              is_count = TRUE, species = "human", input_id_type = "entrezgene_id",
                              use_local = TRUE) {
  
  dir.create(outpath <- paste0(outpath, "/1.DiffAnalysis/"), recursive = TRUE, showWarnings = FALSE)
  
  # ── 样本对齐 ──────────────────────────────────────────────────────────────
  if (!all(colnames(mat) %in% metadata[[1]])) stop("⚠ 样本名与metadata不匹配！")
  metadata <- metadata[match(colnames(mat), metadata[[1]]), ]
  
  modType <- factor(metadata[[2]], c("Control", "Disease"))
  design  <- model.matrix(~ 0 + modType) |> `colnames<-`(c("Control", "Disease"))
  cat(sprintf("Design: %d×%d | 基因: %d | Control=%d Disease=%d\n\n",
              nrow(design), ncol(design), nrow(mat),
              sum(modType == "Control"), sum(modType == "Disease")))
  
  # ── 标准化 & 拟合 ──────────────────────────────────────────────────────────
  if (is_count) {
    cat("count数据 → voom\n")
    dge <- edgeR::DGEList(counts = mat[edgeR::filterByExpr(mat, group = modType), ], group = modType) |>
      edgeR::calcNormFactors()
    pdf(paste0(outpath, "0.Voom_plot.pdf"), 6, 6); v <- limma::voom(dge, design, plot = TRUE); dev.off()
    fit <- limma::lmFit(v, design)
  } else {
    cat("microarray数据 → normalizeBetweenArrays\n")
    rt  <- limma::normalizeBetweenArrays(mat)
    fit <- limma::lmFit(rt, design)
  }
  
  fit2   <- limma::contrasts.fit(fit, limma::makeContrasts(Disease - Control, levels = design)) |> limma::eBayes()
  allDiff <- limma::topTable(fit2, adjust = "fdr", number = Inf)
  allDiff$GeneID <- rownames(allDiff)
  
  # ── 基因ID转换 ─────────────────────────────────────────────────────────────
  db_map      <- c(human = "org.Hs.eg.db", mouse = "org.Mm.eg.db", rat = "org.Rn.eg.db")
  key_map     <- c(entrezgene_id = "ENTREZID", ensembl_gene_id = "ENSEMBL", external_gene_name = "SYMBOL")
  bm_map      <- c(human = "hsapiens_gene_ensembl", mouse = "mmusculus_gene_ensembl", rat = "rnorvegicus_gene_ensembl")
  db_pkg      <- db_map[[species]]
  
  allDiff$SYMBOL <- rownames(allDiff)  # 默认值
  
  if (use_local && key_map[[input_id_type]] != "SYMBOL") {
    if (!requireNamespace(db_pkg, quietly = TRUE)) BiocManager::install(db_pkg)
    tryCatch({
      sym <- AnnotationDbi::mapIds(get(db_pkg), keys = rownames(allDiff),
                                   column = "SYMBOL", keytype = key_map[[input_id_type]], multiVals = "first")
      allDiff$SYMBOL <- ifelse(is.na(sym), rownames(allDiff), as.character(sym))
      cat(sprintf("✓ 本地转换：%d 个基因\n", sum(!is.na(sym))))
    }, error = function(e) cat("⚠ 本地转换失败，尝试在线...\n"))
  }
  
  if (all(allDiff$SYMBOL == rownames(allDiff)) && input_id_type != "external_gene_name") {
    tryCatch({
      mart <- biomaRt::useMart("ensembl", dataset = bm_map[[species]])
      gm   <- biomaRt::getBM(c(input_id_type, "external_gene_name"), input_id_type, rownames(allDiff), mart, uniqueRows = TRUE)
      gm   <- gm[!duplicated(gm[[1]]), ] |> (\(x) setNames(x[[2]], x[[1]]))()
      allDiff$SYMBOL <- ifelse(is.na(gm[rownames(allDiff)]), rownames(allDiff), gm[rownames(allDiff)])
      cat(sprintf("✓ 在线转换：%d 个基因\n", length(gm)))
    }, error = function(e) cat("⚠ 在线转换失败，使用原始ID\n"))
  }
  
  # ── 筛选 & 写出 ────────────────────────────────────────────────────────────
  allDiff  <- allDiff[, c("GeneID", "SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  diffSig  <- subset(allDiff,  abs(logFC) > logFoldChange & P.Value < adjustP)
  diffUp   <- subset(allDiff,      logFC  > logFoldChange & P.Value < adjustP)
  diffDown <- subset(allDiff,      logFC  < -logFoldChange & P.Value < adjustP)
  
  wt <- function(x, f) write.table(x, paste0(outpath, f), sep = "\t", quote = FALSE, row.names = FALSE)
  wt(allDiff, "LimmaTab.txt"); wt(diffSig, "Diff.txt"); wt(diffUp, "Up.txt"); wt(diffDown, "Down.txt")
  
  # ── 火山图 ─────────────────────────────────────────────────────────────────
  pdf(paste0(outpath, "1.Vol.pdf"), 6, 6)
  plot(allDiff$logFC, -log10(allDiff$P.Value), pch = 20, cex = 0.8,
       xlab = "log2FC", ylab = "-log10(P.Value)", main = "Volcano Plot",
       xlim = c(-1, 1) * max(abs(allDiff$logFC)), ylim = c(0, max(-log10(allDiff$P.Value))), yaxs = "i")
  points(diffUp$logFC,   -log10(diffUp$P.Value),   pch = 20, col = "palevioletred1", cex = 0.8)
  points(diffDown$logFC, -log10(diffDown$P.Value),  pch = 20, col = "dodgerblue",     cex = 0.8)
  abline(v = 0, lty = 2, lwd = 3); abline(h = -log10(adjustP), lty = 2)
  legend("topright", c("Up", "Down"), col = c("palevioletred1", "dodgerblue"), pch = 20)
  dev.off()
  
  # ── 热图 ───────────────────────────────────────────────────────────────────
  if (nrow(diffSig) > 0) {
    hmExp <- (if (is_count) v$E else rt)[rownames(diffSig), ]
    hmExp <- t(scale(t(hmExp))); hmExp[hmExp > 5] <- 5; hmExp[hmExp < -5] <- -5
    rownames(hmExp) <- diffSig[rownames(hmExp), "SYMBOL"]
    wt(as.data.frame(hmExp), "Heatmap_data.txt")
    
    ord <- c(diffUp$SYMBOL, diffDown$SYMBOL)
    pdf(paste0(outpath, "2.Heatmap.pdf"), 8, 10)
    pheatmap::pheatmap(hmExp[ord, ],
                       scale = "none", cluster_cols = FALSE, cluster_rows = FALSE,
                       show_rownames = FALSE, show_colnames = FALSE, fontsize = 8,
                       annotation_col = data.frame(Group = metadata[[2]], row.names = metadata[[1]]),
                       annotation_row = data.frame(Direction = rep(c("Up","Down"), c(nrow(diffUp), nrow(diffDown))), row.names = ord),
                       annotation_colors = list(Group = c(Control="blue", Disease="red"), Direction = c(Up="yellow", Down="green")),
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
    dev.off()
  }
  
  # ── 摘要 ───────────────────────────────────────────────────────────────────
  cat(sprintf("\n总基因: %d | 差异: %d | 上调: %d | 下调: %d\n",
              nrow(allDiff), nrow(diffSig), nrow(diffUp), nrow(diffDown)))
  
  invisible(list(allDiff = allDiff, diffSig = diffSig, diffUp = diffUp, diffDown = diffDown))
}
# 单细胞流程化 ------------------------------------------------------------------

require(Seurat)
require(ggplot2)
require(harmony)
require(dplyr)
require(reshape2)
require(scRNAtoolVis)
require(ggpubr)
# require(clusterProfiler)
options(dplyr.print_max = 1e9)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)

ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}

paramSweep_V4 <- function(seu, PCs = 1:10, sct = FALSE, num.cores = 1) {
  require(Seurat)
  require(fields)
  require(parallel)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data) / (1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  # no_cores <- detectCores()-1
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct, mc.cores = num.cores
    )
    stopCluster(cl)
  } else {
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct
    )
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}

doubletFinder_V4 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat)
  require(fields)
  require(KernSmooth)
  
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells / (1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin
      )
      
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
      )
      
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block
      )
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca = orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose = FALSE
      )
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    
    ## Compute PC distance matrix
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells)) / k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets <- neighbors[neighbors > n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            neighbor_types[i, ] <-
              table(doublet_types1[neighbors_that_are_doublets - n_real.cells]) +
              table(doublet_types2[neighbors_that_are_doublets - n_real.cells])
            neighbor_types[i, ] <- neighbor_types[i, ] / sum(neighbor_types[i, ])
          } else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) <- levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, ct]
      }
    }
    return(seu)
  }
}


human2mouse <- function(x) {
  x1 <- substr(x, 1, 1)
  x2 <- substr(x, 2, nchar(x))
  return(paste0(x1, tolower(x2)))
}

getScatterplot <- function(object, gene1, gene2, cor.method = "pearson", jitter.num = 0.15, pos = TRUE) {
  if (!gene1 %in% rownames(object)) {
    print("gene1 was not found")
    if (!gene2 %in% rownames(object)) {
      print("gene2 was not found")
    }
  } else {
    exp.mat <- GetAssayData(object = object, assay = "RNA") %>%
      .[c(gene1, gene2), ] %>%
      as.matrix() %>%
      t() %>%
      as.data.frame()
    if (pos) {
      if (nrow(exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]) > (nrow(exp.mat) * 0.01)) {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 & exp.mat[, 2] > 0), ]
      } else {
        exp.mat <- exp.mat[which(exp.mat[, 1] > 0 | exp.mat[, 2] > 0), ]
      }
    }
    colnames(exp.mat) <- c("Var1", "Var2")
    plots <- ggplot(data = exp.mat, mapping = aes_string(x = "Var1", y = "Var2")) +
      geom_smooth(method = "lm", se = T, color = "red", size = 1) +
      stat_cor(method = cor.method) +
      labs(x = gene1, y = gene2) +
      geom_jitter(width = jitter.num, height = jitter.num, color = "black", size = 1, alpha = 1) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.text = element_text(colour = "black", size = 10),
        axis.text = element_text(colour = "black", size = 10),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(size = 1, linetype = "solid", colour = "black"),
        panel.background = element_rect(fill = "white")
      )
    return(plots)
  }
}

capitalize_first <- function(string) {
  words <- strsplit(string, " ")[[1]]
  capitalized_words <- toupper(substring(words, 1, 1))
  rest_of_words <- substring(words, 2)
  capitalized_string <- paste0(capitalized_words, rest_of_words, collapse = " ")
  return(capitalized_string)
}



#' Main Function
#'
#' @param Path 保存结果的路径
#' @param SeuratObject 创建好seuratV4对象且拥有分组的rds，分组名字为group，肿瘤样本为Normal vs Tumor，非肿瘤样本为Control vs Disease
#' @param Multi 是否为多样本
#' @param ref SingleR参考数据集
#' @param labels 参考数据集的分类，label.main为主要细胞亚群，label.fine为细分亚群
#' @param min.features  细胞所含的最低基因数量
#' @param Idents 注释分辨率选择
#' @param DoubletFinder 是否过滤双细胞
#' @param species 物种，人无需填写，鼠填写"mouse"
#' @param KeyCell 关键细胞，仅在拟时序分析中使用
#' @param KeyGene 关键基因
#' @param logFCfilter logFC阈值
#' @param adjPvalFilter p值阈值
#' @param assay 降维聚类使用的assay，默认RNA
#' @param set.resolutions 分群
#' @param PC 选择主成分个数
#' @param nfeatures 高变基因个数
#' @param npcs 碎石图的主成分个数
#'
#' @return
#' @export
#'
#' @examples
scRNAAutoAnno <- function(Path = ".", SeuratObject = NA, Multi = TRUE, ref, labels,min.features = 200, Idents = "RNA_snn_res.0.2", DoubletFinder = TRUE, species = NA, KeyCell = NA, KeyGene = NA, logFCfilter = 0.585, adjPvalFilter = 0.05,assay = "RNA", set.resolutions = seq(0.2, 1.2, by = 0.1), PC = 20, nfeatures = 2000, npcs = 50) {
  cat("1.质量控制\n")
  cat("2.数据标准化\n")
  cat("3.自动化注释\n")
  cat("4.Cellchat\n")
  cat("5.亚群贡献度\n")
  cat("6.拟时序分析\n")
  cat("7.关键基因表达丰度\n")
  cat("8.疾病基因共表达网络\n")
  cat("9.免疫代谢通路\n")
  cat("10.代谢通路热图\n")
  cat("11.量化热点机制\n")
  cat("12.富集分析\n")
  cat("\033[0;34m请选择要分析的内容，按序号填写，英文逗号分隔")
  AnalysisIndex <- scan(sep = ",", quiet = TRUE)
  
  # AdvanceIndex <- length(which(AnalysisIndex %in% c(1:12) == FALSE)) # 判断是否有其他高级分析
  # if (is.na(KeyCell) & AdvanceIndex > 0) {
  #   message("是否需要在分析过程中挑选关键细胞？yes/no")
  #   isKeyCell <- scan(sep = ",", what = "character", quiet = TRUE)
  #   if (isKeyCell == "yes") {
  #     cat("1.细胞含量变化最大\n")
  #     cat("2.Cellchat\n")
  #     cat("3.亚群贡献度\n")
  #     cat("\033[0;34m从哪一步分析中挑选关键细胞？\n")
  #     KeyCellIndex <- scan(sep = ",", quiet = TRUE)
  #   }
  # }
  
  scRNA <- SeuratObject
  
  if (!is.na(species)) {
    mouse <- ifelse(species == "mouse", T, F)
  } else {
    mouse <- F
  }
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  if (mouse) {
    s.genes <- human2mouse(s.genes)
    g2m.genes <- human2mouse(g2m.genes)
  }
  
  n <- 0
  KeyCells <- KeyCell
  for (step in AnalysisIndex) {
    n <- n + 1
    
    # 质量控制
    if (step == 1) {
      message("*************************质量控制*************************")
      dir <- paste0(Path, "/.QCFilter")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      grep("^[M,m][T,t]-", rownames(scRNA), value = T)
      grep("^R[P,p][SL,sl]", rownames(scRNA), value = T)
      scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, pattern = "^[M,m][T,t]-")
      scRNA[["percent.ribo"]] <- PercentageFeatureSet(scRNA, pattern = "^R[P,p][SL,sl]") # 核糖体基因占比
      
      nomt <- ifelse(max(scRNA[["percent.mt"]]) == 0, T, F)
      if (nomt) {
        MT_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
        if (mouse) {
          MT_genes <- human2mouse(MT_genes)
        }
        MT_genes <- intersect(MT_genes, rownames(scRNA))
        scRNA[["percent.mt"]] <- PercentageFeatureSet(object = scRNA, features = MT_genes)
      }
      
      mt_value <- stats::mad(scRNA$percent.mt)
      lower_mt <- median(scRNA$percent.mt) - 3 * mt_value
      upper_mt <- median(scRNA$percent.mt) + 3 * mt_value
      
      nF_value <- stats::mad(scRNA$nFeature_RNA)
      lower_nF <- median(scRNA$nFeature_RNA) - 3 * nF_value
      upper_nF <- median(scRNA$nFeature_RNA) + 3 * nF_value
      
      nC_value <- stats::mad(scRNA$nCount_RNA)
      lower_nC <- median(scRNA$nCount_RNA) - 3 * nC_value
      upper_nC <- median(scRNA$nCount_RNA) + 3 * nC_value
      
      # nR_value <- stats::mad(scRNA$percent.ribo)
      # lower_nR <- median(scRNA$percent.ribo) - 3 * nR_value
      # upper_nR <- median(scRNA$percent.ribo) + 3 * nR_value
      dir1 <- paste0(Path, "/.QCunFilter")
      if (!dir.exists(dir1)) {
        dir.create(dir1)
      } else {
        print("Dir already exists!")
      }
      p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "Individual_ID", ncol = 4,raster=FALSE)
      # p
      ggsave_fun(filename = paste0(dir1, "/1.VlnPlot"), plot = p, width = 20, height = 4)
      
      plot1 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
        NoLegend()
      plot2 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
        NoLegend()
      plot3 <-
        FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.ribo", raster = F) +
        NoLegend()
      
      plot1 + plot2 + plot3
      ggsave_fun(filename = paste0(dir1, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      
      
      scRNA <- subset(x = scRNA, subset = nFeature_RNA >= min.features & percent.mt <= upper_mt & nFeature_RNA <= upper_nF & nCount_RNA <= upper_nC)
      dim(scRNA)
      
      if (DoubletFinder) {
        require(DoubletFinder)
        scRNA_list <- SplitObject(scRNA, split.by = "orig.ident")
        newscRMA_list <- list()
        for (i in 1:length(scRNA_list)) {
          data <- scRNA_list[[i]]
          data <- NormalizeData(data)
          data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
          data <- ScaleData(data)
          data <- RunPCA(data)
          data <- RunUMAP(data, dims = 1:20)
          
          sweep.res.list <- paramSweep_V4(data, PCs = 1:20, sct = FALSE) # 若使用SCT方法标准化则'sct=T'
          sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
          bcmvn <- find.pK(sweep.stats)
          p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
          
          # 期望doublet数量
          homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) # 可使用注释好的细胞类型
          Doubletrate <- ncol(data) * 8 * 1e-6
          nExp_poi <- round(Doubletrate * ncol(data))
          nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
          
          # 鉴定doublets
          data <- doubletFinder_V4(data, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
          colnames(data@meta.data)[ncol(data@meta.data)] <- "doublet_info"
          newscRMA_list[[i]] <- data
        }
        
        saveRDS(newscRMA_list, paste0(dir, "/doubletFinder_V4.rds"))
        
        plot_list <- list()
        for (i in 1:length(newscRMA_list)) {
          plot_list[[i]] <- DimPlot(newscRMA_list[[i]], group.by = "doublet_info") + ggtitle(newscRMA_list[[i]]$orig.ident[1]) + theme(plot.title = element_text(hjust = 0.5))
        }
        
        
        for (i in 1:length(plot_list)) {
          plot_list[[i]] <- plot_list[[i]] + NoAxes() + NoLegend()
        }
        
        p <- CombinePlots(plot_list, legend = "right")
        
        ggsave_fun(filename = paste0(dir, "/1.DoubletFinder"), plot = p, width = 12, height = 10)
        
        if (Multi) {
          scRNA <- merge(newscRMA_list[[1]], newscRMA_list[2:length(newscRMA_list)])
        } else {
          scRNA <- newscRMA_list[[1]]
        }
        
        scRNA <- subset(scRNA, subset = doublet_info == "Singlet")
        
        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4,raster=F)
        # p
        ggsave_fun(filename = paste0(dir, "/2.VlnPlot"), plot = p, width = 20, height = 4)
        
        plot1 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
          NoLegend()
        plot2 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
          NoLegend()
        plot3 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) +
          NoLegend()
        
        plot1 + plot2 + plot3
        ggsave_fun(filename = paste0(dir, "/3.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      } else {
        p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "Individual_ID", ncol = 4, raster=F)
        # p
        ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 20, height = 4)
        
        plot1 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
          NoLegend()
        plot2 <-
          FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
          NoLegend()
        plot3 <-
          FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.ribo", raster = F) +
          NoLegend()
        
        plot1 + plot2 + plot3
        ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
      }
    }
    
    # 数据标准化
    if (step == 2) {
      message("*************************数据标准化*************************")
      dir <- paste0(Path, "/.Cluster")
      
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      DefaultAssay(scRNA) <- assay
      
      scRNA <- NormalizeData(object = scRNA)
      scRNA <- CellCycleScoring(object = scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
      scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = nfeatures)
      top10 <- head(VariableFeatures(scRNA), 10)
      plot1 <- VariableFeaturePlot(scRNA)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot1 + plot2
      ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
      scRNA <- ScaleData(object = scRNA, vars.to.regress = c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
      scRNA <- RunPCA(scRNA, verbose = T, npcs = npcs, features = VariableFeatures(object = scRNA))
      p <- ElbowPlot(object = scRNA, ndims = npcs)
      ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)
      
      if (Multi) {
        scRNA <- RunHarmony(object = scRNA, group.by.vars = "orig.ident", assay.use = assay, verbose = FALSE)
        scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:PC, verbose = T)
        # scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:PC, verbose = T)
        scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "harmony", verbose = T) # 使用harmony替代PCA
        scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
        p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
      } else {
        scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:PC, verbose = T)
        # scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:PC, verbose = T)
        scRNA <- FindNeighbors(scRNA, dims = 1:PC, reduction = "pca", verbose = T)
        scRNA <- FindClusters(scRNA, resolution = set.resolutions, verbose = T)
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
      }
      
      pdf(paste0(dir, "/7.data.merge.harmony.pdf"))
      # p <- clustree::clustree(scRNA)
      # print(p)
      p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = "orig.ident", raster = F) + NoLegend()
      print(p)
      merge.res <- sapply(set.resolutions, function(x) {
        p <- DimPlot(object = scRNA, reduction = "umap", label = TRUE, group.by = paste0(assay, "_snn_res.", x), raster = F) + NoLegend()
        print(p)
      })
      dev.off()
      saveRDS(scRNA, paste0(dir, "/data.merge.harmony.2000.rds"))
      
      Idents(scRNA) <- Idents
      scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
      scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05,]
      write.csv(scRNA.markers, paste0(dir, "/Findall.markers.RNA_snn_res.0.2.csv"))
    }
    
    # 自动化注释
    if (step == 3) {
      message("*************************自动化注释*************************")
      
      dir <- paste0(Path, "/.CellAnnotate")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      # p <- DimPlot(scRNA,group.by = Idents,label = T,reduction = "tsne")+NoLegend();p
      # ggsave_fun(filename = "2.Cluster/5.Tsne.Cluster",plot = p)
      
      Cluster.dir <- list.files(Path, pattern = "^.Cluster$", all.files = TRUE)
      
      p <- DimPlot(scRNA, group.by = Idents, label = T, reduction = "umap", raster = F) + NoLegend()
      print(p)
      ggsave_fun(filename = paste0(Path, "/", Cluster.dir, "/5.Umap.Cluster"), plot = p)
      
      CELL <- SingleR::SingleR(
        test = as.matrix(scRNA@assays$RNA@data),
        ref = ref,
        labels = labels,
        clusters = scRNA@meta.data[[Idents]]
      )
      
      cells <- CELL[, ncol(CELL):1]
      cells <- data.frame(cells)
      cell_id <- cells$pruned.labels
      names(cell_id) <- levels(scRNA@meta.data[[Idents]])
      Idents(scRNA) <- Idents
      
      if (anyNA(cell_id)) {
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(cell_id)
        p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
        print(p)
        message("某些细胞簇没注释到，是否强行为NA的细胞命名为unknow？yes/no")
        exit <- scan(what = "character", sep = ",", quiet = TRUE)
        
        if (exit == "yes") {
          cell_id[which(is.na(cell_id))] <- "unknow"
        }
        
        if (exit == "no") {
          message("某些细胞簇没注释到，具体请看右侧umap图，建议选择其他分辨率或参考数据集重新运行此函数，跳过1.质量控制和2.数据标准化，选择3重新注释。")
          return(scRNA)
        }
      }
      
      if (!anyNA(cell_id)) {
        p <- DimPlot(scRNA, group.by = Idents, raster = F,label = T)+NoLegend()
        print(p)
        print(paste0(paste0(names(cell_id),":",as.character(cell_id)), collapse = "; "))
        
        # a <- scan(what = "character")
        while (T) {
          message("请输入基因进行检验，跳过该步骤请按回车键：")
          a <- scan(what = "character")
          gene <- intersect(a,rownames(scRNA))
          if(length(gene) > 0){
            p <- FeaturePlot(scRNA,features = gene,raster = F,label = T,cols = c("lightgrey","red"))
            print(p)
          }
          if(length(a) == 0){break}
        }
        
        cat("\033[0;34m由于使用了多个参考数据，参考数据细胞名字不同，请按顺序对下述细胞进行重命名，英文;分隔。\n例：将\033[31mCD8+ T cells,Neutrophils,CD8+ T-cells\033[0;34m更改为\033[31mCD8+ T cells,Neutrophils,CD8+ T cells\033[0;34m以统一细胞名字\n")
        print(paste0(paste0(names(cell_id),":",as.character(cell_id)), collapse = "; "))
        cell_id <- scan(what = "character", sep = ";", quiet = FALSE)
        cell_id <- gsub(".*:","",cell_id)
        names(cell_id) <- levels(scRNA@meta.data[[Idents]])
        Idents(scRNA) <- Idents
        scRNA <- RenameIdents(scRNA, cell_id)
        scRNA$cellType_1 <- Idents(scRNA)
        print(table(scRNA$cellType_1))
        p <- DimPlot(scRNA, group.by = "cellType_1", raster = F)
        print(p)
        saveRDS(scRNA, paste0(Path, "/result.rds"))
        
        # Idents(scRNA) <- "cellType_1"
        # scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
        # scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05,]
        # write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))
        
        scRNA.markers %>%
          group_by(cluster) %>%
          top_n(n = 5, wt = avg_log2FC) -> features
        
        print((features))
        
        features <- features$gene
        cat("\033[0;34m是否要改写基因？yes/no\n")
        RenameGene <- scan(what = "character", sep = ",", quiet = TRUE)
        
        if (RenameGene == "yes") {
          print(paste0(as.character(features), collapse = ","))
          cat("\033[0;34m请改写\n")
          features <- scan(what = "character", sep = ",", quiet = FALSE)
        }
        
        # plot_list <- FeaturePlot(
        #   scRNA,
        #   features = unique(features),
        #   reduction = "umap",
        #   raster = F,
        #   combine = FALSE, cols = c("lightgrey", "red")
        # )
        #
        # for (i in 1:length(plot_list)) {
        #   plot_list[[i]] <- plot_list[[i]] + NoLegend() + NoAxes()
        # }
        #
        # p <- CombinePlots(plot_list)
        # # p
        # ggsave_fun(filename = paste0(dir, "/1.Featureplot.umap"), plot = p, width = 12, height = 12)
        
        # p <- DimPlot(scRNA, group.by = "cellType_1", label = T, label.size = 3, reduction = "umap", raster = F) + NoLegend()
        # p
        # ggsave_fun(filename = paste0(dir, "/2.Umap.plot"), plot = p, width = 7, height = 7)
        my_pal2 <- c(
          "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A"
        )
        
        umap <- scRNA@reductions$umap@cell.embeddings %>%
          as.data.frame() %>%
          cbind(cellType = scRNA@meta.data$cellType_1)
        
        ## 计算标签中心位置区域
        celltypepos <- umap %>%
          group_by(cellType) %>%
          summarise(
            umap_1 = median(UMAP_1),
            umap_2 = median(UMAP_2)
          )
        
        p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) +
          ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType),
                                    fontface = "bold",
                                    data = celltypepos,
                                    box.padding = 0.5, show.legend = FALSE
          )
        print(p)
        ggsave_fun(filename = paste0(dir, "/1.Umap1.plot"), plot = p, width = 10, height = 7)
        
        # add circle
        # p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T, cornerTextSize = 3.5, addCircle = TRUE)
        # print(p)
        # ggsave_fun(filename = paste0(dir, "/3.Umap2.plot"), plot = p, width = 10, height = 7)
        
        p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 10)) + theme(
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.position = "top",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.8, "cm"),
        )
        print(p)
        if (length(unique(features)) >= 50) {
          ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 20, height = 9)
        }
        if (length(unique(features)) < 50) {
          ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 17, height = 9)
        }
        
        tryCatch(
          {
            cellRatio <- cellRatioPlot(
              object = scRNA,
              sample.name = "group",
              celltype.name = "cellType_1",
              flow.curve = 0.5,
              fill.col = my_pal2
            ) + theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1
              )
            )
            print(cellRatio)
            
            ggsave_fun(filename = paste0(dir, "/3.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 7, height = 7)
            
            data_p <- scRNA@meta.data[, c("group", "cellType_1")]
            
            p <- ggstatsplot::ggbarstats(
              data = data_p,
              x = group,
              y = cellType_1,
              palette = "Set3"
            ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
            print(p)
            leg <- length(unique(data_p$cellType_1))
            ggsave_fun(filename = paste0(dir, "/4.Dysfunction"), plot = p, width = 2+1.5*leg, height = 7)
            
            EachSampleClusterDis <- lapply(unique(scRNA@meta.data$orig.ident), function(Sample) {
              seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$orig.ident %in% Sample, ]))
              seur@meta.data$cellType_1 %>%
                table() %>%
                data.frame() %>%
                magrittr::set_colnames(c("CellTypes", "Number")) %>%
                dplyr::mutate(Per = 100 * Number / sum(Number))
            })
            names(EachSampleClusterDis) <- unique(scRNA@meta.data$orig.ident)
            EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
            p1 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
              geom_col(position = "fill", width = 0.8) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Percentage of cell") +
              scale_fill_manual(values = my_pal2, name = "") +
              theme(
                axis.text = element_text(color = "black"),
                panel.background = element_blank(), # panel.grid=element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_text(color = "black"),
                axis.line = element_line(color = "black"),
                axis.title.y = element_blank(), # 轴标题
                legend.position = "right",
                # legend.key.height = unit(0.3, "cm"),
                # legend.key.width = unit(0.3, "cm"),
                legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
                legend.text = element_text(size = 10)
              ) +
              coord_flip()
            
            EachSampleClusterDis <- lapply(unique(scRNA@meta.data$group), function(Sample) {
              seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$group %in% Sample, ]))
              seur@meta.data$cellType_1 %>%
                table() %>%
                data.frame() %>%
                magrittr::set_colnames(c("CellTypes", "Number")) %>%
                dplyr::mutate(Per = 100 * Number / sum(Number))
            })
            
            my_pal1 <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#6387C5")
            
            names(EachSampleClusterDis) <- unique(scRNA@meta.data$group)
            EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
            p2 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
              geom_col(position = "fill", width = 0.8) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Percentage of cell") +
              scale_fill_manual(values = my_pal1, name = "") +
              theme(
                axis.text = element_text(color = "black", ),
                panel.background = element_blank(), panel.grid = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(), axis.title.x = element_blank(), # 删除轴文本
                axis.ticks.y.left = element_blank(),
                axis.title.y = element_blank(), # 轴标题
                axis.line.x = element_line(color = "black"),
                legend.position = "right",
                legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
                # legend.key.height = unit(0.3, "cm"),
                # legend.key.width = unit(0.3, "cm"),
                legend.text = element_text(size = 10)
              ) +
              coord_flip()
            
            # 计算每种细胞类型的细胞数量
            cell_counts <- table(scRNA@meta.data$cellType_1)
            
            # 创建数据框
            cell_counts_df <- data.frame(
              CellType = names(cell_counts),
              Count = as.numeric(cell_counts / 1000)
            )
            
            # 创建柱状图
            cell_counts_df <- cell_counts_df[order(cell_counts_df$Count, decreasing = T), ]
            
            p3 <- ggplot(cell_counts_df, aes(x = reorder(CellType, Count), y = Count)) +
              scale_y_continuous(expand = c(0, 0)) +
              labs(y = "Number of cell(10^3)") +
              geom_bar(stat = "identity", fill = "#6387C5") +
              theme(
                axis.text = element_text(color = "black"),
                panel.background = element_blank(), panel.grid = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(),
                # axis.title.x = element_blank(),#删除轴文本
                # plot.tag=element_blank(),
                axis.ticks.y.left = element_blank(),
                # axis.title.x.top = T,
                axis.title.y = element_blank(), # 轴标题
                axis.line.x = element_line(color = "black"),
                legend.text = element_text(size = 10)
              ) +
              coord_flip()
            print(p1 + p2 + p3)
            ggsave_fun(filename = paste0(dir, "/5.MergeRatio"), plot = p1 + p2 + p3, width = 15, height = 6)
            # if (isKeyCell == "yes" && is.na(KeyCell) && 1 %in% KeyCellIndex) {
            pB2_df <- table(scRNA@meta.data$cellType_1, scRNA@meta.data$group) %>% melt()
            colnames(pB2_df) <- c("Cluster", "Sample", "Number")
            pB2_df$Cluster <- factor(pB2_df$Cluster)
            
            pB2_df <- pB2_df %>%
              group_by(Sample) %>%
              mutate(Percentage = Number / sum(Number))
            
            pB2_df_tumor <- pB2_df %>%
              filter(Sample %in% c("Disease", "Tumor", "Cancer")) %>%
              select(Cluster, Tumor_Number = Number)
            
            pB2_df_control <- pB2_df %>%
              filter(Sample == c("Control", "Normal")) %>%
              select(Cluster, Control_Number = Number)
            
            pB2_df_combined <- merge(pB2_df_tumor, pB2_df_control, by = "Cluster")
            # pB2_df_combined <- pB2_df_combined[pB2_df_combined$Control_Number >= 50 && pB2_df_combined$Tumor_Number >= 50,]
            pB2_df_combined$rate <- (pB2_df_combined$Tumor_Number - pB2_df_combined$Control_Number) / pB2_df_combined$Control_Number
            KeyCell <- pB2_df_combined[abs(pB2_df_combined$rate) == max(abs(pB2_df_combined$rate)), ]$Cluster %>% as.character()
            print(paste0("变化最大的细胞：", KeyCell))
            write.csv(pB2_df_combined, paste0(dir, "/Cell.rate.csv"))
            # }
          },
          error = function(e) {
            message("没有分组信息，跳过此步骤\n", e)
            
            cellRatio <- cellRatioPlot(
              object = scRNA,
              sample.name = "orig.ident",
              celltype.name = "cellType_1",
              flow.curve = 0.5,
              fill.col = my_pal2
            ) + theme(
              axis.text.x = element_text(
                angle = 45,
                hjust = 1,
                vjust = 1
              )
            )
            print(cellRatio)
            
            ggsave_fun(filename = paste0(dir, "/3.Sample.ratio.plot"), plot = cellRatio, width = 7 + 0.25 * (length(unique(scRNA$orig.ident))), height = 7)
          }
        )
      }
    }
    
    # 细胞通讯
    if (step == 4) {
      message("*************************细胞通讯*************************")
      require(CellChat)
      
      dir <- paste0(Path, "/.Cellchat")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      options(stringsAsFactors = FALSE)
      cellchat <- createCellChat(object = scRNA, group.by = "cellType_1")
      cellchat
      groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
      
      if (mouse) {
        CellChatDB <- CellChatDB.mouse
      } else {
        CellChatDB <- CellChatDB.human
      }
      
      CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
      cellchat@DB <- CellChatDB.use
      
      cellchat <- subsetData(cellchat)
      # future::plan("multisession", workers = 4)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- projectData(cellchat, PPI.human)
      
      cellchat <- computeCommunProb(cellchat,seed.use = 123456)
      cellchat <- filterCommunication(cellchat)
      df.net <- subsetCommunication(cellchat)
      write.table(df.net, file = paste0(dir, "/net_lr.txt"), quote = F, sep = "\t", row.names = F)
      
      cellchat <- computeCommunProbPathway(cellchat)
      df.netp <- subsetCommunication(cellchat)
      write.table(df.netp, file = paste0(dir, "/net_pathway.txt"), quote = F, sep = "\t", row.names = F)
      
      ##### -----------show----------#####
      cellchat <- aggregateNet(cellchat)
      saveRDS(cellchat,paste0(dir,"/cellchat.rds"))
      groupSize <- as.numeric(table(cellchat@idents))
      
      # p1 <- netVisual_circle(cellchat@net$count,
      #   vertex.weight = groupSize, weight.scale = T,
      #   label.edge = F, title.name = "Number of interactions"
      # )
      # 
      # p2 <- netVisual_circle(cellchat@net$weight,
      #   vertex.weight = groupSize, weight.scale = T,
      #   label.edge = F, title.name = "Interaction weigths/strength"
      # )
      # # pdf(file="1.Net_number_strength.pdf",width = 10,height = 5)
      # p <- CombinePlots(plots = list(p1, p2), ncol = 2)
      # ggsave_fun(paste0(dir, "/1.Net_number_strength"), plot = p, width = 14)
      pdf(file=paste0(dir, "/1.Net_number_strength.pdf"),width = 10,height = 5)
      par(mfrow = c(1,2), xpd=TRUE)
      netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
      dev.off()
      
      df.net <- read.table(paste0(dir, "/net_lr.txt"), sep = "\t", check.names = F, header = T)
      data <- as.data.frame(table(c(
        df.net$source,
        df.net$target
      )))
      colnames(data) <- c("Cell_Type", "all_sum")
      data <- data[order(data$all_sum, decreasing = T), ]
      data$Cell_Type <- factor(data$Cell_Type, levels = data$Cell_Type)
      head(data)
      sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
      sample_color2 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], ifelse(data$all_sum > summary(data$all_sum)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
      sample_color1 <- ifelse(data$all_sum > summary(data$all_sum)[2], ifelse(data$all_sum > summary(data$all_sum)[3], ifelse(data$all_sum > summary(data$all_sum)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))
      
      pB2 <- ggplot(data = data, aes(x = Cell_Type, y = all_sum, fill = Cell_Type, colour = Cell_Type)) +
        geom_bar(stat = "identity", width = 0.8) +
        scale_fill_manual(values = sample_color1) +
        scale_colour_manual(values = sample_color2) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(x = "Cell Type", y = "Counts") +
        theme(axis.text.y = element_text(size = 12, colour = "black")) +
        theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
      
      ggsave_fun(paste0(dir, "/2.Interaction Count"), plot = pB2, width = 10, height = 8)
      
      # if (isKeyCell == "yes" && is.na(KeyCell) && 2 %in% KeyCellIndex) {
      #   KeyCell <- as.character(data[1, 1])
      # }
    }
    
    # 亚群贡献
    if (step == 5) {
      message("*************************亚群贡献*************************")
      
      dir <- paste0(Path, "/.Subpopulation contribution")
      if (!dir.exists(paste0(dir))) {
        dir.create(paste0(dir))
      } else {
        print("Dir already exists!")
      }
      
      scRNA$group <- factor(scRNA$group, c(grep("Control|Normal", unique(scRNA$group), value = T), grep("Disease|Tumor|Cancer", unique(scRNA$group), value = T)))
      print(names(table(scRNA$group)))
      if (names(table(scRNA$group))[2] %in% c("Disease","Tumor","Cancer") & names(table(scRNA$group))[1]  %in% c("Control","Normal")){
        Bulk.DEGs <- FindMarkers(scRNA,
                                 group.by = "group", min.pct = 0.25, logfc.threshold = logFCfilter,
                                 ident.1 = names(table(scRNA$group))[2], ident.2 = names(table(scRNA$group))[1]
        )
        Bulk.DEGs <- Bulk.DEGs[Bulk.DEGs$p_val_adj < adjPvalFilter, ]
        Bulk.DEGs <- Bulk.DEGs[Bulk.DEGs$avg_log2FC > 0,]
        print(nrow(Bulk.DEGs))
        Bulk.DEGs$symbol <- rownames(Bulk.DEGs)
        # Bulk.DEGs <- arrange(Bulk.DEGs, Bulk.DEGs$p_val_adj, -Bulk.DEGs$avg_log2FC)
        # Bulk.DEGs <- Bulk.DEGs[1:100, ]
        subset.DEGs <- lapply(SplitObject(scRNA, split.by = "cellType_1"), function(subset) {
          deg <- FindMarkers(subset,
                             features = Bulk.DEGs$symbol, min.pct = 0, logfc.threshold = 0,
                             group.by = "group", ident.1 = names(table(scRNA$group))[2], ident.2 = names(table(scRNA$group))[1]
          )
          deg$symbol <- rownames(deg)
          deg$celltype <- unique(subset$cellType_1)
          return(deg)
        })
        subset.DEGs <- do.call(rbind, subset.DEGs)
        subset.DEGs$FCexp <- 2^subset.DEGs$avg_log2FC
        subset.DEGs$FCprop <- subset.DEGs$pct.1 / subset.DEGs$pct.2
        subset.DEGs$FCscore <- sqrt(subset.DEGs$FCexp * subset.DEGs$FCprop)
        subset.DEGs$FCscore[is.infinite(subset.DEGs$FCscore)] <- NA
        FCscore <- dcast(subset.DEGs, symbol ~ celltype, measure.var = "FCscore")
        write.table(FCscore, paste0(dir, "/output_FCscore.txt"),
                    sep = "\t", row.names = F, col.names = T, quote = F
        )
        
        
        plot.data <- data.frame(
          "subtypes" = colnames(FCscore[, -1]),
          "FCscore" = colMeans(FCscore[, -1], na.rm = T)
        )
        plot.data <- arrange(plot.data, plot.data$FCscore)
        plot.data$subtypes <- factor(plot.data$subtypes, levels = plot.data$subtypes)
        my_pal2 <- c(
          "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
          "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00"
        )
        p <- ggplot(plot.data, aes(x = subtypes, y = FCscore, fill = subtypes)) +
          geom_bar(stat = "identity") +
          scale_fill_manual(values = my_pal2[1:nrow(plot.data)]) +
          geom_hline(
            yintercept = median(plot.data$FCscore),
            color = "grey", linetype = "dashed"
          ) +
          coord_polar() +
          theme_classic() +
          theme(
            axis.text = element_blank(), axis.title = element_blank(),
            axis.line = element_blank(), axis.ticks = element_blank()
          )
        # p
        # ggsave_fun(filename = paste0(dir, "/1.ContributionScore.pdf"), plot = p, width = 6, height = 4)
        
        plot.data <- plot.data[order(plot.data$FCscore, decreasing = T), ]
        plot.data$subtypes <- factor(plot.data$subtypes, levels = plot.data$subtypes)
        sample_color <- c("#FB040B", "#F6A717", "#BA06FA", "#172D7A")
        sample_color2 <- ifelse(plot.data$FCscore > summary(plot.data$FCscore)[2], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[3], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[5], sample_color[1], sample_color[2]), sample_color[3]), sample_color[4])
        sample_color1 <- ifelse(plot.data$FCscore > summary(plot.data$FCscore)[2], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[3], ifelse(plot.data$FCscore > summary(plot.data$FCscore)[5], alpha(sample_color[1], 0.9), alpha(sample_color[2], 0.9)), alpha(sample_color[3], 0.9)), alpha(sample_color[4], 0.9))
        
        pB2 <- ggplot(data = plot.data, aes(x = subtypes, y = FCscore, fill = subtypes, colour = subtypes)) +
          geom_bar(stat = "identity", width = 0.8) +
          scale_fill_manual(values = sample_color1) +
          scale_colour_manual(values = sample_color2) +
          theme_bw() +
          theme(panel.grid = element_blank()) +
          labs(x = "Cell Type", y = "Counts") +
          theme(axis.text.y = element_text(size = 12, colour = "black")) +
          theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, colour = "black"))
        ggsave_fun(filename = paste0(dir, "/1.CombinedContribution"), plot = p + pB2, width = 12, height = 6)
        
        write.csv(plot.data, paste0(dir, "/plot.data.csv"), row.names = F)
      }
      # if (isKeyCell == "yes" && is.na(KeyCell) && 3 %in% KeyCellIndex) {
      #   KeyCell <- as.character(plot.data[1, 1])
      # }
    }
    
    # 拟时序分析
    if (step == 6) {
      message("*************************Pseudo-time*************************")
      require(monocle)
      
      dir <- paste0(Path, "/.Pseudo-time")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      scRNASub <- subset(scRNA, subset = cellType_1 == KeyCell)
      
      monocle.matrix <- GetAssayData(object = scRNASub, slot = "data", assay = "RNA")
      monocle.sample <- scRNASub@meta.data
      monocle.geneAnn <- data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
      
      data <- as(as.matrix(monocle.matrix), "sparseMatrix")
      pd <- new("AnnotatedDataFrame", data = monocle.sample)
      fd <- new("AnnotatedDataFrame", data = monocle.geneAnn)
      cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
      if ("group" %in% names(pData(cds))) {
        names(pData(cds))[names(pData(cds)) == "group"] <- "Cluster"
      } else {
        names(pData(cds))[names(pData(cds)) == "orig.ident"] <- "Cluster"
      }
      pData(cds)[, "Cluster"] <- paste0("cluster", pData(cds)[, "Cluster"])
      
      cds <- estimateSizeFactors(cds)
      cds <- estimateDispersions(cds)
      
      disp_table <- dispersionTable(cds)
      disp_table <- arrange(disp_table, -dispersion_empirical)
      track_gene <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
      
      cds <- setOrderingFilter(cds, track_gene)
      cds <- reduceDimension(cds, max_components = 2, reduction_method = "DDRTree")
      cds <- orderCells(cds)
      
      p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
      ggsave(paste0(dir, "/1.Trajectory.Pseudotime.pdf"), plot = p)
      
      p <- plot_cell_trajectory(cds, color_by = "State")
      ggsave(paste0(dir, "/2.Trajectory.State.pdf"), plot = p)
      
      p <- plot_cell_trajectory(cds, color_by = "Cluster")
      ggsave(paste0(dir, "/3.Trajectory.Cluster.pdf"), plot = p)
      
      save(cds, file = paste0(dir, "/cds.rda"))
      
      pseudo_time_diff <- differentialGeneTest(cds, cores = 4, fullModelFormulaStr = "~sm.ns(Pseudotime)")
      pseudo_time_diff <- arrange(pseudo_time_diff, qval)
      write.csv(pseudo_time_diff,file = paste0(dir,'/heatmap_pseudo_time_gene.csv'))
      pseudo_time_genes <- rownames(pseudo_time_diff)
      plot_heatmap_gene <- pseudo_time_genes[1:50]
      len <- length(plot_heatmap_gene)
      p <- monocle::plot_pseudotime_heatmap(cds[plot_heatmap_gene, ], num_clusters = 3, show_rownames = T, return_heatmap = T, hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100))
      ggsave(filename = paste0(dir, "/4.Heatmap.pdf"), plot = ggplotify::as.ggplot(p), width = 5, height = len * 0.12)
      if (!is.null(dev.list())) {
        dev.off()
      }
      
      # clusters <- cutree(p$tree_row, k = 3)
      # clustering <- data.frame(clusters)
      # clustering[,1] <- as.character(clustering[,1])
      # colnames(clustering) <- "Gene_Clusters"
      # clustering <- arrange(clustering,Gene_Clusters)
      # write.csv(clustering,file = paste0(dir,'/','heatmap_cluster_gene.csv'))
      
      if (!class(KeyGene) == "logical") {
        cds_subset <- cds[KeyGene, ]
        p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "Cluster")
        p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
        p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
        p1 | p2 | p3
        ggsave(filename = paste0(dir, "/5.KeyGene.pseudotime.pdf"), plot = p1 | p2 | p3, width = 12, height = length(KeyGene) * 2)
      }
    }
    
    # 关键基因表达丰度
    # 7
    if (step == 7) {
      message("*************************关键基因表达丰度*************************")
      
      dir <- file.path(Path, ".KeyGene expression abundance")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      # 获取在 scRNA 中存在的 KeyGene
      keyGenes_in_data <- intersect(rownames(scRNA), KeyGene)
      leg <- length(unique(keyGenes_in_data))
      
      # 使用 FeaturePlot + patchwork 替代 FeatureCornerAxes
      if (length(keyGenes_in_data) > 0) {
        library(Seurat)
        library(patchwork)
        
        features <- unique(keyGenes_in_data)
        
        # 每个基因一个 FeaturePlot，组合成网格
        plt_list <- lapply(features, function(g) {
          FeaturePlot(scRNA, features = g, reduction = "umap", pt.size = 0.1,raster=FALSE)
        })
        
        p <- wrap_plots(plt_list, ncol = min(3, length(features)))
        
        # 根据基因数量调整宽高，与原逻辑尽量保持一致
        if (length(features) <= 3) {
          w <- length(features) * 5
          h <- 5
        } else {
          w <- 15
          h <- ceiling(length(features) / 3) * 5
        }
        
        ggsave_fun(filename = paste0(dir, "/1.Keygene.FeaturePlot"), plot = p, width = w, height = h)
      } else {
        message("No key genes found in the scRNA object.")
      }
      
      # DotPlot 部分，删除尾随逗号
      p <- DotPlot(scRNA, group.by = "cellType_1", features = unique(keyGenes_in_data), cols = "RdYlBu") +
        scale_size_continuous(range = c(0, 10)) + theme(
          panel.border = element_rect(colour = "black"),
          axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          ),
          legend.position = "right",
          text = element_text(size = 10)
        )
      
      if (leg <= 3) {
        ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p,
                   width = 5 + ceiling(leg / 3), height = 7)
      } else {
        ggsave_fun(filename = paste0(dir, "/2.Keygene.DotPlot"), plot = p,
                   width = 5 + 0.5 * leg, height = 7)
      }
    }
    
    # 疾病基因共表达
    if (step == 8) {
      message("*************************疾病基因共表达网络*************************")
      
      dir <- paste0(Path, "/.Disease gene co-expression")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      showGenes <- intersect(rownames(scRNA), KeyGene)
      cat("\033[0;34m请复制疾病相关基因：\n")
      geneCard <- scan(sep = "\n", what = "character", quiet = FALSE)
      geneCard <- intersect(geneCard,rownames(scRNA))[1:6]
      
      for (i in showGenes) {
        PlotNum <- 0
        path <- paste0(dir, "/", which(showGenes == i), ".", i)
        if (!dir.exists(path)) {
          dir.create(path)
        } else {
          print("Dir already exists!")
        }
        for (j in geneCard) {
          PlotNum <- PlotNum + 1
          p1 <- FeaturePlot(scRNA,
                            features = c(i, j),
                            blend = TRUE, cols = c("gray80", "red", "green"),
                            pt.size = 0.5, raster = F, reduction = "umap"
          ) +
            theme(aspect.ratio = 1)
          p2 <- getScatterplot(scRNA,
                               gene1 = j, gene2 = i,
                               jitter.num = 0.15, pos = TRUE
          ) +
            theme(aspect.ratio = 1)
          p <- CombinePlots(plots = list(p1, p2), ncol = 2, rel_widths = c(4, 1))
          ggsave_fun(paste0(path, "/", PlotNum, ".", j, " ~ ", i), plot = p, width = 15, height = 4)
        }
      }
    }
    
    # 免疫代谢通路
    if (step == 9) {
      message("*************************免疫代谢通路*************************")
      
      dir <- paste0(Path, "/.Immunometabolic pathways")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      ShowGen <- intersect(rownames(scRNA), KeyGene)
      Idents(scRNA) <- "cellType_1"
      geneset <- clusterProfiler::read.gmt(file.path(outputfile,"refdata/Homo/h.all.v7.5.1.symbols.gmt"))
      if (mouse) {
        geneset <- clusterProfiler::read.gmt("refdata/Mus/mh.all.v2023.2.Mm.symbols.gmt")
      }
      geneset <- split(geneset$gene, geneset$term)
      genesetInfo <- read.delim(file.path(outputfile,"refdata/GenesetInfo.txt"), sep = ",")
      genesetInfo <- subset(genesetInfo, Classification != "")
      levels <- c("Immune", "Metabolism", "Signaling", "Proliferation")
      genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
      genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
      genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
      mat <- GetAssayData(object = scRNA, assay = "RNA", slot = "data")#slot = "data"V4使用layerV5使用
      # rm(scRNA.plot2)
      cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 1, plotStats = F)
      score <- AUCell::AUCell_calcAUC(geneset, cells_rankings,
                                      nCores = 1,
                                      aucMaxRank = nrow(cells_rankings) * 0.05
      )
      score <- AUCell::getAUC(score) # score@assays@data$AUC
      
      # compare <<- "Hexp-Lexp"
      # adjust.method <<- "bonferroni"
      DP.list <- base::lapply(ShowGen, function(cellType) {
        # cellType <- ShowGen[1]
        # print(cellType)
        # print(1)
        DatGroup <- FetchData(scRNA, vars = c(cellType), slot = "data")#slot = "data"V4使用
        # print(2)
        group <- setNames(
          object = ifelse(DatGroup[, 1] > median(DatGroup[, 1]), "Hexp", "Lexp"),
          nm = rownames(DatGroup)
        )
        # print(3)
        if (length(unique(group)) > 1) {
          # print(4)
          design <- model.matrix(~ 0 + factor(group))
          # print(5)
          colnames(design) <- levels(factor(group))
          # print(6)
          rownames(design) <- names(group)
          # print(7)
          contrast.matrix <- limma::makeContrasts("Hexp-Lexp", levels = design) # should be Test-Control
          # print(8)
          fit <- limma::lmFit(score[, names(group)], design)
          # print(9)
          fit2 <- limma::contrasts.fit(fit, contrast.matrix)
          # print(10)
          fit2 <- limma::eBayes(fit2)
          # print(11)
          DPs <- limma::topTable(fit2, coef = 1, n = Inf, adjust.method = "bonferroni")
          # print(12)
          DPs$KeyGenes <- cellType
          # print(13)
          DPs$Pathway <- rownames(DPs)
          return(DPs)
        }
      })
      DPs <- do.call(rbind, DP.list)
      write.table(DPs, file = paste0(dir, "/output_hallmark.txt"), sep = "\t", row.names = F, col.names = T, quote = F)
      
      plot.data <- DPs
      plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
      plot.data$KeyGenes <- factor(plot.data$KeyGenes)
      plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset)
      plot.data$FDR <- cut(plot.data$adj.P.Val,
                           breaks = c(0, 1e-125, 1e-75, 1e-25, 1),
                           include.lowest = T
      )
      plot.data$FDR <- factor(as.character(plot.data$FDR),
                              levels = rev(levels(plot.data$FDR))
      )
      levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))
      
      color <- c("#4682B4", "#FFFFFF", "#CD2626")
      class.color <- c(
        "Immune" = "#D58986", "Metabolism" = "#80554C",
        "Signaling" = "#71AC7A", "Proliferation" = "#E8D4B4"
      )
      p1 <- ggplot(plot.data, aes(x = Pathway, y = KeyGenes, color = logFC, size = FDR)) +
        geom_point() +
        scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
        geom_hline(
          yintercept = seq(
            min(as.numeric(plot.data$KeyGenes)) - 0.5,
            max(as.numeric(plot.data$KeyGenes)) + 0.5
          ),
          color = "grey80"
        ) +
        geom_vline(
          xintercept = seq(
            min(as.numeric(plot.data$Pathway)) - 0.5,
            max(as.numeric(plot.data$Pathway)) + 0.5
          ),
          color = "grey80"
        ) +
        theme_classic() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.line = element_blank(),
          legend.position = "top",
          legend.key.height = unit(0.3, "cm"),
          legend.key.width = unit(0.8, "cm"),
        )
      
      p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) +
        geom_tile() +
        theme_classic() +
        scale_fill_manual(values = class.color) +
        theme(
          axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "bottom"
        )
      p <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(10, 2))
      ggsave_fun(filename = paste0(dir, "/1.Hallmark"), width = 15, height = 9, plot = p)
    }
    
    # 代谢通路热图
    if (step == 10) {
      message("*************************代谢通路热图*************************")
      
      require(GSVA)
      require(GSEABase)
      require(limma)
      require(pheatmap)
      
      # === 参数设置 ===
      dir <- paste0(Path, "/.Metabolic pathway")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        cat("\033[0;33mDir already exists!\033[0m\n")
      }
      
      gmtFile <- "refdata/Homo/immune.gmt"
      cat("\033[0;34m请输入分组名字，例：group\n")
      Metagroup <- scan(what = "character", n = 1)
      
      # === 准备表达矩阵 ===
      rt <- scRNA@assays$RNA@data
      group <- data.frame(
        ID = colnames(scRNA), 
        Group = scRNA@meta.data[[Metagroup]],
        stringsAsFactors = FALSE
      )
      
      exp <- as.matrix(rt)
      dimnames <- list(rownames(exp), colnames(exp))
      mat <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
      mat <- limma::avereps(mat)
      mat <- mat[rowMeans(mat) > 0, ]
      
      # === 加载基因集 ===
      geneSet <- getGmt(gmtFile, geneIdType = SymbolIdentifier())
      cat("[1] 'Gene set loaded:', length(geneSet), 'pathways\n")
      
      # === GSVA分析 ===
      param <- ssgseaParam(mat, geneSet, 
                           normalize = TRUE,
                           minSize = 2,
                           alpha = 0.25)
      
      ssgseaScore <- gsva(param, verbose = FALSE)
      cat("[1] 'GSVA calculation completed'\n")
      
      # === 归一化函数 ===
      normalize <- function(x) {
        x <- as.numeric(x)
        min_val <- min(x, na.rm = TRUE)
        max_val <- max(x, na.rm = TRUE)
        
        if (is.infinite(min_val) || is.infinite(max_val) || min_val == max_val) {
          return(rep(0.5, length(x)))
        }
        return((x - min_val) / (max_val - min_val))
      }
      
      # 按行进行归一化
      ssgseaOut_norm <- t(apply(ssgseaScore, 1, normalize))
      
      # 确保行名和列名正确
      rownames(ssgseaOut_norm) <- rownames(ssgseaScore)
      colnames(ssgseaOut_norm) <- colnames(ssgseaScore)
      
      # 保存为标准格式
      write.table(ssgseaOut_norm, 
                  file = paste0(dir, "/ssgseaOut.txt"), 
                  sep = "\t", 
                  quote = FALSE, 
                  row.names = TRUE,      # 保留通路名
                  col.names = NA)      # 保留样本ID
      
      cat("[1] 'Normalization and save completed'\n")
      
      # === 直接使用结果（不从文件重读） ===
      plotdata <- ssgseaOut_norm  # 直接用内存中的矩阵
      
      # 检查列名是否匹配
      cat("[DEBUG] plotdata列数:", ncol(plotdata), "\n")
      cat("[DEBUG] group样本数:", nrow(group), "\n")
      cat("[DEBUG] group$ID前5个:", paste(head(group$ID, 5), collapse=", "), "\n")
      cat("[DEBUG] plotdata列名前5个:", paste(colnames(plotdata)[1:min(5, ncol(plotdata))], collapse=", "), "\n")
      
      # 按分组排序样本
      up <- group[order(group$Group), ]
      
      # 确保列名匹配
      common_samples <- intersect(colnames(plotdata), up$ID)
      if (length(common_samples) == 0) {
        cat("\033[0;31m[ERROR] 样本ID不匹配！\033[0m\n")
        cat("plotdata列名:", colnames(plotdata), "\n")
        cat("group$ID:", up$ID, "\n")
        stop("样本ID无法匹配")
      }
      
      plotdata <- plotdata[, common_samples, drop = FALSE]
      up <- up[up$ID %in% common_samples, ]
      up <- up[order(match(up$ID, colnames(plotdata))), ]
      plotdata <- plotdata[, up$ID, drop = FALSE]
      
      up$group <- up$Group
      
      # 读取通路注释
      down <- read.table("refdata/Homo/type.txt", sep = "\t", 
                         header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
      
      # 取交集
      compath <- intersect(rownames(plotdata), down$pathway)
      down <- down[down$pathway %in% compath, , drop = FALSE]
      plotdata <- plotdata[down$pathway, , drop = FALSE]
      
      cat("[1] 'Pathway count:', nrow(plotdata), '\n")
      cat("[1] 'Sample count:', ncol(plotdata), '\n")
      
      # === 准备注释信息 ===
      annCol <- data.frame(
        Group = up$group,
        row.names = up$ID,
        stringsAsFactors = FALSE
      )
      
      annRow <- data.frame(
        Direct = down$type,
        row.names = down$pathway,
        stringsAsFactors = FALSE
      )
      
      # 获取唯一分组并设置颜色
      unique_groups <- unique(group$Group)
      annGroup <- c("#E7B800", "#2E9FDF", "#FF6B6B", "#4ECDC4")[1:length(unique_groups)]
      names(annGroup) <- unique_groups
      
      # 获取唯一类型并设置颜色
      unique_types <- unique(down$type)
      typeColors <- c(
        "1Amino acid metabolism relevant signatures" = "#247BA0",
        "2lipid metabolism relevant signatures" = "#70C1B3",
        "3Drug metabolism relevant signatures" = "#B2DBBF",
        "4Other metabolism signatures" = "#F3FFBD",
        "5C3 specific metabolism signatures" = "#6CD3A7"
      )
      
      typeColors <- typeColors[names(typeColors) %in% unique_types]
      
      annColors <- list(
        Group = annGroup,
        Direct = typeColors
      )
      
      # === 标准化绘图数据 ===
      plotdata_scaled <- t(scale(t(plotdata)))
      plotdata_scaled[plotdata_scaled > 1] <- 1
      plotdata_scaled[plotdata_scaled < -1] <- -1
      
      # === 绘制热图 ===
      pdf(file = paste0(dir, "/1.Heatmap.pdf"), width = 25, height = 16)
      
      tryCatch({
        # 创建基于分组的列聚类函数
        col_order <- order(annCol$Group)
        
        pheatmap(plotdata_scaled[, col_order],
                 scale = "none",
                 annotation_row = annRow,
                 annotation_col = annCol[col_order, , drop = FALSE],
                 annotation_colors = annColors,
                 color = colorRampPalette(c("#009BC7", "#F3F3F1", "#F15E4C"))(20),
                 fontsize_row = 12,
                 fontsize_col = 8,
                 fontsize = 12,
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 show_colnames = FALSE,
                 border_color = NA)
        cat("\033[0;32m[1] 'Heatmap generated successfully'\033[0m\n")
      }, error = function(e) {
        cat("\033[0;31m[ERROR] Heatmap generation failed:", conditionMessage(e), "\033[0m\n")
      })
      
      dev.off()
      
      # === 两组差异通路分析 ===
      tryCatch({
        cat("\033[0;34m进行两组差异通路分析...\033[0m\n")
        
        unique_groups <- sort(unique(up$Group))
        
        if (length(unique_groups) >= 2) {
          # 获取前两个分组
          group1_name <- unique_groups[1]
          group2_name <- unique_groups[2]
          
          cat("[INFO] 比较分组: ", group1_name, " vs ", group2_name, "\n")
          
          # 获取两个分组的样本
          group1_samples <- up$ID[up$Group == group1_name]
          group2_samples <- up$ID[up$Group == group2_name]
          
          if (length(group1_samples) == 0 || length(group2_samples) == 0) {
            stop("一个或两个分组没有样本")
          }
          
          # 提取两组的数据
          data_group1 <- plotdata[, colnames(plotdata) %in% group1_samples, drop = FALSE]
          data_group2 <- plotdata[, colnames(plotdata) %in% group2_samples, drop = FALSE]
          
          # 计算两组的平均值
          mean_group1 <- rowMeans(data_group1)
          mean_group2 <- rowMeans(data_group2)
          
          # 计算fold change
          fc <- mean_group2 - mean_group1
          
          # 使用t检验计算p值
          pvalue <- sapply(1:nrow(plotdata), function(i) {
            group1_vals <- as.numeric(plotdata[i, group1_samples])
            group2_vals <- as.numeric(plotdata[i, group2_samples])
            
            if (length(group1_vals) > 1 && length(group2_vals) > 1) {
              tryCatch({
                t.test(group2_vals, group1_vals)$p.value
              }, error = function(e) { NA })
            } else {
              NA
            }
          })
          
          # 调整p值
          padj <- p.adjust(pvalue, method = "BH")
          
          # 创建差异通路表格
          diffPath <- data.frame(
            Pathway = rownames(plotdata),
            Type = down$type[match(rownames(plotdata), down$pathway)],
            Mean_Group1 = mean_group1,
            Mean_Group2 = mean_group2,
            FC = fc,
            Pvalue = pvalue,
            Padj = padj,
            stringsAsFactors = FALSE
          )
          
          # 按FC绝对值排序
          diffPath <- diffPath[order(abs(diffPath$FC), decreasing = TRUE), ]
          rownames(diffPath) <- NULL
          
          # 保存差异通路表格
          write.table(diffPath, 
                      file = paste0(dir, "/2.Differential_Pathway.txt"), 
                      sep = "\t", 
                      quote = FALSE, 
                      row.names = FALSE)
          
          cat("[1] 'Differential pathway saved'\n")
          
          # setwd("<PATH>/000project/阿兹海默症/Rscript")
          # dir <- "."
          # diffPath <- read.table( 
          #                   file = paste0(dir, "/.Metabolic pathway/2.Differential_Pathway.txt"), 
          #                  sep = "\t",header=T)
          
          # === 绘制差异通路热图 ===
          sig_pathways <- diffPath$Pathway[!is.na(diffPath$Pvalue) & diffPath$Pvalue < 0.05]#可调P值
          
          if (length(sig_pathways) > 0) {
            valid_pathways <- intersect(sig_pathways, down$pathway)
            
            if (length(valid_pathways) > 0) {
              # 数据准备
              plotdata_diff <- plotdata[valid_pathways, , drop = FALSE]
              plotdata_diff_scaled <- t(scale(t(plotdata_diff)))
              
              # 数据清理
              plotdata_diff_scaled[is.nan(plotdata_diff_scaled)] <- 0
              plotdata_diff_scaled[is.infinite(plotdata_diff_scaled)] <- 0
              plotdata_diff_scaled[plotdata_diff_scaled > 2] <- 2
              plotdata_diff_scaled[plotdata_diff_scaled < -2] <- -2
              
              # 列排序
              col_order_diff <- order(annCol[colnames(plotdata_diff_scaled), 1, drop = TRUE])
              ordered_cols <- colnames(plotdata_diff_scaled)[col_order_diff]
              plotdata_diff_scaled <- plotdata_diff_scaled[, ordered_cols, drop = FALSE]
              
              # 注释准备
              annCol_diff <- data.frame(
                Group = annCol[ordered_cols, 1, drop = TRUE],
                row.names = ordered_cols,
                stringsAsFactors = FALSE
              )
              
              annRow_diff <- data.frame(
                Type = down$type[match(valid_pathways, down$pathway)],
                row.names = valid_pathways,
                stringsAsFactors = FALSE
              )
              
              # 生成色彩
              color_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(20)
              color_breaks <- seq(-2, 2, length.out = 21)
              
              pdf(file = paste0(dir, "/2.Differential_Heatmap.pdf"), width = 22, height = 14)
              
              pheatmap(plotdata_diff_scaled,
                       scale = "none",
                       annotation_row = annRow_diff,
                       annotation_col = annCol_diff,
                       annotation_colors = annColors,
                       color = color_palette,
                       breaks = color_breaks,
                       fontsize_row = 11,
                       fontsize_col = 8,
                       fontsize = 12,
                       cluster_cols = FALSE,
                       cluster_rows = FALSE,
                       show_colnames = FALSE,
                       border_color = NA,
                       main = paste0("Differential Pathways: ", group1_name, " vs ", group2_name))
              
              dev.off()
              
              cat("\033[0;32m[1] 'Differential heatmap generated successfully'\033[0m\n")
              cat("\033[0;32m差异通路热图已保存至:", paste0(dir, "/2.Differential_Heatmap.pdf"), "\033[0m\n")
            } else {
              cat("\033[0;33m[WARNING] 没有有效的通路\033[0m\n")
            }
          } else {
            cat("\033[0;33m[WARNING] 无显著差异通路 (padj < 0.05)\033[0m\n")
          }
        } else {
          cat("\033[0;33m[WARNING] 分组少于2个\033[0m\n")
        }
      }, error = function(e) {
        cat("\033[0;31m[ERROR] 差异分析失败:", conditionMessage(e), "\033[0m\n")
      })
      
      
      cat("\033[0;32m热图已保存至:", paste0(dir, "/1.Heatmap.pdf"), "\033[0m\n")
      
    }  
    
    # 量化热点机制
    if (step == 11) {
      message("*************************量化热点机制*************************")
      
      require(GSVA)
      require(GSEABase)
      require(limma)
      
      dir <- paste0(Path, "/.Hotspot mechanism")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      normalize <- function(x) {
        return((x - min(x)) / (max(x) - min(x)))
      }
      mat <- GetAssayData(object = scRNA, assay = "RNA")
      mat <- as(mat[rowMeans(mat) > 0.3, ], "matrix")
      geneSet <- GSEABase::getGmt("ssGSEA.gmt", geneIdType = GSEABase::SymbolIdentifier())
      
      ssgseaScore <- gsva(mat, geneSet, method = "ssgsea", parallel.sz = 10, kcdf = "Gaussian", abs.ranking = TRUE)
      ssgseaOut <- normalize(ssgseaScore)
      save(ssgseaOut, file = paste0(dir, "/ssgseaOut.rda")) # load("ssgseaOut.rda")
      ssgseaOut <- rbind(id = colnames(ssgseaOut), ssgseaOut)
      write.table(ssgseaOut, file = paste0(dir, "/ssgseaOutAll.txt"), sep = "\t", quote = F, col.names = F)
      
      load(paste0(dir, "/ssgseaOut.rda"))
      ssgseaOut <- as.data.frame(t(ssgseaOut))
      ScoGroup <- ifelse(ssgseaOut[rownames(scRNA[["orig.ident"]]), names(geneSet)] > median(ssgseaOut[rownames(scRNA[["orig.ident"]]), names(geneSet)]), "Hsco", "Lsco")
      scRNA$ScoGroup <- data.frame(row.names = rownames(scRNA[["orig.ident"]]), ScoGroup)
      
      load(paste0(dir, "/ssgseaOut.rda"))
      meta <- as.data.frame(t(ssgseaOut))
      scRNA <- AddMetaData(scRNA, meta)
      saveRDS(scRNA, paste0(dir, "/result.rds"))
      
      tryCatch(
        {
          for (i in rownames(ssgseaOut)) {
            tmp <- scRNA@meta.data[, c("cellType_1", "group", i)]
            colnames(tmp)[3] <- "val"
            p1 <- ggplot(
              tmp,
              aes(
                x = cellType_1, y = val,
                fill = group,
                color = group
              )
            ) +
              geom_boxplot(
                notch = F, alpha = 0.95,
                outlier.shape = 16,
                outlier.size = 0.65
              ) +
              xlab("") +
              ylab("") +
              scale_fill_manual(values = c("#D5EBFB", "#FBEEB7", "#B4FBCD", "#F5B3FC")) +
              scale_color_manual(values = c("#0073C2", "#EFC000", "#00C244", "#C501D7")) +
              ggtitle("") +
              theme_classic() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(vjust = 0.5, size = 12),
                axis.title.y = element_text(angle = 90, size = 15)
              ) +
              theme(legend.position = "top") +
              stat_compare_means(method = "wilcox.test", hjust = 0.5, vjust = 0, hide.ns = T, label = "p.signif", show.legend = F)
            p2 <- FeaturePlot(object = scRNA, ncol = 1, features = i, cols = c("green", "red"))
            p <- CombinePlots(plots = list(p2, p1), rel_widths = c(1, 1.6))
            ggsave_fun(paste0(dir, "/1.", i), width = 12, height = 4, plot = p)
          }
        },
        error = function(e) {
          for (i in rownames(ssgseaOut)) {
            tmp <- scRNA@meta.data[, c("cellType_1", i)]
            colnames(tmp)[2] <- "val"
            p1 <- ggplot(
              tmp,
              aes(
                x = cellType_1, y = val,
                fill = cellType_1
              )
            ) +
              geom_boxplot(
                notch = F, alpha = 0.95,
                outlier.shape = 16,
                outlier.size = 0.65
              ) +
              xlab("") +
              ylab("") +
              ggtitle("") +
              theme_classic() +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                axis.text.y = element_text(vjust = 0.5, size = 12),
                axis.title.y = element_text(angle = 90, size = 15)
              ) +
              theme(legend.position = "top")
            p2 <- FeaturePlot(object = scRNA, ncol = 1, features = i, cols = c("green", "red"))
            p <- CombinePlots(plots = list(p2, p1), rel_widths = c(1, 1.6))
            ggsave_fun(paste0(dir, "/1.", i), width = 12, height = 4, plot = p)
          }
        }
      )
      
      message("是否挑选评分最显著的细胞做火山图？yes/no")
      isSlect <- scan(what = "character", quiet = T)
      if (isSlect == "yes") {
        cat("\033[0;34m请输入细胞名字\n")
        sig.Cell <- scan(what = "character", quiet = T)
        pbmc <- subset(scRNA, subset = cellType_1 %in% sig.Cell)
      }
      if (isSlect == "no") {
        pbmc <- scRNA
      }
      
      PbmcMarkers <- FindMarkers(pbmc,
                                 ident.1 = "Hsco",
                                 min.pct = 0.1,
                                 group.by = "ScoGroup",
                                 assay = "RNA",
                                 logfc.threshold = 0.1
      )
      
      write.csv(PbmcMarkers, file = paste0(dir, "/DiffGene.csv"), row.names = T, quote = F)
      
      PbmcMarkers <- read.csv(paste0(dir, "/DiffGene.csv"), header = T, row.names = 1, check.names = F)
      
      require(ggrepel)
      cols <- c("#DC143C", "#00008B", "#808080")
      names(cols) <- c("Up", "Down", "NoSignifi")
      
      PbmcMarkers <- PbmcMarkers[-log10(PbmcMarkers$p_val_adj) > 0, ]
      PbmcMarkers$gene <- rownames(PbmcMarkers)
      PbmcMarkers$lab <- ""
      PbmcMarkers[order(PbmcMarkers$avg_log2FC), ][c(1:10, (nrow(PbmcMarkers) - 9):nrow(PbmcMarkers)), ]$lab <- PbmcMarkers[order(PbmcMarkers$avg_log2FC), ][c(1:10, (nrow(PbmcMarkers) - 9):nrow(PbmcMarkers)), ]$gene
      PbmcMarkers$logP <- -log10(PbmcMarkers$p_val_adj)
      
      p <- ggplot(PbmcMarkers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = avg_log2FC)) +
        geom_point(aes(size = logP), alpha = 0.9) +
        scale_color_gradient2(
          low = "#0500FF", high = "#FF0000", mid = "#EFEFEF",
          midpoint = 0, space = "Lab", # , limit = c(-1, 1)
          name = "logFC"
        ) +
        scale_size("-log10(p)") +
        geom_text_repel(
          data = PbmcMarkers,
          aes(x = avg_log2FC, y = -log10(p_val_adj), label = lab),
          size = 4, box.padding = unit(0.5, "lines"),
          point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = F
        ) +
        theme_bw() +
        ylab("-log10 (p_val_adj)") +
        xlab("avg_log2FC") +
        geom_vline(xintercept = c(-logFCfilter, logFCfilter), lty = 3, col = "black", lwd = 0.5) +
        geom_hline(yintercept = -log10(adjPvalFilter), lty = 3, col = "black", lwd = 0.5)
      p
      ggsave_fun(filename = paste0(dir, "/2.VlnPlot"), plot = p, width = 10, height = 8)
    }
    
    # 富集分析
    if (step == 12) {
      message("*************************富集分析*************************")
      
      dir <- paste0(Path, "/.Enrichment analysis")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      require(enrichplot)
      require(clusterProfiler)
      # cat("\033[0;34m请复制基因：\n")
      # genes <- scan(sep = "\n", what = "character", quiet = FALSE)
      genes <- KeyGene
      go <<- kk <<- NULL
      if(mouse){
        require(org.Mm.eg.db)
        entrezIDs <- mget(genes, org.Mm.egSYMBOL2EG, ifnotfound = NA)
        entrezIDs <- as.character(entrezIDs)
        out <- cbind(genes,entrezID=entrezIDs)
        write.table(out, file = paste0(dir,"/entrezIDs.txt"), sep = "\t", quote = F, row.names = F)
        
        go <<- enrichGO(
          gene = entrezIDs,
          OrgDb = org.Mm.eg.db,
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.05,
          ont = "all",
          readable = T
        )
        write.csv(go, file = paste0(dir, "/GO.csv"), quote = F, row.names = F)
        
        kk <<- enrichKEGG(gene = entrezIDs, organism = "mmu", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
        write.csv(kk, file = paste0(dir, "/KEGG.csv"), quote = F, row.names = F)
        
      }else{
        require(org.Hs.eg.db)
        entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
        entrezIDs <- as.character(entrezIDs)
        out <- cbind(genes,entrezID=entrezIDs)
        write.table(out, file = paste0(dir,"/entrezIDs.txt"), sep = "\t", quote = F, row.names = F)
        go <<- enrichGO(gene = entrezIDs,
                        OrgDb = org.Hs.eg.db, 
                        pvalueCutoff =0.05, 
                        qvalueCutoff = 0.05,
                        ont="all",
                        readable =T)
        write.csv(go, file = paste0(dir, "/GO.csv"), quote = F, row.names = F)
        
        kk <<- enrichKEGG(gene = entrezIDs, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)
        write.csv(kk, file = paste0(dir, "/KEGG.csv"), quote = F, row.names = F)
      }
      # view GO
      if(nrow(go@result) > 0){
        go.bar <- barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale='free')
        ggsave_fun(filename = paste0(dir,"/1.GOBarplot"),plot = go.bar,width = 12,height = 15)
        
        go.dot <- dotplot(go,showCategory = 10,split="ONTOLOGY",label_format=100) + facet_grid(ONTOLOGY~., scale='free')
        ggsave_fun(filename = paste0(dir,"/2.GODotplot"),plot = go.dot,width = 12,height = 15)
      }
      if(nrow(kk@result) > 0){
        if(mouse){
          kk@result$Description <- gsub(" - Mus musculus \\(house mouse\\)","",kk@result$Description)
        }
        kk.bar <- barplot(kk, drop = TRUE, showCategory = 30,label_format=100)
        ggsave_fun(filename = paste0(dir,"/3.KEGGBarplot"),plot = kk.bar,width = 12,height = 15)
        kk.dot <- dotplot(kk, showCategory = 30,label_format=100)
        ggsave_fun(filename = paste0(dir,"/4.KEGGBarplot"),plot = kk.dot,width = 12,height = 15)
      }
    }
    
  }
  return(scRNA)
}
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
set.seed(123456)

ggsave_fun <- function(filename, plot = last_plot(), plot_Device = c(".pdf"), width = 7, height = 7, dpi = 600, ...) {
  filenames <- paste0(filename, plot_Device)
  for (filename_tmp in filenames) {
    ggsave(filename = filename_tmp, plot = plot, width = width, height = height, dpi = dpi, ...)
  }
}

RePlot <- function(scRNA = scRNA,Ident,Multi = TRUE) {
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(reshape2)
  require(scRNAtoolVis)
  require(ggpubr)
  n <- 0
  for (step in 1:3) {
    n <- n + 1
    
    # 质量控制
    if (step == 1) {
      message("*************************质量控制*************************")
      
      dir <- paste0(n, ".QCFilter")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
      ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 12, height = 4)
      
      plot1 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
        NoLegend()
      plot2 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
        NoLegend()
      plot3 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.ribo", raster = F) +
        NoLegend()
      
      plot1 + plot2 + plot3
      ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
    }
    
    # 数据标准化
    if (step == 2) {
      message("*************************数据标准化*************************")
      
      dir <- paste0(n, ".Cluster")
      
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      assay <- "RNA"
      set.resolutions <- seq(0.2, 1.2, by = 0.1)
      PC <- 20
      nfeatures <- 2000
      npcs <- 50
      
      top10 <- head(VariableFeatures(scRNA), 10)
      plot1 <- VariableFeaturePlot(scRNA)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot1 + plot2
      ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
      p <- ElbowPlot(object = scRNA, ndims = npcs)
      ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)
      
      if (Multi) {
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
        p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
      } else {
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
      }
    }
    
    # 自动化注释
    if (step == 3) {
      message("*************************自动化注释*************************")
      
      dir <- paste0(n, ".CellAnnotate")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      
      # Cluster.dir <- list.files(pattern = ".Cluster$", all.files = TRUE)
      
      p <- DimPlot(scRNA, group.by = Ident, label = T, reduction = "umap") + NoLegend()
      print(p)
      ggsave_fun(filename = paste0(dir, "/1.Umap.Cluster"), plot = p)
      
      
      # Idents(scRNA) <- "cellType_1"
      # scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      # write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))
      
      scRNA.markers <- read.csv(".CellAnnotate/Findall.markers.cellType_1.csv",header = T,row.names = 1)
      
      scRNA.markers %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC) -> features
      
      print((features))
      
      features <- features$gene
      cat("\033[0;34m是否要改写基因？yes/no\n")
      RenameGene <- scan(what = "character", sep = ",", quiet = TRUE)
      
      if (RenameGene == "yes") {
        if (exists("marker", envir = globalenv())) {
          features <- marker
        }else{
          print(paste0(as.character(features), collapse = ","))
          cat("\033[0;34m请改写\n")
          features <- scan(what = "character", sep = ",", quiet = FALSE)
        }
      }
      
      # plot_list <- FeaturePlot(
      #   scRNA,
      #   features = unique(features),
      #   reduction = "umap",
      #   raster = F,
      #   combine = FALSE, cols = c("lightgrey", "red")
      # )
      # 
      # for (i in 1:length(plot_list)) {
      #   plot_list[[i]] <- plot_list[[i]] + NoLegend() + NoAxes()
      # }
      # 
      # p <- CombinePlots(plot_list)
      # ggsave_fun(filename = paste0(dir, "/1.Featureplot.umap"), plot = p, width = 12, height = 12)
      
      my_pal2 <- c(
        "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
        "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
        "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
        "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A"
      )
      
      umap <- scRNA@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        cbind(cellType = scRNA@meta.data$cellType_1)
      
      ## 计算标签中心位置区域
      celltypepos <- umap %>%
        group_by(cellType) %>%
        summarise(
          umap_1 = median(UMAP_1),
          umap_2 = median(UMAP_2)
        )
      
      p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) +
        ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType),
                                  fontface = "bold",
                                  data = celltypepos,
                                  box.padding = 0.5, show.legend = FALSE
        )
      print(p)
      ggsave_fun(filename = paste0(dir, "/2.Umap1.plot"), plot = p, width = 10, height = 7)
      
      # add circle
      # p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T, cornerTextSize = 3.5, addCircle = TRUE)
      # print(p)
      # ggsave_fun(filename = paste0(dir, "/3.Umap2.plot"), plot = p, width = 10, height = 7)
      
      # Idents(scRNA) <- "cellType_1"
      # 
      # library(plot1cell)
      # 
      # # devtools::install_github("TheHumphreysLab/plot1cell")
      # # ## or the development version, devtools::install_github("HaojiaWu/plot1cell")
      # # 
      # # ## You might need to install the dependencies below if they are not available in your R library.
      # # bioc.packages <- c("biomaRt","GenomeInfoDb","EnsDb.Hsapiens.v86","GEOquery","simplifyEnrichment","ComplexHeatmap")
      # # BiocManager::install(bioc.packages)
      # # dev.packages <- c("chris-mcginnis-ucsf/DoubletFinder","Novartis/hdf5r","mojaveazure/loomR")
      # # devtools::install_github(dev.packages)
      # # ## If you can't get the hdf5r package installed, please see the fix here:
      # # ## http<PATH>//github.com/hhoeflin/hdf5r/issues/94
      # 
      # circ_data <- prepare_circlize_data(scRNA, scale = 0.8)
      # # set.seed(1234)
      # cluster_colors<-rand_color(length(levels(scRNA)))
      # group_colors<-rand_color(length(names(table(scRNA$group))))
      # rep_colors<-rand_color(length(names(table(scRNA$orig.ident))))
      # 
      # pdf(paste0(dir,'/2.Circlize_plot.pdf'), width = 10, height = 10)
      # plot_circlize(circ_data,do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.6)
      # add_track(circ_data, group = "group", colors = group_colors, track_num = 2) ## can change it to one of the columns in the meta data of your seurat object
      # add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3) ## can change it to one of the columns in the meta data of your seurat object
      # dev.off()
      
      p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 10)) + theme(
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        ),
        legend.position = "top",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.8, "cm"),
      )
      print(p)
      if (length(unique(features)) >= 50) {
        ggsave_fun(filename = paste0(dir, "/3.Dotplot"), plot = p, width = 20, height = 9)
      }
      if (length(unique(features)) < 50) {
        ggsave_fun(filename = paste0(dir, "/3.Dotplot"), plot = p, width = 17, height = 9)
      }
      
      mynames <-   table(scRNA$cellType_1) %>% names()
      myratio <-  table(scRNA$cellType_1) %>% as.numeric()
      pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")
      
      if(length(pielabel) <= 20){
        cols <-c('#E64A35','#4DBBD4' ,'#01A187','#6BD66B','#3C5588'  ,'#F29F80'  ,
                 '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
                 '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
        )
        pdf(paste0(dir,"/4.CellType_1.ratio.plot.pdf"))
        pie(myratio, labels=pielabel,
            radius = 1.0,clockwise=T,
            main = "celltype",col = cols) 
        dev.off()
      }
      
      tryCatch(
        {
          cellRatio <- cellRatioPlot(
            object = scRNA,
            sample.name = "group",
            celltype.name = "cellType_1",
            flow.curve = 0.5,
            fill.col = my_pal2
          ) + theme(
            axis.text.x = element_text(
              angle = 45,
              hjust = 1,
              vjust = 1
            )
          )
          print(cellRatio)
          
          ggsave_fun(filename = paste0(dir, "/4.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 7, height = 7)
          
          data_p <- scRNA@meta.data[, c("group", "cellType_1")]
          
          p <- ggstatsplot::ggbarstats(
            data = data_p,
            x = group,
            y = cellType_1,
            palette = "Set3"
          ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          print(p)
          leg <- length(unique(data_p$cellType_1))
          ggsave_fun(filename = paste0(dir, "/4.Dysfunction"), plot = p, width = 2+1.5*leg, height = 7)
          
          EachSampleClusterDis <- lapply(unique(scRNA@meta.data$orig.ident), function(Sample) {
            seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$orig.ident %in% Sample, ]))
            seur@meta.data$cellType_1 %>%
              table() %>%
              data.frame() %>%
              magrittr::set_colnames(c("CellTypes", "Number")) %>%
              dplyr::mutate(Per = 100 * Number / sum(Number))
          })
          names(EachSampleClusterDis) <- unique(scRNA@meta.data$orig.ident)
          EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
          p1 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
            geom_col(position = "fill", width = 0.8) +
            scale_y_continuous(expand = c(0, 0)) +
            labs(y = "Percentage of cell") +
            scale_fill_manual(values = my_pal2, name = "") +
            theme(
              axis.text = element_text(color = "black"),
              panel.background = element_blank(), # panel.grid=element_blank(),
              legend.title = element_blank(),
              axis.text.y = element_text(color = "black"),
              axis.line = element_line(color = "black"),
              axis.title.y = element_blank(), # 轴标题
              legend.position = "right",
              # legend.key.height = unit(0.3, "cm"),
              # legend.key.width = unit(0.3, "cm"),
              legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
              legend.text = element_text(size = 10)
            ) +
            coord_flip()
          
          EachSampleClusterDis <- lapply(unique(scRNA@meta.data$group), function(Sample) {
            seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$group %in% Sample, ]))
            seur@meta.data$cellType_1 %>%
              table() %>%
              data.frame() %>%
              magrittr::set_colnames(c("CellTypes", "Number")) %>%
              dplyr::mutate(Per = 100 * Number / sum(Number))
          })
          
          my_pal1 <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#6387C5")
          
          names(EachSampleClusterDis) <- unique(scRNA@meta.data$group)
          EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
          p2 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
            geom_col(position = "fill", width = 0.8) +
            scale_y_continuous(expand = c(0, 0)) +
            labs(y = "Percentage of cell") +
            scale_fill_manual(values = my_pal1, name = "") +
            theme(
              axis.text = element_text(color = "black", ),
              panel.background = element_blank(), panel.grid = element_blank(),
              legend.title = element_blank(),
              axis.text.y = element_blank(), axis.title.x = element_blank(), # 删除轴文本
              axis.ticks.y.left = element_blank(),
              axis.title.y = element_blank(), # 轴标题
              axis.line.x = element_line(color = "black"),
              legend.position = "right",
              legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
              # legend.key.height = unit(0.3, "cm"),
              # legend.key.width = unit(0.3, "cm"),
              legend.text = element_text(size = 10)
            ) +
            coord_flip()
          
          # 计算每种细胞类型的细胞数量
          cell_counts <- table(scRNA@meta.data$cellType_1)
          
          # 创建数据框
          cell_counts_df <- data.frame(
            CellType = names(cell_counts),
            Count = as.numeric(cell_counts / 1000)
          )
          
          # 创建柱状图
          cell_counts_df <- cell_counts_df[order(cell_counts_df$Count, decreasing = T), ]
          
          p3 <- ggplot(cell_counts_df, aes(x = reorder(CellType, Count), y = Count)) +
            scale_y_continuous(expand = c(0, 0)) +
            labs(y = "Number of cell(10^3)") +
            geom_bar(stat = "identity", fill = "#6387C5") +
            theme(
              axis.text = element_text(color = "black"),
              panel.background = element_blank(), panel.grid = element_blank(),
              legend.title = element_blank(),
              axis.text.y = element_blank(),
              # axis.title.x = element_blank(),#删除轴文本
              # plot.tag=element_blank(),
              axis.ticks.y.left = element_blank(),
              # axis.title.x.top = T,
              axis.title.y = element_blank(), # 轴标题
              axis.line.x = element_line(color = "black"),
              legend.text = element_text(size = 10)
            ) +
            coord_flip()
          print(p1 + p2 + p3)
          ggsave_fun(filename = paste0(dir, "/4.MergeRatio"), plot = p1 + p2 + p3, width = 15, height = 6)
          # if (isKeyCell == "yes" && is.na(KeyCell) && 1 %in% KeyCellIndex) {
          pB2_df <- table(scRNA@meta.data$cellType_1, scRNA@meta.data$group) %>% melt()
          colnames(pB2_df) <- c("Cluster", "Sample", "Number")
          pB2_df$Cluster <- factor(pB2_df$Cluster)
          
          pB2_df <- pB2_df %>%
            group_by(Sample) %>%
            mutate(Percentage = Number / sum(Number))
          
          pB2_df_tumor <- pB2_df %>%
            filter(Sample %in% c("Disease", "Tumor", "Cancer")) %>%
            select(Cluster, Tumor_Number = Number)
          
          pB2_df_control <- pB2_df %>%
            filter(Sample == c("Control", "Normal")) %>%
            select(Cluster, Control_Number = Number)
          
          pB2_df_combined <- merge(pB2_df_tumor, pB2_df_control, by = "Cluster")
          # pB2_df_combined <- pB2_df_combined[pB2_df_combined$Control_Number >= 50 && pB2_df_combined$Tumor_Number >= 50,]
          pB2_df_combined$rate <- (pB2_df_combined$Tumor_Number - pB2_df_combined$Control_Number) / pB2_df_combined$Control_Number
          KeyCell <- pB2_df_combined[abs(pB2_df_combined$rate) == max(abs(pB2_df_combined$rate)), ]$Cluster %>% as.character()
          print(paste0("变化最大的细胞：", KeyCell))
          write.csv(pB2_df_combined, paste0(dir, "/Cell.rate.csv"))
          # }
        },
        error = function(e) {
          message("没有分组信息，跳过此步骤\n", e)
          
          cellRatio <- cellRatioPlot(
            object = scRNA,
            sample.name = "orig.ident",
            celltype.name = "cellType_1",
            flow.curve = 0.5,
            fill.col = my_pal2
          ) + theme(
            axis.text.x = element_text(
              angle = 45,
              hjust = 1,
              vjust = 1
            )
          )
          print(cellRatio)
          
          ggsave_fun(filename = paste0(dir, "/4.Sample.ratio.plot"), plot = cellRatio, width = 7 + 0.25 * (length(unique(scRNA$orig.ident))), height = 7)
        }
      )
    }
  }
}


# scRNA <- readRDS("result.rds")
# RePlot(scRNA = scRNA,Ident = "RNA_snn_res.0.2",Multi = TRUE)

# 单细胞审核专用函数 ---------------------------------------------------------------
RePlot <- function(scRNA = scRNA,Ident,Multi = TRUE) {
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(reshape2)
  require(scRNAtoolVis)
  require(ggpubr)
  n <- 0
  for (step in 1:3) {
    n <- n + 1
    
    # 质量控制
    if (step == 1) {
      message("*************************质量控制*************************")
      
      dir <- paste0(n, ".QCFilter")
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      # p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 4)
      # ggsave_fun(filename = paste0(dir, "/1.VlnPlot"), plot = p, width = 20, height = 4)
      
      p <- VlnPlot(scRNA, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), group.by = "orig.ident", ncol = 2)  
      ggsave_fun(filename = paste0(dir, "1.VlnPlot"), plot = p, width = 20, height = 8)
      
      plot1 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F) +
        NoLegend()
      plot2 <-
        FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F) +
        NoLegend()
      plot3 <-
        FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.ribo", raster = F) +
        NoLegend()
      
      plot1 + plot2 + plot3
      ggsave_fun(filename = paste0(dir, "/2.FeatureScatter"), plot = plot1 + plot2 + plot3, height = 7, width = 21)
    }
    
    # 数据标准化
    if (step == 2) {
      message("*************************数据标准化*************************")
      
      dir <- paste0(n, ".Cluster")
      
      if (!dir.exists(dir)) {
        dir.create(dir)
      } else {
        print("Dir already exists!")
      }
      
      assay <- "RNA"
      set.resolutions <- seq(0.2, 1.2, by = 0.1)
      PC <- 20
      nfeatures <- 2000
      npcs <- 50
      
      top10 <- head(VariableFeatures(scRNA), 10)
      plot1 <- VariableFeaturePlot(scRNA)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot1 + plot2
      ggsave_fun(filename = paste0(dir, "/1.VariableFeaturePlot"), plot = plot1 + plot2, width = 14)
      p <- ElbowPlot(object = scRNA, ndims = npcs)
      ggsave_fun(filename = paste0(dir, "/2.ElbowPlot"), plot = p)
      
      if (Multi) {
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
        p <- DimPlot(object = scRNA, reduction = "harmony", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/4.HarmonyPlot"), plot = p)
      } else {
        p <- DimPlot(object = scRNA, reduction = "pca", label = F, group.by = "orig.ident", raster = F) + NoLegend()
        ggsave_fun(filename = paste0(dir, "/3.PcaPlot"), plot = p)
      }
    }
    
    # # 自动化注释
    # if (step == 3) {
    #   message("*************************自动化注释*************************")
    #   
    #   dir <- paste0(n, ".CellAnnotate")
    #   if (!dir.exists(dir)) {
    #     dir.create(dir)
    #   } else {
    #     print("Dir already exists!")
    #   }
    # 
    # 
    #   Cluster.dir <- list.files(pattern = ".Cluster$", all.files = TRUE)
    # 
    #   p <- DimPlot(scRNA, group.by = Ident, label = T, reduction = "umap") + NoLegend()
    #   print(p)
    #   ggsave_fun(filename = paste0(Cluster.dir, "/5.Umap.Cluster"), plot = p)
    # 
    # 
    #   # Idents(scRNA) <- "cellType_1"
    #   # scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    #   # write.csv(scRNA.markers, paste0(dir, "/Findall.markers.cellType_1.csv"))
    #   
    #   scRNA.markers <- read.csv("Findall.markers.cellType_1.csv",header = T,row.names = 1)
    # 
    #   scRNA.markers %>%
    #     group_by(cluster) %>%
    #     top_n(n = 5, wt = avg_log2FC) -> features
    # 
    #   print((features))
    #   
    #   features <- features$gene
    #   cat("\033[0;34m是否要改写基因？yes/no\n")
    #   RenameGene <- scan(what = "character", sep = ",", quiet = TRUE)
    #   
    #   if (RenameGene == "yes") {
    #     print(paste0(as.character(features), collapse = ","))
    #     cat("\033[0;34m请改写\n")
    #     features <- scan(what = "character", sep = ",", quiet = FALSE)
    #   }
    # 
    #   # plot_list <- FeaturePlot(
    #   #   scRNA,
    #   #   features = unique(features),
    #   #   reduction = "umap",
    #   #   raster = F,
    #   #   combine = FALSE, cols = c("lightgrey", "red")
    #   # )
    #   # 
    #   # for (i in 1:length(plot_list)) {
    #   #   plot_list[[i]] <- plot_list[[i]] + NoLegend() + NoAxes()
    #   # }
    #   # 
    #   # p <- CombinePlots(plot_list)
    #   # ggsave_fun(filename = paste0(dir, "/1.Featureplot.umap"), plot = p, width = 12, height = 12)
    # 
    #   my_pal2 <- c(
    #     "#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
    #     "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
    #     "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A","#D4477D", "#D24B27", "#4DBBD5", "#6387C5", "#6E4B9E", "#C10020", "#1E78B4", "#FCBF6E", "#83AD00", "#9ebcda",
    #     "#74a9cf", "#fbdf72", "#FF8E00", "#F37B7D", "#CF4A31", "#F37B7D", "#FF8E00", "#00B3F1", "#00538A"
    #   )
    #   
    #   umap <- scRNA@reductions$umap@cell.embeddings %>%
    #     as.data.frame() %>%
    #     cbind(cellType = scRNA@meta.data$cellType_1)
    #   
    #   ## 计算标签中心位置区域
    #   celltypepos <- umap %>%
    #     group_by(cellType) %>%
    #     summarise(
    #       umap_1 = median(UMAP_1),
    #       umap_2 = median(UMAP_2)
    #     )
    #   
    #   p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T) +
    #     ggrepel::geom_label_repel(aes(x = umap_1, y = umap_2, label = cellType, color = cellType),
    #                               fontface = "bold",
    #                               data = celltypepos,
    #                               box.padding = 0.5, show.legend = FALSE
    #     )
    #   print(p)
    #   ggsave_fun(filename = paste0(dir, "/1.Umap1.plot"), plot = p, width = 10, height = 7)
    #   
    #   # add circle
    #   # p <- clusterCornerAxes(object = scRNA, reduction = "umap", pSize = 0.1, clusterCol = "cellType_1", noSplit = T, cornerTextSize = 3.5, addCircle = TRUE)
    #   # print(p)
    #   # ggsave_fun(filename = paste0(dir, "/3.Umap2.plot"), plot = p, width = 10, height = 7)
    #   
    #   p <- DotPlot(scRNA, features = unique(features), cols = "RdYlBu", group.by = "cellType_1") + scale_size_continuous(range = c(0, 10)) + theme(
    #     panel.border = element_rect(colour = "black"),
    #     axis.text.x = element_text(
    #       angle = 90,
    #       hjust = 1,
    #       vjust = 0.5
    #     ),
    #     legend.position = "top",
    #     legend.key.height = unit(0.3, "cm"),
    #     legend.key.width = unit(0.8, "cm"),
    #   )
    #   print(p)
    #   if (length(unique(features)) >= 50) {
    #     ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 20, height = 9)
    #   }
    #   if (length(unique(features)) < 50) {
    #     ggsave_fun(filename = paste0(dir, "/2.Dotplot"), plot = p, width = 17, height = 9)
    #   }
    #   
    #   tryCatch(
    #     {
    #       cellRatio <- cellRatioPlot(
    #         object = scRNA,
    #         sample.name = "group",
    #         celltype.name = "cellType_1",
    #         flow.curve = 0.5,
    #         fill.col = my_pal2
    #       ) + theme(
    #         axis.text.x = element_text(
    #           angle = 45,
    #           hjust = 1,
    #           vjust = 1
    #         )
    #       )
    #       print(cellRatio)
    #       
    #       ggsave_fun(filename = paste0(dir, "/3.CellAnnoted.cellType.group.ratio.plot"), plot = cellRatio, width = 7, height = 7)
    #       
    #       data_p <- scRNA@meta.data[, c("group", "cellType_1")]
    #       
    #       p <- ggstatsplot::ggbarstats(
    #         data = data_p,
    #         x = group,
    #         y = cellType_1,
    #         palette = "Set3"
    #       ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #       print(p)
    #       leg <- length(unique(data_p$cellType_1))
    #       ggsave_fun(filename = paste0(dir, "/4.Dysfunction"), plot = p, width = 2+1.5*leg, height = 7)
    #       
    #       EachSampleClusterDis <- lapply(unique(scRNA@meta.data$orig.ident), function(Sample) {
    #         seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$orig.ident %in% Sample, ]))
    #         seur@meta.data$cellType_1 %>%
    #           table() %>%
    #           data.frame() %>%
    #           magrittr::set_colnames(c("CellTypes", "Number")) %>%
    #           dplyr::mutate(Per = 100 * Number / sum(Number))
    #       })
    #       names(EachSampleClusterDis) <- unique(scRNA@meta.data$orig.ident)
    #       EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
    #       p1 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
    #         geom_col(position = "fill", width = 0.8) +
    #         scale_y_continuous(expand = c(0, 0)) +
    #         labs(y = "Percentage of cell") +
    #         scale_fill_manual(values = my_pal2, name = "") +
    #         theme(
    #           axis.text = element_text(color = "black"),
    #           panel.background = element_blank(), # panel.grid=element_blank(),
    #           legend.title = element_blank(),
    #           axis.text.y = element_text(color = "black"),
    #           axis.line = element_line(color = "black"),
    #           axis.title.y = element_blank(), # 轴标题
    #           legend.position = "right",
    #           # legend.key.height = unit(0.3, "cm"),
    #           # legend.key.width = unit(0.3, "cm"),
    #           legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
    #           legend.text = element_text(size = 10)
    #         ) +
    #         coord_flip()
    #       
    #       EachSampleClusterDis <- lapply(unique(scRNA@meta.data$group), function(Sample) {
    #         seur <- subset(scRNA, cells = rownames(scRNA@meta.data[scRNA@meta.data$group %in% Sample, ]))
    #         seur@meta.data$cellType_1 %>%
    #           table() %>%
    #           data.frame() %>%
    #           magrittr::set_colnames(c("CellTypes", "Number")) %>%
    #           dplyr::mutate(Per = 100 * Number / sum(Number))
    #       })
    #       
    #       my_pal1 <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#6387C5")
    #       
    #       names(EachSampleClusterDis) <- unique(scRNA@meta.data$group)
    #       EachSampleClusterDisA <- dplyr::bind_rows(EachSampleClusterDis) %>% dplyr::mutate(Sample = rep(names(EachSampleClusterDis), times = unlist(lapply(EachSampleClusterDis, nrow))))
    #       p2 <- ggplot(data = EachSampleClusterDisA, aes(x = reorder(CellTypes, Number), fill = Sample, y = Per)) + # log2(as.numeric(CopyNumber)))
    #         geom_col(position = "fill", width = 0.8) +
    #         scale_y_continuous(expand = c(0, 0)) +
    #         labs(y = "Percentage of cell") +
    #         scale_fill_manual(values = my_pal1, name = "") +
    #         theme(
    #           axis.text = element_text(color = "black", ),
    #           panel.background = element_blank(), panel.grid = element_blank(),
    #           legend.title = element_blank(),
    #           axis.text.y = element_blank(), axis.title.x = element_blank(), # 删除轴文本
    #           axis.ticks.y.left = element_blank(),
    #           axis.title.y = element_blank(), # 轴标题
    #           axis.line.x = element_line(color = "black"),
    #           legend.position = "right",
    #           legend.direction = "vertical", # vertical  legend.direction = 'horizontal'
    #           # legend.key.height = unit(0.3, "cm"),
    #           # legend.key.width = unit(0.3, "cm"),
    #           legend.text = element_text(size = 10)
    #         ) +
    #         coord_flip()
    #       
    #       # 计算每种细胞类型的细胞数量
    #       cell_counts <- table(scRNA@meta.data$cellType_1)
    #       
    #       # 创建数据框
    #       cell_counts_df <- data.frame(
    #         CellType = names(cell_counts),
    #         Count = as.numeric(cell_counts / 1000)
    #       )
    #       
    #       # 创建柱状图
    #       cell_counts_df <- cell_counts_df[order(cell_counts_df$Count, decreasing = T), ]
    #       
    #       p3 <- ggplot(cell_counts_df, aes(x = reorder(CellType, Count), y = Count)) +
    #         scale_y_continuous(expand = c(0, 0)) +
    #         labs(y = "Number of cell(10^3)") +
    #         geom_bar(stat = "identity", fill = "#6387C5") +
    #         theme(
    #           axis.text = element_text(color = "black"),
    #           panel.background = element_blank(), panel.grid = element_blank(),
    #           legend.title = element_blank(),
    #           axis.text.y = element_blank(),
    #           # axis.title.x = element_blank(),#删除轴文本
    #           # plot.tag=element_blank(),
    #           axis.ticks.y.left = element_blank(),
    #           # axis.title.x.top = T,
    #           axis.title.y = element_blank(), # 轴标题
    #           axis.line.x = element_line(color = "black"),
    #           legend.text = element_text(size = 10)
    #         ) +
    #         coord_flip()
    #       print(p1 + p2 + p3)
    #       ggsave_fun(filename = paste0(dir, "/5.MergeRatio"), plot = p1 + p2 + p3, width = 15, height = 6)
    #       # if (isKeyCell == "yes" && is.na(KeyCell) && 1 %in% KeyCellIndex) {
    #       pB2_df <- table(scRNA@meta.data$cellType_1, scRNA@meta.data$group) %>% melt()
    #       colnames(pB2_df) <- c("Cluster", "Sample", "Number")
    #       pB2_df$Cluster <- factor(pB2_df$Cluster)
    #       
    #       pB2_df <- pB2_df %>%
    #         group_by(Sample) %>%
    #         mutate(Percentage = Number / sum(Number))
    #       
    #       pB2_df_tumor <- pB2_df %>%
    #         filter(Sample %in% c("Disease", "Tumor", "Cancer")) %>%
    #         select(Cluster, Tumor_Number = Number)
    #       
    #       pB2_df_control <- pB2_df %>%
    #         filter(Sample == c("Control", "Normal")) %>%
    #         select(Cluster, Control_Number = Number)
    #       
    #       pB2_df_combined <- merge(pB2_df_tumor, pB2_df_control, by = "Cluster")
    #       # pB2_df_combined <- pB2_df_combined[pB2_df_combined$Control_Number >= 50 && pB2_df_combined$Tumor_Number >= 50,]
    #       pB2_df_combined$rate <- (pB2_df_combined$Tumor_Number - pB2_df_combined$Control_Number) / pB2_df_combined$Control_Number
    #       KeyCell <- pB2_df_combined[abs(pB2_df_combined$rate) == max(abs(pB2_df_combined$rate)), ]$Cluster %>% as.character()
    #       print(paste0("变化最大的细胞：", KeyCell))
    #       write.csv(pB2_df_combined, paste0(dir, "/Cell.rate.csv"))
    #       # }
    #     },
    #     error = function(e) {
    #       message("没有分组信息，跳过此步骤\n", e)
    #       
    #       cellRatio <- cellRatioPlot(
    #         object = scRNA,
    #         sample.name = "orig.ident",
    #         celltype.name = "cellType_1",
    #         flow.curve = 0.5,
    #         fill.col = my_pal2
    #       ) + theme(
    #         axis.text.x = element_text(
    #           angle = 45,
    #           hjust = 1,
    #           vjust = 1
    #         )
    #       )
    #       print(cellRatio)
    #       
    #       ggsave_fun(filename = paste0(dir, "/3.Sample.ratio.plot"), plot = cellRatio, width = 7 + 0.25 * (length(unique(scRNA$orig.ident))), height = 7)
    #     }
    #   )
    # }
  }
}


# scRNA <- readRDS("result.rds")
# Idents(scRNA) <- "cellType_1"
# scRNA.markers <- FindAllMarkers(scRNA, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
# scRNA.markers <- scRNA.markers[scRNA.markers$p_val_adj < 0.05,]
#write.csv(scRNA.markers, "./3.CellAnnotate/Findall.markers.cellType_1.csv")

# RePlot(scRNA = scRNA,Ident = "RNA_snn_res.0.2",Multi = TRUE)
# ClusterGvis -------------------------------------------------------------
require(ClusterGVis)
require(dplyr)
require(jjAnno)
require(Seurat)
plan("multicore", workers = 8)
options(future.globals.maxSize = 5000000 * 1024^2)
set.seed(123456)

#' Title
#'
#' @param object Seurat object
#' @param OutPath 设置输出路径
#' @param scRNA.markers FindAllMarker输出对象，默认NULL，将自动运行FindAllMarker，如果有将跳过这一步
#' @param celltype Seurat对象中存放细胞类型信息的表头，例如cellType_1
#' @param species 物种，人无需填写，鼠填写"mouse"
#' @param OnlyPlot 是否只画图，复现填True
#'
#' @return
#' @export
#'
#' @examples
ClusterGVis_module <- function(object,
                               OutPath,
                               scRNA.markers = NULL,
                               celltype = "cellType_1",
                               species = NA,
                               OnlyPlot = F){
  
  if(!dir.exists(paste0(OutPath,"/ClusterGvis")))dir.create(paste0(OutPath,"/ClusterGvis"))
  setwd(OutPath)
  
  if(!OnlyPlot){
    if (!is.na(species)) {
      mouse <- ifelse(species == "mouse", T, F)
    } else {
      mouse <- F
    }
    Idents(object) <- celltype
    if(is.null(scRNA.markers)){
      scRNA.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      write.csv(scRNA.markers,paste0("ClusterGvis/Findall.markers.",celltype,".csv"))
    }else{
      scRNA.markers <- scRNA.markers
    }
    scRNA.top.markers <- scRNA.markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n=5,wt=avg_log2FC)
    write.csv(scRNA.top.markers,"ClusterGvis/TopMarkers.csv")
    
    st.data <- prepareDataFromscRNA(object = object,diffData = scRNA.top.markers,showAverage = T)
    
    if(mouse){
      require(org.Mm.eg.db)
      enrich <- enrichCluster(object = st.data,
                              OrgDb = org.Mm.eg.db,
                              type = "BP",
                              organism = "mmu",
                              pvalueCutoff = 0.05,
                              topn = 5,
                              seed = 123456)
    }else{
      require(org.Hs.eg.db)
      enrich <- enrichCluster(object = st.data,
                              OrgDb = org.Hs.eg.db,
                              type = "BP",
                              organism = "hsa",
                              pvalueCutoff = 0.05,
                              topn = 5,
                              seed = 123456)
    }
    saveRDS(st.data,paste0('ClusterGvis/Heatmap.plotdata.rds'))
    write.csv(enrich,paste0('ClusterGvis/Enrich.plotdata.csv'))
  }
  
  st.data <- readRDS(paste0('ClusterGvis/Heatmap.plotdata.rds'))
  enrich <- read.csv(paste0('ClusterGvis/Enrich.plotdata.csv'),row.names = 1)
  scRNA.top.markers <- read.csv("ClusterGvis/TopMarkers.csv")
  pdf("ClusterGvis/1.Heatmap.pdf",height = 15,width =17,onefile = F)
  visCluster(object = st.data,
             plot.type = "both",
             column_title_rot = 45,
             markGenes = unique(scRNA.top.markers$gene),
             markGenes.side = "left",
             annoTerm.data = enrich,
             # show_column_names = F,
             line.side = "left",
             cluster.order = c(1:length(levels(object@meta.data[[celltype]]))),
             add.bar = T)
  dev.off()
}






# 空间转录组互作分析 ---------------------------------------------------------------

# Copyright (c) [2021] [Ricardo O. Ramirez Flores]
# roramirezf@uni-heidelberg.de

#' Catalog of MISTy utilities
#' 
run_misty_seurat <- function(visium.slide,
                             # Seurat object with spatial transcriptomics data.
                             view.assays,
                             # Named list of assays for each view.
                             view.features = NULL,
                             # Named list of features/markers to use.
                             # Use all by default.
                             view.types,
                             # Named list of the type of view to construct
                             # from the assay.
                             view.params,
                             # Named list with parameters (NULL or value)
                             # for each view.
                             spot.ids = NULL,
                             # spot IDs to use. Use all by default.
                             out.alias = "results"
                             # folder name for output
) {
  
  mistyR::clear_cache()
  
  # Extracting geometry
  geometry <- GetTissueCoordinates(visium.slide,
                                   cols = c("row", "col"), scale = NULL
  )
  
  # Extracting data
  view.data <- map(view.assays,
                   extract_seurat_data,
                   geometry = geometry,
                   visium.slide = visium.slide
  )
  
  # Constructing and running a workflow
  build_misty_pipeline(
    view.data = view.data,
    view.features = view.features,
    view.types = view.types,
    view.params = view.params,
    geometry = geometry,
    spot.ids = spot.ids,
    out.alias = out.alias
  )
}


# Extracts data from an specific assay from a Seurat object
# and aligns the IDs to the geometry
extract_seurat_data <- function(visium.slide,
                                assay,
                                geometry) {
  print(assay)
  data <- GetAssayData(visium.slide, assay = assay) %>%
    as.matrix() %>%
    t() %>%
    as_tibble(rownames = NA)
  
  return(data %>% slice(match(rownames(.), rownames(geometry))))
}

# Filters data to contain only features of interest
filter_data_features <- function(data,
                                 features) {
  if (is.null(features)) features <- colnames(data)
  
  return(data %>% rownames_to_column() %>%
           select(rowname, all_of(features)) %>% rename_with(make.names) %>%
           column_to_rownames())
}

# Builds views depending on the paramaters defined
create_default_views <- function(data,
                                 view.type,
                                 view.param,
                                 view.name,
                                 spot.ids,
                                 geometry) {
  
  mistyR::clear_cache()
  
  view.data.init <- create_initial_view(data)
  
  if (!(view.type %in% c("intra", "para", "juxta"))) {
    view.type <- "intra"
  }
  
  if (view.type == "intra") {
    data.red <- view.data.init[["intraview"]]$data %>%
      rownames_to_column() %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "para") {
    view.data.tmp <- view.data.init %>%
      add_paraview(geometry, l = view.param)
    
    data.ix <- paste0("paraview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  } else if (view.type == "juxta") {
    view.data.tmp <- view.data.init %>%
      add_juxtaview(
        positions = geometry,
        neighbor.thr = view.param
      )
    
    data.ix <- paste0("juxtaview.", view.param)
    data.red <- view.data.tmp[[data.ix]]$data %>%
      mutate(rowname = rownames(data)) %>%
      filter(rowname %in% spot.ids) %>%
      select(-rowname)
  }
  
  if (is.null(view.param) == TRUE) {
    misty.view <- create_view(
      paste0(view.name),
      data.red
    )
  } else {
    misty.view <- create_view(
      paste0(view.name, "_", view.param),
      data.red
    )
  }
  
  return(misty.view)
}

# Builds automatic MISTy workflow and runs it
build_misty_pipeline <- function(view.data,
                                 view.features,
                                 view.types,
                                 view.params,
                                 geometry,
                                 spot.ids = NULL,
                                 out.alias = "default") {
  
  # Adding all spots ids in case they are not defined
  if (is.null(spot.ids)) {
    spot.ids <- rownames(view.data[[1]])
  }
  
  # First filter the features from the data
  view.data.filt <- map2(view.data, view.features, filter_data_features)
  
  # Create initial view
  views.main <- create_initial_view(view.data.filt[[1]] %>%
                                      rownames_to_column() %>%
                                      filter(rowname %in% spot.ids) %>%
                                      select(-rowname))
  
  # Create other views
  view.names <- names(view.data.filt)
  
  all.views <- pmap(list(
    view.data.filt[-1],
    view.types[-1],
    view.params[-1],
    view.names[-1]
  ),
  create_default_views,
  spot.ids = spot.ids,
  geometry = geometry
  )
  
  pline.views <- add_views(
    views.main,
    unlist(all.views, recursive = FALSE)
  )
  
  
  # Run MISTy
  run_misty(pline.views, out.alias, cached = FALSE)
}

#
# Bug in collecting results
#

collect_results_v2 <- function(folders){
  samples <- R.utils::getAbsolutePath(folders)
  message("\nCollecting improvements")
  improvements <- samples %>% furrr::future_map_dfr(function(sample) {
    performance <- readr::read_table2(paste0(sample, .Platform$file.sep, 
                                             "performance.txt"), na = c("", "NA", "NaN"), col_types = readr::cols()) %>% 
      dplyr::distinct()
    performance %>% dplyr::mutate(sample = sample, gain.RMSE = 100 * 
                                    (.data$intra.RMSE - .data$multi.RMSE)/.data$intra.RMSE, 
                                  gain.R2 = 100 * (.data$multi.R2 - .data$intra.R2), 
    )
  }, .progress = TRUE) %>% tidyr::pivot_longer(-c(.data$sample, 
                                                  .data$target), names_to = "measure")
  message("\nCollecting contributions")
  contributions <- samples %>% furrr::future_map_dfr(function(sample) {
    coefficients <- readr::read_table2(paste0(sample, .Platform$file.sep, 
                                              "coefficients.txt"), na = c("", "NA", "NaN"), col_types = readr::cols()) %>% 
      dplyr::distinct()
    coefficients %>% dplyr::mutate(sample = sample, .after = "target") %>% 
      tidyr::pivot_longer(cols = -c(.data$sample, .data$target), 
                          names_to = "view")
  }, .progress = TRUE)
  improvements.stats <- improvements %>% dplyr::filter(!stringr::str_starts(.data$measure, 
                                                                            "p\\.")) %>% dplyr::group_by(.data$target, .data$measure) %>% 
    dplyr::summarise(mean = mean(.data$value), sd = stats::sd(.data$value), 
                     cv = .data$sd/.data$mean, .groups = "drop")
  contributions.stats <- dplyr::inner_join((contributions %>% 
                                              dplyr::filter(!stringr::str_starts(.data$view, "p\\.") & 
                                                              .data$view != "intercept") %>% dplyr::group_by(.data$target, 
                                                                                                             .data$view) %>% dplyr::summarise(mean = mean(.data$value), 
                                                                                                                                              .groups = "drop_last") %>% dplyr::mutate(fraction = abs(.data$mean)/sum(abs(.data$mean))) %>% 
                                              dplyr::ungroup()), (contributions %>% dplyr::filter(stringr::str_starts(.data$view, 
                                                                                                                      "p\\.") & !stringr::str_detect(.data$view, "intercept")) %>% 
                                                                    dplyr::group_by(.data$target, .data$view) %>% dplyr::mutate(view = stringr::str_remove(.data$view, 
                                                                                                                                                           "^p\\.")) %>% dplyr::summarise(p.mean = mean(.data$value), 
                                                                                                                                                                                          p.sd = stats::sd(.data$value), .groups = "drop")), by = c("target", 
                                                                                                                                                                                                                                                    "view"))
  message("\nCollecting importances")
  importances <- samples %>% furrr::future_map(function(sample) {
    targets <- contributions.stats %>% dplyr::pull(.data$target) %>% 
      unique() %>% sort()
    views <- contributions.stats %>% dplyr::pull(.data$view) %>% 
      unique()
    maps <- views %>% furrr::future_map(function(view) {
      all.importances <- targets %>% purrr::map(~readr::read_csv(paste0(sample, 
                                                                        .Platform$file.sep, "importances_", .x, "_", 
                                                                        view, ".txt"), col_types = readr::cols()) %>% 
                                                  dplyr::distinct() %>% dplyr::rename(feature = target))
      features <- all.importances %>% purrr::map(~.x$feature) %>% 
        unlist() %>% unique() %>% sort()
      pvalues <- contributions %>% dplyr::filter(sample == 
                                                   !!sample, view == paste0("p.", !!view)) %>% dplyr::mutate(value = 1 - 
                                                                                                               .data$value)
      all.importances %>% purrr::imap_dfc(~tibble::tibble(feature = features, 
                                                          zero.imp = 0) %>% dplyr::left_join(.x, by = "feature") %>% 
                                            dplyr::arrange(.data$feature) %>% dplyr::mutate(imp = scale(.data$imp)[, 
                                                                                                                   1], `:=`(!!targets[.y], .data$zero.imp + (.data$imp * 
                                                                                                                                                               (pvalues %>% dplyr::filter(target == targets[.y]) %>% 
                                                                                                                                                                  dplyr::pull(.data$value))))) %>% dplyr::select(targets[.y])) %>% 
        dplyr::mutate(Predictor = features)
    }) %>% `names<-`(views)
  }, .progress = TRUE) %>% `names<-`(samples)
  message("\nAggregating")
  importances.aggregated <- importances %>% purrr::reduce(function(acc, 
                                                                   l) {
    acc %>% purrr::map2(l, ~(((.x %>% dplyr::select(-.data$Predictor)) + 
                                (.y %>% dplyr::select(-.data$Predictor))) %>% dplyr::mutate(Predictor = .x %>% 
                                                                                              dplyr::pull(.data$Predictor))))
  }) %>% purrr::map(~.x %>% dplyr::mutate_if(is.numeric, ~./length(samples)))
  return(list(improvements = improvements, improvements.stats = improvements.stats, 
              contributions = contributions, contributions.stats = contributions.stats, 
              importances = importances, importances.aggregated = importances.aggregated))
}



# Pipeline definition:定义函数
run_colocalization <- function(slide, 
                               assay, 
                               useful_features, 
                               out_label, 
                               misty_out_alias = "./results/tissue_structure/misty/cell_map/cm_") {
  
  # Define assay of each view ---------------
  view_assays <- list("main" = assay,
                      "juxta" = assay,
                      "para" = assay)
  # Define features of each view ------------
  view_features <- list("main" = useful_features, 
                        "juxta" = useful_features,
                        "para" = useful_features)
  # Define spatial context of each view -----
  view_types <- list("main" = "intra", 
                     "juxta" = "juxta",
                     "para" = "para")
  # Define additional parameters (l in case of paraview,
  # n of neighbors in case of juxta) --------
  view_params <- list("main" = NULL, 
                      "juxta" = 5,
                      "para" = 15)
  
  misty_out <- paste0(misty_out_alias, 
                      out_label, "_", assay)
  
  run_misty_seurat(visium.slide = slide,
                   view.assays = view_assays,
                   view.features = view_features,
                   view.types = view_types,
                   view.params = view_params,
                   spot.ids = NULL,
                   out.alias = misty_out)
  
  return(misty_out)
}









# 多算法基因量化基因集评分 ------------------------------------------------------------
require(future)
require(GSEABase)
require(UCell)
require(irGSEA)
require(GSVA)
require(Seurat)
require(ggplot2)
require(ggpubr)
require(gghalves)
require(ggrepel)
require(RColorBrewer)
require(patchwork)
require(dplyr)
require(tidyr)


#' signature_module
#'
#' @param object Seurat object
#' @param OutPath 设置输出路径
#' @param signature list类型的基因集
#' @param group_by Seurat对象中存放分组信息的表头，例如group
#' @param celltype Seurat对象中存放细胞类型信息的表头，例如cellType_1
#' @param OnlyPlot 是否只画图，复现填True
#'
#' @return
#' @export
#'
#' @examples
signature_module <- function(object,
                             OutPath,
                             signature,
                             group_by = "all",
                             celltype = "cellType_1",
                             OnlyPlot = F) {
  setwd(OutPath)
  if (!dir.exists("Signature"))
    dir.create("Signature")
  
  if (!OnlyPlot) {
    mat <- GetAssayData(object = object, assay = "RNA")
    cells_rankings <- AUCell::AUCell_buildRankings(mat, nCores = 4, plotStats = F)
    score <- AUCell::AUCell_calcAUC(
      signature,
      cells_rankings,
      nCores = 4,
      aucMaxRank = nrow(cells_rankings) * 0.05
    )
    aucs <- as.numeric(AUCell::getAUC(score)[names(signature), ])
    
    # Ucells and singscore
    object <- irGSEA.score(
      object = object,
      assay = "RNA",
      slot = "data",
      seeds = 5201314,
      ncores = 4,
      msigdb = F,
      minGSSize = 0,
      maxGSSize = 100000,
      ucell.MaxRank = 100000,
      custom = T,
      geneset = signature,
      method = c("UCell", "singscore"),
      kcdf = "Gaussian"
    )
    
    uc <- as.data.frame(object@assays$UCell@counts)
    uc <- as.data.frame(t(uc))
    
    singscore <- as.data.frame(object@assays$singscore@counts)
    singscore <- as.data.frame(t(singscore))
    
    ## 3. GSVA
    exp <- as.matrix(object@assays$RNA@data)
    matrix <- gsva(
      exp,
      signature,
      kcdf = "Gaussian",
      method = "ssgsea",
      abs.ranking = T
    )
    ssgsea <- as.data.frame(t(matrix))
    
    # 5. addmodulescore
    object <- AddModuleScore(object, features = signature, name = "Add")
    
    # 组合
    score <- data.frame(
      AUCell = aucs,
      UCell = uc[[names(signature)]],
      singscore = singscore[[names(signature)]],
      ssGSEA = ssgsea[[names(signature)]],
      AddModuleScore = object$Add1
    )
    # scale 标准化
    score <- scale(score)
    
    # 0-1标准化
    normalize <- function(x) {
      return((x - min(x)) / (max(x) - min(x)))
    }
    
    score <- apply(score, 2, normalize)
    score <- as.data.frame(score)
    score$Scoring <- rowSums(score)
    colnames(object@meta.data)
    
    object$Add1 <- NULL
    object@meta.data <- cbind(object@meta.data, score)
    object$ScoreGroup <- ifelse(object$Scoring > median(object$Scoring), "Hexp", "Lexp")
    
    option <- list(signature = names(signature),
                   size = length(signature[[1]]),
                   group_by = group_by)
    save(option, file = "Signature/option.rdata")
    saveRDS(object, "Signature/seurat_signature.rds")
  }
  
  object <- readRDS("Signature/seurat_signature.rds")
  
  p1 <- DotPlot(
    object,
    cols = "RdYlBu",
    group.by = celltype,
    features = c(
      "AUCell",
      "UCell",
      "singscore",
      "ssGSEA",
      "AddModuleScore",
      "Scoring"
    )
  ) +
    theme_bw() +
    scale_size_continuous(range = c(0, 8)) + theme(
      axis.text.x = element_text(
        angle = 45,
        size = 12,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15)
    ) + geom_point(
      shape = 21,
      color = "black",
      stroke = .5,
      mapping = aes_string(size = "pct.exp")
    ) + xlab("Methods") + ylab("Cell Type") + coord_flip()
  ggsave(
    paste0("Signature/1.Dotplot of ", names(signature), ".pdf"),
    plot = p1,
    width = 5.5,
    height = 8
  )
  
  p2 <- FeaturePlot(
    object = object,
    ncol = 1,
    features = "Scoring",
    pt.size = 0.1 * 20000 / ncol(object)
  ) +
    viridis::scale_color_viridis() +
    labs(title = names(signature)) +
    tidydr::theme_dr() + theme(axis.title.y = element_text(size = 15),
                               axis.title.x = element_text(size = 15))
  
  
  if (group_by == "all") {
    tmp <- object@meta.data[, c(celltype, "Scoring")]
    colnames(tmp) <-  c("cellType_1", "val")
    p3 <- ggplot(data = tmp,
                 aes(
                   x = cellType_1,
                   y = val,
                   fill = cellType_1,
                   color = cellType_1
                 )) +
      geom_boxplot(
        width = 0.2,
        fill  = "transparent",
        size  = 0.4,
        outlier.shape = NA
      ) +
      geom_violin(color = NA, alpha = 0.35) +
      # geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
      # geom_half_point_panel(side = "l", shape=21, size=3, color="white") +
      theme_classic() +
      scale_fill_manual(values = c(
        brewer.pal(12, 'Paired'),
        brewer.pal(9, 'Pastel1'),
        brewer.pal(8, 'Pastel2')
      )) +
      scale_color_manual(values = c(
        brewer.pal(12, 'Paired'),
        brewer.pal(9, 'Pastel1'),
        brewer.pal(8, 'Pastel2')
      )) +
      theme(
        legend.position = "right",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 12
        ),
        axis.text.y = element_text(vjust = 0.5, size = 12),
        axis.title.y = element_text(angle = 90, size = 15),
        axis.title.x = element_text(size = 15)
      ) +
      xlab("Cell Type") +
      ylab(paste0(names(signature), " Score"))
    
  } else if (group_by != "all") {
    tmp <- object@meta.data[, c(celltype, "Scoring", group_by)]
    colnames(tmp) <- c("cellType_1", "val", "group")
    tmp$cellType_1 <- as.factor(tmp$cellType_1)
    p3 <- ggplot(data = tmp, aes(
      x = cellType_1,
      y = val,
      fill = group,
      color = group
    )) +
      geom_boxplot(
        show.legend = T,
        width = .6,
        # position = position_dodge(0.9),
        alpha = 0.5,
        outlier.color = 'grey50'
      ) +
      # geom_violin(alpha = 0.5,width=.6,
      #             position = position_dodge(.6),
      #             trim = T,
      #             color = NA) +
      scale_fill_manual(values = c("#398AB9", "red")) +
      scale_color_manual(values = c("#398AB9", "red")) +
      theme_classic() +
      stat_compare_means(
        method = "wilcox.test",
        hjust = 0.5,
        vjust = 0,
        hide.ns = T,
        label = "p.signif",
        show.legend = F
      ) +
      theme(
        legend.position = "top",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 12
        ),
        axis.text.y = element_text(vjust = 0.5, size = 12),
        axis.title.y = element_text(angle = 90, size = 15),
        axis.title.x = element_text(size = 15)
      ) +
      xlab("Cell Type") +
      ylab(paste0(names(signature), " Score"))
  }
  
  p4 <- CombinePlots(plots = list(p2, p3), rel_widths = c(1, 1.5))
  ggsave(
    paste0("Signature/2.Boxplot of ", names(signature), ".pdf"),
    plot = p4,
    width = 12.5,
    height = 5
  )
  
  save(p1, p2, p3, file = "Signature/plot.rdata")
  return(object)
}

#' signature_volcano
#'
#' @param object Seurat object
#' @param ShowCell 挑选要做差异分析的细胞
#' @param group_by Seurat对象中存放分组信息的表头，例如group、ScoreGroup
#' @param celltype Seurat对象中存放细胞类型信息的表头，例如cellType_1
#' @param idents 分组名字，例如c("Hexp","Lexp")、c("Control","Disease")
#' @param logfc.threshold 差异倍数阈值，默认0.25
#' @param OnlyPlot 是否只画图，复现填True
#'
#' @return
#' @export
#'
#' @examples
signature_volcano <- function(object,
                              ShowCell = NA,
                              group_by = "ScoreGroup",
                              celltype = "cellType_1",
                              idents = c("Hexp", "Lexp"),
                              logfc.threshold = 0.25,
                              OnlyPlot = F) {
  if (!OnlyPlot) {
    object <- object[, object@meta.data[[celltype]] %in% ShowCell]
    
    if (group_by == "ScoreGroup") {
      object$ScoreGroup_2 <- ifelse(object$Scoring > median(object$Scoring), "Hexp", "Lexp")
    }
    
    Idents(object) <- group_by
    DEGs <- FindMarkers(
      object,
      assay = "RNA",
      ident.1 = idents[1],
      ident.2 = idents[2],
      logfc.threshold = 0,
      min.pct = 0.1
    )
    sig.DEGs <- DEGs[DEGs$p_val_adj < 0.05 &
                       abs(DEGs$avg_log2FC) > logfc.threshold, ]
    write.csv(
      DEGs,
      file = paste0("Signature/", idents[1], "~", idents[2], "-AllDEGs.csv"),
      row.names = T,
      quote = F
    )
    write.csv(
      sig.DEGs,
      file = paste0(
        "Signature/",
        idents[1],
        "~",
        idents[2],
        "-LogFC",
        logfc.threshold,
        ".csv"
      ),
      row.names = T,
      quote = F
    )
    
    load("Signature/option.rdata")
    option$ShowCell <- ShowCell
    option$idents <- idents
    option$LogFC <- logfc.threshold
    option$numDEGs <- nrow(sig.DEGs)
    save(option, file = "Signature/option.rdata")
  }
  
  DEGs <- read.csv(
    paste0("Signature/", idents[1], "~", idents[2], "-AllDEGs.csv"),
    header = T,
    row.names = 1,
    check.names = F
  )
  
  logFCfilter <- logfc.threshold
  p_valFilter <- 0.05
  
  DEGs <- DEGs[-log10(DEGs$p_val_adj) > 0, ]
  DEGs$gene <- rownames(DEGs)
  DEGs$lab <- ""
  DEGs[order(DEGs$avg_log2FC), ][c(1:10, (nrow(DEGs) - 9):nrow(DEGs)), ]$lab <- DEGs[order(DEGs$avg_log2FC), ][c(1:10, (nrow(DEGs) - 9):nrow(DEGs)), ]$gene
  DEGs$logP <- -log10(DEGs$p_val_adj)
  
  p5 <-
    ggplot(DEGs, aes(
      x = avg_log2FC,
      y = -log10(p_val_adj),
      color = avg_log2FC
    )) +
    geom_point(aes(size = logP), alpha = 0.9) +
    scale_color_gradientn(
      colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
      name = "logFC",
      values = seq(0, 1, 0.2)
    ) +
    scale_fill_gradientn(
      colours = c("#3288bd", "#66c2a5", "#ffffbf", "#f46d43", "#9e0142"),
      name = "logFC",
      values = seq(0, 1, 0.2)
    ) +
    scale_size("-log10(p)") +
    geom_text_repel(
      data = DEGs,
      aes(
        x = avg_log2FC,
        y = -log10(p_val_adj),
        label = lab
      ),
      size = 4,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.8, "lines"),
      segment.color = "black",
      show.legend = F
    ) +
    theme_bw() +
    ylab("-log10 (p_val_adj)") +
    xlab("avg_log2FC") +
    geom_vline(
      xintercept = c(-logFCfilter, logFCfilter),
      lty = 3,
      col = "black",
      lwd = 0.5
    ) +
    geom_hline(
      yintercept = -log10(p_valFilter),
      lty = 3,
      col = "black",
      lwd = 0.5
    ) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    )
  ggsave(
    filename = paste0("Signature/3.Volcano diagram of ", option$signature, ".pdf"),
    plot = p5,
    width = 10,
    height = 8
  )
  
  
  load("Signature/plot.rdata")
  #   layout <- "
  # 11122
  # 11122
  # 33344
  # 33344
  # "
  # ggsave(filename = paste0("Signature/4.Merge.Plot.pdf"), plot = p1 + p2 + p3 + p5 + plot_layout(design = layout), width = 15, height = 10)
  
  ggsave(
    filename = paste0("Signature/4.Merge.Plot.pdf"),
    plot = p1 + p2 + p3 + p5 + plot_layout(ncol = 2, widths = c(2, 1)),
    width = 15,
    height = 10
  )
}

# 单细胞单基因伪样本分析流程函数 --------------------------------------------------------------

create_pseudobulk <- function(
    seurat_obj,
    key_gene,
    min_cells_per_ps = 50,
    n_pseudobulk_target = 100,
    seed = 42,
    chunk_size = 5000  # 每次处理5000个基因
) {
  cat("\n\033[0;34m====== 开始伪样本化分析======\033[0m\n")
  
  # 提取表达矩阵（保持稀疏格式）
  mat <- seurat_obj@assays$RNA@data
  metadata <- seurat_obj@meta.data
  
  cat(sprintf("原始数据: %d genes × %d cells\n", nrow(mat), ncol(mat)))
  cat(sprintf("矩阵格式: %s\n", class(mat)[1]))
  cat(sprintf("稀疏度: %.2f%%\n\n", (1 - length(mat@x) / (nrow(mat) * ncol(mat))) * 100))
  
  # 检查基因存在
  if (!key_gene %in% rownames(mat)) {
    stop(paste0("基因 '", key_gene, "' 未在数据中找到。"))
  }
  
  # === 伪样本化步骤 ===
  cat("\033[0;34m[伪批量化步骤]\033[0m\n")
  set.seed(seed)
  
  n_cells <- ncol(mat)
  max_ps <- max(1, n_cells %/% min_cells_per_ps)
  n_pseudobulk <- min(n_pseudobulk_target, max_ps)
  
  # 随机分配细胞到伪样本
  sample_labels <- sample(1:n_pseudobulk, n_cells, replace = TRUE)
  
  cat(sprintf("✓ 将 %d 个细胞分配到 %d 个伪样本\n", n_cells, n_pseudobulk))
  
  # === 分块处理：避免一次性分配大内存 ===
  cat("✓ 采用分块处理（chunk size: %d 基因）\n\n", chunk_size)
  
  n_genes <- nrow(mat)
  n_chunks <- ceiling(n_genes / chunk_size)
  
  # 预初始化结果列表（不预分配内存）
  pseudobulk_chunks <- list()
  
  for (chunk_idx in 1:n_chunks) {
    # 计算当前chunk的基因范围
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_genes)
    gene_indices <- start_idx:end_idx
    
    cat(sprintf("  处理 Chunk %d/%d (基因 %d-%d)...", 
                chunk_idx, n_chunks, start_idx, end_idx))
    
    # 提取当前chunk的数据（稀疏矩阵）
    mat_chunk <- mat[gene_indices, ]
    
    # 创建当前chunk的伪样本矩阵
    pseudobulk_chunk <- Matrix::Matrix(0, 
                                       nrow = nrow(mat_chunk), 
                                       ncol = n_pseudobulk,
                                       sparse = TRUE)
    
    colnames(pseudobulk_chunk) <- paste0("PS_", 1:n_pseudobulk)
    rownames(pseudobulk_chunk) <- rownames(mat_chunk)
    
    # 逐个伪样本汇聚
    for (i in 1:n_pseudobulk) {
      cells_idx <- which(sample_labels == i)
      if (length(cells_idx) > 0) {
        pseudobulk_chunk[, i] <- Matrix::rowSums(mat_chunk[, cells_idx, drop = FALSE])
      }
    }
    
    # 移除该chunk中的空行
    non_zero_rows <- Matrix::rowSums(pseudobulk_chunk) > 0
    pseudobulk_chunk <- pseudobulk_chunk[non_zero_rows, ]
    
    pseudobulk_chunks[[chunk_idx]] <- pseudobulk_chunk
    
    cat(sprintf(" ✓ (%.1f MB)\n", object.size(pseudobulk_chunk) / 1024^2))
    
    # 主动进行垃圾回收
    gc()
  }
  
  cat("\n✓ 合并所有chunks...\n")
  
  # 合并所有chunks
  pseudobulk_mat <- do.call(rbind, pseudobulk_chunks)
  
  cat(sprintf("✓ 合并后: %d genes × %d samples\n", 
              nrow(pseudobulk_mat), ncol(pseudobulk_mat)))
  
  # 移除全空列
  non_empty_cols <- which(Matrix::colSums(pseudobulk_mat) > 0)
  pseudobulk_mat <- pseudobulk_mat[, non_empty_cols]
  
  cat(sprintf("✓ 随机分组: %d 个伪样本 (平均每个 %.0f 细胞)\n", 
              ncol(pseudobulk_mat), n_cells / ncol(pseudobulk_mat)))
  cat(sprintf("✓ 伪批量化后: %d genes × %d samples\n\n", 
              nrow(pseudobulk_mat), ncol(pseudobulk_mat)))
  
  # === 规范化步骤 ===
  cat("\033[0;34m[规范化步骤]\033[0m\n")
  
  # CPM规范化（保持稀疏格式）
  cat("✓ CPM规范化处理...\n")
  col_sums <- Matrix::colSums(pseudobulk_mat)
  pseudobulk_cpm <- Matrix::t(Matrix::t(pseudobulk_mat) / col_sums) * 1e6
  
  # Log转换
  cat("✓ Log2转换处理...\n")
  pseudobulk_log <- log2(as.matrix(pseudobulk_cpm) + 1)
  
  cat(sprintf("✓ 规范化矩阵维度: %d genes × %d samples\n\n", 
              nrow(pseudobulk_log), ncol(pseudobulk_log)))
  
  # === 在伪样本上基于key_gene分组 ===
  pb_gene_expr <- as.numeric(pseudobulk_log[key_gene, ])
  pb_median_expr <- median(pb_gene_expr)
  
  pb_group_labels <- ifelse(pb_gene_expr > pb_median_expr,
                            paste0(key_gene, "_High"),
                            paste0(key_gene, "_Low"))
  
  cat(sprintf("[1] 在伪样本中检测到 %s 基因表达\n", key_gene))
  cat(sprintf("    中位数表达值: %.2f\n", pb_median_expr))
  cat(sprintf("    High组: %d 个伪样本\n", sum(pb_gene_expr > pb_median_expr)))
  cat(sprintf("    Low组: %d 个伪样本\n", sum(pb_gene_expr <= pb_median_expr)))
  
  # === 创建元数据 ===
  sample_ids <- colnames(pseudobulk_log)
  
  pb_metadata <- data.frame(
    Group = pb_group_labels,
    GeneExpr = pb_gene_expr,
    stringsAsFactors = FALSE,
    row.names = sample_ids
  )
  
  cat("\n✓ 元数据验证:\n")
  cat("  行名与列名是否一致:", all(rownames(pb_metadata) == colnames(pseudobulk_log)), "\n")
  cat("  分组统计:\n")
  print(table(pb_metadata$Group))
  
  # 清理内存
  rm(pseudobulk_chunks, pseudobulk_cpm, pseudobulk_mat)
  gc()
  
  return(list(
    mat_log = pseudobulk_log,        # 规范化后的表达矩阵
    metadata = pb_metadata,
    gene_name = key_gene,
    n_samples = ncol(pseudobulk_log),
    median_expr = pb_median_expr
  ))
}




Step1_GSEA_scRNA_pseudobulk <- function(
    seurat_obj,
    KeyGene = NULL,
    organism = "human",
    Path = "./",
    all_pathway = FALSE,
    mouse = FALSE,
    group_by_var = NULL,
    min_cells_per_ps = 5,
    isKEGG = TRUE,
    top_pathway = 10
) {
  
  library(Seurat)
  library(limma)
  library(clusterProfiler)
  library(enrichplot)
  library(GOplot)
  library(tidygraph)
  library(ggraph)
  library(cowplot)
  library(tidyverse)
  library(Matrix)
  library(ggrepel)
  library(igraph)
  
  # ========== 初始化 ==========
  if (class(seurat_obj)[1] == "Seurat") {
    mat <- GetAssayData(seurat_obj, slot = "data")
    metadata <- seurat_obj@meta.data
  } else {
    mat <- seurat_obj
    metadata <- NULL
  }
  
  # 设置数据库
  if (mouse) {
    library(org.Mm.eg.db)
    OrgDb <- org.Mm.eg.db
    kegg_org <- "mmu"
  } else {
    library(org.Hs.eg.db)
    OrgDb <- org.Hs.eg.db
    kegg_org <- "hsa"
  }
  
  if (is.null(KeyGene)) KeyGene <- rownames(mat)[1:5]
  
  dir <- paste0(Path, "/.GSEA")
  dir.create(dir, showWarnings = F, recursive = T)
  
  cat(sprintf("✓ 原始数据: %d genes × %d cells\n", nrow(mat), ncol(mat)))
  
  # ========== 颜色方案 ==========
  mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", 
             "#088247", "#58CDD9", "#7A142C", "#5D90BA", "#431A3D", 
             "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
  
  # ========== 伪批量化数据加载 ==========
  pseudobulk_mat <- pseudobulk_data$mat_log
  metadata <- pseudobulk_data$metadata
  key_gene <- pseudobulk_data$gene_name
  
  
  # ========== 循环分析基因 ==========
  for (GeneIndex in 1:length(KeyGene)) {
    gene <- KeyGene[GeneIndex]
    cat(sprintf("\n[%d/%d] 分析基因: %s\n", GeneIndex, length(KeyGene), gene))
    
    if (!gene %in% rownames(pseudobulk_mat)) {
      cat(sprintf("  ✗ 基因 %s 不存在\n", gene))
      next
    }
    
    # ========== 分组 ==========
    expr <- as.numeric(pseudobulk_mat[gene, ])
    names(expr) <- colnames(pseudobulk_mat)
    
    threshold <- median(expr)
    low_idx <- expr <= threshold
    high_idx <- !low_idx
    
    n_low <- sum(low_idx)
    n_high <- sum(high_idx)
    
    if (n_low < 2 || n_high < 2) {
      cat(sprintf("  ✗ 分组样本数不足: Low=%d, High=%d (需要≥2)\n", n_low, n_high))
      next
    }
    
    mat_temp <- pseudobulk_mat[, c(which(low_idx), which(high_idx))]
    Type <- c(rep("con", n_low), rep("treat", n_high))
    
    cat(sprintf("  分组: Low(%d) vs High(%d)\n", n_low, n_high))
    
    # ========== 差异分析 ==========
    design <- model.matrix(~ 0 + factor(Type))
    colnames(design) <- c("con", "treat")
    
    mat_temp_dense <- as.matrix(mat_temp)
    keep_genes <- Matrix::rowSums(mat_temp_dense > 0) >= 2
    mat_temp_dense <- mat_temp_dense[keep_genes, ]
    
    fit <- limma::lmFit(mat_temp_dense, design)
    cont.matrix <- limma::makeContrasts(treat - con, levels = design)
    fit2 <- limma::contrasts.fit(fit, cont.matrix)
    fit2 <- limma::eBayes(fit2)
    
    allDiff <- limma::topTable(fit2, adjust = "fdr", number = Inf, 
                               confint = TRUE)
    
    if (nrow(allDiff) == 0) {
      cat("  ✗ 差异分析失败\n")
      next
    }
    
    gsym.fc <- data.frame(
      SYMBOL = rownames(allDiff), 
      logFC = allDiff$logFC,
      padj = allDiff$adj.P.Val,
      stringsAsFactors = FALSE
    )
    
    gene_dir <- file.path(dir, paste0(GeneIndex, ".", gene))
    dir.create(gene_dir, showWarnings = F, recursive = T)
    write.table(allDiff, file.path(gene_dir, paste0(gene, ".xls")), 
                sep = "\t", quote = F)
    
    # ========== ID转换 ==========
    suppressMessages({
      gsym.id <- clusterProfiler::bitr(gsym.fc$SYMBOL, 
                                       fromType = "SYMBOL", 
                                       toType = "ENTREZID", 
                                       OrgDb = OrgDb)
    })
    
    if (nrow(gsym.id) == 0) {
      cat("  ✗ 无有效的ENTREZID映射\n")
      next
    }
    
    gsym.id <- gsym.id[!duplicated(gsym.id$SYMBOL), ]
    gsym.fc.id <- gsym.fc[gsym.fc$SYMBOL %in% gsym.id$SYMBOL, ]
    gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T), ]
    
    id.fc <- setNames(gsym.fc.id.sorted$logFC, 
                      gsym.id$ENTREZID[match(gsym.fc.id.sorted$SYMBOL, gsym.id$SYMBOL)])
    
    cat(sprintf("  ✓ 转换 %d 个ENTREZID\n", length(id.fc)))
    
    # ========== 富集分析 ==========
    cat("  进行富集分析...\n")
    
    tryCatch({
      if (isKEGG) {
        kk <- clusterProfiler::gseKEGG(id.fc, 
                                       organism = kegg_org, 
                                       pvalueCutoff = 1,
                                       eps = 1e-10)
      } else {
        kk <- clusterProfiler::gseGO(id.fc, 
                                     OrgDb = OrgDb, 
                                     pvalueCutoff = 1,
                                     eps = 1e-10)
      }
      
      if (is.null(kk) || nrow(kk) == 0) {
        cat("  ✗ 无显著通路\n")
        next
      }
      
      kk <- clusterProfiler::setReadable(kk, OrgDb, "ENTREZID")
      sortkk <- kk@result[order(kk@result$enrichmentScore, decreasing = T), ]
      
      
      # ✅ 明确选择需要的列（按顺序）
      export_cols <- c(
        "ID",                # 通路ID
        "Description",       # 通路名称
        "setSize",          # 基因集大小
        "enrichmentScore",  # 富集分数（可以是负数）
        "NES",              # 标准化富集分数（可以是负数）
        "pvalue",           # P值（0-1）✅ 这个不能是负数
        "p.adjust",         # 调整后P值（0-1）
        "qvalue",           # Q值（0-1）
        "rank",             # 排名
        "leading_edge",     # 主要边缘基因
        "core_enrichment"   # 核心富集基因
      )
      
      # ✅ 只导出存在的列
      available_cols <- export_cols[export_cols %in% colnames(sortkk)]
      export_df <- sortkk[, available_cols, drop = FALSE]
      
      write.csv(export_df, 
                file.path(gene_dir, "gsea_output.csv"), 
                quote = FALSE, 
                row.names = FALSE,
                na = "")
      
      cat(sprintf("  ✓ %d 条通路\n", nrow(sortkk)))
      
    }, error = function(e) {
      cat(sprintf("  ✗ 富集分析异常: %s\n", e$message))
    })
    
    if (!exists("sortkk") || nrow(sortkk) < 1) {
      next
    }
    
    # ========== 可视化1：Dotplot ==========
    tryCatch({
      library(enrichplot)
      library(ggplot2)
      
      n_top <- min(top_pathway, nrow(sortkk))
      
      # ✅ 方式1：使用enrichplot::dotplot（推荐）
      p_dot <- enrichplot::dotplot(
        kk, 
        showCategory = n_top,
        title = paste0("GSEA: ", gene),
        split = NULL,
        font.size = 12
      ) +
        # ✅ 手动添加颜色映射
        ggplot2::scale_color_gradient(
          low = "#3182BD",      # 蓝色（p值小）
          high = "#E6550D",     # 橙色（p值大）
          name = "pvalue"
        ) +
        # ✅ 调整主题
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.y = ggplot2::element_text(size = 10),
          legend.position = "right"
        )
      
      ggplot2::ggsave(
        file.path(gene_dir, "01_dotplot.pdf"), 
        p_dot, 
        width = 10, 
        height = 6, 
        dpi = 300
      )
      
      cat("  ✓ 点图（有颜色）\n")
    }, error = function(e) {
      cat(sprintf("  ⚠ 点图失败: %s\n", e$message))
    })
    
    # ========== 可视化2：自定义GSEA图（参考原代码风格） ==========
    # ✓ 注意：需要用户选择通路或按 all_pathway 参数处理
    
    cat("  进行自定义GSEA绘图...\n")
    
    if (all_pathway) {
      # ✓ 方案A：绘制所有通路
      geneSetID <- sortkk$ID
      
      Annex_dir <- file.path(gene_dir, "1.Annex")
      dir.create(Annex_dir, showWarnings = F)
      
      # 绘制所有通路的详细图
      for (i in 1:min(nrow(sortkk), 20)) {  # 限制最多20个以避免文件过多
        tryCatch({
          id <- sortkk$ID[i]
          names <- sortkk$Description[i]
          
          p_gsea <- enrichplot::gseaplot2(kk, 
                                          geneSetID = id,
                                          pvalue_table = TRUE,
                                          color = "#228B22",
                                          title = names)
          
          # 清理文件名
          clean_name <- gsub("[/\\\\:|*?\"<>]", "_", names)
          clean_name <- substr(clean_name, 1, 50)
          
          ggplot2::ggsave(file.path(Annex_dir, 
                                    paste0(sprintf("%02d", i), ".", clean_name, ".pdf")),
                          p_gsea, width = 8, height = 6, dpi = 300)
          
        }, error = function(e) {
          cat(sprintf("    ⚠ 通路 %d 绘图失败\n", i))
        })
      }
      
      cat(sprintf("  ✓ 绘制 %d 条通路详细图\n", min(nrow(sortkk), 20)))
      
    } else {
      # ✓ 方案B：用户选择通路或选择前3个
      cat("  \033[0;34m请输入通路ID (逗号分隔，留空则自动选择前3个): \033[0m\n")
      geneSetID_input <- tryCatch({
        readLines(n = 1)
      }, error = function(e) {
        cat("  (自动跳过交互式输入)\n")
        return("")
      })
      
      if (geneSetID_input == "" || is.na(geneSetID_input)) {
        # 自动选择前3个高显著的通路
        geneSetID <- sortkk$ID[1:min(3, nrow(sortkk))]
      } else {
        geneSetID <- trimws(strsplit(geneSetID_input, ",")[[1]])
        geneSetID <- geneSetID[geneSetID %in% sortkk$ID]
      }
      
      if (length(geneSetID) == 0) {
        cat("  ✗ 无有效通路\n")
        next
      }
      
      # ========== 绘制多通路GSEA对比图（参考原代码） ==========
      tryCatch({
        
        # 提取GSEA数据
        geneList <- position <- NULL
        gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = kk))
        gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL, length(geneSetID))
        
        # 获取Description信息
        description.grep <- sortkk[sortkk$ID %in% geneSetID, ]$Description
        
        # ========== 第1层：运行分数曲线 ==========
        p.res <- ggplot(gsdata, aes_(x = ~x)) +
          xlab(NULL) +
          geom_line(aes_(y = ~runningScore, color = ~Description), size = 1.2) +
          scale_color_manual(values = mycol[1:length(unique(gsdata$Description))]) +
          geom_hline(yintercept = 0, lty = "longdash", lwd = 0.3) +
          ylab("Enrichment Score") +
          theme_bw() +
          theme(
            panel.grid = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            legend.background = element_rect(fill = "transparent"),
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank(),
            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, unit = "cm")
          )
        
        # ========== 第2层：基因hits ==========
        p2 <- ggplot(gsdata, aes_(x = ~x)) +
          geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description),
                         size = 1) +
          xlab(NULL) +
          ylab(NULL) +
          scale_color_manual(values = mycol[1:length(unique(gsdata$Description))]) +
          theme_bw() +
          theme(
            panel.grid = element_blank(),
            legend.position = "none",
            plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_blank()
          ) +
          scale_y_continuous(expand = c(0, 0))
        
        # ========== 第3层：秩次指标柱状图 ==========
        df2 <- p.res$data
        df2$y <- p.res$data$geneList[df2$x]
        df2$gsym <- p.res$data$gsym[df2$x]
        
        p.pos <- ggplot(df2, aes(x, y, fill = Description, color = Description)) +
          geom_segment(
            data = df2, aes_(x = ~x, xend = ~x, y = ~y, yend = 0),
            color = "grey70"
          ) +
          geom_bar(position = "dodge", stat = "identity", alpha = 0.8) +
          scale_fill_manual(values = mycol[1:length(unique(df2$Description))],
                            guide = "none") +
          scale_color_manual(values = mycol[1:length(unique(df2$Description))],
                             guide = "none") +
          geom_hline(yintercept = 0, lty = 2, lwd = 0.3) +
          ylab("Ranked List Metric") +
          xlab("Rank in Ordered Dataset") +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = 11, face = "bold"),
            axis.text.x = element_text(size = 11, face = "bold"),
            panel.grid = element_blank(),
            plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, l = 0.2, unit = "cm")
          )
        
        # ========== 组合三层图 ==========
        rel_heights <- c(1.5, 0.5, 1.5)
        plotlist <- list(p.res, p2, p.pos)
        plotlist[[3]] <- plotlist[[3]] +
          theme(
            axis.line.x = element_line(),
            axis.ticks.x = element_line()
          )
        
        p_combined <- cowplot::plot_grid(
          plotlist = plotlist,
          ncol = 1,
          align = "v",
          rel_heights = rel_heights
        )
        
        ggplot2::ggsave(file.path(gene_dir, "02_GSEA_multi_pathways.pdf"),
                        p_combined, width = 10, height = 10, dpi = 300)
        cat("  ✓ 多通路对比GSEA图\n")
        
      }, error = function(e) {
        cat(sprintf("  ⚠ 多通路GSEA图失败: %s\n", e$message))
      })
      
      # ✓ 可视化3：通路-基因互作弦图（修复版）
      if (length(geneSetID) > 1) {
        tryCatch({
          
          sortkk_filtered <- sortkk[sortkk$Description %in% description.grep[1:min(3, length(description.grep))], ]
          
          go <- data.frame(
            Category = if (isKEGG) "KEGG" else "GO",
            ID = sortkk_filtered$ID,
            Term = sortkk_filtered$Description,
            Genes = gsub("/", ", ", sortkk_filtered$core_enrichment),
            adj_pval = sortkk_filtered$p.adjust,
            stringsAsFactors = FALSE
          )
          
          genelist <- data.frame(ID = gsym.fc.id.sorted$SYMBOL, 
                                 logFC = gsym.fc.id.sorted$logFC,
                                 stringsAsFactors = FALSE)
          
          df <- GOplot::circle_dat(go, genelist)[, c(3, 5, 6)]
          
          # ✓ 手动构建nodes（简洁版）
          terms <- unique(df$term)
          genes <- unique(df$genes)
          
          nodes <- data.frame(
            node = c("all", terms, genes),
            node.level = c("root", rep("term", length(terms)), rep("gene", length(genes))),
            node.branch = c("0", 1:length(terms), rep(1:length(terms), length.out = length(genes))),
            node.size = c(1, rep(3, length(terms)), rep(1, length(genes))),
            node.short_name = c("all", substr(terms, 1, 20), substr(genes, 1, 15)),
            stringsAsFactors = FALSE
          )
          nodes$node.level <- as.character(nodes$node.level)
          nodes$node.branch <- as.character(nodes$node.branch)
          rownames(nodes) <- nodes$node
          
          # ✓ 手动构建edges（修复版）
          edges <- df[, c("term", "genes", "logFC")]
          colnames(edges) <- c("from", "to", "logFC")
          
          # 加上root到term的连接
          root_edges <- data.frame(
            from = "all",
            to = terms,
            logFC = 0
          )
          edges <- rbind(root_edges, edges)
          
          graph <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
          
          width <- 600 / (nrow(nodes) + 100)
          
          gc1 <- ggraph::ggraph(graph, layout = "dendrogram", circular = TRUE) +
            ggraph::geom_edge_diagonal(
              aes(
                color = node1.node.branch,
                filter = node1.node.level != "root"
              ),
              alpha = 0.5,
              edge_width = 2.5
            ) +
            scale_edge_color_manual(
              values = c("#61C3ED", "red", "purple", "darkgreen", "orange", "brown")
            ) +
            ggraph::geom_node_point(
              aes(size = node.size, filter = node.level != "root"),
              color = "#61C3ED"
            ) +
            scale_size(range = c(width, 10)) +
            ggraph::geom_node_text(
              aes(
                x = 1.05 * x, y = 1.05 * y,
                label = node.short_name,
                angle = ggraph::node_angle(x, y),
                filter = node.level == "gene"
              ),
              color = "black", size = 2.5, hjust = "outward", vjust = 0.5
            ) +
            ggraph::geom_node_text(
              aes(label = node.short_name, filter = node.level == "term"),
              color = "black", fontface = "bold", size = 6
            ) +
            theme_void() +
            theme(legend.position = "none", panel.background = element_rect(fill = NA)) +
            coord_cartesian(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4))
          
          ggplot2::ggsave(file.path(gene_dir, "02_Ccgraph.pdf"),
                          gc1, width = 18, height = 18, dpi = 300)
          cat("  ✓ 通路-基因互作弦图\n")
          
        }, error = function(e) {
          cat(sprintf("  ⚠ 弦图失败: %s\n", e$message))
        })
      }
    }
  }
  
  cat("\n✅ 分析完成！\n")
  cat(sprintf("📁 结果保存在: %s\n", dir))
  
  invisible(list(pseudobulk_mat = pseudobulk_mat, dir = dir))
}



Step2_GSVA_scRNA <- function(
    seurat_obj,
    KeyGene = NULL,
    organism = "human",
    gene_sets = NULL,
    grouping_var = NULL,
    group_by_var = NULL,
    min_cells_per_ps = 10,
    score_cutoff = 1,
    output_dir = "./GSVA_results"
) {
  
  library(Seurat)
  library(GSVA)
  library(limma)
  library(tidyverse)
  library(stringr)
  library(Matrix)
  library(ComplexHeatmap)
  library(BiocParallel)
  
  # ========== 数据准备 ==========
  if (class(seurat_obj)[1] == "Seurat") {
    mat <- GetAssayData(seurat_obj, slot = "data")
    metadata <- seurat_obj@meta.data
  } else {
    mat <- seurat_obj
    metadata <- NULL
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ========== 基因集加载 ==========
  if (is.null(gene_sets)) {
    message("自动加载基因集...")
    
    if (organism == "mouse") {
      hallmark_file <- "refdata/Mus/Mm.symbols.RData"
      if (!file.exists(hallmark_file)) {
        warning("小鼠基因集文件不存在，请手动加载或修改路径")
        return(invisible(NULL))
      }
      load(hallmark_file)
      
    } else {
      hallmark_file <- "refdata/Homo/hallmark.gs.RData"
      if (!file.exists(hallmark_file)) {
        warning("人类基因集文件不存在，请手动加载或修改路径")
        return(invisible(NULL))
      }
      load(hallmark_file)
    }
  } else {
    gs <- gene_sets
  }
  
  message(paste0("✓ 加载的基因集数量: ", length(gs)))
  
  
  # ========== 加载伪样本数据 ==========
  pseudobulk_log <- pseudobulk_data$mat_log
  metadata <- pseudobulk_data$metadata
  key_gene <- pseudobulk_data$gene_name
  
  cat("矩阵维度:", nrow(pseudobulk_log), "genes ×", ncol(pseudobulk_log), "samples\n")
  cat("分组统计:\n")
  print(table(metadata$Group))
  cat("\n")
  
  # ========== GSVA计算 ==========
  message("\n计算GSVA评分...")
  
  tryCatch({
    message("使用新 GSVA API (gsvaParam)...")
    
    if (!is.list(gs)) {
      gs <- as.list(gs)
    }
    
    gsva_param <- gsvaParam(
      expr = as.matrix(pseudobulk_log),
      gsets = gs
    )
    
    gsva_es <- gsva(
      gsva_param,
      BPPARAM = SerialParam(progressbar = TRUE)
    )
    
    message(paste0("✓ GSVA评分矩阵维度: ", nrow(gsva_es), " 通路 × ", ncol(gsva_es), " 伪批量"))
    
  }, error = function(e) {
    message("新 API 计算失败，尝试兼容旧API...")
    message(paste("错误信息:", e$message))
    
    tryCatch({
      message("使用基础 GSVA 方法...")
      
      gsva_es <<- gsva(
        expr = as.matrix(pseudobulk_log),
        gset.idx.list = gs,
        method = "gsva",
        verbose = FALSE,
        parallel.sz = 1
      )
      
      message(paste0("✓ GSVA评分矩阵维度: ", nrow(gsva_es), " 通路 × ", ncol(gsva_es), " 伪批量"))
      
    }, error = function(e2) {
      message("基础方法也失败，请检查GSVA包版本...")
      message(paste("错误信息:", e2$message))
      stop("无法计算GSVA评分")
    })
  })
  
  # ========== 基因准备 ==========
  if (is.null(KeyGene)) {
    KeyGene <- rownames(mat)[1:10]
  }
  
  # ========== 循环分析每个关键基因 ==========
  for (GeneIndex in 1:length(KeyGene)) {
    i <- KeyGene[GeneIndex]  # 使用变量名 i 保持一致
    
    if (!i %in% rownames(pseudobulk_log)) {
      warning(paste0(i, " 不在表达矩阵中，跳过！"))
      next
    }
    
    message(paste0("\n========== 正在分析基因: ", i, " [", GeneIndex, "/", length(KeyGene), "] =========="))
    
    # ========== 分组方案 ==========
    subexpr <- as.numeric(pseudobulk_log[i, ])
    names(subexpr) <- colnames(pseudobulk_log)
    
    if (!is.null(grouping_var) && !is.null(metadata)) {
      group_info <- metadata[[grouping_var]]
      unique_groups <- unique(group_info)
      
      if (length(unique_groups) < 2) {
        warning("分组数少于2，无法比较！")
        next
      }
      
      lsam <- names(subexpr)[1:length(lsam)]
      hsam <- names(subexpr)[-(1:length(lsam))]
      group_name1 <- unique_groups[1]
      group_name2 <- unique_groups[2]
      
    } else {
      # 按基因表达中位数分组
      lsam <- names(subexpr[subexpr < median(subexpr)])
      hsam <- names(subexpr[subexpr >= median(subexpr)])
      group_name1 <- "Lexp"
      group_name2 <- "Hexp"
    }
    
    message(paste0("按基因表达水平分组: ", group_name1, "(n=", length(lsam), 
                   ") vs ", group_name2, "(n=", length(hsam), ")"))
    
    # 验证样本数
    if (length(lsam) < 2 || length(hsam) < 2) {
      warning(paste0("样本数不足，跳过基因 ", i))
      next
    }
    
    # ========== 设计矩阵构建 ==========
    group_factor <- c(rep(group_name1, length(lsam)), 
                      rep(group_name2, length(hsam)))
    all_samples <- c(lsam, hsam)
    
    design <- model.matrix(~ 0 + factor(group_factor))
    colnames(design) <- c(group_name1, group_name2)
    rownames(design) <- all_samples
    
    message(paste0("设计矩阵维度: ", nrow(design), " × ", ncol(design)))
    
    # ========== 对比矩阵构建 ==========
    contrast_formula <- paste0(group_name2, "-", group_name1)
    message(paste0("对比设置: ", contrast_formula))
    
    contrast.matrix <- makeContrasts(
      contrasts = contrast_formula,
      levels = colnames(design)
    )
    
    # ========== 差异分析 ==========
    gsva_subset <- gsva_es[, all_samples, drop = FALSE]
    
    fit <- lmFit(gsva_subset, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
    
    # ========== 整理结果 ==========
    pathway <- stringr::str_replace(row.names(x), "HALLMARK_", "")
    df <- data.frame(
      ID = pathway,
      score = x$t,
      logFC = x$logFC,
      pvalue = x$P.Value,
      adj.pvalue = x$adj.P.Val,
      stringsAsFactors = FALSE
    )
    
    df$group <- cut(df$score, 
                    breaks = c(-Inf, -score_cutoff, score_cutoff, Inf), 
                    labels = c("Downregulated", "Neutral", "Upregulated"))
    
    # 保存原始 df（用于文本标签）
    sortdf_char <- df[order(df$score), ]
    
    # ========== 保存结果 ==========
    write.table(sortdf_char, 
                file = file.path(output_dir, 
                                 paste0(sprintf("%02d", GeneIndex), ".", i, ".GSVA_result.txt")),
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    # ========== 高级绘图 (集成原始风格) ==========
    
    # 设置 ID 为有序因子用于绘图
    sortdf_plot <- sortdf_char
    sortdf_plot$ID <- factor(sortdf_plot$ID, levels = sortdf_plot$ID)
    
    # 创建绘图数据（包括下面和上面的标签）
    df_lower_text <- subset(sortdf_plot, score < 0)
    df_upper_text <- subset(sortdf_plot, score > 0)
    
    # 主图
    p <- ggplot(sortdf_plot, aes(x = ID, y = score, fill = group)) +
      # 条形图
      geom_bar(stat = "identity", width = 0.7) +
      coord_flip() +
      # 填充颜色
      scale_fill_manual(
        values = c(
          "Downregulated" = "palegreen3", 
          "Neutral" = "snow3", 
          "Upregulated" = "dodgerblue4"
        ), 
        guide = "none"  # 不显示图例
      ) +
      # 阈值线
      geom_hline(
        yintercept = c(-score_cutoff, score_cutoff),
        color = "white",
        linetype = 2,
        linewidth = 0.3  # ✅ 改为 linewidth
      ) +
      # 下方标签（得分 < 0）
      geom_text(
        data = df_lower_text,
        aes(x = ID, y = 0, label = paste0(" ", ID), color = group),
        size = 3,
        hjust = "inward",
        inherit.aes = FALSE
      ) +
      # 上方标签（得分 > 0）
      geom_text(
        data = df_upper_text,
        aes(x = ID, y = -0.05, label = ID, color = group),
        size = 3, 
        hjust = "outward",
        inherit.aes = FALSE
      ) +
      # 文本颜色
      scale_colour_manual(
        values = c(
          "Downregulated" = "black", 
          "Neutral" = "snow3", 
          "Upregulated" = "black"
        ), 
        guide = "none"  # 不显示图例
      ) +
      # 坐标轴标签
      labs(
        x = "",
        y = paste0("t value of GSVA score\n", group_name2, " vs ", group_name1, " group of ", i)
      ) +
      # 主题
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.6),  # ✅ 改为 linewidth
        axis.line.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
      )
    
    # 保存图形
    ggsave(
      file.path(output_dir, paste0(sprintf("%02d", GeneIndex), ".GSVA_plot_", i, ".pdf")),
      plot = p, 
      width = 8, 
      height = 7,
      device = "pdf",
      dpi = 300
    )
    
    message(paste0("✓ 绘图已保存: ", i))
    
    # ========== 热图可视化 ==========
    top_n <- min(5, nrow(sortdf_char))
    top_pathways <- unique(c(
      as.character(head(sortdf_char$ID, top_n)),
      as.character(tail(sortdf_char$ID, top_n))
    ))
    
    valid_pathways <- intersect(top_pathways, rownames(gsva_subset))
    
    if (length(valid_pathways) > 0) {
      gsva_hm <- gsva_subset[valid_pathways, , drop = FALSE]
      
      ha <- HeatmapAnnotation(
        Group = factor(group_factor, levels = c(group_name1, group_name2)),
        col = list(Group = c(
          group_name1 = "lightblue", 
          group_name2 = "lightcoral"
        )),
        show_legend = TRUE
      )
      
      p_heatmap <- Heatmap(
        as.matrix(gsva_hm),
        name = "GSVA\nScore",
        top_annotation = ha,
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        show_column_names = FALSE
      )
      
      pdf(
        file.path(output_dir, paste0(sprintf("%02d", GeneIndex), ".GSVA_heatmap_", i, ".pdf")),
        width = 10, 
        height = 6
      )
      print(p_heatmap)
      dev.off()
      
      message(paste0("✓ 热图已保存: ", i))
    }
    
    message(paste0("✓ 基因 ", i, " 分析完成！"))
  }
  
  message("\n========== GSVA分析全部完成！==========")
  message(paste0("结果保存在: ", output_dir))
  
  return(invisible(list(
    gsva_es = gsva_es,
    pseudobulk_log = pseudobulk_log
  )))
}

Step4_Immune_scRNA_Pseudobulk <- function(
    seurat_obj, 
    KeyGene = NULL, 
    organism = "human", 
    immune_gene_file = NULL, 
    cor_method = "spearman", 
    pvalue_cutoff = 0.05, 
    min_cells_expr = 10, 
    group_by_var = NULL,
    min_cells_per_ps = 30,
    n_pseudobulk = 50,
    output_dir = "./Immune_association") {
  
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  
  start_time <- Sys.time()
  
  # ========== 辅助函数 ==========
  capitalize_first <- function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  }
  
  ggsave_fun <- function(filename, width, height, plot) {
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      width = width,
      height = height,
      device = "pdf",
      dpi = 300
    )
    cat(sprintf("  [OK] Saved: %s.pdf\n", basename(filename)))
  }
  
  # 创建分隔线函数
  print_sep <- function(char = "=", length = 80) {
    cat(paste0(rep(char, length), collapse = ""), "\n")
  }
  
  # ========== 数据准备 ==========
  cat("\n")
  print_sep("=", 80)
  cat("Step4: Single-Cell Immune Regulatory Factor Analysis\n")
  print_sep("=", 80)
  cat("\n")
  
  if (inherits(seurat_obj, "Seurat")) {
    mat <- GetAssayData(seurat_obj, slot = "data")
    metadata <- seurat_obj@meta.data
  } else {
    mat <- seurat_obj
    metadata <- data.frame(sample = "All", row.names = colnames(mat))
  }
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("Matrix class: %s\n", class(mat)))
  cat(sprintf("Original size: %d genes x %d cells\n\n", nrow(mat), ncol(mat)))
  
  # ========== 伪样本化数据加载 ==========
  pseudobulk_mat <- pseudobulk_data$mat_log
  metadata <- pseudobulk_data$metadata
  key_gene <- pseudobulk_data$gene_name
  
  non_empty_cols <- which(Matrix::colSums(pseudobulk_mat) > 0)
  pseudobulk_mat <- pseudobulk_mat[, non_empty_cols]
  
  a_2 <- pseudobulk_mat
  cat(sprintf("[OK] Pseudobulk matrix: %d genes x %d samples\n\n", 
              nrow(a_2), ncol(a_2)))
  
  # ========== 免疫基因库加载 ==========
  cat("[2/6] Loading immune gene database...\n")
  
  if (is.null(immune_gene_file)) {
    immune_file <- if (organism == "mouse") 
      "refdata/Mus/immune.txt" else "refdata/Homo/Immunomodulator_and_chemokines.txt"
    if (!file.exists(immune_file)) {
      warning(sprintf("Immune gene file not found: %s", immune_file))
      return(invisible(NULL))
    }
    b_1 <- read.table(immune_file, header = TRUE, sep = "\t", 
                      fill = TRUE, stringsAsFactors = FALSE)
  } else {
    b_1 <- read.table(immune_gene_file, header = TRUE, sep = "\t", 
                      fill = TRUE, stringsAsFactors = FALSE)
  }
  
  cat(sprintf("[OK] Immune genes loaded: %d\n", nrow(b_1)))
  cat(sprintf("     Types: %s\n\n", paste(unique(b_1$type), collapse = ", ")))
  
  # ========== 关键基因准备 ==========
  cat("[3/6] Preparing key genes...\n")
  
  if (is.null(KeyGene)) {
    KeyGene <- rownames(a_2)[1:10]
    cat("     WARNING: No KeyGene specified, using first 10 genes\n")
  }
  KeyGene <- KeyGene[KeyGene %in% rownames(a_2)]
  
  if (length(KeyGene) == 0) {
    warning("No key genes found in expression matrix")
    return(invisible(NULL))
  }
  
  cat(sprintf("[OK] Key genes: %d\n", length(KeyGene)))
  cat(sprintf("     Genes: %s\n\n", paste(KeyGene, collapse = ", ")))
  
  # ========== 基因过滤 ==========
  cat("[4/6] Filtering genes...\n")
  
  genes_keep <- rownames(a_2)[Matrix::rowSums(a_2 > 0) >= min_cells_expr]
  a_2 <- a_2[genes_keep, ]
  b_1 <- b_1[b_1$Id %in% genes_keep, ]
  a_3 <- data.frame(Id = rownames(a_2), stringsAsFactors = FALSE)
  
  cat(sprintf("[OK] Genes after filtering: %d (min %.0f%% samples expressed)\n", 
              nrow(a_2), min_cells_expr / ncol(a_2) * 100))
  cat(sprintf("     Immune genes available: %d\n\n", nrow(b_1)))
  
  # ========== 内存优化的相关性计算函数 ==========
  compute_sparse_cor <- function(mat, KeyGene_vec, ImmuneGene_vec, cor_method) {
    result <- data.frame(
      row = character(),
      column = character(),
      cor = numeric(),
      pvalue = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (kg in KeyGene_vec) {
      if (!kg %in% rownames(mat)) next
      kg_expr <- as.numeric(mat[kg, ])
      
      for (ig in ImmuneGene_vec) {
        if (!ig %in% rownames(mat)) next
        ig_expr <- as.numeric(mat[ig, ])
        
        valid_idx <- !is.na(kg_expr) & !is.na(ig_expr) & 
          !is.infinite(kg_expr) & !is.infinite(ig_expr)
        
        if (sum(valid_idx) < 3) next
        
        tryCatch({
          ct <- cor.test(kg_expr[valid_idx], ig_expr[valid_idx], 
                         method = cor_method)
          result <- rbind(result, data.frame(
            row = kg,
            column = ig,
            cor = as.numeric(ct$estimate),
            pvalue = ct$p.value,
            stringsAsFactors = FALSE
          ))
        }, error = function(e) NULL)
      }
    }
    
    return(result)
  }
  
  # ========== 参考代码风格的循环分析 ==========
  cat("[5/6] Computing correlations...\n\n")
  
  plotlist <- list()
  result_summary <- list()
  immune_types <- unique(b_1$type)
  
  for (i in seq_along(immune_types)) {
    immune_type <- immune_types[i]
    
    # 提取免疫类型基因
    b_2 <- b_1[b_1$type == immune_type, ]
    data1 <- dplyr::inner_join(b_2, a_3, by = "Id")
    
    if (nrow(data1) == 0) {
      msg1 <- sprintf("[%d/%d] %s: SKIP - no genes found\n", 
                      i, length(immune_types), immune_type)
      cat(msg1)
      next
    }
    
    immune_genes_vec <- data1$Id
    
    # 计算相关性
    result_1 <- compute_sparse_cor(
      mat = a_2,
      KeyGene_vec = KeyGene,
      ImmuneGene_vec = immune_genes_vec,
      cor_method = cor_method
    )
    
    if (nrow(result_1) == 0) {
      msg2 <- sprintf("[%d/%d] %s: SKIP - no correlations\n", 
                      i, length(immune_types), immune_type)
      cat(msg2)
      next
    }
    
    # 整理结果
    result_2 <- result_1[result_1$row %in% KeyGene, ]
    b_2$column <- b_2$Id
    result_3 <- dplyr::inner_join(result_2, b_2, by = "column")
    result1 <- result_3[, 1:4]
    
    # 添加调控方向
    result1$Regulation <- ifelse(result1$cor > 0, "positive", "negative")
    colnames(result1) <- c("gene", "immuneGene", "cor", "pvalue", "Regulation")
    
    # 保存全部结果
    out_all <- sprintf("%s_correlation_all.txt", immune_type)
    out_all_path <- file.path(output_dir, out_all)
    write.table(result1, file = out_all_path,
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    # 添加显著性标记
    data2 <- result1
    data2$pvalue_mark <- ifelse(data2$pvalue < 0.01, "**",
                                ifelse(data2$pvalue < 0.05, "*", ""))
    
    # 排序
    data3 <- data2[order(data2$immuneGene, data2$cor), ]
    
    # 筛选显著结果
    data4 <- data3[data3$pvalue < pvalue_cutoff, ]
    
    if (nrow(data4) == 0) {
      msg3 <- sprintf("[%d/%d] %s: SKIP - no significant results\n", 
                      i, length(immune_types), immune_type)
      cat(msg3)
      next
    }
    
    # 保存显著结果
    out_sig <- sprintf("%s_correlation_sig.txt", immune_type)
    out_sig_path <- file.path(output_dir, out_sig)
    write.table(data4, file = out_sig_path,
                row.names = FALSE, sep = "\t", quote = FALSE)
    
    # 动态高度调整
    n_rows <- nrow(data4)
    if (n_rows > 40) {
      dotheight <- n_rows * 0.18
    } else if (n_rows < 20) {
      dotheight <- 10
    } else {
      dotheight <- 15
    }
    
    data4 <- na.omit(data4)
    
    # 绘制气泡图
    p <- ggplot(data4, aes(x = gene, y = immuneGene)) +
      geom_point(aes(colour = cor, size = -log10(pvalue)), alpha = 0.8) +
      labs(x = "", y = immune_type, 
           title = sprintf("%s Correlation", immune_type)) +
      scale_colour_gradient2(
        low = "blue", 
        high = "red", 
        mid = "white",
        midpoint = 0, 
        limit = c(-1, 1), 
        space = "Lab",
        name = sprintf("%s\nCorrelation", toupper(cor_method))
      ) +
      scale_size_continuous(name = "-log10(p-value)", range = c(2, 8)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.text = element_text(size = 15, colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 15),
            axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 15),
            axis.title = element_text(size = 20, face = "bold"),
            plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
            text = element_text(size = 15),
            legend.position = "right")
    
    # 保存图片
    pdf_base <- sprintf("%s_%s", i, capitalize_first(immune_type))
    pdf_path <- file.path(output_dir, pdf_base)
    ggsave_fun(pdf_path, width = 7, height = dotheight, plot = p)
    
    # 保存到列表
    plotlist[[i]] <- p
    result_summary[[i]] <- data.frame(
      immune_type = immune_type,
      total_genes = length(immune_genes_vec),
      total_correlations = nrow(result1),
      significant = nrow(data4),
      stringsAsFactors = FALSE
    )
    
    # 打印进度
    progress_line <- sprintf("[%d/%d] %s: DONE - %d/%d significant (p < %g)",
                             i, length(immune_types), immune_type,
                             nrow(data4), nrow(result1), pvalue_cutoff)
    cat(progress_line)
    cat("\n")
  }
  
  # ========== 汇总统计 ==========
  cat("[6/6] Summary Statistics\n")
  print_sep("=", 80)
  
  if (length(result_summary) > 0) {
    summary_df <- do.call(rbind, result_summary)
    rownames(summary_df) <- NULL
    
    cat("\nImmune Type Analysis Summar<PATH>/n")
    print(summary_df)
    
    cat("\n")
    cat(sprintf("Total immune types analyzed: %d\n", nrow(summary_df)))
    cat(sprintf("Total key genes: %d\n", length(KeyGene)))
    cat(sprintf("Total significant correlations: %d\n", sum(summary_df$significant)))
    
    if (sum(summary_df$total_correlations) > 0) {
      pct_sig <- sum(summary_df$significant) / sum(summary_df$total_correlations) * 100
      cat(sprintf("Percentage of significant: %.2f%%\n", pct_sig))
    }
  } else {
    cat("WARNING: No significant results found!\n")
  }
  
  # ========== 完成统计 ==========
  elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
  cat("\n")
  print_sep("=", 80)
  cat("[DONE] Analysis Complete!\n")
  cat(sprintf("Results saved to: %s\n", normalizePath(output_dir)))
  cat(sprintf("Time elapsed: %.2f minutes\n", as.numeric(elapsed_time)))
  print_sep("=", 80)
  cat("\n")
  
  # ========== 返回结果 ==========
  return_list <- list(
    pseudobulk_mat = a_2,
    plots = plotlist,
    summary = if (length(result_summary) > 0) do.call(rbind, result_summary) else NULL,
    key_genes = KeyGene
  )
  
  return(invisible(return_list))
}


# 机器学习 --------------------------------------------------------------------
#SVM
# Copyright (C) 2011  John Colby
# htt<PATH>//github.com/johncolby/SVM-RFE

# svmRFE.wrap <- function(test.fold, X, ...) {
#   # Wrapper to run svmRFE function while omitting a given test fold
#   train.data = X[-test.fold, ]
#   test.data  = X[test.fold, ]
#   
#   # Rank the features
#   features.ranked = svmRFE(train.data, ...)
#   
#   return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
# }
# 
# svmRFE <- function(X, k=1, halve.above=5000) {
#   # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
#   n = ncol(X) - 1
#   
#   # Scale data up front so it doesn't have to be redone each pass
#   cat('Scaling data...')
#   X[, -1] = scale(X[, -1])
#   cat('Done!\n')
#   flush.console()
#   
#   pb = txtProgressBar(1, n, 1, style=3)
#   
#   i.surviving = 1:n
#   i.ranked    = n
#   ranked.list = vector(length=n)
#   
#   # Recurse through all the features
#   while(length(i.surviving) > 0) {
#     if(k > 1) {
#       # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
#       folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
#       folds = lapply(1:k, function(x) which(folds == x))
#       
#       # Obtain weights for each training set
#       w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
#       w = do.call(rbind, w)
#       
#       # Normalize each weights vector
#       w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
#       
#       # Compute ranking criteria
#       v    = w * w
#       vbar = apply(v, 2, mean)
#       vsd  = apply(v, 2, sd)
#       c    = vbar / vsd
#     } else {
#       # Only do 1 pass (i.e. regular SVM-RFE)
#       w = getWeights(NULL, X[, c(1, 1+i.surviving)])
#       c = w * w
#     }
#     
#     # Rank the features
#     ranking = sort(c, index.return=T)$ix
#     if(length(i.surviving) == 1) {
#       ranking = 1
#     }
#     
#     if(length(i.surviving) > halve.above) {
#       # Cut features in half until less than halve.above
#       nfeat = length(i.surviving)
#       ncut  = round(nfeat / 2)
#       n     = nfeat - ncut
#       
#       cat('Features halved from', nfeat, 'to', n, '\n')
#       flush.console()
#       
#       pb = txtProgressBar(1, n, 1, style=3)
#       
#     } else ncut = 1
#     
#     # Update feature list
#     ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
#     i.ranked    = i.ranked - ncut
#     i.surviving = i.surviving[-ranking[1:ncut]]
#     
#     setTxtProgressBar(pb, n-length(i.surviving))
#     flush.console()
#   }
#   
#   close(pb)
#   
#   return (ranked.list)
# }
# 
# getWeights <- function(test.fold, X) {
#   # Fit a linear SVM model and obtain feature weights
#   train.data = X
#   if(!is.null(test.fold)) train.data = X[-test.fold, ]
#   
#   svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
#                  scale=F, type="C-classification", kernel="linear")
#   
#   t(svmModel$coefs) %*% svmModel$SV
# }
# 
# WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
#   # Compile feature rankings across multiple folds
#   featureID = sort(apply(sapply(results, function(x) sort(x$feature.ids, index.return=T)$ix), 1, mean), 
#                    index.return=T)$ix
#   avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
#   feature.name = colnames(input[, -1])[featureID]
#   features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
#   if(save==T) {
#     write.table(features.ranked, file=file, quote=F, row.names=F)
#   } else {
#     features.ranked
#   }
# }
# 
# FeatSweep.wrap <- function(i, results, input) {
#   # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
#   svm.list = lapply(results, function(x) tune(svm,
#                                               train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
#                                               train.y      = input[x$train.data.ids, 1],
#                                               validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
#                                               validation.y = input[x$test.data.ids, 1],
#                                               # Optimize SVM hyperparamters
#                                               ranges       = tune(svm,
#                                                                   train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
#                                                                   train.y = input[x$train.data.ids, 1],
#                                                                   ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
#                                               tunecontrol  = tune.control(sampling='fix'))$perf)
#   
#   error = mean(sapply(svm.list, function(x) x$error))
#   return(list(svm.list=svm.list, error=error))
# }


svmRFE.wrap <- function(test.fold, X, ...) {
  train.data = X[-test.fold, ]
  test.data  = X[test.fold, ]
  features.ranked = svmRFE(train.data, ...)
  
  return(list(feature.ids=features.ranked, 
              train.data.ids=row.names(train.data), 
              test.data.ids=row.names(test.data)))
}

svmRFE <- function(X, halve.above=5000) {
  # 简化版：不做内部多折，直接使用单个SVM
  n = ncol(X) - 1
  
  cat('Scaling data...')
  X[, -1] = scale(X[, -1])
  cat('Done!\n')
  flush.console()
  
  pb = txtProgressBar(1, n, 1, style=3)
  
  i.surviving = 1:n
  i.ranked    = n
  ranked.list = vector(length=n)
  
  while(length(i.surviving) > 0) {
    # 单个SVM权重计算
    tryCatch({
      w = getWeights(NULL, X[, c(1, 1+i.surviving)])
      c = w * w
    }, error=function(e) {
      cat("Warning in SVM:", e$message, "\n")
      c <<- rep(0, length(i.surviving))
    })
    
    ranking = sort(c, index.return=T)$ix
    if(length(i.surviving) == 1) {
      ranking = 1
    }
    
    if(length(i.surviving) > halve.above) {
      nfeat = length(i.surviving)
      ncut  = round(nfeat / 2)
      n     = nfeat - ncut
      
      cat('\nFeatures halved from', nfeat, 'to', n, '\n')
      flush.console()
      
      pb = txtProgressBar(1, n, 1, style=3)
      
    } else ncut = 1
    
    ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
    i.ranked    = i.ranked - ncut
    i.surviving = i.surviving[-ranking[1:ncut]]
    
    setTxtProgressBar(pb, n-length(i.surviving))
    flush.console()
  }
  
  close(pb)
  return(ranked.list)
}

getWeights <- function(test.fold, X) {
  train.data = X
  if(!is.null(test.fold)) train.data = X[-test.fold, ]
  
  tryCatch({
    svmModel = svm(train.data[, -1], train.data[, 1], 
                   cost=1, gamma=1/(ncol(train.data)-1),
                   scale=F, type="C-classification", kernel="linear")
    
    if(length(svmModel$coefs) == 0) {
      return(rep(0, ncol(train.data)-1))
    }
    
    weights = t(svmModel$coefs) %*% svmModel$SV
    return(as.vector(weights))
    
  }, error=function(e) {
    cat("Error in getWeights:", e$message, "\n")
    return(rep(0, ncol(train.data)-1))
  })
}

WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
  # 汇总5折结果
  feature.ranks <- sapply(results, function(x) {
    ranked_order <- match(1:(ncol(input)-1), x$feature.ids)
    return(ranked_order)
  })
  
  avg.rank <- rowMeans(feature.ranks, na.rm=T)
  featureID <- order(avg.rank)
  avg.rank <- sort(avg.rank)
  
  feature.name = colnames(input[, -1])[featureID]
  features.ranked = data.frame(FeatureName=feature.name, 
                               FeatureID=featureID, 
                               AvgRank=avg.rank)
  
  if(save==T) {
    write.table(features.ranked, file=file, quote=F, row.names=F)
  } else {
    return(features.ranked)
  }
}

FeatSweep.wrap <- function(i, results, input) {
  svm.list = lapply(results, function(x) {
    tryCatch({
      idx_cols = c(1, 1+x$feature.ids[1:min(i, length(x$feature.ids))])
      
      train.x = input[as.character(x$train.data.ids), idx_cols]
      train.y = input[as.character(x$train.data.ids), 1]
      test.x = input[as.character(x$test.data.ids), idx_cols]
      test.y = input[as.character(x$test.data.ids), 1]
      
      if(nrow(train.x) < 5 | nrow(test.x) < 2 | ncol(train.x) < 2) {
        return(NULL)
      }
      
      svmModel <- svm(train.x[, -1], train.y,
                      cost=1, gamma=1/(ncol(train.x)-1),
                      scale=F, type="C-classification", kernel="linear")
      
      pred <- predict(svmModel, test.x[, -1])
      error <- mean(pred != test.y)
      
      return(list(error=error))
      
    }, error=function(e) {
      return(NULL)
    })
  })
  
  errors <- sapply(svm.list, function(x) ifelse(is.null(x), NA, x$error))
  error = mean(errors, na.rm=T)
  
  return(list(error=error))
}

PlotErrors <- function(errors, errors2=NULL, no.info=0.5, 
                       ylim=range(c(errors, errors2), na.rm=T), 
                       xlab='Number of Features',  ylab='5 x CV Error') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='dodgerblue') {
    lines(which(!is.na(errors)), na.omit(x), col=col,lwd=3)
    points(which.min(x), min(x, na.rm=T), col='firebrick3')
    text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', 
                                              format(min(x, na.rm=T), dig=3)), pos=2, col='red', cex=1.15)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=2)
}


Plotaccuracy <- function(errors, errors2=NULL, no.info=0.5, 
                         ylim=range(c(errors, errors2), na.rm=T), 
                         xlab='Number of Features',  ylab='5 x CV Accuracy') {
  # Makes a plot of average generalization error vs. number of top features
  AddLine <- function(x, col='dodgerblue') {
    lines(which(!is.na(errors)), na.omit(x), col=col,lwd=3)
    points(which.max(x), max(x, na.rm=T), col='firebrick3')
    text(which.max(x), max(x, na.rm=T), paste(which.max(x), '-', 
                                              format(max(x, na.rm=T), dig=3)), pos=2, col='red', cex=1.15)
  }
  
  plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
  AddLine(errors)
  if(!is.null(errors2)) AddLine(errors2, 'gray30')
  abline(h=no.info, lty=2)
}





