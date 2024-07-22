library(monocle3)
library(ggplot2)
library(dplyr)
library(Seurat)

dir_work <- '/home/xuj2/beegfs/xuj2/TREM2/code'
setwd(dir_work)


nn_micro_so <- readRDS(file = file.path('./output/micro.rds'))

oc_nuclei <- colnames(nn_micro_so)[which(nn_micro_so$region == "OC")]

nn_micro_oc_so <- nn_micro_so[, oc_nuclei]

gene_annotation <- data.frame(gene_short_name = rownames(nn_micro_oc_so))
rownames(gene_annotation) <- rownames(nn_micro_oc_so)

cell_metadata <- data.frame(barcode = nn_micro_oc_so@assays[["RNA"]]@counts@Dimnames[[2]], batch = nn_micro_oc_so$orig.ident, celltype = nn_micro_oc_so$celltype)
rownames(cell_metadata) <- nn_micro_oc_so@assays[["RNA"]]@counts@Dimnames[[2]]


expression_matrix <- nn_micro_oc_so@assays[["RNA"]]@counts

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

print(str(recreate.partition))
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


list_cluster <- nn_micro_oc_so@meta.data[[sprintf("seurat_clusters")]]
names(list_cluster) <- nn_micro_oc_so@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

### Could be a space-holder, but essentially fills out louvain parameters

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

#### Assign UMAP coordinate

cds_from_seurat@int_colData$reducedDims@listData[["UMAP"]] <- nn_micro_oc_so@reductions[["umap"]]@cell.embeddings



get_earliest_principal_node <- function(cds_from_seurat, celltype="Homeostasis"){
  cell_ids <- which(colData(cds_from_seurat)[, "celltype"] == celltype)

  closest_vertex <- cds_from_seurat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_from_seurat), ])
  root_pr_nodes <- igraph::V(principal_graph(cds_from_seurat)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]

  return(root_pr_nodes)
}


# find all possible partitions
all_partitions <- unique(cds_from_seurat@clusters$UMAP$partitions)
all_partitions <- all_partitions[all_partitions != "1"]
# set all partitions to 1
cds_from_seurat@clusters$UMAP$partitions[cds_from_seurat@clusters$UMAP$partitions %in% all_partitions] <- "1"

cds_from_seurat <- learn_graph(cds_from_seurat)
cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = "UMAP",  root_pr_nodes=get_earliest_principal_node(cds_from_seurat))

celltype <- c(cds_from_seurat@colData@listData$celltype)
pseudotime <- c(cds_from_seurat@principal_graph_aux@listData$UMAP$pseudotime)
cells <- c(colnames(cds_from_seurat))

pdf(file.path('../output_05012023', "Trajectory_157827_pc_micro_new_2023_rasterize.pdf"), width = 6, height = 6)
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           norm_method = "size_only",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.0, rasterize = TRUE)

dev.off()
