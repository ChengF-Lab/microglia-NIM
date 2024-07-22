
library(Seurat)
library(ggplot2)

dir_work <- '/mnt/beegfs/chengflab/xuj2/TL_test/code'
setwd(dir_work)

PATH_I <- '../data'
PATH_O <- '../output/'

tl_so <- readRDS(file = file.path(PATH_I, "integrated_data.rds"))
#
nuclei <- c(read.table(file.path(PATH_O, 'nuclei.txt'), sep = '\t', header = F)$V1)
tl_so <- tl_so[, nuclei]

latent <- read.table(file.path(PATH_O, 'scVI_hvg_3000_correction_all_True_latent_30_lr_0.0005_epoch_100_latent.txt'), sep = '\t', header = F)

latent_mat <- as.matrix(latent)

print(dim(latent_mat))

rownames(latent_mat) <- colnames(tl_so)


tl_so[["scvi"]] <- CreateDimReducObject(embeddings = latent_mat, key = "scvi_", assay = DefaultAssay(tl_so))

tl_so <- FindNeighbors(tl_so, dims = 1:dim(latent_mat)[2], reduction = "scvi")
tl_so <- FindClusters(tl_so, resolution = 0.2)
tl_so <- RunUMAP(tl_so, dims = 1:dim(latent_mat)[2], reduction = "scvi", n.components = 2)


DefaultAssay(tl_so) <- "RNA"
#
pdf(file.path(PATH_O, "3_DimPlot_Cluster#_res02_package.pdf"))
print(DimPlot(tl_so, label = TRUE))
dev.off()
#
BATCH <- "orig.ident"
if (BATCH != FALSE) {
 pdf(file.path(PATH_O, "3_DimPlot_Batch_res02_package.pdf"), width = 15)
 print(DimPlot(tl_so, group.by = BATCH))
 dev.off()
}
#
tl_so <- NormalizeData(tl_so)



pdf(file.path(PATH_O, "4_DotPlot_Marker#_res02_package.pdf"), width = 10, height=20)
print(DotPlot(tl_so, features = MARKERS) + RotatedAxis())
dev.off()




tl_so <- RenameIdents( # Assign cell type. Note: in order to avoid typo, we predefine the cell type names by variables
  tl_so,
  `0` = 'Microglia',
  `1` = 'Oligodendrocytes',
  `2` = 'Astrocytes',
  `3` = 'Astrocytes',
  `4` = 'Excitatory Neuron',
  `5` = 'Endothelial',
  `6` = 'OPC',
  `7` = 'Fibroblasts',
  `8` = 'Inhibitory Neuron',
  `9` = 'Mural',
  `10` = 'Inhibitory Neuron',
  `11` = 'Unknown',
  `12` = 'Lymphocytes',
  `13` = 'Macrophage',
  `14` = 'Excitatory Neuron',
  `15` = 'Excitatory Neuron',
  `16` = 'Astrocytes',
  `17` = 'Excitatory Neuron',
  `18` = 'Excitatory Neuron',
  `19` = 'Excitatory Neuron',
  `20` = 'Excitatory Neuron',
  `21` = 'Excitatory Neuron',
  `22` = 'Oligodendrocytes',
  `23` = 'Astrocytes',
  `24` = 'Unknown',
  `25` = 'Excitatory Neuron',
  `26` = 'Inhibitory Neuron',
  `27` = 'Unknown',
  `28` = 'Unknown',
  `29` = 'Microglia',
  `30` = 'Microglia'
)

tl_so@meta.data$celltype <- tl_so@active.ident



MARKER_HUMAN_Microglia <- c("CSF1R", "C3", "CIITA", "P2RY12", "CX3CR1")
MARKER_HUMAN_Astrocyte <- c("SLC1A2", "ADGRV1", "GPC5", "RYR3", "GFAP")
MARKER_HUMAN_Endothelial <- c("CLDN5", "FLT1", "ABCB1", "EBF1")
MARKER_HUMAN_Excit_neuron <- c("RALYL", "KCNIP4", "CBLN2", "LDB2", "KCNQ5")
MARKER_HUMAN_Inhib_neuron <- c("NXPH1", "LHFPL3", "PCDH15", "GRIK1", "ADARB2")
MARKER_HUMAN_Oligodendrocyte <- c("ST18", "PLP1", "CTPNASA3", "MBP","PIP4K2A")
MARKER_HUMAN_OPC <- c("DSCAM", "PCDH15", "MEGF11", "PDGFRA", "SOX10", "VCAN", "CSPG4", "OLIG1")
MARKER_HUMAN_Macrophage <- c("F13A1", "CD163", "MRC1")
MARKER_HUMAN_Fibroblasts <- c("DCN")
MARKER_HUMAN_Mural <- c("DCN", "PDGFRB")
MARKER_HUMAN_Lymphocytes <- c("IL7R", "SLIT2")
MARKER_HUMAN_Neurons <- c("MAP1B", "RBFOX3")

MARKERS <- c(MARKER_HUMAN_Microglia, MARKER_HUMAN_Astrocyte, MARKER_HUMAN_Endothelial,
            MARKER_HUMAN_Excit_neuron, MARKER_HUMAN_Inhib_neuron,
            MARKER_HUMAN_Oligodendrocyte, MARKER_HUMAN_OPC, MARKER_HUMAN_Macrophage,
           MARKER_HUMAN_Fibroblasts, MARKER_HUMAN_Mural, MARKER_HUMAN_Lymphocytes,
           MARKER_HUMAN_Neurons) # change to markers
MARKERS <- c(unique(MARKERS))



# pdf(file.path(PATH_O, "5_DotPlot_Marker_res02.pdf"), width = 10)
png(file.path(PATH_O, "5_DotPlot_Marker_res02.png"), height = 650, width =400)
print(DotPlot(tl_so, features = MARKERS) + coord_flip() + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)))
dev.off()



png(file.path(PATH_O, "5_DimPlot_all_no_label.png"))
print(DimPlot(tl_so, label = FALSE, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
    'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) +  NoLegend() + theme(
plot.title = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text.x = element_text(size=20),
axis.text.y = element_text(size=20)
))
dev.off()







png(file.path(PATH_O, "5_DimPlot_all_no_label_no_axes.png"))
print(DimPlot(tl_so, label = FALSE, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
    'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Neuron' = '#F4D03F', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) +  NoLegend() + NoAxes()
    +theme(plot.title = element_blank())
)
dev.off()


# png(file.path(PATH_O, "5_DimPlot_all.png"))
# print(DimPlot(tl_so, label = TRUE, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
#     'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Neuron' = '#F4D03F', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) + theme(
# plot.title = element_blank(),
# axis.title.x = element_blank(),
# axis.title.y = element_blank(),
# axis.text.x = element_text(size=20),
# axis.text.y = element_text(size=20)
# ))
# dev.off()


png(file.path(PATH_O, "5_DimPlot_all.png"), width=625, height=625)
print(DimPlot(tl_so, label = TRUE, label.size =7.5, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
    'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) +  NoLegend() + theme(
plot.title = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text.x = element_text(size=20),
axis.text.y = element_text(size=20)
))
dev.off()



pc_so <- tl_so[, which(tl_so$region == "PC")]

# pdf(file.path(PATH_O, "5_DimPlot_CellType_updated_others_08302023_GSE157827.pdf"))
# print(DimPlot(pc_so, label = TRUE, group.by = "celltype", raster = FALSE))
# dev.off()

png(file.path(PATH_O, "5_DimPlot_all_no_label_gse157827.png"))
print(DimPlot(pc_so, label = FALSE, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
    'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) +  NoLegend() + theme(
plot.title = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text.x = element_text(size=20),
axis.text.y = element_text(size=20)
))
dev.off()


ec_nuclei <- colnames(tl_so)[which(tl_so$region == "EC")]
sfg_nuclei <- colnames(tl_so)[which(tl_so$region == "SFG")]
gse147528_nuclei <- c(ec_nuclei, sfg_nuclei)

gse147528_so <- tl_so[, gse147528_nuclei]


png(file.path(PATH_O, "5_DimPlot_all_no_label_gse147528.png"))
print(DimPlot(gse147528_so, label = FALSE, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
    'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) +  NoLegend() + theme(
plot.title = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text.x = element_text(size=20),
axis.text.y = element_text(size=20)
))
dev.off()



oc_nuclei <- colnames(tl_so)[which(tl_so$region == "OC")]
otc_nuclei <- colnames(tl_so)[which(tl_so$region == "OTC")]
gse148822_nuclei <- c(oc_nuclei, otc_nuclei)

gse148822_so <- tl_so[, gse148822_nuclei]



png(file.path(PATH_O, "5_DimPlot_all_no_label_gse148822.png"))
print(DimPlot(gse148822_so, label = FALSE, group.by = "celltype", raster = FALSE,  cols = c('Macrophage'='#C0392B','Microglia'='#F1948A','Astrocytes'='#52BE80','Fibroblasts'='#AAB7B8','Mural'='#808000',
    'Endothelial'='#A569BD','OPC'='#16A085','Oligodendrocytes'='#5499C7','Lymphocytes'='#85C1E9','Excitatory Neuron'='#DC7633', 'Inhibitory Neuron'='#F5B041', 'Unknown' = '#566573')) +  NoLegend() + theme(
plot.title = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
axis.text.x = element_text(size=20),
axis.text.y = element_text(size=20)
))
dev.off()


all_nuclei <- colnames(tl_so)

macrophage <- all_nuclei[which(tl_so$celltype == "Macrophage")] # 8353
astrocytes <- all_nuclei[which(tl_so$celltype == "Astrocytes")] # 127571
fibroblasts <- all_nuclei[which(tl_so$celltype == "Fibroblasts")] # 17571
microglia <- all_nuclei[which(tl_so$celltype == "Microglia")] # 174682
mural <- all_nuclei[which(tl_so$celltype == "Mural")] # 11231
endothelial <- all_nuclei[which(tl_so$celltype == "Endothelial")] # 25588
unknown <- all_nuclei[which(tl_so$celltype == "Unknown")] # 11841
opc <- all_nuclei[which(tl_so$celltype == "OPC")] # 18836
oligo <- all_nuclei[which(tl_so$celltype == "Oligodendrocytes")] # 112769
lymphocytes <- all_nuclei[which(tl_so$celltype == "Lymphocytes")] # 9314
excit_neuron <- all_nuclei[which(tl_so$celltype == "Excitatory Neuron")] # 61296
inhib_neuron <- all_nuclei[which(tl_so$celltype == "Inhibitory Neuron")] # 25827
