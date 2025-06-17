
library(Seurat)
library(ggplot2)

#-------------------------------------------------------------------------------------------
# IMMUNE CLUSTERS

integrated <- readRDS("/data/Totaro/merging2/complete_atlas(not_ann).rds")

print("subsetting...")

imm <- subset(integrated, subset = RNA_snn_res.1.3 %in% c(33, 29, 23))

rm(integrated)

gc()

imm <- JoinLayers(imm)

print("scaling...")

imm <- ScaleData(imm, features = rownames(imm), vars.to.regress = "nCount_RNA")

imm <- FindNeighbors(imm, dims = 1:10, reduction = "harmony")

print("clustering...")

imm <- FindClusters(object = imm, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023, 1.2, 1.3, 1.4, 1.5), verbose = TRUE)

imm <- RunUMAP(imm, dims = 1:10)

print("saving immune...")

saveRDS(imm, "/data/Totaro/merging2/immune_atlas_not_annotated.rds")

for(i in c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023, 1.2, 1.3, 1.4, 1.5)){
  
  resolution <- paste0("RNA_snn_res.", i)
  
  p <- DimPlot(object = imm, reduction = "umap", label = TRUE, raster = F, group.by = resolution, shuffle = T) +  NoAxes() + ggtitle(paste0("immune clusters with resolution ", i))
  
  ggplot2::ggsave(paste0("/data/Totaro/merging2/plots/immune/",resolution,".png"), plot = p, width = 15, height = 10, dpi = 300)
  
}

markers <- list(
  
  neutrophil =  c("CXCR2", "IFITM2", "S100A8", "S100A9", "CSF3R", "FCGR3B", "AQP9", "SMCHD1", "CEACAM8", "FUT4"),
  
  t_cell = c("CD3D", "CD3G"),
  
  b_cell = c("MS4A1", "CD79A"),
  
  plasma = c("IGKC", "JCHAIN"),
  
  macrophage = c("FCER1G", "C1QB", "C1QA", "C1QC"),
  
  macrophage_lipid = c("APOE", "GPNMB", "CHIT1", "TM4SF19"),
  
  dendtritic = c("CPVL", "ID2", "CCL17", "CCL22")
  
)

for (i in 1:length(markers)) {
  
  p <- DotPlot(imm2,
               features = rev(markers[[i]]),
               cols = c("lightgray", "blue"),
               dot.scale = 6,  # Controls the size of dots
               group.by = "RNA_snn_res.1.3") +  # Make sure to use the correct column name for cell types
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels for readability
          axis.text.y = element_text(size = 10),  # Adjust y-axis label size (gene names)
          axis.title.x = element_text(size = 12),  # Adjust x-axis title size
          axis.title.y = element_text(size = 12),  # Adjust y-axis title size
          plot.title = element_text(size = 14, face = "bold"),  # Adjust plot title size
          legend.title = element_text(size = 12),  # Adjust legend title size
          legend.text = element_text(size = 10),  # Adjust legend text size
          panel.grid = element_blank()) +  # Optional: remove grid lines
    coord_flip() +  # Flip the coordinates for better space utilization
    labs(title = paste0("Gene Expression of ", names(markers)[i], " Markers Among the Clusters"))
  
  ggplot2::ggsave(paste0("/data/Totaro/merging2/plots/immuneold/res 1.3/dotplot ", names(markers)[i],".png"), plot = p, width = 15, height = 10, dpi = 300)
  
}

clusters <- list(
  
  macro = c(6),
  
  dend = c(13),
  
  t = c(5,3,2,0,1,11),
  
  b = c(8,4,14),
  
  plasma = c(9,15)
  
)

cell_type = c("macrophage", "dendritic cell", "t cell", "b cell", "plasma cell")

imm$general_cell_type = "unknown"

for (i in seq_along(clusters)) {
  
  cluster <- subset(imm, subset = RNA_snn_res.1.2 %in% clusters[[i]])
  
  barcodes <- colnames(cluster)
  
  imm$general_cell_type[barcodes] <- cell_type[i]
  
}

