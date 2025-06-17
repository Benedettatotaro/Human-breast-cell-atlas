
library(Seurat)

integrated <- readRDS("/data/Totaro/merging2/complete_atlas(not_ann).rds")

markers <- list(
  
  epithelial = c("EPCAM", "CDH1", "KRT17"),
  lactocyte = c("CSN2", "CSN3", "LALBA"),
  LASP = c("KIT", "ALDH1A3", "SLPI", "SERPINB3", "MUC16"),
  lsh = c("AREG", "ANKRD30A", "FOXA1"),
  bmyo = c("ACTG2", "ACTA2", "TAGLN"),
  
  fb=c("DCN", "LUM"),
  pv=c("PDGFRB", "NOTCH3", "MYH11"),
  ve=c("PLVAP", "PECAM1", "CLDN5"),
  le=c("MMRN1", "PDPN", "CCL21"),
  immune = c("PTPRC", "CD3D", "CD3E", "CD3G", "CD19", "MS4A1", "CD79A", "CEACAM8")
  
)


marker_genes <- unique(unlist(markers))

for (i in 1:length(markers)) {
  
  p <- DotPlot(integrated, 
               features = rev(marker_genes), 
               cols = c("lightgray", "#AD0327"), 
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
  
  ggplot2::ggsave(paste0("/data/Totaro/merging2/plots/res 1.3/dotplot ", "allmarkers",".png"), plot = p, width = 15, height = 7, dpi = 300)
  
}


clusters <- list(
  
  basal = c(17,37,32,8),
  
  lhs = c(36,16,4,0,35,20),
  
  lasp = c(3,6,27),
  
  pv = c(10,18,5),
  
  fb = c(7,14,12,13,9),
  
  stroma = c(21,25,28),
  
  le = 24,
  
  ve = c(31, 30, 11, 15, 26, 1),
  
  bcell = 29,
  
  myeloid = 33,
  
  tcell = 23
  
)

cell_type = c("basal-myoepithelial cell of mammary gland", "luminal hormone-sensing cell of mammary gland",
              "luminal adaptive secretory precursor cell of mammary gland", "perivascular cell",
              "fibroblast of mammary gland", "stromal cell", "endothelial cell of lymphatic vessel",
              "vascular endothelial cell",
              "B cell", "myeloid cell", "T cell")

integrated$cell_type = "unknown"

for (i in seq_along(clusters)) {
  
  cluster <- subset(integrated, subset = RNA_snn_res.1.3 %in% clusters[[i]])
  
  barcodes <- colnames(cluster)
  
  integrated$cell_type[barcodes] <- cell_type[i]
  
}

p <- DimPlot(object = integrated, reduction = "umap", label = TRUE, raster = F, group.by = "cell_type", shuffle = T) +  NoAxes() + ggtitle(paste0("Breast Atlas"))

ggplot2::ggsave(paste0("/data/Totaro/merging2/plots/","cell_type",".png"), plot = p, width = 20, height = 10, dpi = 300)

saveRDS(integrated, "/data/Totaro/merging2/complete_atlas.rds")
