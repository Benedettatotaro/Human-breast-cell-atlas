
library(Seurat)
library(readxl)
library(writexl)
library(feature)
library(presto)
library(dplyr)
library(DoubletFinder)

path <- "/data/Totaro/fastq_files_cellranger_stroma/donor_"

path2 <- "/outs/per_sample_outs/donor_"

path3 <- "/count/sample_filtered_feature_bc_matrix"

# Define filtering parameters
filter_params <- list(
  min_UMIs = 500,        # Minimum UMI count for high-quality cells (single-cell)
  min_genes_sc = 200,    # Minimum genes for high-quality cells (single-cell)
  max_UMIs = 20000,      # Maximum UMI count threshold for detecting multiplets
  max_genes = 5000,      # Maximum gene count threshold for detecting multiplets
  max_mito = 10,         # Maximum mitochondrial content (%)
  max_ribo = 50          # Maximum ribosomal content (%)
)

info_stroma <- read_excel("/data/Totaro/merging2/info_stroma.xlsx")

dirs_list_stroma <- c(
  
  donor_1 = paste0(path,1,path2,1,path3),
  donor_2 = paste0(path,2,path2,2,path3),
  donor_3 = paste0(path,3,path2,3,path3),
  donor_4 = paste0(path,4,path2,4,path3),
  donor_5 = paste0(path,5,path2,5,path3),
  donor_6 = paste0(path,6,path2,6,path3),
  donor_7 = paste0(path,7,path2,7,path3),
  donor_8 = paste0(path,8,path2,8,path3),
  donor_9 = paste0(path,9,path2,9,path3),
  donor_10 = paste0(path,10,path2,10,path3),
  
  donor_11 = paste0(path,11,path2,11,path3),
  donor_12 = paste0(path,12,path2,12,path3),
  donor_13 = paste0(path,13,path2,13,path3),
  donor_14 = paste0(path,14,path2,14,path3),
  donor_15 = paste0(path,15,path2,15,path3),
  donor_16 = paste0(path,16,path2,16,path3),
  donor_16_2 = "",
  donor_17 = paste0(path,17,path2,17,path3),
  donor_18 = paste0(path,18,path2,18,path3),
  donor_18_2 = "",
  donor_19 = paste0(path,19,path2,19,path3),
  donor_19_2 = "",
  donor_20 = paste0(path,20,path2,20,path3),
  
  donor_21 = paste0(path,21,path2,21,path3),
  donor_22 = paste0(path,22,"/outs/filtered_feature_bc_matrix"),
  donor_23 = paste0(path,23,path2,23,path3),
  donor_24 = paste0(path,24,path2,24,path3),
  donor_25 = paste0(path,25,path2,25,path3),
  donor_26 = paste0(path,26,path2,26,path3),
  donor_27 = paste0(path,27,path2,27,path3),
  donor_28 = paste0(path,28,"/outs/filtered_feature_bc_matrix"),
  donor_29 = paste0(path,29,path2,29,path3),
  donor_29_2 = "",
  donor_30 = paste0(path,30,path2,30,path3),
  
  donor_31 = paste0(path,31,path2,31,path3),
  donor_32 = paste0(path,32,path2,32,path3),
  donor_32_2 = "",
  donor_33 = paste0(path,33,path2,33,path3),
  donor_33_2 = "",
  donor_34 = paste0(path,34,path2,34,path3),
  donor_35 = paste0(path,35,path2,35,path3),
  donor_36 = paste0(path,36,path2,36,path3),
  donor_37 = paste0(path,37,path2,37,path3),
  donor_38 = paste0(path,38,path2,38,path3),
  donor_38_2 = "",
  donor_39 = paste0(path,39,path2,39,path3),
  donor_40 = paste0(path,40,path2,40,path3),
  
  donor_41 = paste0(path,41,path2,41,path3),
  donor_42 = paste0(path,42,path2,42,path3),
  donor_43 = paste0(path,43,"_1",path2,43,"_1",path3),
  donor_43_2 = paste0(path,43,"_2",path2,43,"_2",path3),
  donor_44 = paste0(path,44,path2,44,path3),
  donor_44_2 = "",
  donor_45 = paste0(path,45,path2,45,path3),
  donor_46 = paste0(path,46,path2,46,path3),
  donor_47 = paste0(path,47,path2,47,path3),
  donor_48 = paste0(path,48,path2,48,path3),
  donor_49 = paste0(path,49,path2,49,path3),
  
  donor_50 = "",
  donor_51 = paste0(path,51,path2,51,path3),
  donor_52 = paste0(path,52,path2,52,path3),
  donor_53 = paste0(path,53,path2,53,path3),
  donor_54 = paste0(path,54,path2,54,path3),
  donor_54_2 = "",
  donor_55 = paste0(path,55,path2,55,path3)
  
)

path <- "/data/Totaro/fastq_files_cellranger_epithelial/donor_"

info_epithelial <- read_excel("/data/Totaro/merging2/info_epithelial.xlsx")

dirs_list_epithelial <- c(
  
  donor_1 = paste0(path,1,path2,1,path3),
  donor_2 = paste0(path,2,path2,2,path3),
  donor_3 = paste0(path,3,path2,3,path3),
  donor_4 = paste0(path,4,path2,4,path3),
  donor_5 = paste0(path,5,path2,5,path3),
  donor_6 = paste0(path,6,path2,6,path3),
  donor_7 = paste0(path,7,path2,7,path3),
  donor_8 = paste0(path,8,path2,8,path3),
  donor_9 = paste0(path,9,path2,9,path3),
  donor_10 = paste0(path,10,path2,10,path3),
  
  donor_11 = paste0(path,11,path2,11,path3),
  donor_12 = paste0(path,12,path2,12,path3),
  donor_13 = paste0(path,13,path2,13,path3),
  donor_14 = paste0(path,14,path2,14,path3),
  donor_15 = paste0(path,15,path2,15,path3),
  donor_16 = paste0(path,16,path2,16,path3),
  donor_16_2 = paste0(path,16,"_2",path2,16,"_2",path3),
  donor_17 = paste0(path,17,path2,17,path3),
  donor_18 = paste0(path,18,path2,18,path3),
  donor_18_2 = paste0(path,18,"_2",path2,18,"_2",path3),
  donor_19 = paste0(path,19,path2,19,path3),
  donor_19_2 = paste0(path,19,"_2",path2,19,"_2",path3),
  donor_20 = paste0(path,20,path2,20,path3),
  
  donor_21 = paste0(path,21,path2,21,path3),
  donor_22 = paste0(path,22,path2,22,path3),
  donor_23 = paste0(path,23,path2,23,path3),
  donor_24 = "",
  donor_25 = paste0(path,25,path2,25,path3),
  donor_26 = paste0(path,26,path2,26,path3),
  donor_27 = paste0(path,27,path2,27,path3),
  donor_28 = paste0(path,28,path2,28,path3),
  donor_29 = paste0(path,29,path2,29,path3),
  donor_29_2 = paste0(path,29,"_2",path2,29,"_2",path3),
  donor_30 = paste0(path,30,path2,30,path3),
  
  donor_31 = paste0(path,31,path2,31,path3),
  donor_32 = paste0(path,32,path2,32,path3),
  donor_32_2 = paste0(path,32,"_2",path2,32,"_2",path3),
  donor_33 = paste0(path,33,path2,33,path3),
  donor_33_2 = paste0(path,33,"_2",path2,33,"_2",path3),
  donor_34 = paste0(path,34,path2,34,path3),
  donor_35 = paste0(path,35,path2,35,path3),
  donor_36 = paste0(path,36,path2,36,path3),
  donor_37 = paste0(path,37,path2,37,path3),
  donor_38 = paste0(path,38,path2,38,path3),
  donor_38_2 = paste0(path,38,"_2",path2,38,"_2",path3),
  donor_39 = paste0(path,39,path2,39,path3),
  donor_40 = paste0(path,40,path2,40,path3),
  
  donor_41 = paste0(path,41,path2,41,path3),
  donor_42 = paste0(path,42,path2,42,path3),
  donor_43 = paste0(path,43,path2,43,path3),
  donor_43_2 = "",
  donor_44 = paste0(path,44,path2,44,path3),
  donor_44_2 = paste0(path,44,"_2",path2,44,"_2",path3),
  donor_45 = paste0(path,45,path2,45,path3),
  donor_46 = paste0(path,46,path2,46,path3),
  donor_47 = paste0(path,47,path2,47,path3),
  donor_48 = paste0(path,48,path2,48,path3),
  donor_49 = paste0(path,49,path2,49,path3),
  
  donor_50 = "",
  donor_51 = "",
  donor_52 = "",
  donor_53 = "",
  donor_54 = paste0(path,54,path2,54,path3),
  donor_54_2 = paste0(path,54,"_2",path2,54,"_2",path3),
  donor_55 = paste0(path,55,path2,55,path3)
  
)

path <- "/data/Totaro/fastq_files_cellranger_lasp/donor_"

info_LASP <- read_excel("/data/Totaro/merging2/info_LASP.xlsx")

dirs_list_LASP <- c(
  
  donor_1 = paste0(path,1,path2,1,path3),
  donor_2 = paste0(path,2,path2,2,path3),
  donor_3 = paste0(path,3,path2,3,path3),
  donor_4 = paste0(path,4,path2,4,path3),
  donor_5 = paste0(path,5,path2,5,path3),
  donor_6 = paste0(path,6,path2,6,path3),
  donor_7 = paste0(path,7,path2,7,path3),
  donor_8 = paste0(path,8,path2,8,path3),
  donor_9 = paste0(path,9,path2,9,path3),
  donor_10 = paste0(path,10,path2,10,path3),
  
  donor_11 = paste0(path,11,path2,11,path3),
  donor_12 = paste0(path,12,path2,12,path3),
  donor_13 = "",
  donor_14 = paste0(path,14,path2,14,path3),
  donor_15 = "",
  donor_16 = paste0(path,16,path2,16,path3),
  donor_16_2 = "",
  donor_17 = paste0(path,17,path2,17,path3),
  donor_18 = paste0(path,18,path2,18,path3),
  donor_18_2 = "",
  donor_19 = "",
  donor_19_2 = "",
  donor_20 = paste0(path,20,path2,20,path3),
  
  donor_21 = "",
  donor_22 = paste0(path,22,path2,22,path3),
  donor_23 = "",
  donor_24 = "",
  donor_25 = paste0(path,25,path2,25,path3),
  donor_26 = paste0(path,26,path2,26,path3),
  donor_27 = paste0(path,27,path2,27,path3),
  donor_28 = "",
  donor_29 = paste0(path,29,path2,29,path3),
  donor_29_2 = "",
  donor_30 = "",
  
  donor_31 = "",
  donor_32 = paste0(path,32,path2,32,path3),
  donor_32_2 = "",
  donor_33 = paste0(path,33,path2,33,path3),
  donor_33_2 = "",
  donor_34 = paste0(path,34,path2,34,path3),
  donor_35 = "",
  donor_36 = "",
  donor_37 = paste0(path,37,path2,37,path3),
  donor_38 = paste0(path,38,path2,38,path3),
  donor_38_2 = paste0(path,38,"_2",path2,38,"_2",path3),
  donor_39 = paste0(path,39,path2,39,path3),
  donor_40 = paste0(path,40,path2,40,path3),
  
  donor_41 = paste0(path,41,path2,41,path3),
  donor_42 = "",
  donor_43 = paste0(path,43,path2,43,path3),
  donor_43_2 = "",
  donor_44 = paste0(path,44,path2,44,path3),
  donor_44_2 = "",
  donor_45 = paste0(path,45,path2,45,path3),
  donor_46 = paste0(path,46,path2,46,path3),
  donor_47 = paste0(path,47,path2,47,path3),
  donor_48 = paste0(path,48,path2,48,path3),
  donor_49 = paste0(path,49,path2,49,path3),
  donor_50 = "",
  
  donor_51 = "",
  donor_52 = "",
  donor_53 = "",
  donor_54 = "",
  donor_54_2 = "",
  donor_55 = ""
  
)

seu_list <- list(
  
  donor_1 = "",
  donor_2 = "",
  donor_3 = "",
  donor_4 = "",
  donor_5 = "",
  donor_6 = "",
  donor_7 = "",
  donor_8 = "",
  donor_9 = "",
  donor_10 = "",
  
  donor_11 = "",
  donor_12 = "",
  donor_13 = "",
  donor_14 = "",
  donor_15 = "",
  donor_16 = "",
  donor_16_2 = "",
  donor_17 = "",
  donor_18 = "",
  donor_18_2 = "",
  donor_19 = "",
  donor_19_2 = "",
  donor_20 = "",
  
  donor_21 = "",
  donor_22 = "",
  donor_23 = "",
  donor_24 = "",
  donor_25 = "",
  donor_26 = "",
  donor_27 = "",
  donor_28 = "",
  donor_29 = "",
  donor_29_2 = "",
  donor_30 = "",
  
  donor_31 = "",
  donor_32 = "",
  donor_32_2 = "",
  donor_33 = "",
  donor_33_2 = "",
  donor_34 = "",
  donor_35 = "",
  donor_36 = "",
  donor_37 = "",
  donor_38 = "",
  donor_38_2 = "",
  donor_39 = "",
  donor_40 = "",
  
  donor_41 = "",
  donor_42 = "",
  donor_43 = "",
  donor_43_2 = "",
  donor_44 = "",
  donor_44_2 = "",
  donor_45 = "",
  donor_46 = "",
  donor_47 = "",
  donor_48 = "",
  donor_49 = "",
  donor_50 = "",
  
  donor_51 = "",
  donor_52 = "",
  donor_53 = "",
  donor_54 = "",
  donor_54_2 = "",
  donor_55 = ""
  
)

no_good <- ""

barcodes_to_keep <- ""

var_regex = "TR[ABDG][VDJ]|^IG[HJKL]|PLCG2"

stress_sig = list(c("JUNB", "JUN", "FOS","FOSB", "DNAJB1", "HSPA1A", "HSPA1B", "CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                    "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                    "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5"))

for (i in 1:65) {
  
  print(paste0("donor ", i))
  
  if(i!=59){
    
    if(dirs_list_stroma[[i]]!=""){
      
      sparse_matrix_stroma <- Seurat::Read10X(data.dir = dirs_list_stroma[[i]])
      colnames(sparse_matrix_stroma) <- gsub("-1$", "", colnames(sparse_matrix_stroma))
      colnames(sparse_matrix_stroma) <- paste0(colnames(sparse_matrix_stroma), "-", info_stroma$sampleID[i])
      
    }
    
    if(dirs_list_epithelial[[i]]!=""){
      
      sparse_matrix_epithelial <- Seurat::Read10X(data.dir = dirs_list_epithelial[[i]])
      colnames(sparse_matrix_epithelial) <- gsub("-1$", "", colnames(sparse_matrix_epithelial))
      colnames(sparse_matrix_epithelial) <- paste0(colnames(sparse_matrix_epithelial), "-", info_epithelial$sampleID[i])
      
    }
    
    if(dirs_list_LASP[[i]]!=""){
      
      sparse_matrix_lasp <- Seurat::Read10X(data.dir = dirs_list_LASP[[i]])
      colnames(sparse_matrix_lasp) <- gsub("-1$", "", colnames(sparse_matrix_lasp))
      colnames(sparse_matrix_lasp) <- paste0(colnames(sparse_matrix_lasp), "-", info_LASP$sampleID[i])
      
    }
    
    stromal_empty <- dirs_list_stroma[[i]] == ""
    epithelial_empty <- dirs_list_epithelial[[i]] == ""
    lasp_empty <- dirs_list_LASP[[i]] == ""
    
    if (!stromal_empty & !epithelial_empty & !lasp_empty) {
      # Case 1: All non-empty      
      if(identical(rownames(sparse_matrix_epithelial), rownames(sparse_matrix_lasp)) && identical(rownames(sparse_matrix_lasp), rownames(sparse_matrix_stroma))){
        
        sparse_matrix <- cbind(sparse_matrix_epithelial, sparse_matrix_lasp, sparse_matrix_stroma)
        
      } else {
        
        no_good <- c(no_good, paste0("donor_", i))
        
      }
      
    } else if (!stromal_empty & !epithelial_empty & lasp_empty) {
      # Case 2: stromal + epithelial non-empty, lasp empty
      if(identical(rownames(sparse_matrix_epithelial), rownames(sparse_matrix_stroma))){
        
        sparse_matrix <- cbind(sparse_matrix_epithelial, sparse_matrix_stroma)
        
      } else {
        
        no_good <- c(no_good, paste0("donor_", i))
        
      }
      
    } else if (!stromal_empty & epithelial_empty & !lasp_empty) {
      # Case 3: stromal + lasp non-empty, epithelial empty
      if(identical(rownames(sparse_matrix_stroma), rownames(sparse_matrix_lasp))){
        
        sparse_matrix <- cbind(sparse_matrix_lasp, sparse_matrix_stroma)
        
      } else {
        
        no_good <- c(no_good, paste0("donor_", i))
        
      }
      
    } else if (stromal_empty & !epithelial_empty & !lasp_empty) {
      # Case 4: epithelial + lasp non-empty, stromal empty
      if(identical(rownames(sparse_matrix_epithelial), rownames(sparse_matrix_lasp))){
        
        sparse_matrix <- cbind(sparse_matrix_epithelial, sparse_matrix_lasp)
        
      } else {
        
        no_good <- c(no_good, paste0("donor_", i))
        
      }
    } else if (!stromal_empty & epithelial_empty & lasp_empty) {
      # Case 5: only stromal non-empty
      
      sparse_matrix <- sparse_matrix_stroma
      
    } else if (stromal_empty & !epithelial_empty & lasp_empty) {
      # Case 6: only epithelial non-empty
      
      sparse_matrix <- sparse_matrix_epithelial
      
    } else if (stromal_empty & epithelial_empty & !lasp_empty) {
      # Case 7: only lasp non-empty
      
      sparse_matrix <- sparse_matrix_lasp
      
    }
    
    seu <- Seurat::CreateSeuratObject(counts = sparse_matrix) 
    
    seu <- Seurat::PercentageFeatureSet(seu, 
                                        pattern = "^MT-", 
                                        col.name = "percent.mito")
    # ribosomal genes
    seu <- Seurat::PercentageFeatureSet(seu, 
                                        pattern = "^RP[SL]",
                                        col.name = "percent.ribo")
    
    seu <- subset(
      seu,
      nFeature_RNA >= filter_params$min_genes_sc &
        nCount_RNA >= filter_params$min_UMIs &
        nFeature_RNA <= filter_params$max_genes &
        nCount_RNA <= filter_params$max_UMIs &
        percent.mito <= filter_params$max_mito &
        percent.ribo <= filter_params$max_ribo
    )
    
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu,
                                selection.method = "vst",
                                nfeatures = 5000)
    
    # Get current variable features
    var_feats <- VariableFeatures(seu, assay = "RNA")
    
    # Remove features matching a regex
    var_feats <- var_feats[!stringr::str_detect(var_feats, var_regex)]
    
    # Remove features in stress_sig (flatten if it's a list of character vectors)
    stress_genes <- unlist(stress_sig)
    var_feats <- var_feats[!var_feats %in% stress_genes]
    
    # Set updated variable features
    VariableFeatures(seu, assay = "RNA") <- var_feats
    
    seu <- ScaleData(seu, features = VariableFeatures(seu))
    
    seu <- RunPCA(seu)
    
    seurat_filtered <- seu
    
    sweep.res.list <- paramSweep(seurat_filtered, PCs = 1:10, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    
    
    nExp_poi <- round(0.02*nrow(seurat_filtered@meta.data))  ## Assuming 2% doublet formation rate
    
    pk <- as.numeric(as.vector(bcmvn$pK)[which.max(bcmvn$BCmetric)])
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    seurat_filtered <- doubletFinder(seurat_filtered, PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi, sct = FALSE)
    
    # Replace with your actual column name from the metadata
    doublet_col <- grep("DF.classifications", colnames(seurat_filtered@meta.data), value = TRUE)
    keep.cells = rownames(seurat_filtered@meta.data[seurat_filtered@meta.data[[doublet_col]]=="Singlet",])
    data.final = seurat_filtered[,keep.cells]
    seurat_filtered_f = subset(seurat_filtered, cells = colnames(data.final))
    
    barcodes_to_keep <- c(barcodes_to_keep, colnames(seurat_filtered_f))
    
    rm(seurat_filtered_f)
    
    seu_list[[i]] <- seu
    
  }
  
}

if(no_good == ""){
  
  print("analysis working...")

  seu_list <- seu_list[-59]
  
  merged <- merge(x = seu_list[[1]], y = seu_list[2:64], collapse = FALSE)
  
  rm(seu_list)
  
  print("scaling data...")
  
  merged <- ScaleData(merged, features = VariableFeatures(merged))
  
  merged <- RunPCA(merged)
  
  print("integrating data...")
  
  integrated <- IntegrateLayers(
    object = merged, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE
  )
  
  rm(merged)
  
  p <- ElbowPlot(integrated)
  
  ggplot2::ggsave("/data/Totaro/merging2/elbowplot_harmony.png", plot = p, width = 15, height = 10, dpi = 300)
  
  integrated <- RunUMAP(object = integrated, reduction = "harmony", dims = 1:10, metric = 'manhattan') 
  
  options(future.globals.maxSize = 600 * 1024^2)
  integrated <- FindNeighbors(object = integrated, dims = 1:10, reduction = "harmony", k.param = 30)
  integrated <- FindClusters(object = integrated, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
  
  integrated$donor_id <- NA
  integrated$sampleID <- NA
  integrated$sample_type <- NA
  integrated$brca_status <- NA
  
  info_complete <- read_excel("/data/Totaro/merging2/info.xlsx")
  
  for (i in 1:nrow(info_complete)) {
    
    if(info_complete$sampleID[i]!="0"){
      
      # Identify cells with the specified suffix
      cells_donor <- grepl(paste0("-", info_complete$sampleID[i],"$"), colnames(integrated))
      
      integrated$donor_id[cells_donor] <- info_complete$donor_id[i]
      integrated$sampleID[cells_donor] <- info_complete$sampleID[i]
      integrated$brca_status[cells_donor] <- info_complete$brca_status[i]
      integrated$sample_type[cells_donor] <- info_complete$sample_type[i]
      
    }
    
  }
  
  barcodes_to_keep <- barcodes_to_keep[-1]
  
  integrated <- integrated[, colnames(integrated) %in% barcodes_to_keep]
  
  cluster_metrics <- integrated@meta.data %>%
    group_by(RNA_snn_res.1.3) %>%
    summarise(
      mean_nCount = mean(nCount_RNA),
      sd_nCount = sd(nCount_RNA),
      mean_mito = mean(percent.mito),
      sd_mito = sd(percent.mito)
    )
  
  # Identify outlier clusters (greater or less than 2 SD from the mean)
  outlier_clusters <- cluster_metrics %>%
    filter(
      mean_nCount > mean(mean_nCount) + 2 * sd(mean_nCount) |
        mean_nCount < mean(mean_nCount) - 2 * sd(mean_nCount) |
        mean_mito > mean(mean_mito) + 2 * sd(mean_mito) |
        mean_mito < mean(mean_mito) - 2 * sd(mean_mito)
    ) %>%
    pull(RNA_snn_res.1.3)
  
  seu <- subset(seu, RNA_snn_res.1.3 %in% outlier_clusters)
  
  print("joining layers of data...")
  
  integrated <- JoinLayers(integrated)
  
  print("find varible genes...")
  
  deg_list <- lapply(unique(seu$RNA_snn_res.1.3), function(cluster) {
    FindMarkers(
      seu,
      ident.1 = cluster,
      only.pos = TRUE  # Focus on upregulated genes in each cluster
    ) %>%
      top_n(15, avg_log2FC) %>%
      mutate(cluster = cluster)
  })
  
  deg_df <- bind_rows(deg_list)
  
  deg_df$gene <- rownames(deg_df)
  
  saveRDS(deg_df, "/data/Totaro/merging2/dge all clusters.rds")
  
  print("saving data...")
  
  saveRDS(integrated, paste0("/data/Totaro/merging2/complete_atlas(not_ann).rds"))
    
} else {
  
  print("the analysis can't be done")
  
  write_xlsx(as.data.frame(no_good), "/data/Totaro/merging2/patients not suitable for merging.xlsx")
  
}

