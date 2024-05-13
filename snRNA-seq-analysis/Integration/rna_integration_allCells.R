library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)
library(edgeR)
`%notin%` = Negate('%in%')

## rpca ####
seurat.rna = readRDS('Seurat_Objects/Global_seurat_object.rds')
mtx = ceiling(seurat.rna@assays$RNA@counts)
mdata = seurat.rna@meta.data
rcells = names(which(mtx['HBB', ] < 2))
mtx = mtx[, rcells]
mdata = mdata[rcells, ]
mdata$nCount_RNA <- NULL
mdata$nFeature_RNA <- NULL
mdata$nCount_SCT <- NULL
mdata$nFeature_SCT <- NULL

seurat.rna = CreateSeuratObject(counts = mtx, meta.data = mdata)
seurat.rna[['percent.mt']] <- PercentageFeatureSet(seurat.rna, pattern = 'MT-')

seurat.rna = subset(seurat.rna, percent.mt < 10 & nCount_RNA > 1000 &
                      nCount_RNA < 40000)

freq_pt = sapply(unique(seurat.rna$Patient_No), function(x){
       return(rowSums(mtx[, seurat.rna$Patient_No == x]))
})
freq_pt = edgeR::cpm(freq_pt, log = T, prior.count = 1)
gini_pt = apply(freq_pt, 1, gini)
pt_specific_gene = names(which(gini_pt > 0.8))

seurat.list <- SplitObject(seurat.rna, split.by = "Patient_No")
rm(seurat.rna)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2000)
})

sele.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 2000)
rare_genes1 = names(which(rowMeans(mtx > 0) < 0.002))

rare_genes = union(rare_genes1, pt_specific_gene)
sele.features = setdiff(sele.features, rare_genes)

seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, verbose = FALSE, features = sele.features)
  x <- RunPCA(x, verbose = FALSE, npcs = 50, features = sele.features)
})

seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                         anchor.features = sele.features, 
                                         reduction = "rpca", k.anchor = 10, dims = 1:50)
seurat.rna <- IntegrateData(anchorset = seurat.anchors)
DefaultAssay(seurat.rna) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.rna <- ScaleData(seurat.rna, verbose = FALSE)
seurat.rna <- RunPCA(seurat.rna, npcs = 50, verbose = FALSE)
seurat.rna <- RunUMAP(seurat.rna, reduction = "pca", dims = 1:50, 
                      min.dist = 0.2, n.neighbors = 50)
seurat.rna <- FindNeighbors(seurat.rna, reduction = "pca", dims = 1:50)
seurat.rna <- FindClusters(seurat.rna, resolution = 0.2)

p1 <- DimPlot(seurat.rna, group.by = 'Patient_No', raster = T) + NoLegend()
p2 <- DimPlot(seurat.rna, label = T) + NoLegend()
p3 <- DimPlot(seurat.rna, group.by = 'cell_type', raster = T, label = T)
p4 <- DimPlot(seurat.rna, group.by = 'cell_state', raster = T, label = T)

ggsave(p1 + p2 + p3 + p4 + plot_layout(ncol = 2, nrow = 2), 
       filename = 'rna_all_cells_rpca0.pdf', device = 'pdf', 
       width = 14, height = 10)

saveRDS(seurat.rna, file = 'Seurat_Objects/seurat_rna_allCells_rpca.rds')
