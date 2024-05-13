source('scDataAnalysis_Utilities_simp.R')
library(Signac)
library(harmony)
library(patchwork)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')


integrate_by = 'signac'


seuratOutFile0 = "Seurat_Objects/seurat_atac_non_binary_vap20000.rds" # pooled version (see integrate_non_binary.R)
seurat.atac = readRDS(seuratOutFile0)
seurat.atac = UpdateSeuratObject(seurat.atac)
seurat.list = SplitObject(seurat.atac, split.by = 'Patient_No')
for(sample0 in names(seurat.list)){
  seurat.list[[sample0]] = RunTFIDF(seurat.list[[sample0]])
  seurat.list[[sample0]] = FindTopFeatures(seurat.list[[sample0]],
                                           min.cutoff = as.integer(0.01 * ncol(seurat.list[[sample0]])))
  seurat.list[[sample0]] = RunSVD(seurat.list[[sample0]])
  
}
#seurat.list[['3463']] <- NULL # too few cells

### merge ####
seurat.merged = merge(seurat.list[[1]], seurat.list[-1])

seurat.merged = FindTopFeatures(seurat.merged, min.cutoff = 100)

## select anchor features
mtx = seurat.merged@assays$ATAC@counts
sIDs = seurat.merged$sampleID
mtx_pbulk <- sapply(unique(sIDs), function(x){
  rowMeans(mtx[, sIDs == x] >0)
})
rmaxs = apply(mtx_pbulk, 1, max)
filtered_peaks1 <- names(which(rmaxs < 0.03))

sele_features = setdiff(VariableFeatures(seurat.merged), filtered_peaks1)
message(paste(length(sele_features), 'selected anchor features!'))

VariableFeatures(seurat.merged) <- sele_features
seurat.merged = RunTFIDF(seurat.merged)
seurat.merged = RunSVD(seurat.merged)
seurat.merged = RunUMAP(seurat.merged, reduction = 'lsi', dims = 2:50)

### integration ####

integration.anchors <- FindIntegrationAnchors(
  object.list = seurat.list,
  anchor.features = 10000,
  reduction = "rlsi", k.anchor = 30,
  dims = 2:50)

# integrate LSI embeddings
seurat.integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = seurat.merged[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50,
  k.weight = 50
)
# create a new UMAP using the integrated embeddings
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "integrated_lsi", 
                             dims = 2:50)

seuratPath = paste0('Seurat_Objects/seurat_atac_signac_byPatient.rds')
plotOutFile = paste0('Seurat_Objects/seurat_atac_signac_byPatient.pdf')

mdata = readRDS('MetaData/metadata_cellTypeAnnotated_seurat_atac_non_binary_vap20000.rds')
seurat.integrated$cell_type0 = mdata[colnames(seurat.integrated),]$cell_type

seurat.integrated = FindNeighbors(seurat.integrated, reduction = 'integrated_lsi', dims = 2:50)
seurat.integrated = FindClusters(seurat.integrated, resolution = 0.2)
seurat.integrated = FindClusters(seurat.integrated, resolution = 0.25)
seurat.integrated = FindClusters(seurat.integrated, resolution = 0.4)


p0 <- DimPlot(seurat.integrated, group.by = 'cell_type0', label = T) + ggtitle('')
p1 <- DimPlot(seurat.integrated, group.by = 'Patient_No', label = T) + ggtitle('')
p2 <- DimPlot(seurat.integrated, group.by = 'Stage_Code', label = T) + ggtitle('')
p3 <- DimPlot(seurat.integrated, group.by = 'seurat_clusters', label = T) + ggtitle('')

ggsave(p0 + p1 + p2 + p3 + plot_layout(ncol = 2, heights = c(2, 2)), 
       filename = plotOutFile, device = 'pdf', 
       width = 14, height = 10)

saveRDS(seurat.integrated, file = seuratPath)



