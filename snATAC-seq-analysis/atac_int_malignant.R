source('scDataAnalysis_Utilities_simp.R')
library(Signac)
library(harmony)
`%notin%` = Negate('%in%')

nvap = 10000

## start from the pooled version
seuratPath = 'Seurat_Objects/seurat_atac_non_binary_vap20000.rds'
seurat.atac = readRDS(seuratPath)
mdata1 = readRDS('MetaData/metadata_seurat_atac_signac_byPatient.rds') # the annotation based on latest global integration

seurat.atac = AddMetaData(seurat.atac, metadata = subset(mdata1, select = c(cell_type, cell_type0)))

cl_pt_freq = diag(1/table(seurat.atac$seurat_clusters)) %*% table(seurat.atac$seurat_clusters, seurat.atac$Patient_No)
rownames(cl_pt_freq) = as.character(0:27)
apply(cl_pt_freq, 1, max)


seurat.atac = subset(seurat.atac, cell_type %in% c('Adrenergic', 'Schwann', 'Fibroblast') & cell_type0 == 'Adrenergic')
seurat.atac = subset(seurat.atac, Stage_Code %in% c('IDX', 'PTI'))

seurat.list = SplitObject(seurat.atac, split.by = 'Patient_No')
seurat.list <- lapply(seurat.list, function(X){
  #X = RunTFIDF(X)
  X = NormalizeData(X)
  X = FindTopFeatures(X, min.cutoff = as.integer(0.01 * ncol(X)))
  X = RunSVD(X)
  return(X)
})

### merge 
seurat.merged = merge(seurat.list[[1]], seurat.list[-1])

seurat.merged = FindTopFeatures(seurat.merged, min.cutoff = 200)

## select anchor features
mtx = seurat.merged@assays$ATAC@counts
frac_pt <- sapply(unique(seurat.merged$Patient_No), function(x){
  mtx0 = mtx[, seurat.merged$Patient_No == x]
  return(rowMeans(mtx0 > 0))
})

frac_nmax = rowSums(frac_pt > 0.01)
sele.pks1 = names(which(frac_nmax >= 2)) 

sele_features = intersect(VariableFeatures(seurat.merged), sele.pks1)
message(paste(length(sele_features), 'selected anchor features!'))

VariableFeatures(seurat.merged) <- sele_features
seurat.merged = NormalizeData(seurat.merged)
seurat.merged = RunSVD(seurat.merged)

### integration ####

integration.anchors <- FindIntegrationAnchors(
  object.list = seurat.list,
  anchor.features = nvap,
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
DepthCor(seurat.integrated, reduction = 'integrated_lsi')
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "integrated_lsi", 
                             dims = 2:50, n.neighbors = 50, min.dist = 0.1)
seurat.integrated = FindNeighbors(seurat.integrated, dims = 2:30, 
                                  reduction = 'integrated_lsi')

seurat.integrated = FindClusters(seurat.integrated, resolution = 0.3)
seurat.integrated = FindClusters(seurat.integrated, resolution = 0.2)

p0 <- DimPlot(seurat.integrated, group.by = 'sampleID', label = T) + ggtitle('')
p1 <- DimPlot(seurat.integrated, group.by = 'Patient_No', label = T) + ggtitle('')
p2 <- DimPlot(seurat.integrated, group.by = 'Stage_Code', label = T) + ggtitle('')
p3 <- DimPlot(seurat.integrated, group.by = 'seurat_clusters', label = T) + ggtitle('')


p0 <- DimPlot(seurat.integrated, group.by = 'sampleID', label = T) + ggtitle('')
p1 <- DimPlot(seurat.integrated, group.by = 'Patient_No', label = T) + ggtitle('')
p2 <- DimPlot(seurat.integrated, group.by = 'Stage_Code', label = T) + ggtitle('')
p3 <- DimPlot(seurat.integrated, group.by = 'seurat_clusters', label = T) + ggtitle('')

seuratPath = paste0('Seurat_Objects/seurat_atac_adrenergic_signac_modified_Integrated_byPatient.rds')
plotOutFile = paste0('Seurat_Objects/seurat_atac_adrenergic_signac_modified_Integrated_byPatient.pdf')
ggsave(p0 + p1 + p2 + p3 + plot_layout(ncol = 2, heights = c(2, 2)), 
       filename = plotOutFile, device = 'pdf', 
       width = 14, height = 10)
saveRDS(seurat.integrated, file = seuratPath)



