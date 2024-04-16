library(Seurat)
library(Signac)
library(patchwork)
library(Matrix)

# do normalization using log, tf-idf, or none, regress out confounds on pca or not 
doBasicSeurat_atac_updated <- function(mtx, npc = 30, top.variable = 5000, 
                                       norm_by = 'tf-idf',
                                       doScale = T, doCenter = T, assay = 'ATAC',
                                       reg.var = 'nCount_ATAC', regressOnPca = T,
                                       project = 'scATAC', 
                                       meta.data = NULL, excludePks.fromVAP = NULL){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-',
                                  meta.data = meta.data)
  
  if(norm_by == 'log') seurat.obj@assays[[assay]]@data <- log1p(mtx) / log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays[[assay]]@data <- TF_IDF(mtx, verbose = F)
  if(norm_by == 'logNormalize') seurat.obj <- NormalizeData(
    object = seurat.obj,
    assay = assay,
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.obj[[paste0('nCount_', assay)]][, 1])
  ) ## treat the data as non-binary
  
  nvap = ifelse(top.variable > 1, top.variable, floor(top.variable * ncol(mtx)))
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  vaps = VariableFeatures(seurat.obj)
  
  #peak.frac = rowMeans(mtx > 0)
  #excludePks.fromVAP = names(which(peak.frac < vap.min.frac))
  
  if(!is.null(excludePks.fromVAP)) vaps = setdiff(vaps, excludePks.fromVAP)
  
  if(length(vaps) < 10) stop('Top few VAPs left!')
  ## redo normalization using vap
  if(norm_by == 'tf-idf'){
    mtx.norm = TF_IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    rm(data0, tmp, mtx.norm)
  }
  
  
  if(regressOnPca){
    reg.var0 = NULL
  }else{
    reg.var0 = reg.var
  }
  VariableFeatures(seurat.obj) <- vaps
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = vaps,
                          vars.to.regress = reg.var0, do.scale = doScale,
                          do.center = doCenter)
  
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = vaps,
                       verbose = FALSE, npc = npc)
  if(length(reg.var) > 0 & regressOnPca) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  return(seurat.obj)
}

# pool and filter variable peaks:  treat data non-binary ####
seurat.atac = readRDS(file = 'Seurat_Objects/seurat_atac_pooled.rds')

## recreate seurat object using the sample only
mtx = seurat.atac@assays$ATAC@counts

rnames = rownames(mtx)
pks.nochrm = rnames[!grepl(rnames, pattern = '^chrM')]
mtx = mtx[pks.nochrm, ]

rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.001))
sampleIDs = seurat.atac$sampleID
freq_by_sample <- sapply(unique(sampleIDs), function(x){
  return(Matrix::rowMeans(mtx[, sampleIDs == x] > 0))
})
freq.max = apply(freq_by_sample, 1, max)
names(freq.max) = rownames(freq_by_sample)
filtered.pks0 = names(which(freq.max < 0.01))

filtered.pks = union(filtered.pks, filtered.pks0)

nvap = 20000
inherited.mdata <- seurat.atac@meta.data
npc = 30
seurat.atac = doBasicSeurat_atac_updated(mtx, npc = npc,
                                         norm_by = 'logNormalize',
                                         top.variable = nvap,
                                         regressOnPca = TRUE,
                                         reg.var = 'nCount_ATAC',
                                         excludePks.fromVAP = filtered.pks,
                                         meta.data = inherited.mdata)

vaps1 = VariableFeatures(seurat.atac)
length(vaps1)

seurat.atac <- RunUMAP(seurat.atac, dims = 1:npc)
p0 <- DimPlot(seurat.atac, group.by = 'sampleID', label = T) + ggtitle('')
p1 <- DimPlot(seurat.atac, group.by = 'Patient_No', label = T) + ggtitle('')
p2 <- DimPlot(seurat.atac, group.by = 'Stage_Code', label = T) + ggtitle('')

seurat.atac <- FindNeighbors(seurat.atac, reduction = 'pca', 
                             dims = 1:npc)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.3)
p3 <- DimPlot(seurat.atac, group.by = 'seurat_clusters', label = T) + ggtitle('')


plotOutFile0 = paste0('Seurat_Objects/seurat_atac_non_binary_vap', nvap, '.pdf')
seuratOutFile0 = paste0('Seurat_Objects/seurat_atac_non_binary_vap', nvap, '.rds')
ggsave(p0 + p1 + p2 + p3 + plot_layout(ncol = 2, heights = c(2, 2)), filename = plotOutFile0, device = 'pdf', 
       width = 14, height = 6)
saveRDS(seurat.atac, seuratOutFile0)

