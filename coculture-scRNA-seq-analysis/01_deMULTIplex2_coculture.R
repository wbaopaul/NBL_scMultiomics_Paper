library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
library(DoubletFinder)
library(deMULTIplex2)
library(harmony)
library(RColorBrewer)
library(enrichR)
library(viridis)

rmDoublets <- function(seu.obj, expt.rate = 0.04){
  
  seu.obj <- NormalizeData(seu.obj, scale.factor = median(seu.obj$nCount_RNA)) %>% FindVariableFeatures() %>% 
    ScaleData() %>% RunPCA(npcs = 20, verbose = F) %>% RunUMAP(dims = 1:20, verbose = F)
  seu.obj <- FindNeighbors(seu.obj, dims = 1:20) %>% FindClusters(resolution = 0.2)       
  ## remove doublets
  ## pK Identification (no ground-truth) ##
  sweep.res <- paramSweep(seu.obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  ## Homotypic Doublet Proportion Estimate ##
  annotations = seu.obj$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(expt.rate*nrow(seu.obj@meta.data))  ## Assuming 2.5% doublet formation rate - since we had filtered out many cells already
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ##
  seu.obj <- doubletFinder(seu.obj, PCs = 1:10, pN = 0.25, pK = 0.09, 
                            nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  seu.obj <- doubletFinder(seu.obj, PCs = 1:10, pN = 0.25, pK = 0.09, 
                            nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), 
                            sct = FALSE)
  doublet_var = paste0('DF.classifications_0.25_0.09_', nExp_poi.adj)
  seu.obj[['Doublet_Singlet']] = seu.obj[[doublet_var]]
  
  mdata = seu.obj@meta.data
  mdata[, grep(rownames(mdata), pattern = '0.25_0.09', value = T)] <- NULL
  seu.obj@meta.data = mdata
  DimPlot(seu.obj, group.by = 'Doublet_Singlet')
  seu.obj = subset(seu.obj, Doublet_Singlet == 'Singlet')
  return(seu.obj)
  
}


stateColor = brewer.pal(n = 6, name = "Set1")
names(stateColor) = c('ADRN-Calcium', 'ADRN-Baseline', 'Interm-OxPhos',
                      'ADRN-Dopaminergic', 'ADRN-Proliferating', 'MES')


dir0 = '/mnt/isilon/tan_lab/thadia/Scripts/Neuroblastoma/mkcount/NBL_THP1_Coculture/'
poolNames = dir(dir0)

## demultiplex by deMULTIplex2 ####
for(pool_name in poolNames){
  dir1 = paste0(dir0, pool_name, '/', pool_name, '/outs/')
  
  mtx = Read10X(paste0(dir1, 'filtered_feature_bc_matrix/'))
  
  mtx.umi = mtx[[1]]
  mtx.hto = mtx[[2]]
  
  res <- demultiplexTags(t(mtx.hto), # Required, the tag count matrix from your experiment, can be either dense or sparse
                         plot.path = "MetaData/coculture/tmp", # Where to output a summary plot
                         plot.name = pool_name, # text append to the name of the summary plot file
                         plot.diagnostics = T) # Whether to output diagnostics plots for each tag
  # set plot.umap = "none" to not produce any plots
  table(res$final_assign)
  
  ## update seurat obj
  saveRDS(res$final_assign, file = paste0('MetaData/coculture/deMULTIplex2_demultiplexed_', 
                                          pool_name, '.rds'))
  
  ## create seurat 
  seurat.obj = CreateSeuratObject(mtx.umi)
  seurat.obj$bc = colnames(seurat.obj)
  seurat.obj = subset(seurat.obj, bc %in% names(res$final_assign))
  all(names(res$final_assign) == colnames(seurat.obj))
  
  seurat.obj$deMULTIplex2 = res$final_assign
  
  hto <- CreateAssayObject(mtx.hto[, colnames(seurat.obj)], assay = 'HTO')
  seurat.obj[['HTO']] = hto
  seurat.obj <- NormalizeData(seurat.obj, assay = 'HTO', normalization.method = 'CLR')
  RidgePlot(seurat.obj, assay = 'HTO', rownames(seurat.obj[["HTO"]])[1:2], ncol = 2,
            group.by = 'deMULTIplex2')
  
  seurat.obj[['HTO']] = NULL
  
  seurat.obj <- subset(seurat.obj, deMULTIplex2 %in% c('negative', 'multiplet'), invert = T)
  
  table(seurat.obj$deMULTIplex2)
  
  ## remove doublets
  seurat.obj <- rmDoublets(seurat.obj, expt.rate = 0.04)
  table(seurat.obj$deMULTIplex2)
  
  ## filtering low-quality cells
  seurat.obj[['percent.mt']] = PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  Idents(seurat.obj) = seurat.obj$orig.ident
  VlnPlot(seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3)
  seurat.obj = subset(seurat.obj, percent.mt < 5 & 
                      nFeature_RNA > 500 & nFeature_RNA < 10000 &
                      nCount_RNA > 2000)
  
  seurat.obj <- seurat.obj %>% NormalizeData(scale.factor = median(seurat.obj$nCount_RNA)) %>%
    FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) %>%
    RunUMAP(verbose = F, dims = 1:20)
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:20) %>% FindClusters(resolution = 0.1)
  
  FeaturePlot(seurat.obj, features = c('percent.mt', 'nCount_RNA'))
  FeaturePlot(seurat.obj, features = c('PHOX2B', 'CD68'))
  DimPlot(seurat.obj, group.by = 'deMULTIplex2')
  saveRDS(seurat.obj, file = paste0('Seurat_Objects/coculture/seurat_deMULTIplex2_', 
                                          pool_name, '_Singlets.rds'))
  
}
