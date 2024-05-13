library(data.table)
library(magrittr)
library(Seurat)
library(limma)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(GSEABase)
library(stringr)
library(RColorBrewer)
library(viridis)
library(openxlsx)
library(edgeR)
library(Signac)
`%notin%` = Negate(`%in%`)

## within each condition, aggregate data by sample
#filter sample by its fraction and ncell within a condition
aggrBySamplePerCondition <- function(counts, cell_condition, sampleIDs, min_frac = 0.05, 
                                     min_ncell = 20){
  mtx_cl_sample = sample_frac_cl = list()
  for(cl0 in sort(unique(cell_condition))){
    mtx = counts[, cell_condition == cl0]
    sampleIDs0 = sampleIDs[cell_condition == cl0]
    mtx_by_sample <- sapply(sort(unique(sampleIDs0)), function(x) {
      cl_data <- mtx[, sampleIDs0 == x, drop = F]
      return(Matrix::rowSums(cl_data))
    })
    sample_frac_cl <- sapply(sort(unique(sampleIDs0)), function(x) {
      cl_data <- mtx[, sampleIDs0 == x, drop = F]
      return(ncol(cl_data)/ncol(mtx))
    })
    names(sample_frac_cl) = colnames(mtx_by_sample) = paste0(cl0, '_sample', 
                                                             sort(unique(sampleIDs0)))
    sele.cls = which(sample_frac_cl > min_frac & (sample_frac_cl * ncol(mtx)) > min_ncell)
    mtx_cl_sample[[cl0]] = mtx_by_sample[, sele.cls, drop = F]
  }
  
  mtx_cl_sample = do.call('cbind', mtx_cl_sample)
  return(mtx_cl_sample)
}

diff_edgeR <- function(mtx_cl_sample, design.cls = NULL, design.samples = NULL,
                       min_fc = 0.5, max_pvalue = 0.05){
  
  if(is.null(design.cls)) design.cls = sapply(colnames(mtx_cl_sample), function(x) unlist(strsplit(x, '_'))[1])
  if(is.null(design.samples)) design.samples = sapply(colnames(mtx_cl_sample), function(x) unlist(strsplit(x, '_'))[2])
  dap_bulk <- list()
  for(cl0 in unique(design.cls)){
    design.cls0 = ifelse(design.cls == cl0, cl0, 'Others')
    #filter peaks that are not open in half of the samples
    mtx_cl0 = mtx_cl_sample[, design.cls0 == cl0]
    
    dds <- DGEList(counts = mtx_cl_sample, group = factor(design.cls0))
    
    keep <- filterByExpr(y = dds, min.count = 5,
                         min.prop = 0.1, min.total.count = 50)
    dds <- dds[keep, , keep.lib.sizes=FALSE]
    dds <- calcNormFactors(object = dds)
    dds <- estimateDisp(y = dds)
    design <- model.matrix(~ 0+group, data=dds$samples)
    colnames(design) <- levels(dds$samples$group)
    
    fit = glmQLFit(dds, design)
    test = glmQLFTest(fit, contrast=c(1, -1))
    
    res = topTags(object = test, n = "Inf")
    res = res$table
    res$cluster = cl0
    res = data.table(res, keep.rownames = T)
    names(res)[1] = 'peak'
    dap_bulk[[cl0]] = res[PValue < 0.05 & logFC > 0]
  }
  dap_bulk = do.call('rbind', dap_bulk)
  dap_bulk = dap_bulk[logFC > min_fc]
  return(dap_bulk)
}

# pseudo bulk dap btw predicted cell states ####
seurat.atac <- readRDS("Seurat_Objects/signac_obj_adrenergic_signac_modified_integration_v3.rds")
labelTransRes = readRDS('MetaData/labelTransferResult_cca_signac_obj_adrenergic_signac_modified_integration_v2.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = labelTransRes)
seurat.atac$ATAC_snn_res.0.4 = as.character(seurat.atac$ATAC_snn_res.0.4)
seurat.tmp = subset(seurat.atac, prediction.score.max > 0.6)

counts = seurat.tmp@assays$ATAC@counts
cell_clusters = seurat.tmp$predicted.id

mtx_cl_sample <- aggrBySamplePerCondition(counts, cell_clusters, 
                                          sampleIDs = seurat.tmp$sampleID,
                                          min_frac = 0.05, min_ncell = 20)
dir.create('data/intermediate', showWarnings = F, recursive = T)
saveRDS(mtx_cl_sample, "data/intermediate/new_atac_pseudo_bulk_counts_byStates_per_sample.rds")

mtx_cl_sample = readRDS('data/intermediate/new_atac_pseudo_bulk_counts_byStates_per_sample.rds')
design.samples = paste0('s', stringr::str_split_i(colnames(mtx_cl_sample), '_s', 2))
design.cls = stringr::str_split_i(colnames(mtx_cl_sample), '_s', 1)

## filter some peaks by global analysis
peaks.mean.sample <- sapply(unique(seurat.tmp$sampleID), function(x){
  rowMeans(seurat.tmp@assays$ATAC@counts[, seurat.tmp$sampleID == x] > 0)
})

peaks.mean.cl <- sapply(unique(seurat.tmp$cell_state), function(x){
  rowMeans(seurat.tmp@assays$ATAC@counts[, seurat.tmp$cell_state == x] > 0)
})

rmax = apply(peaks.mean.sample, 1,  max)
rmax1 = apply(peaks.mean.cl, 1, max )
sele.peaks = names(which(rmax > 0.01 & rmax1 > 0.001))

mtx_cl_sample = mtx_cl_sample[sele.peaks, ]

## filter cluster that mainly from one sample
#mtx_cl_sample = mtx_cl_sample[, design.cls %notin% c('11', '14', '18')]
dap_bulk = diff_edgeR(mtx_cl_sample, design.cls, design.samples,
                      min_fc = 0.5, max_pvalue = 0.05)
names(dap_bulk)[1] = 'peak'
saveRDS(dap_bulk, file = 'data/intermediate/new_atac_daps_byStates_bulk_edgeR.rds')



