library(Seurat)
library(ggplot2)
library(ggridges)
`%notin%` = Negate('%in%')
library(pheatmap)
library(stringr)
library(harmony)
library(edgeR)

## use robust zscore, (x-median)/mad
signatureScore_zscore_robust <- function(seurat.obj, pfeatures, nfeatures = NULL,
                                  vars.to.regress = NULL,
                                  score.name = 'quiescent'){
  
  pfeatures = pfeatures[pfeatures %in% rownames(seurat.obj)]
  if(!is.null(nfeatures)) nfeatures = nfeatures[nfeatures %in% rownames(seurat.obj)]
  
  ## regress out
   seurat.obj <- ScaleData(seurat.obj, features = c(pfeatures, nfeatures),
                          do.scale = F, do.center = F, 
                          vars.to.regress = vars.to.regress)
  
  mtx = GetAssayData(seurat.obj, slot = 'scale.data')
  
  ## calculated robust zscore
  rzscore <- function(x){
    sd.robust = sd(x)/2 + mad(x)/2
    x.robust = (x - median(x))/sd.robust
    if(sd.robust <= 0) x.robust = x - median(x)
    return(x.robust)
  }
  
  pscore = mtx[pfeatures, ,drop = F]
  pscore = t(apply(pscore, 1, rzscore))
  pscore = Matrix::colSums(pscore)
  
  if(length(nfeatures) > 0) {
    nscore = mtx[nfeatures, ,drop = F]
    nscore = t(apply(nscore, 1, rzscore))
    nscore = Matrix::colSums(nscore)
    scores = (pscore - nscore)/length(c(pfeatures, nfeatures))
  }else{
    scores = pscore/length(pfeatures)
  }
  seurat.obj <- AddMetaData(seurat.obj, metadata = scores, 
                            col.name = score.name)
  return(seurat.obj)
}

## pool ####
seurat.rna = readRDS('Seurat_Objects/seurat_rna_allCells_clean.rds')

seurat.rna[['RNA']] <- as(seurat.rna[['RNA']], Class = 'Assay5')
seurat.rna = subset(seurat.rna, malignancy == 'Malignant' & 
                      cell_type %in% c('Fibroblast', 'Schwann', 'Adrenergic'))

DefaultAssay(seurat.rna) <- 'RNA'

## gene expr freq by sample
mtx = seurat.rna@assays$RNA@layers$counts
rownames(mtx) = rownames(seurat.rna@assays$RNA@features)
freq_pt = sapply(unique(seurat.rna$Patient_No), function(x){
  return(rowSums(mtx[, seurat.rna$Patient_No == x]))
})
freq_pt = edgeR::cpm(freq_pt, log = F, prior.count = 1)
gini_pt = apply(freq_pt, 1, gini)
pt_specific_genes = names(which(gini_pt > 0.8))

seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = 2500)

sele.features <- VariableFeatures(seurat.rna)

sele.features = setdiff(sele.features, pt_specific_genes)

VariableFeatures(seurat.rna) = sele.features


# Run the standard workflow for visualization and clustering
seurat.rna <- ScaleData(seurat.rna, verbose = FALSE)
seurat.rna <- RunPCA(seurat.rna, npcs = 30, verbose = FALSE)
seurat.rna <- RunUMAP(seurat.rna, reduction = "pca", dims = 1:30)
seurat.rna <- FindNeighbors(seurat.rna, reduction = "pca", dims = 1:30)
seurat.rna <- FindClusters(seurat.rna, resolution = 0.2)


p1 <- DimPlot(seurat.rna, group.by = 'Patient_No', raster = T) + NoLegend()
p2 <- DimPlot(seurat.rna, label = T) + NoLegend()
p3 <- DimPlot(seurat.rna, group.by = 'cell_type', raster = T, label = F)
p4 <- DimPlot(seurat.rna, group.by = 'cell_state', raster = T, label = F)

cl_pt_freq = diag(1/table(seurat.rna$seurat_clusters)) %*% table(seurat.rna$seurat_clusters, seurat.rna$Patient_No)
rownames(cl_pt_freq) = as.character(0:41)
apply(cl_pt_freq, 1, max)

## + harmony ####
if(T){
  seurat.rna = subset(seurat.rna, seurat_clusters %notin% c(33))
  seurat.rna <- RunHarmony(seurat.rna, group.by.vars = 'Patient_No')
  seurat.rna <- RunUMAP(seurat.rna, reduction = 'harmony', dims = 1:20)
  
  seurat.rna <- FindNeighbors(seurat.rna, reduction = "harmony", dims = 1:20)
  seurat.rna <- FindClusters(seurat.rna, resolution = 0.2)
  
  p1 <- DimPlot(seurat.rna, group.by = 'Patient_No', raster = T) + 
    NoLegend()
  p2 <- DimPlot(seurat.rna, label = T) + NoLegend()
  p3 <- DimPlot(seurat.rna, group.by = 'cell_type', raster = T, label = T)
  p4 <- DimPlot(seurat.rna, group.by = 'cell_state', raster = T, label = T)
  
  
  saveRDS(seurat.rna, file = 'Seurat_Objects/seurat_rna_malignant_harmony.rds')
  
}

## further explore ####
seurat.rna = readRDS(file = 'Seurat_Objects/seurat_rna_malignant_harmony.rds')

## remove tiny clusters
seurat.rna = subset(seurat.rna, RNA_snn_res.0.2 %in% c(0:5))
p0 <- DimPlot(seurat.rna, label = T) + scale_color_brewer(palette = 'Set1')
p00 <- DimPlot(seurat.rna, group.by = 'Patient_No') 


adrn_tfs = c('PHOX2B', 'HAND2', 'GATA3', 'ISL1')
mes_tfs = c('SMAD3', 'SNAI2', 'FOSL1', 'FOSL2', 'JUN')
adrn_genes = c( 'TH', 'PHOX2A',  'DBH', 'KLF7', 'ZNF536',
                'DLK1',   'EYA1')
mes_genes = c( 'FN1', 'VIM', 'YAP1', 'CD44',
               'ELK4','CREG1','DCAF6', 'ID1', 'SMAD3', 'SIX4', 'PRRX1', 'WWTR1',
               'COL1A1')
key_tfs = c('PHOX2B', 'HAND2', 'GATA3', 'ISL1', 
            'SIX4', 'SOX11',  'RARA',
            'RXRA','SMAD3', 'SNAI2', 'FOSL1', 'FOSL2', 'JUN')

seurat.rna = signatureScore_zscore_robust(seurat.rna, 
                                          pfeatures = unique(c(mes_tfs, mes_genes)),
                                          score.name = 'MES_Score')

seurat.rna = AddModuleScore(seurat.rna, features = list(unique(c(mes_tfs, mes_genes))),
                            name = 'MES_Module_Score')

seurat.rna = signatureScore_zscore_robust(seurat.rna, 
                                          pfeatures = unique(c(adrn_tfs, adrn_genes)),
                                          score.name = 'ADRN_Score')
seurat.rna = AddModuleScore(seurat.rna, features = list(unique(c(adrn_tfs, adrn_genes))),
                            name = 'ADRN_Module_Score')
p1 <- FeaturePlot(seurat.rna, features = 'MES_Score', 
                  max.cutoff = 'q95') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

p2 <- FeaturePlot(seurat.rna, features = 'ADRN_Score', 
                  max.cutoff = 'q95') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

p3 <- FeaturePlot(seurat.rna, features = 'MES_Module_Score1', 
                  max.cutoff = 'q95') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

p4 <- FeaturePlot(seurat.rna, features = 'ADRN_Module_Score1', 
                  max.cutoff = 'q95') +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

seurat.rna$MES_ADRN_diff = seurat.rna$MES_Score - seurat.rna$ADRN_Score 
seurat.rna$MES_ADRN_diff_module = seurat.rna$MES_Module_Score1 - seurat.rna$ADRN_Module_Score1

p5 <- FeaturePlot(seurat.rna, features = 'MES_ADRN_diff', 
                  max.cutoff = 'q95', raster = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])


p6 <- FeaturePlot(seurat.rna, features = 'MES_ADRN_diff_module', 
                  max.cutoff = 'q95', min.cutoff = -0.5, raster = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])


p7 <- VlnPlot(seurat.rna, features = 'MES_ADRN_diff', pt.size = 0, idents = 0:5, y.max = 4) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')


p8 <- VlnPlot(seurat.rna, features = 'MES_ADRN_diff_module', idents = 0:5, pt.size = 0) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')

p9 <- VlnPlot(seurat.rna, features = 'ADRN_Score', pt.size = 0, idents = 0:5) + 
  geom_hline(yintercept=median(seurat.rna$ADRN_Score), linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')
p10 <- VlnPlot(seurat.rna, features = 'ADRN_Module_Score1', pt.size = 0, idents = 0:5) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')

p11 <- VlnPlot(seurat.rna, features = 'MES_Score', pt.size = 0, idents = 0:5) + 
  geom_hline(yintercept=median(seurat.rna$MES_Score), linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')
p12 <- VlnPlot(seurat.rna, features = 'MES_Module_Score1', pt.size = 0, idents = 0:5) + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')




## cell cycle ####
# Read in the expression matrix The first row is a header row, the first column is rownames
exp.mat <- read.table(file = "MetaData/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurat.rna <- CellCycleScoring(seurat.rna, s.features = s.genes, 
                               g2m.features = g2m.genes, set.ident = F)

p0 <- VlnPlot(seurat.rna, features = 'G2M.Score', pt.size = 0, group.by = 'seurat_clusters') + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')

p1 <- VlnPlot(seurat.rna, features = 'S.Score', pt.size = 0, group.by = 'seurat_clusters') + 
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_brewer(palette = 'Set1')

p2 <- FeaturePlot(seurat.rna, features = 'G2M.Score', 
            max.cutoff = 'q95',  raster = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

p3 <- FeaturePlot(seurat.rna, features = 'S.Score', 
            max.cutoff = 'q95',  raster = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

saveRDS(seurat.rna@meta.data, file = 'MetaData/metadata_seurat_rna_malignant_harmony_clean.rds')
saveRDS(seurat.rna, file = 'Seurat_Objects/seurat_rna_malignant_harmony_clean.rds')
