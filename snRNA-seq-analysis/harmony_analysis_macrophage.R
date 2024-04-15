library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(Matrix)
library(matrixStats)
library(RColorBrewer)
library(uwot)
library(harmony)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')
`%notin%` = Negate(`%in%`)

seurat.clean <- readRDS('Seurat_Objects/seurat_rna_allCells_clean.rds')





## Macrophage ####
seurat.sele <- subset(seurat.clean, cell_type == 'Macrophage')

## routine seurat pipeline
DefaultAssay(seurat.sele) <- 'RNA'

seurat.sele = FindVariableFeatures(seurat.sele)
seurat.sele = ScaleData(seurat.sele)
seurat.sele = RunPCA(seurat.sele, verbose = F, npcs = 20)
seurat.sele = RunUMAP(seurat.sele, dims = 1:20)
DimPlot(seurat.sele, group.by = 'malignancy')
DimPlot(seurat.sele, group.by = 'Patient_No')

## harmony
seurat.sele = RunHarmony(seurat.sele, group.by.vars = 'Patient_No')
seurat.sele = RunUMAP(seurat.sele, dims = 1:20, reduction = 'harmony', reduction.name = 'harmonyUMAP')
seurat.sele = FindNeighbors(seurat.sele, reduction = 'harmony', dims = 1:20)
seurat.sele = FindClusters(seurat.sele, resolution = 0.2)
seurat.sele = FindClusters(seurat.sele, resolution = 0.4)

### < remove contaminated clusters and redo ####
## 2,7-malignant, 6,9,11-mes, 10-T_cells
seurat.sele = subset(seurat.sele, seurat_clusters %in% c(0, 1, 3, 4, 5,  8, 12) & 
                       malignancy == 'Nonmalignant')

seurat.sele = FindVariableFeatures(seurat.sele, nfeatures = 1000)
seurat.sele = ScaleData(seurat.sele)
seurat.sele = RunPCA(seurat.sele, verbose = F, npcs = 30)
seurat.sele = RunUMAP(seurat.sele, dims = 1:30)
DimPlot(seurat.sele, group.by = 'Patient_No')

## harmony
seurat.sele = RunHarmony(seurat.sele, group.by.vars = 'Patient_No')
seurat.sele = RunUMAP(seurat.sele, dims = 1:30, reduction = 'harmony', 
                      n.neighbors = 50, min.dist = 0.2)
seurat.sele = FindNeighbors(seurat.sele, reduction = 'harmony', dims = 1:30)
seurat.sele = FindClusters(seurat.sele, resolution = 0.1)
seurat.sele = FindClusters(seurat.sele, resolution = 0.2)


p1 <- DimPlot(seurat.sele, group.by = 'Patient_No')
p2 <- DimPlot(seurat.sele,  label = T) + scale_color_brewer(palette = 'Paired')

p3 <- StackedVlnPlot(seurat.sele, features = c('SPP1', 'C1QC', 'HBEGF', 'CD68'), 
                     group.by = 'seurat_clusters', 
                     myColors = brewer.pal(7, "Paired")) 

p5 <-StackedVlnPlot(seurat.sele, features = c('PHOX2B', 'ISL1', 'DBH'), 
              group.by = 'seurat_clusters',  myColors = brewer.pal(7, "Paired")) 
p6 <- StackedVlnPlot(seurat.sele, features = c('WWTR1', 'YAP1', 'DCN'), 
              group.by = 'seurat_clusters',  myColors = brewer.pal(7, "Paired")) 
p7 <- StackedVlnPlot(seurat.sele, features = c('CD68', 'CD163', 'ITGAM', 'CD86', 'MRC1', 'MSR1'), 
              group.by = 'seurat_clusters',  myColors = brewer.pal(7, "Paired")) 

degs_macro = FindAllMarkers(seurat.sele, test.use = 'LR', 
                      latent.vars = c('nCount_RNA', 'percent.mt'),
                      max.cells.per.ident = 500, only.pos = T)
degs_macro = data.table(degs_macro)


### < further remove ALB+ cluster6  and redo ####
seurat.sele = subset(seurat.sele, seurat_clusters %in% c(0:5))

seurat.sele = FindVariableFeatures(seurat.sele, nfeatures = 1000)
seurat.sele = ScaleData(seurat.sele)
seurat.sele = RunPCA(seurat.sele, verbose = F, npcs = 20)
seurat.sele = RunUMAP(seurat.sele, dims = 1:20)
DimPlot(seurat.sele, group.by = 'Patient_No')

## harmony
seurat.sele = RunHarmony(seurat.sele, group.by.vars = 'Patient_No', dims.use = 1:20)
seurat.sele = RunUMAP(seurat.sele, dims = 1:20, reduction = 'harmony', 
                      n.neighbors = 50)
seurat.sele = FindNeighbors(seurat.sele, reduction = 'harmony', dims = 1:20)
seurat.sele = FindClusters(seurat.sele, resolution = 0.1)
seurat.sele = FindClusters(seurat.sele, resolution = 0.2)

degs_macro = FindAllMarkers(seurat.sele, test.use = 'LR', 
                            latent.vars = c('nCount_RNA', 'percent.mt'),
                            max.cells.per.ident = 500, only.pos = T)
degs_macro = data.table(degs_macro)

saveRDS(degs_macro, file = 'data/intermediate/RNA/new_degs/degs_by_clusters_macrophage_byHarmony_clean2.rds')


p1 <- DimPlot(seurat.sele, group.by = 'Patient_No')
p2 <- DimPlot(seurat.sele,  label = T) + scale_color_brewer(palette = 'Paired')

p3 <- StackedVlnPlot(seurat.sele, features = c('SPP1', 'C1QC', 'HBEGF', 'CD68'), 
                     group.by = 'seurat_clusters', 
                     myColors = brewer.pal(7, "Paired")) 

p4 <- StackedVlnPlot(seurat.sele, features = c('PHOX2B', 'ISL1', 'DBH'), 
                    group.by = 'seurat_clusters',  myColors = brewer.pal(7, "Paired")) 
p5 <- StackedVlnPlot(seurat.sele, features = c('WWTR1', 'YAP1', 'DCN'), 
                     group.by = 'seurat_clusters',  myColors = brewer.pal(7, "Paired")) 
p6 <- StackedVlnPlot(seurat.sele, features = c('CD68', 'CD163', 'ITGAM', 'CD86', 'MRC1', 'MSR1'), 
                     group.by = 'seurat_clusters',  myColors = brewer.pal(7, "Paired")) 


saveRDS(seurat.sele, file = 'Seurat_Objects/byCellType/seurat_macrophage_byHarmony_clean2.rds')

### < increase cl resolution ####
seurat.macro = readRDS(file = 'Seurat_Objects/byCellType/seurat_macrophage_byHarmony_clean2.rds')
seurat.macro = FindClusters(seurat.macro, resolution = 0.3)
seurat.macro = FindClusters(seurat.macro, resolution = 0.4)

pp <-  DimPlot(seurat.macro, label = T, group.by = 'RNA_snn_res.0.1') +
  DimPlot(seurat.macro, label = T, group.by = 'RNA_snn_res.0.2') + 
  DimPlot(seurat.macro, label = T, group.by = 'RNA_snn_res.0.3') +
  DimPlot(seurat.macro, label = T, group.by = 'RNA_snn_res.0.4') 

degs_macro = FindAllMarkers(seurat.macro, test.use = 'LR', 
                            latent.vars = c('nCount_RNA', 'percent.mt'),
                            logfc.threshold = 0.25, min.pct = 0.1, min.diff.pct = 0.05,
                            max.cells.per.ident = 1000, only.pos = T)
degs_macro = data.table(degs_macro)

saveRDS(degs_macro, file = 'data/intermediate/RNA/new_degs/degs_by_clusters_macrophage_byHarmony_clean2_res0.4.rds')
saveRDS(seurat.macro, file = 'Seurat_Objects/seurat_macrophage_byHarmony_clean2_res0.4.rds')




