library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
library(harmony)
library(RColorBrewer)
library(enrichR)
library(viridis)

seuratProj <- function(seurat.query, seurat.ref, queryName4plot = 'WT',
                       ref_ctype = 'cell_type', 
                       dims = 1:30, ref_reduction = 'pca'){
  anchor <- FindTransferAnchors(
    reference = seurat.ref,
    query = seurat.query,
    reference.reduction = ref_reduction,
    dims = dims
  )
  
  seurat.query <- MapQuery(
    anchorset = anchor,
    query = seurat.query,
    reference = seurat.ref,
    refdata = list(cell_type = ref_ctype,
                   ref_cluster = 'seurat_clusters'),
    reference.reduction = ref_reduction,
    reference.dims = dims,
    reduction.model = "umap"
  )
  
  umap.query = seurat.query@reductions$ref.umap@cell.embeddings
  colnames(umap.query) = c('umap_1', 'umap_2')
  umap.ref = seurat.ref@reductions$umap@cell.embeddings
  umap.query.ref = rbind(umap.ref, umap.query)
  umap.query.ref = data.table(umap.query.ref)
  umap.query.ref$ctype = rep(c('Ref', queryName4plot), c(nrow(umap.ref), nrow(umap.query)))
  
  myColors = c('red', '#cccccc')
  names(myColors) = c(queryName4plot, 'Ref')
  
  pp <- ggplot(umap.query.ref, aes(x = umap_1, y = umap_2, col = ctype)) + 
    geom_point(size = 0.4) + 
    scale_colour_manual(name = "annotation", values = myColors) + theme_classic()
  
  return(list('plot' = pp, 'seurat.query' = seurat.query))
}


stateColor = brewer.pal(n = 6, name = "Set1")
names(stateColor) = c('ADRN-Calcium', 'ADRN-Baseline', 'Interm-OxPhs',
                      'ADRN-Dopaminergic', 'ADRN-Proliferating', 'MES')


## CHLA15 ####
seu.ch15.cr = readRDS('Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15_THP1_Monoculture_Control_Singlets.rds')
seu.ch15.tr = readRDS('Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15_THP1_CRM_Afatinib_Singlets.rds')

seu.ch15 = merge(seu.ch15.cr, seu.ch15.tr,
                 add.cell.ids = c('control', 'treated'))
seu.ch15 <- JoinLayers(seu.ch15)
seu.ch15$deMULTIplex2[seu.ch15$deMULTIplex2== 'CHLA15_THP1_Control'] = 'CHLA15_THP1_Coculture'

seu.ch15 <- seu.ch15 %>% NormalizeData(scale.factor = median(seu.ch15$nCount_RNA)) %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) 

seu.ch15 <- FindNeighbors(seu.ch15, dims = 1:20) %>% FindClusters(resolution = 0.1) %>%
  RunUMAP(verbose = F, dims = 1:20)


DimPlot(seu.ch15, label = T)
DimPlot(seu.ch15, group.by = 'deMULTIplex2')
FeaturePlot(seu.ch15, features = c('PHOX2B', 'CD68'))
FeaturePlot(seu.ch15, features = c('MKI67', 'EZH2'))
FeaturePlot(seu.ch15, features = c('percent.mt', 'nFeature_RNA'))

## plot umaps
p0 <- FeaturePlot(seu.ch15, features = c('PHOX2B', 'CD68'),
                  pt.size = 0.1)
ggsave(p0, file = 'Figures/RNA/coculture/CHLA15_umap_PHOX2B_CD68.pdf',
       width = 6, height = 3, device = 'pdf')

#seu.ch15 <- RunHarmony(seu.ch15, 'deMULTIplex2') %>% 
#  RunUMAP(dims = 1:20, reduction = 'harmony')

## < work on malignant cells ###
seu.ch15.malignant = subset(seu.ch15, seurat_clusters %in% c(1, 2, 3))

seu.ch15.malignant <- seu.ch15.malignant %>% NormalizeData(scale.factor = median(seu.ch15.malignant$nCount_RNA)) %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) 

seu.ch15.malignant <- FindNeighbors(seu.ch15.malignant, dims = 1:20) %>% FindClusters(resolution = 0.1) %>%
  RunUMAP(verbose = F, dims = 1:20)

DimPlot(seu.ch15.malignant, label = T)
DimPlot(seu.ch15.malignant, group.by = 'deMULTIplex2')
FeaturePlot(seu.ch15.malignant, features = c('PHOX2B', 'CD68'))
FeaturePlot(seu.ch15.malignant, features = c('MKI67', 'EZH2'))
saveRDS(seu.ch15, file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15.rds')

## projection 
seurat.malignant = readRDS('Seurat_Objects/seurat_rna_malignant_harmony_cellStateAnnotated.rds')

seurat.malignant = RunUMAP(seurat.malignant, reduction = "harmony", dims = 1:20, 
                           return.model = T)
projRes = seuratProj(seu.ch15.malignant, seurat.malignant, 'coculture', dims = 1:20,
                     ref_reduction = 'harmony', ref_ctype = 'cell_state')
seu.ch15.malignant = projRes$seurat.query

p1 <- DimPlot(seu.ch15.malignant, group.by = 'predicted.cell_type') +
  scale_color_manual(values = stateColor)
p2 <- FeaturePlot(seu.ch15.malignant, features = 'predicted.cell_type.score',
                  max.cutoff = 'q99')

saveRDS(seu.ch15.malignant, file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15_malignant.rds')


### < proj% changes in malignant cells ####
seu.ch15.malignant = readRDS(file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15_malignant.rds')

mdata = subset(seu.ch15.malignant@meta.data, 
               select = c('deMULTIplex2', 'predicted.cell_type'))
mdata = data.table(mdata)
mdata[, 'N' := .N, by = deMULTIplex2]
mdata[, 'n' := .N, by = list(deMULTIplex2, predicted.cell_type)]
mdata[, 'frac' := n/N]
mdata <- unique(mdata)
mdata$experiment = stringr::str_split_i(mdata$deMULTIplex2, '_', 3)

library(dplyr)

mdata1 <- mdata %>% arrange(experiment, predicted.cell_type) %>%
  group_by(experiment) %>%
  mutate(label_y = cumsum(frac))

mdata1$experiment = factor(mdata1$experiment, levels = c('Monoculture', 'Coculture',
                                                         'Afatinib', 'CRM'))

p1 <- ggplot(data = mdata1, aes(x = experiment, y = frac, fill = predicted.cell_type)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = stateColor) + 
  theme_classic() + xlab('') + ylab('Fraction') + ggtitle('CHLA15') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(
    aes(label = round(frac, 2), y = 1-label_y + 0.04),
    colour = "black", size = 3,
    vjust = 1.5, position = position_dodge(.5)
  )
ggsave(p1, file = 'Figures/RNA/coculture/CHLA15_malignant_cell_state_frac.pdf',
       width = 4, height = 4, device = 'pdf')

## compare some changes by proportion test

experiment1 = 'Coculture'
for(cell_state in unique(mdata$predicted.cell_type)){
  for(experiment2 in c('Monoculture', 'Afatinib', 'CRM')){
    aa = mdata[predicted.cell_type == cell_state & experiment == experiment1]
    bb = mdata[predicted.cell_type == cell_state & experiment == experiment2]
    cc = rbind(aa, bb)
    message(paste('CHLA15:', cell_state, 'change in', experiment1, 'vs', experiment2, 'p-value =',
                  prop.test(cc$n, cc$N)$p.val/2))
  }
}



### < proj% changes in macrophage subsets ####
seu.ch15 = readRDS(file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15.rds')
seu.ch15.macro = subset(seu.ch15, seurat_clusters == 0)

seu.ch15.macro <- seu.ch15.macro %>% NormalizeData(scale.factor = median(seu.ch15.macro$nCount_RNA)) %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) 

seu.ch15.macro <- FindNeighbors(seu.ch15.macro, dims = 1:20) %>% FindClusters(resolution = 0.1) %>%
  RunUMAP(verbose = F, dims = 1:20)

DimPlot(seu.ch15.macro, label = T)
DimPlot(seu.ch15.macro, group.by = 'deMULTIplex2')
FeaturePlot(seu.ch15.macro, features = c('CD68', 'PTPRC'))
FeaturePlot(seu.ch15.macro, features = c('MKI67', 'SPP1', 'C1QC', 
                                         'VCAN', 'VEGFA',
                                         'IL18', 'CCL4', 'CD163', 'MRC1',
                                         'F13A1', 'HS3ST2', 'THY1'))
seurat.macro = readRDS(file = 'seurat_macrophage_byHarmony_clean2_cDCremoved.rds')
seurat.macro$cell_state = gsub('Macro', '', seurat.macro$cell_state)

seurat.macro = RunUMAP(seurat.macro, reduction = "harmony", dims = 1:20, 
                       return.model = T, n.neighbors = 50)
projRes = seuratProj(seu.ch15.macro, seurat.macro, 'coculture', dims = 1:20,
                     ref_reduction = 'harmony', ref_ctype = 'cell_state')
seu.ch15.macro = projRes$seurat.query

macroColors = brewer.pal(8, 'Dark2')
names(macroColors) = c('THY1+', 'HS3ST2+',
                       'F13A1+', 'CCL4+', 
                       'IL18+','VCAN+',
                       'C1QC+SPP1+', 'Proliferating')


p1 <- DimPlot(seu.ch15.macro, group.by = 'predicted.cell_type') +
  scale_color_manual(values = macroColors)
p2 <- FeaturePlot(seu.ch15.macro, features = 'predicted.cell_type.score',
                  max.cutoff = 'q99')

mdata = subset(seu.ch15.macro@meta.data, 
               select = c('deMULTIplex2', 'predicted.cell_type'))
mdata = data.table(mdata)
mdata[, 'N' := .N, by = deMULTIplex2]
mdata[, 'n' := .N, by = list(deMULTIplex2, predicted.cell_type)]
mdata[, 'frac' := n/N]
mdata <- unique(mdata)
mdata$experiment = stringr::str_split_i(mdata$deMULTIplex2, '_', 3)


mdata1 <- mdata %>% arrange(experiment, predicted.cell_type) %>%
  group_by(experiment) %>%
  mutate(label_y = cumsum(frac))

mdata1$experiment = factor(mdata1$experiment, levels = c('Monoculture', 'Coculture',
                                                         'CRM', 'Afatinib'))

p1 <- ggplot(data = mdata1, aes(x = experiment, y = frac, fill = predicted.cell_type)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = macroColors) + 
  theme_classic() + xlab('') + ylab('Fraction') + ggtitle('CHLA15') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(
    aes(label = round(frac, 2), y = 1-label_y + 0.04),
    colour = "black", size = 3,
    vjust = 1.5, position = position_dodge(.5)
  )
ggsave(p1, file = 'Figures/RNA/coculture/CHLA15_macro_cell_state_frac.pdf',
       width = 4, height = 4, device = 'pdf')


## CHLA20 ####
seu.ch20.cr = readRDS('Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA20_THP1_Monoculture_Control_Singlets.rds')
seu.ch20.tr = readRDS('Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA20_THP1_CRM_Afatinib_Singlets.rds')

seu.ch20 = merge(seu.ch20.cr, seu.ch20.tr,
                 add.cell.ids = c('control', 'treated'))
seu.ch20 <- JoinLayers(seu.ch20)
seu.ch20$deMULTIplex2[seu.ch20$deMULTIplex2== 'CHLA20_THP1_Control'] = 'CHLA20_THP1_Coculture'

seu.ch20 <- seu.ch20 %>% NormalizeData(scale.factor = median(seu.ch20$nCount_RNA)) %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) 

seu.ch20 <- FindNeighbors(seu.ch20, dims = 1:20) %>% FindClusters(resolution = 0.1) %>%
  RunUMAP(verbose = F, dims = 1:20)


DimPlot(seu.ch20, label = T)
DimPlot(seu.ch20, group.by = 'deMULTIplex2')
FeaturePlot(seu.ch20, features = c('PHOX2B', 'CD68'))
FeaturePlot(seu.ch20, features = c('MKI67', 'EZH2'))
FeaturePlot(seu.ch20, features = c('percent.mt', 'nFeature_RNA'))

p0 <- FeaturePlot(seu.ch20, features = c('PHOX2B', 'CD68'),
                  pt.size = 0.1)
ggsave(p0, file = 'Figures/RNA/coculture/CHLA20_umap_PHOX2B_CD68.pdf',
       width = 6, height = 3, device = 'pdf')

## < work on malignant cells ###
seu.ch20.malignant = subset(seu.ch20, seurat_clusters %in% c(1))

seu.ch20.malignant <- seu.ch20.malignant %>% NormalizeData(scale.factor = median(seu.ch20.malignant$nCount_RNA)) %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) 

seu.ch20.malignant <- FindNeighbors(seu.ch20.malignant, dims = 1:20) %>% 
  FindClusters(resolution = 0.2) %>%
  RunUMAP(verbose = F, dims = 1:20)

DimPlot(seu.ch20.malignant, label = T)
DimPlot(seu.ch20.malignant, group.by = 'deMULTIplex2')
FeaturePlot(seu.ch20.malignant, features = c('PHOX2B', 'CD68'))
FeaturePlot(seu.ch20.malignant, features = c('MKI67', 'EZH2'))

saveRDS(seu.ch20, file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA20.rds')

projRes = seuratProj(seu.ch20.malignant, seurat.malignant, 'coculture', dims = 1:20,
                     ref_reduction = 'harmony', ref_ctype = 'cell_state')
seu.ch20.malignant = projRes$seurat.query

DimPlot(seu.ch20.malignant, group.by = 'predicted.cell_type')
FeaturePlot(seu.ch20.malignant, features = 'predicted.cell_type.score')

p1 <- DimPlot(seu.ch20.malignant, group.by = 'predicted.cell_type') +
  scale_color_manual(values = stateColor)
p2 <- FeaturePlot(seu.ch20.malignant, features = 'predicted.cell_type.score',
                  max.cutoff = 'q99')

saveRDS(seu.ch20.malignant, file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA20_malignant.rds')

### < proj% changes in malignant cells ####
mdata = subset(seu.ch20.malignant@meta.data, 
               select = c('deMULTIplex2', 'predicted.cell_type'))
mdata = data.table(mdata)
mdata[, 'N' := .N, by = deMULTIplex2]
mdata[, 'n' := .N, by = list(deMULTIplex2, predicted.cell_type)]
mdata[, 'frac' := n/N]
mdata <- unique(mdata)
mdata$experiment = stringr::str_split_i(mdata$deMULTIplex2, '_', 3)

library(dplyr)

mdata1 <- mdata %>% arrange(experiment, predicted.cell_type) %>%
  group_by(experiment) %>%
  mutate(label_y = cumsum(frac))

mdata1$experiment = factor(mdata1$experiment, levels = c('Monoculture', 'Coculture',
                                                         'Afatinib', 'CRM'))

p2 <- ggplot(data = mdata1, aes(x = experiment, y = frac, fill = predicted.cell_type)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = stateColor) + 
  theme_classic() + xlab('') + ylab('Fraction') + ggtitle('CHLA20') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(
    aes(label = round(frac, 2), y = 1-label_y + 0.04),
    colour = "black", size = 3,
    vjust = 1.5, position = position_dodge(.5)
  )
ggsave(p2, file = 'Figures/RNA/coculture/CHLA20_malignant_cell_state_frac.pdf',
       width = 4, height = 4, device = 'pdf')

## comparison by proportion test
experiment1 = 'Coculture'
for(cell_state in unique(mdata$predicted.cell_type)){
  for(experiment2 in c('Monoculture', 'Afatinib', 'CRM')){
    aa = mdata[predicted.cell_type == cell_state & experiment == experiment1]
    bb = mdata[predicted.cell_type == cell_state & experiment == experiment2]
    cc = rbind(aa, bb)
    message(paste('CHLA20:', cell_state, 'change in', experiment1, 'vs', experiment2, 'p-value =',
                  prop.test(cc$n, cc$N)$p.val/2))
  }
}

### < proj% changes in macrophage subsets ####
seu.ch20 = readRDS(file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA20.rds')
seu.ch20.macro = subset(seu.ch20, seurat_clusters %in% c(0, 2)) 

seu.ch20.macro <- seu.ch20.macro %>% NormalizeData(scale.factor = median(seu.ch20.macro$nCount_RNA)) %>%
  FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20, verbose = F) 

seu.ch20.macro <- FindNeighbors(seu.ch20.macro, dims = 1:20) %>% FindClusters(resolution = 0.1) %>%
  RunUMAP(verbose = F, dims = 1:20)

DimPlot(seu.ch20.macro, label = T)
DimPlot(seu.ch20.macro, group.by = 'deMULTIplex2')
FeaturePlot(seu.ch20.macro, features = c('CD68', 'PTPRC'))
FeaturePlot(seu.ch20.macro, features = c('MKI67', 'SPP1', 'C1QC', 
                                         'VCAN', 'VEGFA',
                                         'IL18', 'CCL4', 'CD163', 'MRC1',
                                         'F13A1', 'HS3ST2', 'THY1'))
seurat.macro = readRDS(file = 'Seurat_Objects/seurat_macrophage_byHarmony_clean2_cDCremoved.rds')
seurat.macro$cell_state = gsub('Macro', '', seurat.macro$cell_state)

seurat.macro = RunUMAP(seurat.macro, reduction = "harmony", dims = 1:20, 
                       return.model = T, n.neighbors = 50)
projRes = seuratProj(seu.ch20.macro, seurat.macro, 'coculture', dims = 1:20,
                     ref_reduction = 'harmony', ref_ctype = 'cell_state')
seu.ch20.macro = projRes$seurat.query

macroColors = brewer.pal(8, 'Dark2')
names(macroColors) = c('THY1+', 'HS3ST2+',
                       'F13A1+', 'CCL4+', 
                       'IL18+','VCAN+',
                       'C1QC+SPP1+', 'Proliferating')


p1 <- DimPlot(seu.ch20.macro, group.by = 'predicted.cell_type') +
  scale_color_manual(values = macroColors)
p2 <- FeaturePlot(seu.ch20.macro, features = 'predicted.cell_type.score',
                  max.cutoff = 'q99')

mdata = subset(seu.ch20.macro@meta.data, 
               select = c('deMULTIplex2', 'predicted.cell_type'))
mdata = data.table(mdata)
mdata[, 'N' := .N, by = deMULTIplex2]
mdata[, 'n' := .N, by = list(deMULTIplex2, predicted.cell_type)]
mdata[, 'frac' := n/N]
mdata <- unique(mdata)
mdata$experiment = stringr::str_split_i(mdata$deMULTIplex2, '_', 3)


mdata1 <- mdata %>% arrange(experiment, predicted.cell_type) %>%
  group_by(experiment) %>%
  mutate(label_y = cumsum(frac))

mdata1$experiment = factor(mdata1$experiment, levels = c('Monoculture', 'Coculture',
                                                         'CRM', 'Afatinib'))

p1 <- ggplot(data = mdata1, aes(x = experiment, y = frac, fill = predicted.cell_type)) +
  geom_bar(stat = 'identity') + scale_fill_manual(values = macroColors) + 
  theme_classic() + xlab('') + ylab('Fraction') + ggtitle('CHLA20') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_text(
    aes(label = round(frac, 2), y = 1-label_y + 0.04),
    colour = "black", size = 3,
    vjust = 1.5, position = position_dodge(.5)
  )
ggsave(p1, file = 'Figures/RNA/coculture/CHLA20_macro_cell_state_frac.pdf',
       width = 4, height = 4, device = 'pdf')

