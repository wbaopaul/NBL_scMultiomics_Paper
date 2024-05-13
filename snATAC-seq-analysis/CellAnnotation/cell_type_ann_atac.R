library(Seurat)
library(ggplot2)
library(data.table)
library(harmony)
library(Signac)
`%notin%` = Negate(`%in%`)

getPalette = colorRampPalette(brewer.pal(9, "Paired"))
myColors = c(getPalette(11), '#d285f2')
names(myColors) = c("Fibroblast",  "Neuroblast", "Endothelial", "T_Cell", "Macrophage",
                    "Dendritic" , 'Schwann', "B_Cell", "Kidney_Cell",
                    "Hepatocyte", "Adrn_Cortex", "doublets")

## load integrated seurat object ####
seurat.atac = readRDS('Seurat_Objects/seurat_atac_signac_byPatient.rds')


DimPlot(seurat.atac)
## annotated manually based on pooled data
DimPlot(seurat.atac, group.by = 'cell_type0') + scale_color_manual(values = myColors) 

activity.matrix = readRDS('Seurat_Objects/mtx_gene_activity_signac.rds')
activity.matrix = activity.matrix[, colnames(seurat.atac)]
seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

seurat.atac <- NormalizeData(
  object = seurat.atac,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat.atac$nCount_ACTIVITY)
)

DefaultAssay(seurat.atac) = 'ACTIVITY'

## add label tranfer result from newly annotated rna cells
labelTransRes = readRDS('MetaData/labelTransferResult_seurat_atac_signac_byPatient.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = labelTransRes)
p0 <- DimPlot(seurat.atac, group.by = 'predicted.id', label = T) + 
  scale_color_manual(values = myColors) 
p1 <- DimPlot(seurat.atac, label = T) 



## further refine cell type annotation ####
seurat.atac = readRDS('Seurat_Objects/seurat_atac_signac_byPatient.rds')
labelTransRes = readRDS('MetaData/labelTransferResult_seurat_atac_signac_byPatient.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = labelTransRes)

seurat.atac <- FindClusters(seurat.atac, resolution = 0.5)
DimPlot(seurat.atac, label = T)

seurat.atac$cell_type = 'Neuroblast'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(7)] = 'Adrn_Cortex'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(12, 17)] = 'T_Cell'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(15)] = 'B_Cell'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(6)] = 'Macrophage'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(8, 19)] = 'Fibroblast'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(14)] = 'Schwann'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(13)] = 'Endothelial'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(16)] = 'Hepatocyte'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(20)] = 'Dendritic'
seurat.atac$cell_type[seurat.atac$seurat_clusters %in% c(13) & seurat.atac$predicted.id == 'Kidney_Cell'] = 'Kidney_Cell'

p2 <- DimPlot(seurat.atac, group.by = 'cell_type', label = T) + 
  scale_color_manual(values = myColors) 


saveRDS(seurat.atac@meta.data, file = 'MetaData/metadata_seurat_atac_signac_byPatient.rds')
saveRDS(seurat.atac, file = 'Seurat_Objects/seurat_atac_signac_byPatient_cellTypeAnnotated.rds')


gene_markers = c('PHOX2B', 'ISL1', 
                 'PDGFRB', 'DCN', 
                 'CD247', 'CD96',
                 'PECAM1', 'PTPRB',
                 'CD163', 'CD86',
                 'IRF8', 'FLT3',
                 'CYP11A1', 'CYP11B1',
                 'PLP1', 'CDH19',
                 'PAX5', 'MS4A1',
                 'ALB', 'DCDC2', 'PKHD1')
seurat.atac$cell_type = factor(seurat.atac$cell_type, levels =  rev(c("Neuroblast", "Fibroblast", "T_Cell", "Endothelial", "Macrophage",
                                                                    "Dendritic"  , "Adrn_Cortex", "Schwann",   "B_Cell",  
                                                                    "Hepatocyte", "Kidney_Cell")))

activity.matrix = readRDS('Seurat_Objects/mtx_gene_activity_signac_new.rds')
activity.matrix = activity.matrix[, colnames(seurat.atac)]
seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

seurat.atac <- NormalizeData(
  object = seurat.atac,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat.atac$nCount_ACTIVITY)
)
DefaultAssay(seurat.atac) = 'ACTIVITY'
seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))

p3 <- DotPlot(seurat.atac, assay = 'ACTIVITY', features = gene_markers, 
              group.by = 'cell_type', cols = c('grey100', 'red3')) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab('') + ylab('')

## cell proportion changes ####
mdata = readRDS('MetaData/metadata_seurat_atac_signac_byPatient.rds')
mdata = subset(mdata, cell_type %notin% c('Adrn_Cortex', 'Kidney_Cell', 'Hepatocyte'))

mdata = data.table(mdata, select = c('cell_type', 'Stage_Code'))
mdata = mdata[Stage_Code %in% c('IDX', 'PTI')]
mdata[, 'n' := .N, by = list(cell_type, Stage_Code)]
mdata[, 'N' := .N, by = Stage_Code]
mdata[, 'percentage' := n/N]
mdata = subset(mdata, select = c( 'cell_type', 'Stage_Code', 'percentage')) %>% unique()

p3 <- ggplot(data = mdata, aes(x = Stage_Code, y = percentage, fill = cell_type)) +
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  scale_fill_manual(values = myColors) +
  xlab('') + ylab('Fraction') + guides(fill=guide_legend(title=""))


mdata = readRDS('MetaData/metadata_seurat_atac_signac_byPatient.rds')
mdata = subset(mdata, cell_type %notin% c('Adrn_Cortex', 'Kidney_Cell', 'Hepatocyte'))
mdata = data.table(mdata)
idx_pts = unique(mdata[Stage_Code == 'IDX']$Patient_No)
mdata = mdata[Patient_No %in% idx_pts]

mdata = mdata %>% subset(select = c(Patient_No, Stage_Code, cell_type))
mdata = mdata[Stage_Code %in% c('IDX', 'PTI')]
mdata[, 'N' := .N, by = list(Patient_No, Stage_Code)]

mdata1 = mdata[Stage_Code == 'IDX']
mdata2 = mdata[Stage_Code == 'PTI']
mdata1[, 'n' := .N, by = list(Patient_No, cell_type)]
mdata2[, 'n' := .N, by = list(Patient_No, cell_type)]

mdata = rbind(mdata1, mdata2)
mdata[, 'frac' := n/N]
mdata = subset(mdata, select = c(Stage_Code, cell_type, frac, n, N, Patient_No)) %>% unique()

for(cl0 in unique(mdata$cell_type)){
  mdata0 = mdata[cell_type == cl0]
  setkey(mdata0, Patient_No)
  mdata0[, nSample := .N, by = Patient_No]
  mdata0 = subset(mdata0, select = c(Stage_Code, frac, Patient_No))
  mdata01 = mdata0[ Stage_Code == 'IDX']
  mdata02 = mdata0[ Stage_Code == 'PTI']
  
  pts = union(mdata01$Patient_No, mdata02$Patient_No)
  y = mdata02[pts]$frac
  x = mdata01[pts]$frac
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  pv_cl = wilcox.test(y, x, paired = T)$p.value/2
  
  
  pp <- ggplot(data = mdata0, aes(x = Stage_Code, y = frac, colour = Stage_Code)) +
    geom_boxplot() + ggtitle(paste0(cl0)) + scale_color_brewer(palette = 'Set1', direction = -1) +
    geom_point(aes(fill = Stage_Code), size = 2) + theme_classic() + NoLegend() +
    xlab('') + ylab('Fraction') + scale_fill_brewer(palette = 'Set1', direction = -1) +
    geom_line(aes(group = Patient_No), color='gray', alpha=0.5, size = 1) +
    annotate('text', x = 1.5, y = max(mdata0$frac), label = paste0('p = ', format(pv_cl, digits = 2)))
  
  ggsave(pp, filename = paste0('Figures/ATAC/', cl0, '_frac_byTimepoint_atac.pdf'), device = 'pdf',
         width = 3, height = 4)
  
}


