library(Seurat)
library(ggplot2)
library(data.table)
library(harmony)
library(VisCello)
source('scDataAnalysis_Utilities_simp.R')
`%notin%` = Negate(`%in%`)

getPalette = colorRampPalette(brewer.pal(9, "Paired"))
myColors = c(getPalette(11), '#d285f2')
names(myColors) = c("Fibroblast",  "Adrenergic", "Endothelial", "T_Cell", "Macrophage",
                    "Dendritic" , 'Schwann', "B_Cell", "Kidney_Cell",
                    "Hepatocyte", "Adrn_Cortex", "doublets")

# start from rpca integration ####
seurat.rna = readRDS('Seurat_Objects/seurat_rna_allCells_rpca.rds')
DimPlot(seurat.rna, group.by = 'cell_type', label = T)

DimPlot(seurat.rna, group.by = 'Patient_No', label = T)

## degs between clusters
DimPlot(seurat.rna, label = T)
DefaultAssay(seurat.rna) = 'RNA'
degs = FindAllMarkers(seurat.rna, slot = 'data', test.use = 'LR', 
                      latent.vars = c('percent.mt', 'nCount_RNA'), 
                      max.cells.per.ident = 1000, 
                      only.pos = T)
degs = data.table(degs)

## compare too endo clusters
degs_endo = FindMarkers(seurat.rna, slot = 'data', test.use = 'LR',
                        latent.vars = c('percent.mt', 'nCount_RNA'), 
                        max.cells.per.ident = 500, 
                        ident.1 = 17, ident.2 = 10)
## compare myeloid clusters 
degs_myeloid = FindMarkers(seurat.rna, slot = 'data', test.use = 'LR',
                        latent.vars = c('percent.mt', 'nCount_RNA'), 
                        max.cells.per.ident = 500, 
                        ident.1 = 18, ident.2 = 4)
degs_myeloid = data.table(degs_myeloid, keep.rownames = T)
names(degs_myeloid)[1] = 'gene'
degs_myeloid[avg_log2FC > 1]

## assign cell type ####
seurat.rna$cell_type = 'Neuroblast'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(2, 12, 13, 16)] = 'Fibroblast'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(10, 17)] = 'Endothelial'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(4)] = 'Macrophage'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(5)] = 'Adrn_Cortex'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(7)] = 'T_Cell'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(9)] = 'Schwann'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(10, 17)] = 'Endothelial'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(11)] = 'Hepatocyte'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(14)] = 'Kidney_Cell'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(15)] = 'B_Cell'
seurat.rna$cell_type[seurat.rna$seurat_clusters %in% c(18)] = 'Dendritic'


seurat.rna <- FindClusters(seurat.rna, resolution = 0.4)
seurat.rna <- FindClusters(seurat.rna, resolution = 0.1)
seurat.rna <- FindClusters(seurat.rna, resolution = 0.2)

## remove possible doublets 
seurat.rna <- FindNeighbors(seurat.rna, reduction = 'umap', graph.name = c('umap_nn', 'umap_snn'), dims = 1:2)
seurat.rna <- FindClusters(seurat.rna, resolution = 0.1, 
                           graph.name = 'umap_snn', cluster.name = 'umap_clusters')
DimPlot(seurat.rna, label = T, group.by = 'umap_clusters')

seurat.rna <- FindClusters(seurat.rna, resolution = 0.2)

seurat.rna$cell_type[seurat.rna$umap_clusters %in% c(2, 8) & 
                       seurat.rna$integrated_snn_res.0.2 %in% c(7, 4, 15)] = 'doublets'
seurat.rna$cell_type[seurat.rna$umap_clusters %in% 15 & 
                       seurat.rna$integrated_snn_res.0.2 %in%  4] = 'doublets'

seurat.clean = subset(seurat.rna, cell_type != 'doublets')

p0 <- DimPlot(seurat.clean, group.by = 'cell_type', label = T) + 
  scale_color_manual(values = myColors)
p1 <- DimPlot(seurat.clean,  label = T) 

saveRDS(seurat.clean, file = 'Seurat_Objects/seurat_rna_allCells_clean.rds')
saveRDS(seurat.clean@meta.data, file = 'MetaData/metadata_seurat_rna_allCells_clean.rds')

# cell composition ####

## %pti in each cell type
#seurat.rna = readRDS('Seurat_Objects/seurat_rna_allCells_clean.rds')
#mdata = seurat.rna@meta.data
## cell composition stacked for each timepoint
mdata = readRDS('MetaData/metadata_seurat_rna_allCells_clean.rds')
mdata = subset(mdata, cell_type %notin% c('Adrn_Cortex', 'Kidney_Cell', 'Hepatocyte'))

mdata = data.table(mdata, select = c('cell_type', 'Stage_Code'))
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]
mdata[, 'n' := .N, by = list(cell_type, Stage_Code)]
mdata[, 'N' := .N, by = Stage_Code]
mdata[, 'percentage' := n/N]
mdata = subset(mdata, select = c( 'cell_type', 'Stage_Code', 'percentage')) %>% unique()

p2 <- ggplot(data = mdata, aes(x = Stage_Code, y = percentage, fill = cell_type)) +
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  scale_fill_manual(values = myColors) +
  xlab('') + ylab('Fraction') + guides(fill=guide_legend(title=""))

ggsave(p2, filename = 'Figures/RNA/cell_type_composition_byTimepoint_rna.pdf', device = 'pdf',
       width = 4, height = 5)

## test %change for a cell_type
## test %change for a specific state
mdata = readRDS('MetaData/metadata_seurat_rna_allCells_clean.rds')
mdata = subset(mdata, cell_type %notin% c('Adrn_Cortex', 'Kidney_Cell', 'Hepatocyte'))
mdata = data.table(mdata, keep.rownames = T)

table(mdata$malignancy, mdata$seurat_clusters)
tmp = mdata[seurat_clusters == 12 & malignancy == 'Malignant'] ##most cluster12 were predicted to be malignant (so probably not Fibroblasts)
mdata = mdata[rn %notin% tmp$rn]

mdata = mdata %>% subset(select = c(Patient_No, Stage_Code, cell_type))
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]

idx_pts = unique(mdata[Stage_Code == 'DX']$Patient_No)
mdata = mdata[Patient_No %in% idx_pts]

mdata[, 'N' := .N, by = list(Patient_No, Stage_Code)]

mdata1 = mdata[Stage_Code == 'DX']
mdata2 = mdata[Stage_Code == 'PTX']
mdata1[, 'n' := .N, by = list(Patient_No, cell_type)]
mdata2[, 'n' := .N, by = list(Patient_No, cell_type)]

mdata = rbind(mdata1, mdata2)
mdata[, 'frac' := n/N]
mdata = subset(mdata, select = c(Stage_Code, cell_type, frac, n, N, Patient_No)) %>% unique()

plist = list()
for(cl0 in unique(mdata$cell_type)){
  mdata0 = mdata[cell_type == cl0]
  setkey(mdata0, Patient_No)
  mdata0[, nSample := .N, by = Patient_No]
  mdata0 = subset(mdata0, select = c(Stage_Code, frac, Patient_No))
  mdata01 = mdata0[ Stage_Code == 'DX']
  mdata02 = mdata0[ Stage_Code == 'PTX']

  pts = union(mdata01$Patient_No, mdata02$Patient_No)
  y = mdata02[pts]$frac
  x = mdata01[pts]$frac
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  pv_cl = wilcox.test(y, x, paired = T)$p.value/2
  
  pdata = data.table('Patient_No' = c(pts, pts), 'frac' = c(x, y), 
                     'Stage_Code' = rep(c('DX', 'PTX'), each = length(pts)))
  pp <- ggplot(data = pdata, aes(x = Stage_Code, y = frac, colour = Stage_Code)) +
    geom_boxplot() + ggtitle(paste0(cl0)) + scale_color_brewer(palette = 'Set1', direction = -1) +
    geom_point(aes(fill = Stage_Code), size = 2) + theme_classic() + NoLegend() +
    xlab('') + ylab('Fraction') + scale_fill_brewer(palette = 'Set1', direction = -1) +
    geom_line(aes(group = Patient_No), color='gray', alpha=0.5, size = 1) +
    annotate('text', x = 1.5, y = max(pdata$frac), label = paste0('p = ', format(pv_cl, digits = 2)))
  
  ggsave(pp, filename = paste0('Figures/RNA/', cl0, '_frac_byTimepoint_rna.pdf'), device = 'pdf',
         width = 3, height = 4)
  plist[[cl0]] = pp
}

ggsave(grid.arrange(grobs = plist, ncol = 4),
       filename = 'Figures/RNA/cellType_frac_byTimepoint_rna.pdf', device = 'pdf',
       width = 8, height = 6)

# feauture plots ####
## dot plots ####
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
seurat.rna$cell_type = factor(seurat.rna$cell_type, levels =  rev(c("Adrenergic", "Fibroblast", "T_Cell", "Endothelial", "Macrophage",
                                                                "Dendritic"  , "Adrn_Cortex", "Schwann",   "B_Cell",  
                                                                "Hepatocyte", "Kidney_Cell")))
p5 <- DotPlot(seurat.rna, assay = 'RNA', features = gene_markers, 
        group.by = 'cell_type') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab('') + ylab('')
ggsave(p5, filename = 'Figures/RNA/marker_gene_expression.pdf', device = 'pdf',
       width = 15, height = 4)


## degs btw cell types ####
seurat.rna = readRDS('Seurat_Objects/seurat_rna_allCells_clean.rds')

DefaultAssay(seurat.rna) = 'RNA'
Idents(seurat.rna) = seurat.rna$cell_type
degs <- FindAllMarkers(seurat.rna, 
                       assay = 'RNA', test.use = 'LR', logfc.threshold = 0.0, 
                       latent.vars = c('nCount_RNA', 'percent.mt'),
                       slot = 'data', max.cells.per.ident = 1000)
degs = data.table(degs)
saveRDS(degs, file = 'data/intermediate/RNA/new_degs/degs_cell_type.rds')

degs[, 'pct.diff' := pct.1 - pct.2]
degs = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]

split_degs = split(degs, by = 'cluster')
writexl::write_xlsx(split_degs, path = 'data/intermediate/RNA/new_degs/degs_between_cell_types.xlsx')



