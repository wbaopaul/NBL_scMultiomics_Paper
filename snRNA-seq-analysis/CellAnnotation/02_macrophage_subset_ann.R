source('scDataAnalysis_Utilities_simp.R')
library(openxlsx)
library(ggplot2)
library(ggridges)
`%notin%` = Negate('%in%')
library(pheatmap)
library(stringr)
library(harmony)

seurat.macro = readRDS('Seurat_Objects/seurat_macrophage_byHarmony_clean2_res0.4.rds')
degs = readRDS('data/intermediate/RNA/new_degs/degs_by_clusters_macrophage_byHarmony_clean2_res0.4.rds')
degs = degs[p_val < 0.001 & avg_log2FC > 0.5]
degs = degs[(pct.1 - pct.2) > 0.1]

## immunosuppressive/proinflammatory/phagocytosis/agiogenesis ####
## remove cDC cluster from downstream analysis
seurat.macro = subset(seurat.macro, cell_state != 'CD1C+cDC_6')
Immunosuppressive_markers <- c('CD84','LGALS3','SPP1','A2M','ENTPD1',
                               'TGFB1','NECTIN2','CCL2','CCL17','CCL22',
                               'CD274','PDCD1LG2','PTGER4','AHR',
                               'HLA-E','LGALS9','VSIR',
                               'NECTIN3','IL10',"ARG1","IDO1")
ProInflammatory_markers <- c('IL18','IL12A','IL12B',
                             'IL1B', 'TNF', 'IL6',
                             'IFNG','IFNA1','CXCL10','CXCL9',
                             'CD80','CD86',
                             'NOS2')
Phagocytosis_markers <- c('MRC1', 'CD163', 'MERTK', 'C1QB')
Agiogenesis_markers <- c('CCND2','CCNE1','CD44','CXCR4','E2F3','EDN1',
                         'EZH2','FGF18','FGFR1','FYN','HEY1','ITGAV',
                         'JAG1','JAG2','MMP9','NOTCH1','PDGFA','PTK2',
                         'SPP1','STC1','TNFAIP6','TYMP','VAV2','VCAN','VEGFA')

seurat.macro <- AddModuleScore(seurat.macro, features = list(Immunosuppressive_markers),
                               name = 'immunosuppr_score')
seurat.macro <- AddModuleScore(seurat.macro, features = list(ProInflammatory_markers),
                               name = 'proinflamm_score')
seurat.macro <- AddModuleScore(seurat.macro, features = list(Agiogenesis_markers),
                               name = 'agiogenesis_score')
seurat.macro <- AddModuleScore(seurat.macro, features = list(Phagocytosis_markers),
                               name = 'phagocytosis_score')



p00 <- VlnPlot(seurat.macro, features = 'immunosuppr_score1', pt.size = 0, y.max = 1, group.by = 'cell_state') +
  geom_hline(yintercept=0, linetype="dashed", color = "red", linewidth=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_manual(values = maColors)
p01 <- VlnPlot(seurat.macro, features = 'proinflamm_score1', pt.size = 0, y.max = 1, group.by = 'cell_state') +
  geom_hline(yintercept=0, linetype="dashed", color = "red", linewidth=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_manual(values = maColors)

p20 <- VlnPlot(seurat.macro, features = 'agiogenesis_score1', pt.size = 0, y.max = 1, 
               group.by = 'cell_state') +
  geom_hline(yintercept=0, linetype="dashed", color = "red", linewidth=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_manual(values = maColors)

p21 <- VlnPlot(seurat.macro, features = 'phagocytosis_score1', pt.size = 0, y.max = 1, 
               group.by = 'cell_state') +
  geom_hline(yintercept=0, linetype="dashed", color = "red", linewidth=1) +
  stat_summary(fun.y = median, geom='point', size = 10, 
               colour = "black", shape = 95) + NoLegend() +
  scale_fill_manual(values = maColors)


for(sele.gene in Immunosuppressive_markers){
  pp <- VlnPlot(seurat.macro, features = sele.gene,  pt.size = 0, group.by = 'cell_state') + 
    NoLegend() + scale_fill_manual(values = maColors)
  ggsave(pp, filename = paste0('Figures/RNA/OtherCellTypes/Macrophage/ImmunoSuppr/', sele.gene, '_macro.pdf'),
         width = 4, height = 4, device = 'pdf')
  
}

for(sele.gene in ProInflammatory_markers){
  pp <-  VlnPlot(seurat.macro, features = sele.gene,  pt.size = 0, 
                 group.by = 'cell_state') + 
    NoLegend() + scale_fill_manual(values = maColors)  
  ggsave(pp, filename = paste0('Figures/RNA/OtherCellTypes/Macrophage/ProInflam/', sele.gene, '_macro.pdf'),
         width = 4, height = 4, device = 'pdf')
  
}

for(sele.gene in Phagocytosis_markers){
  pp <-  VlnPlot(seurat.macro, features = sele.gene,  pt.size = 0, group.by = 'cell_state') + 
    NoLegend() + scale_fill_manual(values = maColors)  
  ggsave(pp, filename = paste0('Figures/RNA/OtherCellTypes/Macrophage/Phagocytosis/', sele.gene, '_macro.pdf'),
         width = 4, height = 4, device = 'pdf')
  
}

for(sele.gene in Agiogenesis_markers){
  pp <-  VlnPlot(seurat.macro, features = sele.gene,  pt.size = 0, group.by = 'cell_state') + 
    NoLegend() + scale_fill_manual(values = maColors)  
  ggsave(pp, filename = paste0('Figures/RNA/OtherCellTypes/Macrophage/Agiogenesis/', sele.gene, '_macro.pdf'),
         width = 4, height = 4, device = 'pdf')
}

saveRDS(seurat.macro, file = 'Seurat_Objects/seurat_macrophage_byHarmony_clean2_cDCremoved.rds')

## replot umap & marker gene dot plot ####
seurat.macro = readRDS(file = 'Seurat_Objects/seurat_macrophage_byHarmony_clean2_cDCremoved.rds')
sele.genes = c('THY1', 'CYP27A1', 'HS3ST2', 'LYVE1','F13A1', 
               'CCL4','IL18', 'VEGFA',  'VCAN', 'SPP1', 'C1QC', 
               'MKI67', 'TOP2A')

seurat.macro$cell_state = factor(seurat.macro$cell_state, 
                                 levels = c('THY1+Macro_0', 'HS3ST2+Macro_1',
                                            'F13A1+Macro_2', 'CCL4+Macro_3', 
                                            'IL18+Macro_4','VCAN+Macro_5',
                                            'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8'))  
maColors = brewer.pal(8, 'Dark2')
names(maColors) = c('THY1+Macro_0', 'HS3ST2+Macro_1',
                    'F13A1+Macro_2', 'CCL4+Macro_3', 
                    'IL18+Macro_4','VCAN+Macro_5',
                    'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8')
maColors['CD1C+cDC_6'] = '#cccccc'

mcColors = maColors[c(1:6, 9, 7:8)]
names(mcColors) = 0:8
pm <- DimPlot(seurat.macro, group.by = 'cell_state', cols = maColors, label = T)
pc <- DimPlot(seurat.macro, label = T, cols = c(mcColors))
ggsave(pm, filename = 'Figures/RNA/OtherCellTypes/Macrophage/umap_byState_cDCRemoved.pdf',
       width = 8, height = 6, device = 'pdf')
ggsave(pc, filename = 'Figures/RNA/OtherCellTypes/Macrophage/umap_byCluster_cDCRemoved.pdf',
       width = 7, height = 6, device = 'pdf')

p0 <- DotPlot(seurat.macro, features = rev(sele.genes), 
              cols = c('grey100', 'red3'), group.by = 'cell_state') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab('') + ylab('')


sele.genes1 = c('CD68', 'CD163', 'MRC1', 
                'THY1', 'CYP27A1', 'HS3ST2', 'LYVE1','F13A1', 
                'CCL4','IL18', 'VEGFA',  'VCAN', 'SPP1', 'C1QC', 
                'MKI67', 'TOP2A')
p1 <- DotPlot(seurat.macro, features = rev(sele.genes1), 
              cols = c('grey100', 'red3'), group.by = 'cell_state') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab('') + ylab('')


sele.genes2 = c('S100A8', 'S100A9', 'FCGR3A', 'CD14', 'CD68', 'CD163', 'MRC1', 
                'THY1', 'CYP27A1', 'HS3ST2', 'LYVE1','F13A1', 
                'CCL4','IL18', 'VEGFA',  'VCAN', 'SPP1', 'C1QC', 
                'MKI67', 'TOP2A')
p2 <- DotPlot(seurat.macro, features = rev(sele.genes1), 
              cols = c('grey100', 'red3'), group.by = 'cell_state') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab('') + ylab('')

## composition across time ####
### < overall composition ####
mdata = seurat.macro@meta.data
mdata = data.table(subset(mdata, select = c('cell_state', 'Stage_Code')), 
                   keep.rownames = T)
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]
mdata = mdata[cell_state != 'CD1C+cDC_6']
PTX_expt = length(which(mdata$Stage_Code == 'PTX'))/nrow(mdata)
mdata$cell_state = as.character(mdata$cell_state)
mdata[, 'n' := .N, by = list(cell_state, Stage_Code)]
mdata[, 'N' := .N, by = cell_state]
mdata[, 'percentage' := n/N]
mdata = subset(mdata, select = c( 'cell_state', 'Stage_Code', 'percentage')) %>% unique()

p2 <- ggplot(data = mdata, aes(x = cell_state, y = percentage, fill = Stage_Code)) +
  geom_bar(stat = 'identity') + theme_classic() + 
  scale_fill_brewer(palette = 'Set1', direction = -1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  xlab('') + ylab('Fraction') + guides(fill=guide_legend(title="")) +
  geom_hline(yintercept = PTX_expt,
             linetype = 'dashed')

## cell composition stacked for each timepoint
mdata = data.table(seurat.macro@meta.data, select = c('cell_state', 'Stage_Code'))
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]
mdata = mdata[cell_state != 'CD1C+cDC_6']
mdata[, 'n' := .N, by = list(cell_state, Stage_Code)]
mdata[, 'N' := .N, by = Stage_Code]
mdata[, 'percentage' := n/N]
mdata = subset(mdata, select = c( 'cell_state', 'Stage_Code', 'percentage')) %>% unique()

p3 <- ggplot(data = mdata, aes(x = Stage_Code, y = percentage, fill = cell_state)) +
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  scale_fill_manual(values = maColors) +
  xlab('') + ylab('Fraction') + guides(fill=guide_legend(title=""))


### < patient level composition ####
mdata = seurat.macro@meta.data
mdata = data.table(mdata)

mdata = mdata %>% subset(select = c(Patient_No, Stage_Code, cell_state))
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]
mdata = mdata[cell_state != 'CD1C+cDC_6']
mdata[, 'N' := .N, by = list(Patient_No, Stage_Code)]

mdata1 = mdata[Stage_Code == 'DX']
mdata2 = mdata[Stage_Code == 'PTX']
mdata1[, 'n' := .N, by = list(Patient_No, cell_state)]
mdata2[, 'n' := .N, by = list(Patient_No, cell_state)]

mdata = rbind(mdata1, mdata2)
mdata[, 'frac' := n/N]
mdata = subset(mdata, select = c(Stage_Code, cell_state, frac, n, N, Patient_No)) %>% unique()
mdata = mdata[N > 20]

plist = list()
for(cl0 in unique(mdata$cell_state)){
  mdata0 = mdata[cell_state == cl0]
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
    annotate('text', x = 1.5, y = 0.15, label = paste0('p = ', format(pv_cl, digits = 2)))
  cl1 = gsub(' ', '_', cl0)
  cl1 = gsub('+', '_', cl1, fixed = T)
  ggsave(pp, filename = paste0('Figures/RNA/OtherCellTypes/Macrophage/', cl1, '_frac_byTimepoint_rna.pdf'), device = 'pdf',
         width = 4, height = 4)
  plist[[cl0]] = pp
}


## misc plots ####
p00 <- DimPlot(seurat.macro, label = T) 

p02 <- VlnPlot(subset(seurat.macro, Stage_Code %in% c('DX', 'PTX')), features = 'immunosuppr_score1', 
               pt.size = 0, split.by = 'Stage_Code') +
  geom_hline(yintercept=0, linetype="dashed", color = "red", size=1) +
  scale_fill_brewer(palette = 'Set1', direction = -1)

p03 <- VlnPlot(subset(seurat.macro, Stage_Code %in% c('DX', 'PTX')), features = 'proinflamm_score1', 
               pt.size = 0, split.by = 'Stage_Code') +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=1) +
  scale_fill_brewer(palette = 'Set1', direction = -1)

p04 <- VlnPlot(subset(seurat.macro, Stage_Code %in% c('DX', 'PTX')), features = 'IL18', 
               pt.size = 0.1, split.by = 'Stage_Code')  +
  scale_fill_brewer(palette = 'Set1', direction = -1)

p05 <- VlnPlot(subset(seurat.macro, Stage_Code %in% c('DX', 'PTX')), features = 'MERTK', 
               pt.size = 0.1, split.by = 'Stage_Code')  +
  scale_fill_brewer(palette = 'Set1', direction = -1)


pp <- VlnPlot(subset(seurat.macro, Stage_Code %in% c('DX', 'PTX')), features = 'VCAN', 
               pt.size = 0.1, split.by = 'Stage_Code', group.by = 'cell_state')  +
  scale_fill_brewer(palette = 'Set1', direction = -1)


## degs btw DX and PTX within cluster ####

degs = NULL
for(cl0 in unique(seurat.macro$cell_state)){
  seurat0 = subset(seurat.macro, cell_state == cl0 & Stage_Code %in% c('DX', 'PTX'))
  Idents(seurat0) = 'Stage_Code'
  degs0 = FindAllMarkers(seurat0,  max.cells.per.ident = 500, test.use = 'LR',
                         latent.vars = c('nCount_RNA', 'percent.mt'), only.pos = T)
  degs0 = data.table(degs0)
  degs0$cell_state = cl0
  degs = rbind(degs, degs0)
}
degs = degs[p_val < 0.001 & avg_log2FC > 0.25 & (pct.1 - pct.2) > 0]
saveRDS(degs, file = 'data/intermediate/RNA/new_degs/degs_PTXvsDX_withinNewClusters_macro.rds')
