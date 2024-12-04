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


stateColor = brewer.pal(n = 6, name = "Set1")
names(stateColor) = c('ADRN-Calcium', 'ADRN-Baseline', 'Interm-OxPhs',
                      'ADRN-Dopaminergic', 'ADRN-Proliferating', 'MES')

# CHLA15 ####
## DEG treating cell state as covariate ####
seu.ch15.malignant = readRDS(file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA15_malignant.rds')
seu.ch15.malignant$experiment = stringr::str_split_i(seu.ch15.malignant$deMULTIplex2, '_', 3)
seu.ch15.malignant = subset(seu.ch15.malignant, predicted.cell_type.score > 0.5)
Idents(seu.ch15.malignant) = seu.ch15.malignant$experiment

degs1 = FindMarkers(seu.ch15.malignant, max.cells.per.ident = 200,
                    min.pct = 0.05, test.use = 'LR', 
                    logfc.threshold = 0.1,
                    ident.1 = 'Coculture', ident.2 = 'Monoculture',
                    latent.vars = 'predicted.cell_type')
degs1 = data.table(degs1, keep.rownames = T)
names(degs1)[1] = 'gene'
degs1$comparison = 'Coculture_vs_Monoculture'
degs1 = degs1[p_val < 0.05]
saveRDS(degs1, file = 'data/intermediate/coculture/CHLA15_degs_Coculture_vs_Monoculture_cellStateAsCovariate.rds')


degs2 = FindMarkers(seu.ch15.malignant, max.cells.per.ident = 200,
                    min.pct = 0.05, test.use = 'LR', 
                    logfc.threshold = 0.1,
                    ident.1 = 'Afatinib', ident.2 = 'Coculture',
                    latent.vars = 'predicted.cell_type')
degs2 = data.table(degs2, keep.rownames = T)
names(degs2)[1] = 'gene'
degs2$comparison = 'Afatinib_vs_Coculture'
degs2 = degs2[p_val < 0.05]
saveRDS(degs2, file = 'data/intermediate/coculture/CHLA15_degs_Afatinib_vs_Coculture_cellStateAsCovariate.rds')


degs3 = FindMarkers(seu.ch15.malignant, max.cells.per.ident = 200,
                    min.pct = 0.05, test.use = 'LR', 
                    logfc.threshold = 0.1,
                    ident.1 = 'CRM', ident.2 = 'Coculture',
                    latent.vars = 'predicted.cell_type')
degs3 = data.table(degs3, keep.rownames = T)
names(degs3)[1] = 'gene'
degs3$comparison = 'CRM_vs_Coculture'
degs3 = degs3[p_val < 0.05]
saveRDS(degs3, file = 'data/intermediate/coculture/CHLA15_degs_CRM_vs_Coculture_cellStateAsCovariate.rds')


# CHLA20 ####
## DEG treating cell state as covariate ####
seu.ch20.malignant = readRDS(file = 'Seurat_Objects/coculture/seurat_deMULTIplex2_CHLA20_malignant.rds')
seu.ch20.malignant$experiment = stringr::str_split_i(seu.ch20.malignant$deMULTIplex2, '_', 3)
seu.ch20.malignant = subset(seu.ch20.malignant, predicted.cell_type.score > 0.5)
Idents(seu.ch20.malignant) = seu.ch20.malignant$experiment

degs1 = FindMarkers(seu.ch20.malignant, max.cells.per.ident = 200,
                    min.pct = 0.05, test.use = 'LR', 
                    logfc.threshold = 0.1,
                    ident.1 = 'Coculture', ident.2 = 'Monoculture',
                    latent.vars = 'predicted.cell_type')
degs1 = data.table(degs1, keep.rownames = T)
names(degs1)[1] = 'gene'
degs1$comparison = 'Coculture_vs_Monoculture'
degs1 = degs1[p_val < 0.05]
saveRDS(degs1, file = 'data/intermediate/coculture/CHLA20_degs_Coculture_vs_Monoculture_cellStateAsCovariate.rds')


degs2 = FindMarkers(seu.ch20.malignant, max.cells.per.ident = 200,
                    min.pct = 0.05, test.use = 'LR', 
                    logfc.threshold = 0.1,
                    ident.1 = 'Afatinib', ident.2 = 'Coculture',
                    latent.vars = 'predicted.cell_type')
degs2 = data.table(degs2, keep.rownames = T)
names(degs2)[1] = 'gene'
degs2$comparison = 'Afatinib_vs_Coculture'
degs2 = degs2[p_val < 0.05]
saveRDS(degs2, file = 'data/intermediate/coculture/CHLA20_degs_Afatinib_vs_Coculture_cellStateAsCovariate.rds')


degs3 = FindMarkers(seu.ch20.malignant, max.cells.per.ident = 200,
                    min.pct = 0.05, test.use = 'LR', 
                    logfc.threshold = 0.1,
                    ident.1 = 'CRM', ident.2 = 'Coculture',
                    latent.vars = 'predicted.cell_type')
degs3 = data.table(degs3, keep.rownames = T)
names(degs3)[1] = 'gene'
degs3$comparison = 'CRM_vs_Coculture'
degs3 = degs3[p_val < 0.05]
saveRDS(degs3, file = 'data/intermediate/coculture/CHLA20_degs_CRM_vs_Coculture_cellStateAsCovariate.rds')



