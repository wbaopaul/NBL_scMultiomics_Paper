library(Seurat)
library(data.table)
library(liana)
library(magrittr)
library(ggplot2)

## prepare input: just need to run once ####
dir_input = 'data/intermediate/liana/'
dir.create(dir_input, showWarnings = F, recursive = T)
set.seed(2024)

seurat.macro = readRDS(file = 'seurat_macrophage_byHarmony_clean2_cDCremoved.rds')
DefaultAssay(seurat.macro) = 'RNA'

seurat.rna = readRDS('seurat_rna_allCells_clean.rds')
DefaultAssay(seurat.rna) = 'RNA'

seurat.malignant = readRDS('seurat_rna_malignant_harmony_cellStateAnnotated.rds')
mdata1 = seurat.malignant@meta.data
mdata1$cell_state = as.character(mdata1$cell_state)
rm(seurat.malignant)

seurat.rna <- AddMetaData(seurat.rna, metadata = subset(mdata1, select = 'cell_state' ))
seurat.rna$cell_state = ifelse(!is.na(seurat.rna$cell_state), seurat.rna$cell_state, 
                               seurat.rna$cell_type)
seurat.rna = subset(seurat.rna, cell_state %notin% c('Neuroblast', 'Kidney_Cell', 'Hepatocyte', 
                                                     'Adrn_Cortex'))  #rm non-malignant Neuroblast

seurat.rna$fpn = ifelse(seurat.rna$cell_state %notin% c('ADRN-Baseline', 'ADRN-Calcium', 'ADRN-Dopaminergic',
                                                        'ADRN-Proliferating', 'Interm-OxPhos', 'MES') & 
                          seurat.rna$malignancy == "Malignant", 1, 0)
seurat.rna = subset(seurat.rna, fpn == 0)

seurat.rna$cell_bc = colnames(seurat.rna)
seurat.rna$cell_state = ifelse(seurat.rna$cell_bc %in% colnames(seurat.macro),
                               seurat.macro$cell_state[seurat.rna$cell_bc],
                               seurat.rna$cell_state)
seurat.rna = subset(seurat.rna, cell_state %notin% c('Macrophage')) ## subsets included
seurat.rna$cell_state = gsub('Macro', '', seurat.rna$cell_state) # shorter name

Idents(seurat.rna) = seurat.rna$cell_state
seurat.rna <- subset(seurat.rna, downsample = 5000)

saveRDS(seurat.rna, file = paste0(dir_input, 'seurat_tumor_macro4liana.rds'))


## Run liana ####
##Obtain all purely curated resources
seurat.rna = readRDS(file = paste0(dir_input, 'seurat_tumor_macro4liana.rds'))
liana_res <- liana_wrap(seurat.rna, idents_col = 'cell_state',
                        resource = 'all')

saveRDS(liana_res, file = 'data/intermediate/liana/liana_res_allResource.rds')

## post analysis -- plot top10 LR for VCAN+macro ####
liana_res0 = readRDS('data/intermediate/liana/liana_res_allResource.rds')
liana_res0 %>% dplyr::glimpse()

for(rs0 in names(liana_res0[[1]])[-c(1:2)]){
  liana_res <- liana_res0 %>% liana_aggregate(resource = rs0) 
  
  liana_res %>% dplyr::glimpse()
  
  sele.res <- liana_res %>% subset(target %in% c("ADRN-Baseline", "ADRN-Calcium", 'ADRN-Dopaminergic',
                                                 'ADRN-Proliferating', 'Interm-OxPhos', 'MES'))
  sele.res = data.table(sele.res)
  
  p0 <- sele.res %>%
    liana_dotplot(target_groups = c("ADRN-Baseline", "ADRN-Calcium", 'ADRN-Dopaminergic',
                                    'ADRN-Proliferating', 'Interm-OxPhs', 'MES'),
                  source_groups = c("VCAN+"),
                  ntop = 10) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(p0, file = paste0('Figures/RNA/LIANA/liana_top10_', rs0, '.pdf'), 
         device = 'pdf', width = 9, height = 7)
}

