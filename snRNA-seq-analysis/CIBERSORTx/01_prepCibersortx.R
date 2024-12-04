library(data.table)
library(Seurat)
`%notin%` = Negate(`%in%`)


## < Cangelosi bulk data####
R2_data = fread('data/raw/R2/ps_avgpres_dgc2102a786_dgc2102_box1695308502-datagrabber-.txt')
R2_clinical = t(R2_data[1:10, -c(1,2)])
R2_clinical = data.frame(R2_clinical)
names(R2_clinical) = R2_data[1:10, ]$probeset
R2_clinical = data.table(R2_clinical, keep.rownames = T)
names(R2_clinical)[1] = 'sample'
R2_clinical = R2_clinical[platform != 'rnaseq']
R2_rna = R2_data[-c(1:10), -(1:2)]
R2_rna = apply(R2_rna, 2, as.numeric)
rownames(R2_rna) = R2_data$`#H:hugo`[-(1:10)]
R2_rna = R2_rna[, R2_clinical$sample]/2

cang_rna0 = data.table(data.frame(R2_rna), keep.rownames = T)
names(cang_rna0)[1] = 'Gene'

fwrite(cang_rna0, file = 'data/intermediate/RNA/CIBERSORTx/cangelosi_mtx.txt', 
       sep = '\t', quote = F)



## < Seqc bulk data ####
seqc_data = fread('data/raw/R2/ps_avgpres_gse62564geo498_seqcnb1_box1701183946-datagrabber-rpm.txt')
seqc_clinical = t(seqc_data[1:16, -c(1,2)])
seqc_clinical = data.frame(seqc_clinical)
names(seqc_clinical) = seqc_data[1:16, ]$probeset
seqc_clinical = data.table(seqc_clinical, keep.rownames = T)
names(seqc_clinical)[1] = 'sample'
seqc_rna = seqc_data[-c(1:16), -(1:2)]
seqc_rna = apply(seqc_rna, 2, as.numeric)
rownames(seqc_rna) = seqc_data$`#H:hugo`[-(1:16)]

seqc_rna0 = data.table(data.frame(seqc_rna), keep.rownames = T)
names(seqc_rna0)[1] = 'Gene'

fwrite(seqc_rna0, file = 'data/intermediate/RNA/CIBERSORTx/seqc_mtx.txt', 
       sep = '\t')


## < single-cell ####

### < include normal cell types, distinct macrophage and malignant subsets  ####
seurat.macro = readRDS('seurat_macrophage_byHarmony_clean2_cDCremoved.rds')

seurat.rna = readRDS('seurat_rna_allCells_clean.rds')
DefaultAssay(seurat.rna) = 'RNA'
seurat.rna[['integrated']] <- NULL

seurat.malignant = readRDS('seurat_rna_malignant_harmony_cellStateAnnotated.rds')
mdata1 = seurat.malignant@meta.data
mdata1$cell_state = as.character(mdata1$cell_state)

seurat.rna <- AddMetaData(seurat.rna, metadata = subset(mdata1, select = 'cell_state' ))
seurat.rna$cell_state = ifelse(!is.na(seurat.rna$cell_state), seurat.rna$cell_state, 
                               seurat.rna$cell_type)

seurat.rna = subset(seurat.rna, cell_state %notin% c('Neuroblast', 'Kidney_Cell', 'Hepatocyte', 
                                                     'Adrenal_Cortex'))  #rm neuroblast because those maybe normal

## remove potential false positives
seurat.rna$fpn = ifelse(seurat.rna$cell_state %notin% c('ADRN-Baseline', 'ADRN-Calcium', 'ADRN-Dopaminergic',
                                                        'ADRN-Proliferating', 'Interm-OxPhos', 'MES') & 
                          seurat.rna$malignancy == "Malignant", 1, 0) 
seurat.rna = subset(seurat.rna, fpn == 0)

Idents(seurat.rna) = seurat.rna$cell_state

seurat.rna$cell_bc = colnames(seurat.rna)
seurat.rna$cell_state = ifelse(seurat.rna$cell_bc %in% colnames(seurat.macro),
                               seurat.macro$cell_state[seurat.rna$cell_bc],
                               seurat.rna$cell_state)
seurat.rna = subset(seurat.rna, cell_state %notin% c('Macrophage')) ## subsets alreday included

seurat.rna <- subset(seurat.rna, downsample = 10000)
seurat.rna <- NormalizeData(seurat.rna, normalization.method = 'RC',
                            scale.factor = 1e6)
saveRDS(seurat.rna, file = 'seurat_rna4Cibersortx_eachGroupDownTo10kCells_Data2.rds')





### < to adjust format to be compatible with CIBERSORTx ####
sc_file = 'seurat_rna4Cibersortx_eachGroupDownTo10kCells_combData2.rds'
sc_file_out = paste0('data/intermediate/RNA/CIBERSORTx/', basename(sc_file))
sc_file_out = gsub('.rds', '.txt', sc_file_out, fixed = T)

mdata <- seurat.rna@meta.data
normalized.mtx <- seurat.rna[['RNA']]@data
table(mdata$cell_state)
rm(seurat.rna)
colnames(normalized.mtx) <- mdata$cell_state
normalized.data <- as.data.frame(normalized.mtx)
normalized.data <- data.table(normalized.data, keep.rownames = T)
names(normalized.data) <- c('GeneSymbol', names(normalized.data)[-1])
fwrite(normalized.data, file = sc_file_out, sep = '\t', quote = F)


