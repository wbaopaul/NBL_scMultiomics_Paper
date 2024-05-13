## transfer label from scRNA to scATAC using Seurat 

library(Signac)
library(Seurat)

seuratAtacPath = 'Seurat_Objects/seurat_atac_signac_byPatient.rds'

reduction2use = 'integrated_lsi'

seurat.rna = readRDS('Seurat_Objects/seurat_rna_allCells_clean.rds')
seurat.rna[['integrated']] = as(seurat.rna[['integrated']], Class = 'Assay5') ## redundant

## downsample rna cells
set.seed(2022)
seurat.rna$bc = colnames(seurat.rna)
sele.bcs = sample(colnames(seurat.rna), 100000)
seurat.rna <- subset(seurat.rna, bc %in% sele.bcs)

seurat.atac = readRDS(seuratAtacPath)
## use GAS = promote + gene body accessibility ####
activity.matrix = readRDS('Seurat_Objects/mtx_gene_activity_signac.rds')
seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
DefaultAssay(seurat.atac) = 'ACTIVITY'
seurat.atac <- NormalizeData(
  object = seurat.atac,
  assay = 'ACTIVITY',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat.atac$nCount_ACTIVITY)
)
seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))

seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

## transfer label 
DefaultAssay(seurat.rna) = 'integrated'
genes4anchors = VariableFeatures(object = seurat.rna)
genes4anchors = intersect(genes4anchors, rownames(seurat.atac))
message(paste('Final # of variable genes to use:', length(genes4anchors)))

transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = "integrated",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 20)


celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.rna$cell_type,
                                     weight.reduction = seurat.atac[[reduction2use]],
                                     dims = 2:ncol(seurat.atac[[reduction2use]]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, select = c('predicted.id', 'prediction.score.max'))
seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)


p1 <- DimPlot(seurat.atac, group.by = "predicted.id") 

ggsave(p1, filename = paste0(seuratAtacPath, '_labelTransfer_umap.pdf'), device = 'pdf',
       width = 9, height = 7)
message('Label trasfer done! Now saving data...')

## save the label transfer results
mdata = subset(seurat.atac@meta.data, select = c('predicted.id', 'prediction.score.max'))
saveRDS(mdata, file = 'MetaData/labelTransferResult_seurat_atac_signac_byPatient.rds')
message('Saving data done!')


