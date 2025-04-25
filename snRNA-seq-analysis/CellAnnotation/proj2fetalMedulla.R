library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
library(patchwork)
library(Matrix)
library(matrixStats)
library(RColorBrewer)
library(uwot)
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = 'v5')

## project NCC onto fetal medulla ref ####
ref_type = 'adrenal_medulla' # adrenal_gland or adrenal_medulla

seurat.query <- readRDS('Seurat_Objects/seurat_rna_allCells_clean.rds')


seurat.ref = readRDS('Seurat_Objects/adrenal_medulla_Seurat.RDS')

seurat.ref = UpdateSeuratObject(seurat.ref)
seurat.ref$cell_type = Idents(seurat.ref)
DimPlot(seurat.ref)

seurat.ref = RunUMAP(seurat.ref, reduction = 'pca', dims = 1:20, 
                     return.model = T, seed.use = 42, umap.method = 'uwot',
                     n.neighbors = 30, metric = 'cosine', learning.rate = 1,
                     min.dist = 0.3, spread = 1, set.op.mix.ratio = 1,
                     local.connectivity = 1, negative.sample.rate = 5,
                     uwot.sgd = F, angular.rp.forest = F, verbose = T,
                     reduction.name = 'umap', reduction.key = 'UMAP_')
DimPlot(seurat.ref)
#seurat.ref@reductions$umap@cell.embeddings = -seurat.ref@reductions$umap@cell.embeddings

### using seurat MapQuery ####

DefaultAssay(seurat.query) = 'RNA'
seurat.query = NormalizeData(seurat.query, verbose = F)

anchor <- FindTransferAnchors(
  reference = seurat.ref,
  query = seurat.query,
  reference.reduction = "pca",
  dims = 1:50
)
seurat.query <- MapQuery(
  anchorset = anchor,
  query = seurat.query,
  reference = seurat.ref,
  refdata = list(ctype_fetal_adrenal = "cell_type"),
  reference.reduction = "pca",
  reference.dims = 1:50,
  reduction.model = "umap"
)

DimPlot(seurat.query, group.by = 'predicted.ctype_fetal_adrenal', reduction = 'ref.umap') 

saveRDS(seurat.query@meta.data, 
        file = paste0('MetaData/metadata_allCells_projected_ctype_fetal_',
                      ref_type, '.rds'))

saveRDS(seurat.query, file = 'Seurat_Objects/seurat_rna_projected_ctype_fetal_medulla.rds')

## plot projected umap ####
DimPlot(seurat.query, group.by = 'predicted.ctype_fetal_adrenal', reduction = 'ref.umap') 
seurat.query@reductions$ref.umap@cell.embeddings = -seurat.query@reductions$ref.umap@cell.embeddings
DimPlot(seurat.query, group.by = 'predicted.ctype_fetal_adrenal', reduction = 'ref.umap') 

mdata <- readRDS('MetaData/metadata_seurat_rna_malignant_harmony_cellStateAnnotated.rds')
seurat.query$barcode = colnames(seurat.query)
seurat.query = subset(seurat.query, barcode %in% rownames(mdata))
seurat.query = AddMetaData(seurat.query, metadata = subset(mdata, select = cell_state))

p0 <- DimPlot(seurat.query, pt.size = 1.5,
        group.by = 'predicted.ctype_fetal_adrenal', 
        reduction = 'ref.umap', split.by = 'cell_state') + 
  scale_color_brewer(palette = 'Set3') + NoLegend()


seurat.ref@reductions$umap@cell.embeddings = -seurat.ref@reductions$umap@cell.embeddings
seurat.ref$cell_type = as.character(seurat.ref$cell_type)

p1 <- DimPlot(seurat.ref, group.by = 'cell_type') + 
  scale_color_brewer(palette = 'Set3') + NoLegend()

ggsave(p0, filename = 'Figures/RNA/Malignant/btw_clusters/umap_projected_cell_states.pdf',
       device = 'pdf', width = 15, height = 5)

ggsave(p1, filename = 'Figures/RNA/Malignant/btw_clusters/umap_ref.pdf',
       device = 'pdf', width = 3, height = 5)


