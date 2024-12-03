library(Seurat)
library(data.table)
library(magrittr)
library(anndata)
library(homologene)
library(ggplot2)
library(harmony)
library(RColorBrewer)
library(Matrix)
library(patchwork)



## for malignant cell annotation ####

### prepare human reference data ####
seurat.malignant = readRDS(file = 'seurat_rna_malignant_harmony_cellStateAnnotated.rds')
counts <- seurat.malignant[['RNA']]$data
mdata = seurat.malignant@meta.data
rm(seurat.malignant)

## only keep human homologues xenium available genes
hm = homologene(rownames(counts), outTax = 10090, inTax = 9606)
names(hm)[1:2] = c('human_gene', 'mouse_gene')
hm = hm[, -c(3:4)]
hm = data.table(hm)
hm = hm[!duplicated(mouse_gene),] ## just keep one for duplicates

miss_genes = setdiff(rownames(counts), hm$human_gene)

hm = hm[!duplicated(human_gene), ] ## no duplicate human gene as well
gene_map = rbind(hm, data.table('human_gene' = miss_genes,
                                'mouse_gene' = stringr::str_to_title(miss_genes) ))
gene_map = gene_map[!duplicated(mouse_gene)]
setkey(gene_map, human_gene)

counts = counts[gene_map$human_gene, ]
all(rownames(counts) == gene_map$human_gene)
rownames(counts) = gene_map$mouse_gene

## xenium data
xadata = read_h5ad('nbl_xenium_anndata_final.h5ad')
xcounts = xadata$layers[['counts']]
xmdata = xadata$obs
xcounts = xcounts[xmdata$cell_type == 'Neuroblast', ]
xmdata = xmdata[xmdata$cell_type == 'Neuroblast', ]

rm(adata, xadata)

shared.genes = intersect(rownames(counts), colnames(xcounts))
counts = (counts[shared.genes, ])
xcounts = t(xcounts[, shared.genes])

### cluster and integration with xenium available genes ####

seu_obj <- CreateSeuratObject(counts, meta.data = mdata)
seu_obj_query <- CreateSeuratObject(xcounts, meta.data = xmdata)

seu_obj[['RNA']]$data = seu_obj[['RNA']]$counts ## since it's normalized

VariableFeatures(seu_obj) <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj) %>% RunPCA(npcs = 30, verbose = F) %>% 
  FindNeighbors(dims = 1:30)
seu_obj <- RunUMAP(seu_obj, verbose = F, dims = 1:30)

DimPlot(seu_obj, group.by = 'Patient_No')
DimPlot(seu_obj, group.by = 'cell_state')

seu_obj <- RunHarmony(seu_obj, 'Patient_No') %>% 
  FindNeighbors(dims = 1:30, reduction = 'harmony') 

seu_obj <- RunUMAP(seu_obj, verbose = F, dims = 1:30, 
                   reduction = 'harmony', return.model = T)

## copy harmony as pca
seu_obj@reductions$pca@cell.embeddings = seu_obj@reductions$harmony@cell.embeddings
seu_obj@reductions$pca@feature.loadings.projected = seu_obj@reductions$harmony@feature.loadings.projected
seu_obj@reductions$pca@feature.loadings = seu_obj@reductions$harmony@feature.loadings
colnames(seu_obj@reductions$pca@cell.embeddings) = paste0('PC_', 1:30)
colnames(seu_obj@reductions$pca@feature.loadings) = paste0('PC_', 1:30)
colnames(seu_obj@reductions$pca@feature.loadings.projected) = paste0('PC_', 1:30)


seu_obj_query = NormalizeData(seu_obj_query, 
                              scale.factor = median(seu_obj_query$nCount_RNA))
VariableFeatures(seu_obj_query) <- rownames(seu_obj_query)
seu_obj_query <- ScaleData(seu_obj_query) %>% RunPCA(npcs = 30, verbose = F)


save(seu_obj, seu_obj_query, file = 'seurat_tmp.RData')  



## for macrophage subset annotation ####
### prepare human reference data ####
seurat.macro = readRDS(file = 'seurat_macrophage_byHarmony_clean2_cDCremoved.rds')
counts <- seurat.macro[['RNA']]$datq
mdata = seurat.macro@meta.data
rm(seurat.macro)

## only keep human homologues xenium available genes
hm = homologene(rownames(counts), outTax = 10090, inTax = 9606)
names(hm)[1:2] = c('human_gene', 'mouse_gene')
hm = hm[, -c(3:4)]
hm = data.table(hm)
hm = hm[!duplicated(mouse_gene),] ## just keep one for duplicates

miss_genes = setdiff(rownames(counts), hm$human_gene)

hm = hm[!duplicated(human_gene), ] ## no duplicate human gene as well
gene_map = rbind(hm, data.table('human_gene' = miss_genes,
                                'mouse_gene' = stringr::str_to_title(miss_genes) ))
gene_map = gene_map[!duplicated(mouse_gene)]
setkey(gene_map, human_gene)

counts = counts[gene_map$human_gene, ]
all(rownames(counts) == gene_map$human_gene)
rownames(counts) = gene_map$mouse_gene


## xenium data
xadata = read_h5ad('Revision/nbl_xenium_anndata_final.h5ad')
xcounts = xadata$layers[['counts']]
xmdata = xadata$obs
xcounts = xcounts[xmdata$cell_type == 'Macrophage', ]
xmdata = xmdata[xmdata$cell_type == 'Macrophage', ]
rm(xadata)

shared.genes = intersect(rownames(counts), colnames(xcounts))
counts = counts[shared.genes, ]
xcounts = t(xcounts[, shared.genes])


seu_obj <- CreateSeuratObject(counts, meta.data = mdata)
seu_obj_query <- CreateSeuratObject(xcounts, meta.data = xmdata)

### cluster and integration with xenium available genes ####

seu_obj[['RNA']]$data = seu_obj[['RNA']]$counts ## since it's normalized

VariableFeatures(seu_obj) <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj) %>% RunPCA(npcs = 30, verbose = F) %>% 
  FindNeighbors(dims = 1:30)
seu_obj <- RunUMAP(seu_obj, verbose = F, dims = 1:30)

DimPlot(seu_obj, group.by = 'Patient_No')
DimPlot(seu_obj, group.by = 'cell_state')

seu_obj <- RunHarmony(seu_obj, 'Patient_No') %>% 
  FindNeighbors(dims = 1:30, reduction = 'harmony') 

seu_obj <- RunUMAP(seu_obj, verbose = F, dims = 1:30, 
                   reduction = 'harmony', return.model = T)

DimPlot(seu_obj, group.by = 'Patient_No')
DimPlot(seu_obj, group.by = 'cell_state')

## copy harmony as pca
seu_obj@reductions$pca@cell.embeddings = seu_obj@reductions$harmony@cell.embeddings
seu_obj@reductions$pca@feature.loadings.projected = seu_obj@reductions$harmony@feature.loadings.projected
seu_obj@reductions$pca@feature.loadings = seu_obj@reductions$harmony@feature.loadings
colnames(seu_obj@reductions$pca@cell.embeddings) = paste0('PC_', 1:30)
colnames(seu_obj@reductions$pca@feature.loadings) = paste0('PC_', 1:30)
colnames(seu_obj@reductions$pca@feature.loadings.projected) = paste0('PC_', 1:30)


seu_obj_query = NormalizeData(seu_obj_query, 
                              scale.factor = median(seu_obj_query$nCount_RNA))
VariableFeatures(seu_obj_query) <- rownames(seu_obj_query)
seu_obj_query <- ScaleData(seu_obj_query) %>% RunPCA(npcs = 30, verbose = F)
seu_obj_query <- RunHarmony(seu_obj_query, 'mouse_id') %>% 
  RunUMAP(dims = 1:30, reduction = 'harmony') 


save(seu_obj, seu_obj_query, file = 'seurat_macros.RData')  


