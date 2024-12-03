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

args = commandArgs(T)
mid = args[1] ## mouse_id
celltype = args[2] ## Neurobalst or Macrophage
mdata = NULL

if(celltype == 'Neuroblast') load('seurat_tmp.RData')
if(celltype == 'Macrophage') load('seurat_macros.RData')


seu_obj_query0 = subset(seu_obj_query, mouse_id == mid)
seu_obj_query0 = ScaleData(seu_obj_query0)
seu_obj_query0 = RunPCA(seu_obj_query0, verbose = F, npcs = 30)

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = seu_obj, query = seu_obj_query0, 
                                        features = VariableFeatures(object = seu_obj),
                                        reference.assay = "RNA", query.assay = "RNA",
                                        reduction = "rpca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = seu_obj$cell_state,
                                     weight.reduction = seu_obj_query0[["pca"]],
                                     dims = 1:30)

seu_obj_query0 <- AddMetaData(seu_obj_query0, metadata = celltype.predictions)
mdata0 = seu_obj_query0@meta.data
saveRDS(mdata0, 
        file = paste0('MetaData/Xenium/seurat_xenium_labelTransfer_mouse', mid, '_rpca_', celltype, '.rds'))


