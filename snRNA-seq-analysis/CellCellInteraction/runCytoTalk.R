library(Seurat)
library(data.table)
library(CytoTalk)

## prepare input: just need to run once
dir_input = 'data/intermediate/cytotalk/secRound/inputs/'
dir.create(dir_input, showWarnings = F, recursive = T)
set.seed(2023)
if(T){
  seurat.adrn = readRDS('Seurat_Objects/seurat_rna_malignant_harmony_clean.rds')
   
  seurat.macro = readRDS('Seurat_Objects/seurat_macrophage_byHarmony_clean2_res0.4.rds')
  
  DefaultAssay(seurat.adrn) <- 'RNA'
  DefaultAssay(seurat.macro) <- 'RNA'
  
  seurat.adrn = subset(seurat.adrn, downsample = 5000)
  seurat.macro = subset(seurat.macro, downsample = 2000) 
  
  adrn_mtx = seurat.adrn@assays$RNA@layers$data
  rownames(adrn_mtx) = rownames(seurat.adrn@assays$RNA@features)
  colnames(adrn_mtx) = rownames(seurat.adrn@assays$RNA@cells)
  
  adrn_mdata = seurat.adrn$seurat_clusters
  
  macro_mtx = seurat.macro@assays$RNA@data
  macro_mdata = seurat.macro$seurat_clusters
  all(rownames(adrn_mtx) == rownames(macro_mtx))
  
  ## write matrix per cluster
  for(adrn_cl in 0:5){
    tmp_mtx = adrn_mtx[, seurat.adrn$seurat_clusters == adrn_cl]
    write.csv(tmp_mtx, file = paste0(dir_input, 'scRNAseq_Malignant', adrn_cl, '.csv'),
              quote = F, row.names = T)
  }
  for(macro_cl in c(0:5,7,8)){
    tmp1_mtx = macro_mtx[, seurat.macro$seurat_clusters == macro_cl]
    write.csv(tmp1_mtx, file = paste0(dir_input, 'scRNAseq_Macro', macro_cl, '.csv'),
              quote = F, row.names = T)
  }
}


args = commandArgs(T)
adrn_cl = args[1] ## from 0-5
macro_cl = args[2] ## from 0-6
expr_cutoff = 0.1 
lst_scrna <- CytoTalk:::read_matrix_folder(dir_input)

# set required parameters
adrn_cl = paste0('Malignant', adrn_cl)
macro_cl = paste0('Macro', macro_cl)

# run CytoTalk process
dir_res = paste0('data/intermediate/cytotalk/secRound/exprCutoff', expr_cutoff, '/result_', adrn_cl, '_', macro_cl)
dir.create(dir_res, showWarnings = F)
results <- CytoTalk::run_cytotalk(lst_scrna, adrn_cl, macro_cl,
                                  cutoff_a = expr_cutoff, cutoff_b = expr_cutoff, beta_max = 100,
                                  depth = 5, 
                                  dir_out = dir_res)

saveRDS(results, 
        file = paste0(dir_res, '/resultObj_', adrn_cl, '_', macro_cl, '.rds'))

message(paste0(adrn_cl, '-', macro_cl, 'pair Done!'))


