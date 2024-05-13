library(data.table)
library(magrittr)
library(Seurat)
library(limma)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(GSEABase)
library(stringr)
library(RColorBrewer)
library(viridis)
library(openxlsx)
library(edgeR)



## DEG btw clusters (0-5) ####
### single cell level ####
seurat.rna = readRDS('Seurat_Objects/seurat_rna_malignant_harmony_clean.rds')
## *** output all for GSEA analysis*** #
degs <- FindAllMarkers(seurat.rna, 
                       assay = 'RNA', test.use = 'LR', logfc.threshold = 0.0, 
                       latent.vars = c('nCount_RNA', 'percent.mt'),
                       slot = 'data', max.cells.per.ident = 1000)
degs = data.table(degs)
saveRDS(degs, file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')


### EnrichR ####
deg_type = 'sc' 

if(deg_type == 'sc'){
  degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
  degs[, 'pct.diff' := pct.1 - pct.2]
  degs = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]
}



library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
sele.dbs <- c("GO_Biological_Process_2023", 'KEGG_2021_Human', "GO_Biological_Process_2023",
              "Reactome_2022", 'MSigDB_Hallmark_2020', "GO_Molecular_Function_2023")


## should do it state by state ##
cell_cluster1 = 0

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
degs0 = degs[cluster == cell_cluster1]
#if(cell_cluster1 %notin% c(1)) degs0 = degs0[pct.diff > 0.1]

genes0 = unique(degs0$gene)

message(paste0(cell_cluster1, ': ', length(genes0), ' genes'))
enriched_res = enrichr(genes0, sele.dbs)
if(is.null(enriched_res)){
  message('No enriched terms found!')
}else{
  for(db0 in names(enriched_res)){
    res.state0 = data.table(enriched_res[[db0]])
    res.state0[, 'Count' := as.numeric(unlist(strsplit(Overlap, '/', fixed = T))[1]), 
               by = Overlap]
    res.state0 = res.state0[order(P.value)]
    
    res.state0[, 'score' := -log10(P.value)]
    tmp = res.state0[Count >= 3 & P.value < 0.05]
    if(nrow(tmp) == 0) {
      message(paste0(cell_cluster1, ': no sig terms!')) 
      next
    }
    tmp = tmp[order(P.value)]
    nmax = min(15, nrow(tmp))
    pic.height = 8
    if(nmax < 5) pic.height = 2
    if(nmax >= 5 & nmax < 10) pic.height = 5
    
    
    tmp = tmp[1:nmax, ]
    tmp$Term = factor(tmp$Term, levels = rev(tmp$Term))
    p0 <- ggplot(tmp, aes(Term, score, fill = score)) +
      geom_bar(stat = 'identity') +
      coord_flip() + scale_fill_viridis(oPTXon = 'D') +
      labs(x="", title = paste0(db0, ': ', cell_cluster1),
           y = '-log10(pvalue)') + guides(fill=guide_legend(title="-log10(pvalue)")) +
      theme_classic() + theme(axis.text.y = element_text(size = 12))
    #p0 <- plotEnrich(enriched_res, showTerms = 20, numChar = 40, y = "Count", 
    #                 orderBy = "P.value")
    
    dir.create(paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, '/enrichR_', db0),
               showWarnings = F, recursive = T)
    filekey = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, '/enrichR_', db0, '/', 
                     db0, '_', cell_cluster1, '_', deg_type)
    ggsave(p0, filename = paste0(filekey, '.pdf'),
           device = 'pdf', width = 12, height = pic.height)
    fwrite(res.state0[Count >= 3 & P.value < 0.05], sep = '\t',
           file = paste0(filekey, '.tsv'))
  }
  
}

### Plot EnrichR result ####
deg_type = 'sc' # sc/bulk_limma_voom/bulk_edger

if(deg_type == 'sc'){
  degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
  degs[, 'pct.diff' := pct.1 - pct.2]
  degs = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]
}


#### plot dot plot -- KEGG ####
enrichr_res_dir = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                         '/enrichR_KEGG_2021_Human/')
### collect all top terms
sele_terms = NULL
for(cell_state0 in unique(degs$cluster)){
  tmp = fread(paste0(enrichr_res_dir, 'KEGG_2021_Human_', cell_state0, '_', deg_type, '.tsv'))
  tmp$cell_state = cell_state0
  if(cell_state0 == 1){
    tmp = tmp[P.value < 0.05]
  }else{
    tmp = tmp[Adjusted.P.value < 0.05]
  }
  tmp = tmp[1:min(nrow(tmp), 10), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(cell_state0 in unique(degs$cluster)){
  tmp = fread(paste0(enrichr_res_dir, 'KEGG_2021_Human_', cell_state0, '_', deg_type, '.tsv'))
  tmp$cell_state = cell_state0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  if(cell_state0 == 1){
    tmp = tmp[P.value < 0.05]
  }else{
    tmp = tmp[Adjusted.P.value < 0.05]
  }
  tmp$rk = order(-tmp$score)
  tmp[rk>10]$rk=10
  enrichr_res[[cell_state0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'cell_state'))

dt = dcast(dt, Term ~ cell_state, value.var = 'score')

mat = dt[, -1]
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap(mat, na_col = 'white', border_color = NA,
              cluster_rows = T, cluster_cols = T, angle_col = 45,
              main = 'enrichR: KEGG')
ggsave(aa, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type,
                             '/enrichR_KEGG_2021_Human/heatmap_enrichr_kegg_allClusters_', 
                             deg_type,'.pdf'),
       device = 'pdf', width = 10, height = 8)

### dotplot

enrichr_res$Term = factor(enrichr_res$Term, 
                          levels = rev(rownames(mat)[aa$tree_row$order]))
enrichr_res$cell_state = factor(enrichr_res$cell_state, 
                                levels = colnames(mat)[aa$tree_col$order])
enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$size = as.integer(stringr::str_split_i(enrichr_res$Overlap, '/', 2))
enrichr_res$frac = enrichr_res$Count/enrichr_res$size
enrichr_res$frac = 11-enrichr_res$rk
pk <- ggplot(data = enrichr_res, aes(x = cell_state, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(oPTXon = "D") + theme_bw() +
  ggtitle('enrichr: KEGG') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(padj)"), 
         size=guide_legend(title="Order -log10(paj)")) 

ggsave(pk, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                             '/enrichR_KEGG_2021_Human/dotplot_enrichr_kegg_allCellStates_', 
                             deg_type, '.pdf'),
       device = 'pdf', width = 8, height = 10)



#### plot dot plot -- HALLMARK ####

enrichr_res_dir = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, '/enrichR_MSigDB_Hallmark_2020/')
### collect all top terms
sele_terms = NULL
for(cell_state0 in unique(degs$cluster)){
  efile = paste0(enrichr_res_dir, 'MSigDB_Hallmark_2020_', cell_state0, '_', deg_type, '.tsv')
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$cell_state = cell_state0
  if(cell_state0 %in% c(1, 3)){
    tmp = tmp[P.value < 0.05]
  }else{
    tmp = tmp[Adjusted.P.value < 0.05]
  }
  #tmp = tmp[1:min(nrow(tmp), 10), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(cell_state0 in unique(degs$cluster)){
  efile = paste0(enrichr_res_dir, 'MSigDB_Hallmark_2020_', cell_state0, '_', deg_type, '.tsv')
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$cell_state = cell_state0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  if(cell_state0 %in% c(1, 3)){
    tmp = tmp[P.value < 0.05]
  }else{
    tmp = tmp[Adjusted.P.value < 0.05]
  }
  tmp$rk = order(-tmp$score)
  tmp[rk>10]$rk=10
  enrichr_res[[cell_state0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'cell_state'))

dt = dcast(dt, Term ~ cell_state, value.var = 'score')

mat = dt[, -1]
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap(mat, na_col = 'white', border_color = NA,
              cluster_rows = T, cluster_cols = T, angle_col = 45)
ggsave(aa, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type,
                             '/enrichR_MSigDB_Hallmark_2020/heatmap_enrichr_hallmark_allCellStates_',
                             deg_type, '.pdf'),
       device = 'pdf', width = 6, height = 6)
### dotplot

enrichr_res$Term = factor(enrichr_res$Term, levels = rev(rownames(mat)[aa$tree_row$order]))
enrichr_res$cell_state = factor(enrichr_res$cell_state, 
                                levels = colnames(mat)[aa$tree_col$order])

enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$size = as.integer(stringr::str_split_i(enrichr_res$Overlap, '/', 2))
enrichr_res$frac = enrichr_res$Count/enrichr_res$size
enrichr_res$frac = 15 - enrichr_res$rk
ph <- ggplot(data = enrichr_res, aes(x = cell_state, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(oPTXon = "D") + theme_bw() +
  ggtitle('enrichr: HALLMARK') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(padj)"), 
         size=guide_legend(title="Order -log10(padj)")) 

ggsave(ph, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                             '/enrichR_MSigDB_Hallmark_2020/dotplot_enrichr_hallmark_allCellStates_',
                             deg_type, '.pdf'),
       device = 'pdf', width = 6, height = 7)



#### plot dot plot -- GO ####
enrichr_res_dir = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, '/enrichR_GO_Biological_Process_2023/')
### collect all top terms
sele_terms = NULL
for(cell_state0 in unique(degs$cluster)){
  tmp = fread(paste0(enrichr_res_dir, 'GO_Biological_Process_2023_', cell_state0, '_', deg_type, '.tsv'))
  tmp$cell_state = cell_state0
  if(cell_state0 %in% c(1, 3)){
    tmp = tmp[P.value < 0.05]
  }else{
    tmp = tmp[Adjusted.P.value < 0.05]
  }
  tmp = tmp[1:min(nrow(tmp), 8), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(cell_state0 in unique(degs$cluster)){
  tmp = fread(paste0(enrichr_res_dir, 'GO_Biological_Process_2023_', cell_state0, '_', deg_type, '.tsv'))
  tmp$cell_state = cell_state0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  if(cell_state0 %in% c(1, 3)){
    tmp = tmp[P.value < 0.05]
  }else{
    tmp = tmp[Adjusted.P.value < 0.05]
  }
  enrichr_res[[cell_state0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
enrichr_res$Term = stringr::str_split_i(enrichr_res$Term, '\\(', 1)
enrichr_res$Term = str_to_title(enrichr_res$Term)
dt = subset(enrichr_res, select = c('Term', 'score', 'cell_state'))

dt = dcast(dt, Term ~ cell_state, value.var = 'score')

mat = dt[, -1]
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap(mat, na_col = 'white', border_color = NA, main = 'enrichR:GO_BP',
              cluster_rows = T, cluster_cols = T, angle_col = 45)
ggsave(aa, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                             '/enrichR_GO_Biological_Process_2023/heatmap_enrichr_hallmark_allCellStates_',
                             deg_type, '.pdf'),
       device = 'pdf', width = 12, height = 10)
### dotplot

enrichr_res$Term = factor(enrichr_res$Term, 
                          levels = rev(rownames(mat)[aa$tree_row$order]))
enrichr_res$cell_state = factor(enrichr_res$cell_state, 
                                levels = colnames(mat)[aa$tree_col$order])

enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$size = as.integer(stringr::str_split_i(enrichr_res$Overlap, '/', 2))
enrichr_res$frac = enrichr_res$Count/enrichr_res$size
pg <- ggplot(data = enrichr_res, aes(x = cell_state, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(oPTXon = "D") + theme_bw() +
  ggtitle('enrichr: GO_BP') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(padj)"), 
         size=guide_legend(title="Overlap")) 

ggsave(pg, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                             '/enrichR_GO_Biological_Process_2023/dotplot_enrichr_go_allCellStates_',
                             deg_type, '.pdf'),
       device = 'pdf', width = 7, height = 8)

#### plot dot plot -- reactome ####
enrichr_res_dir = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, '/enrichR_Reactome_2022/')
### collect all top terms
sele_terms = NULL
for(cell_state0 in unique(degs$cluster)){
  tmp = fread(paste0(enrichr_res_dir, 'Reactome_2022_', cell_state0, '_', deg_type, '.tsv'))
  tmp$cell_state = cell_state0
  tmp = tmp[Adjusted.P.value < 0.05]
  tmp = tmp[1:min(nrow(tmp), 8), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(cell_state0 in unique(degs$cluster)){
  tmp = fread(paste0(enrichr_res_dir, 'Reactome_2022_', cell_state0, '_', deg_type, '.tsv'))
  tmp$cell_state = cell_state0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  tmp = tmp[tmp$Adjusted.P.value < 0.05, ]
  enrichr_res[[cell_state0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'cell_state'))

dt = dcast(dt, Term ~ cell_state, value.var = 'score')

mat = dt[, -1]
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap(mat, na_col = 'white', border_color = NA, main = 'enrichR: Reactome',
              cluster_rows = T, cluster_cols = T, angle_col = 45)
ggsave(aa, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                             '/enrichR_Reactome_2022/heatmap_enrichr_reactome_allCellStates_',
                             deg_type, '.pdf'),
       device = 'pdf', width = 12, height = 10)
### dotplot

enrichr_res$Term = factor(enrichr_res$Term, levels = rownames(mat)[rev(aa$tree_row$order)])
enrichr_res$cell_state = factor(enrichr_res$cell_state, 
                                levels = colnames(mat)[aa$tree_col$order])

enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$size = as.integer(stringr::str_split_i(enrichr_res$Overlap, '/', 2))
enrichr_res$frac = enrichr_res$Count/enrichr_res$size
pr <- ggplot(data = enrichr_res, aes(x = cell_state, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(oPTXon = "D") + theme_bw() +
  ggtitle('enrichr: REACTOME') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '-log10(pvalue)') + 
  guides(color=guide_legend(title="-log10(padj)"), 
         size=guide_legend(title="Overlap")) 

ggsave(pr, filename = paste0('Figures/RNA/Malignant/btw_clusters/', deg_type, 
                             '/enrichR_Reactome_2022/dotplot_enrichr_reactome_allCellStates_',
                             deg_type, '.pdf'),
       device = 'pdf', width = 10, height = 10)

#### save degs into excel file ####
## sc 
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
degs[, 'pct.diff' := pct.1 - pct.2]
degs = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]

degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'
split_degs = split(degs, by = 'cell_state')
writexl::write_xlsx(split_degs, path = 'data/intermediate/RNA/new_degs/degs_between_malignant_cell_states.xlsx')

## bulk edger
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_clusters_bulk_edger.rds')
degs = degs[logFC > 0.5 & PValue < 0.05]

degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'
split_degs = split(degs, by = 'cell_state')
writexl::write_xlsx(split_degs, path = 'data/intermediate/RNA/new_degs/degs_between_malignant_cell_states_bulk_edger.xlsx')

## bulk limma-vom
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_clusters_bulk_limma_voom.rds')
degs = degs[logFC > 0.5 & P.Value < 0.05]

degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'
split_degs = split(degs, by = 'cell_state')
writexl::write_xlsx(split_degs, path = 'data/intermediate/RNA/new_degs/degs_between_malignant_cell_states_bulk_limma_voom.xlsx')

## Plot top degs ####
degs_sc = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
degs_sc = degs_sc[avg_log2FC > 0.5 & pct.1 > pct.2]

sele_sc = NULL
for(cl0 in 0:5){
  degs0 = degs_sc[cluster == cl0]
  degs1 = degs_limma[cluster == cl0]
  sele_sc = unique(c(sele_sc, degs0$gene[1:15]))
}
seurat.rna = ScaleData(seurat.rna, features = csele_sc)
seurat.down = subset(seurat.rna, down = 100)
pd <- DoHeatmap(seurat.down, features = sele_sc, disp.min = -2, disp.max = 2) +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])


## DX vs PTX within malignant subcluster ####
seurat.rna = readRDS(file = 'Seurat_Objects/seurat_rna_malignant_harmony_clean.rds')
DefaultAssay(seurat.rna) = 'RNA'
mdata = readRDS('MetaData/metadata_seurat_rna_malignant_harmony_cellStateAnnotated.rds')
seurat.rna <- AddMetaData(seurat.rna, metadata = mdata)

degs = NULL
for(cl0 in 0:5){
  seurat0 = subset(seurat.rna, seurat_clusters == cl0 & Stage_Code %in% c('DX', 'PTX'))
  Idents(seurat0) = seurat0$Stage_Code
  degs0 = FindMarkers(seurat0,  max.cells.per.ident = 500, ident.1 = 'PTX',
                      test.use = 'LR', logfc.threshold = 0.25,
                      min.pct = 0.1, latent.vars = c('nCount_RNA', 'percent.mt'), only.pos = F)
  degs0 = data.table(degs0, keep.rownames = T)
  names(degs0)[1] = 'gene'
  degs0$seurat_clusters = cl0
  degs = rbind(degs, degs0)
}
degs$cell_state[degs$seurat_clusters == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$seurat_clusters == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$seurat_clusters == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$seurat_clusters == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$seurat_clusters == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$seurat_clusters == 5] = 'MES_5'
degs[, 'pct.diff' := pct.1 - pct.2]
degs = degs[abs(avg_log2FC) > 0.25 & p_val < 0.001 & abs(pct.diff) > 0.05]


split_degs = split(degs, by = 'cell_state')
writexl::write_xlsx(split_degs, 
                    path = 'data/intermediate/RNA/new_degs/degs_between_timepoints_within_cell_states.xlsx')


saveRDS(degs, file = 'data/intermediate/RNA/new_degs/degs_PTXvsDX_withinNewClusters.rds')


## Stage specific DEGs with Stage-timepoint specific DEGs ####
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
degs_time = readRDS('data/intermediate/RNA/new_degs/degs_PTXvsDX_withinNewClusters.rds')
degs_time = degs_time[avg_log2FC > 0.25 & p_val < 0.001 & pct.diff > 0]

degs = degs[avg_log2FC > 0.25 & p_val < 0.001 & pct.1 > pct.2]
degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'
degs[, 'pct.diff' := pct.1 - pct.2]

shared_degs = NULL
for(state0 in unique(degs$cell_state)){
  degs0 = degs[cell_state == state0]
  degs_time0 = degs_time[cell_state == state0]
  shared_genes0 = intersect(degs0$gene, degs_time0$gene)
  degs0 = subset(degs0, gene %in% shared_genes0, 
                 select = c(gene, p_val, avg_log2FC, pct.diff, p_val_adj, cell_state))
  degs_time0 = subset(degs_time0, gene %in% shared_genes0, 
                 select = c(gene, p_val, avg_log2FC, pct.diff, p_val_adj, cell_state))
  names(degs_time0)[2:5] = c('p_val_time', 'avg_log2FC_time', 'pct.diff_time', 'p_val_adj_time')
  tmp = merge(degs0, degs_time0, by = c('gene', 'cell_state'))
  shared_degs = rbind(shared_degs, tmp)
}

saveRDS(shared_degs, file = 'data/intermediate/RNA/new_degs/degs_upInstate_upInPTXState.rds')




