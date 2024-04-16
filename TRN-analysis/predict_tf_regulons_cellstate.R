## predict tf regulon per tumor state

library(data.table)
library(Seurat)
library(magrittr)
library(ggplot2)
library(matrixStats)
library(viridis)
library(openxlsx)
library(stringr)
library(enrichR)
`%notin%` = Negate(`%in%`)

## prepare data ####
## diff tfs
diff.tf = readRDS(file = 'data/intermediate/ATAC/seurat_diff_chromVAR_dscore_by_atac_cell_states_allTFsKept.rds')

## diff degs
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
degs[, 'pct.diff' := pct.1 - pct.2]
degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'

degs_down = degs[avg_log2FC < -0.25 & p_val < 0.001] 
degs_up = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]


## degs btw timepoints
degs_time = readRDS('data/intermediate/RNA/new_degs/degs_PTXvsDX_withinNewClusters.rds')


## daps
#daps = readRDS('data/intermediate/new_atac_daps_byClusters_bulk_edgeR_resl0.4.rds')
daps = readRDS('data/intermediate/new_atac_daps_byStates_bulk_edgeR.rds')
fc_peaks = readRDS('data/intermediate/new_atac_peak_fc_byStates.rds')
# weak peaks
if(T){
  seurat.atac <- readRDS("Seurat_Objects/signac_obj_adrenergic_signac_modified_integration_v3.rds")
  labelTransRes = readRDS('MetaData/labelTransferResult_cca_signac_obj_adrenergic_signac_modified_integration_v2.rds')
  seurat.atac = AddMetaData(seurat.atac, metadata = labelTransRes)
  seurat.atac$ATAC_snn_res.0.4 = as.character(seurat.atac$ATAC_snn_res.0.4)
  seurat.tmp = subset(seurat.atac, prediction.score.max > 0.6)
  
  mtx = seurat.tmp@assays$ATAC@counts
  open.frac.pks = sapply(unique(seurat.tmp$predicted.id), function(x) {
    return(rowMeans(mtx[, seurat.tmp$predicted.id==x] > 0))
  })
  pm = rowMaxs(open.frac.pks)
  weak.pks = rownames(mtx)[pm < 0.025]
  saveRDS(weak.pks, file = 'MetaData/weak_peaks.rds')
  saveRDS(open.frac.pks, file = 'MetaData/peak_openFrac_byState.rds')
}

weak.pks = readRDS('MetaData/weak_peaks.rds')
open.frac.pks = readRDS('MetaData/peak_openFrac_byState.rds')

## filter/add TF by gene expression
if(T){
  ## just run once
  seurat.rna = readRDS('Seurat_Objects/seurat_rna_malignant_harmony_cellStateAnnotated.rds')
  mtx = seurat.rna[['RNA']]$data
  expr.frac.state = sapply(unique(seurat.rna$cell_state), function(x){
    return(rowMeans(mtx[, seurat.rna$cell_state == x] > 0))
  })
  colnames(expr.frac.state) = unique(seurat.rna$cell_state)
  saveRDS(expr.frac.state, file = 'MetaData/expr_frac_state.rds')
  
  expr.avg.state = sapply(unique(seurat.rna$cell_state), function(x){
    return(rowMeans(mtx[, seurat.rna$cell_state == x]))
  })
  colnames(expr.avg.state) = unique(seurat.rna$cell_state)
  saveRDS(expr.avg.state, file = 'MetaData/avg_expr_state.rds')
  
  expr.median.state = sapply(unique(seurat.rna$cell_state), function(x){
    return(sparseMatrixStats::rowMedians(mtx[, seurat.rna$cell_state == x]))
  })
  colnames(expr.median.state) = unique(seurat.rna$cell_state)
  rownames(expr.median.state) = rownames(mtx)
  saveRDS(expr.median.state, file = 'MetaData/median_expr_state.rds')
  
  ## diff expressed TFs
  tfs = readRDS('MetaData/snATAC/TF/motif_id_motif_name.rds')
  names(tfs) = NULL
  tfs.deg = degs_up[gene %in% tfs]
  saveRDS(tfs.deg, file = 'MetaData/tfs_deg_state.rds')
}
expr.frac.state = readRDS('MetaData/expr_frac_state.rds')
expr.median.state = readRDS('MetaData/median_expr_state.rds')
expr.avg.state = readRDS('MetaData/avg_expr_state.rds')
tfs.deg = readRDS('MetaData/tfs_deg_state.rds')

##motif scan
motif_pk_match = readRDS('data/intermediate/TRN/peak_motif_match_pv5e-05.rds')
motif_pk_match = as.matrix(motif_pk_match)

## controlling parameters
filtered_links = NULL
filterWeakPeak = TRUE
useDAP = TRUE ## filter peaks by daps/foldChange or by peak open fraction(>0.05)
useAllDEGs = TRUE ## use all DEGs or filter
useEnhancerPeakOnly = TRUE ##use promoter enhancer interaction only or include PPIs
peakFracCutoff = 0.1 ## used to filter peaks per state (default 0.1)
tfexprfracCutoff = 0.2

predicted_links = fread('data/intermediate/TRN/EP_Prediction/new_regrRes4_gene_peak_links_metacell.tsv')

if(useEnhancerPeakOnly){
  predicted_links = fread('data/intermediate/TRN/EP_Prediction/new_regrRes4_gene_peak_ep_links_metacell.tsv')
}
if(filterWeakPeak) predicted_links = predicted_links[peak_name %notin% weak.pks]


## construct TRN per state ####
## provide or read your interested gene and tf lists
edges_list = nodes_list = list()
for(cell_state0 in colnames(expr.frac.state)){
  
  degs0 = degs_up[cell_state == cell_state0]
  degs_down0 = degs_down[cell_state == cell_state0]
  open.frac.pks0 = open.frac.pks[, cell_state0]
  
  interested_genes = unique(degs0$gene)
  
  ## as TF strength
  avg.zscore.state = subset(diff.tf, select = c('gene', 'avg_diff', 'cluster'))
  avg.zscore.state0 = avg.zscore.state[cluster == cell_state0]
  setkey(avg.zscore.state0, gene)
  
  interested_tfs = union(diff.tf[cluster == cell_state0 & avg_diff > 0.5 & p_val_adj < 0.05]$gene, 
                              tfs.deg[cell_state == cell_state0]$gene)
  interested_tfs = interested_tfs[interested_tfs %notin% degs_down0$gene]  
  
  
  ## filtering peaks
  filtered.peaks = daps[cluster == cell_state0]$peak
  if(!useDAP) filtered.peaks = names(which(open.frac.pks0 > peakFracCutoff))
  
  ## further filtering TF by expression
  expr.frac.state0 <- expr.frac.state[, cell_state0]
  expr.median.state0 <- expr.median.state[, cell_state0]
  expr.avg.state0 <- expr.avg.state[, cell_state0]
  
  tf.expr.frac = expr.frac.state0[interested_tfs]
  tf.expr.frac = tf.expr.frac[!is.na(tf.expr.frac)]
  interested_tfs = names(which(tf.expr.frac > tfexprfracCutoff))
  
  if(length(interested_tfs) <= 5) {
    interested_tfs = intersect(names(sort(tf.expr.frac, decreasing = T))[1:10],
                               names(which(tf.expr.frac > 0.05)))
  }
  predicted_links0 = predicted_links[gene_name %in% interested_genes]
  sele.tf.mat = motif_pk_match[, interested_tfs] 
  sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
  sele.peaks = intersect(sele.peaks, filtered.peaks)
  predicted_links0 = predicted_links0[peak_name %in% sele.peaks]
  
 
  ##save filtered links
  predicted_links0 = predicted_links0[Estimate > 0.2 & fdr < 0.001]
  predicted_links0$cell_state = cell_state0
  filtered_links = rbind(filtered_links, predicted_links0)
  
  ## split links by tf
  tf_regulons = list()
  open.frac.pks0 = open.frac.pks[, cell_state0]
  for(TF0 in interested_tfs){
    peaks0 = names(which(sele.tf.mat[, TF0] > 0))
    ep0 = predicted_links0[peak_name %in% peaks0]
    ep0 = subset(ep0, select = c('gene_name', 'Estimate',  'fdr', 'peak_name'))
    ep0$peak_strength = open.frac.pks0[ep0$peak_name]
    #ep0 = ep0[peak_strength > 0.05]
    ep0[, 'score' := -log10(fdr)]
    ep0$TF = TF0
    ep0[, 'edge_strength' := sum(Estimate * peak_strength), by = gene_name]
    ep0[, 'Estimate' := sum(Estimate), by = gene_name]
    ep0[, 'score' := mean(score), by = gene_name]
    ep0[, 'peaks' := paste0(peak_name, collapse = ','), by = gene_name]
    ep0[, 'peak_strengths' := paste0(peak_strength, collapse = ','), by = gene_name]
    ep0[, 'npeaks' := .N, by = gene_name]
    ep0 = subset(ep0, select = c('TF', 'gene_name', 'edge_strength', 'peak_strengths',
                                 'Estimate', 'score', 'peaks', 'npeaks'))
    ## node inf
    ep0$TF_strength = avg.zscore.state0[ep0$TF]$avg_diff
    ep0$gene_strength = expr.avg.state0[ep0$gene_name]
    
    tf_regulons[[TF0]] = ep0[!duplicated(ep0)]
  }
  
  
  tf_regulons.comb = do.call('rbind', tf_regulons)
  
  tf_regulons.comb = tf_regulons.comb[order(-edge_strength)]
  message(paste0(cell_state0, ': ', nrow(tf_regulons.comb), ' regulons'))
  
  nodes1 = subset(tf_regulons.comb, select = c('TF', 'TF_strength'))
  nodes2 = subset(tf_regulons.comb, select = c('gene_name', 'gene_strength'))
  names(nodes1) = names(nodes2) = c('gene', 'node_strength')
  nodes1$node_strength = nodes1$node_strength/max(nodes1$node_strength)
  nodes2$node_strength = nodes2$node_strength/max(nodes2$node_strength)
  nodes1$node_type = 'TF'
  nodes2$node_type = 'gene'
  nodes2 = nodes2[gene %notin% nodes1$gene]
  nodes = rbind(nodes1, nodes2) %>% unique()
  
  #add deg information btw timepoint
  degs_time0 = degs_time[cell_state == cell_state0] ## include up and down information
  setkey(degs_time0, gene)
  nodes[, 'avg_log2FC_time' := degs_time0[nodes$gene]$avg_log2FC]
  nodes[is.na(avg_log2FC_time)]$avg_log2FC_time = 0
  
  if(useDAP) {
    tf_regulon_fname = paste0('data/intermediate/TRN/networks/edges_', cell_state0, '_useDAPs.tsv')
    nodes_fname = paste0('data/intermediate/TRN/networks/nodes_', cell_state0, '_useDAPs.tsv')
    tf_rank_fname = paste0('data/intermediate/TRN/networks/ConnectRank_TF_regulons_', cell_state0, '_useDAPs.pdf')
  }else{
    tf_regulon_fname = paste0('data/intermediate/TRN/networks/edges_', cell_state0, '_peakFracCutoff', peakFracCutoff,  '.tsv')
    nodes_fname = paste0('data/intermediate/TRN/networks/nodes_', cell_state0, '_peakFracCutoff', peakFracCutoff, '.tsv')
    tf_rank_fname = paste0('data/intermediate/TRN/networks/ConnectRank_TF_regulons_', cell_state0, '_peakFracCutoff', peakFracCutoff, '.pdf')
    
  }
  
  if(!useEnhancerPeakOnly){
    tf_regulon_fname = gsub('.tsv', '_IncludePPI.tsv', tf_regulon_fname)
    nodes_fname = gsub('.tsv', '_IncludePPI.tsv', nodes_fname)
    tf_rank_fname = gsub('.pdf', '_IncludePPI.pdf', tf_rank_fname)
  }
  if(useAllDEGs)  {
    tf_regulon_fname = gsub('.tsv', '_useAllDEGs.tsv', tf_regulon_fname)
    nodes_fname = gsub('.tsv', '_useAllDEGs.tsv', nodes_fname)
    tf_rank_fname = gsub('.pdf', '_useAllDEGs.pdf', tf_rank_fname)
  }
  tf_regulons.comb$cell_state = cell_state0
  nodes$cell_state = cell_state0
  write.table(tf_regulons.comb, file = tf_regulon_fname,
              sep = '\t', row.names = F, quote = F)
  write.table(nodes, file = nodes_fname,
              sep = '\t', row.names = F, quote = F)
  edges_list[[cell_state0]] = tf_regulons.comb
  nodes_list[[cell_state0]] = nodes
  #plot tf ranking
  pdata = data.table('TF' = names(sort(table(tf_regulons.comb$TF))),
                     'nTargets' = as.numeric(sort(table(tf_regulons.comb$TF))))
  pdata$TF = factor(pdata$TF, levels = names(sort(table(tf_regulons.comb$TF))))
  p1 <- ggplot(data = pdata, aes(x = TF, y = nTargets, fill = nTargets)) + geom_bar(stat = 'identity') + coord_flip() +
    theme_classic() + scale_fill_viridis_c(oPTXon = "magma") + ggtitle(cell_state0) + xlab('')
  
  ggsave(p1, filename = tf_rank_fname, device = 'pdf', width = 5, height = 2 * ceiling(nrow(pdata)/20))
}

if(useDAP) {
  filtered_ep_file = 'data/intermediate/TRN/EP_Prediction/filtered_regrRes4_gene_peak_links_metacell_useDAPs_ep.tsv'
}else{
  filtered_ep_file = paste0('data/intermediate/TRN/EP_Prediction/filtered_regrRes4_gene_peak_links_metacell_peakFracCutoff',
                            peakFracCutoff, '_ep.tsv')
}
if(!useEnhancerPeakOnly) {
  filtered_ep_file = gsub('ep.tsv', 'IncludePPI.tsv', filtered_ep_file, fixed = T)
}
if(useAllDEGs)  filtered_ep_file = gsub('.tsv', '_useAllDEGs.tsv', filtered_ep_file)

fwrite(filtered_links, file = filtered_ep_file, sep = '\t')

edge_file = 'data/intermediate/TRN/networks/TRN_edges.xlsx'
node_file = 'data/intermediate/TRN/networks/TRN_nodes.xlsx'
if(useDAP){
  edge_file = 'data/intermediate/TRN/networks/TRN_edges_useDAPs.xlsx'
  node_file = 'data/intermediate/TRN/networks/TRN_nodes_useDAPs.xlsx'
}
writexl::write_xlsx(edges_list,
                    path = edge_file)

writexl::write_xlsx(nodes_list,
                    path = node_file)

edges_combined = do.call('rbind', edges_list)
nodes_combined = do.call('rbind', nodes_list)

writexl::write_xlsx(edges_combined,
                    path = 'data/intermediate/TRN/networks/TRN_edges_combined.xlsx')

writexl::write_xlsx(nodes_combined,
                    path = 'data/intermediate/TRN/networks/TRN_nodes_combined.xlsx')
writexl::write_xlsx(daps,
                    path = 'data/intermediate/TRN/networks/DAPs_Between_Neoplastic_Cell_States.xlsx')


## dot plot of all top tfs based on #targets ####
pdata.comb = NULL
sele.tfs = NULL
for(cell_state1 in c('ADRN_Calcium_0', 'MES_5', 'ADRN_Dopaminergic_3',
                      'ADRN_Proliferating_4',  'ADRN_Ribosome_1',
                      'Interm_OxPhos_2')){
  
  tf_regulon_file = paste0('data/intermediate/TRN/networks/edges_',
                           cell_state1, '_useDAPs_useAllDEGs', '.tsv')
  if(cell_state1 == 'Interm_OxPhos_2') tf_regulon_file = paste0('data/intermediate/TRN/networks/edges_',
                                                                  cell_state1, '_peakFracCutoff0.1_useAllDEGs', '.tsv')
  if(!file.exists(tf_regulon_file)) next
  tf_regulons = fread(tf_regulon_file)
  
  pdata = subset(tf_regulons, select = c(TF, gene_name, TF_strength)) %>% unique()
  names(pdata)[3] = 'zscore_avg_diff'
  pdata[, 'nTarget' := .N, by = TF]
  pdata$N = length(unique(pdata$gene_name))
  pdata[, 'frac_target' := nTarget/N]
  pdata = subset(pdata, select = c(TF, frac_target, zscore_avg_diff, N)) %>% unique()
  pdata = pdata[zscore_avg_diff > 0]
  pdata = pdata[order(-frac_target)]
  pdata$state = cell_state1
  pdata.comb = rbind(pdata.comb, pdata)
  
  #pdata = pdata[1:min(15, nrow(pdata))]
  pdata = pdata[frac_target > 0.2]
  sele.tfs = union(sele.tfs, pdata$TF)
  
}
pdata.comb = pdata.comb[TF %in% sele.tfs]
pdata.comb[zscore_avg_diff > 4]$zscore_avg_diff = 4

pdata.comb = pdata.comb[order(state)]

pdata.comb$state = factor(pdata.comb$state, levels = c('ADRN_Calcium_0', 'MES_5', 'ADRN_Dopaminergic_3',
                                                       'ADRN_Proliferating_4',  'ADRN_Ribosome_1',
                                                       'Interm_OxPhos_2'))
pdata.comb$TF = factor(pdata.comb$TF, levels = rev(unique(pdata.comb$TF)))

ph <- ggplot(data = pdata.comb, aes(x = state, y = TF, color = zscore_avg_diff)) +
  geom_point(aes(size = frac_target)) +
  scale_color_viridis_c(oPTXon = "C") + theme_bw() +
  ggtitle('') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') +
  guides(color=guide_legend(title="diff.chromvar.zscore"),
         size=guide_legend(title="targets.coverage"))

ggsave(ph, filename = 'data/intermediate/TRN/networks/dotplot_TFs_allStates.pdf',
       device = 'pdf', width = 4, height = 10)



## sub-network for presentation ####
## load data
diff.tf = readRDS(file = 'data/intermediate/ATAC/seurat_diff_chromVAR_dscore_by_atac_cell_states.rds')
diff.tf = diff.tf[avg_diff > 0.5 & p_val_adj < 0.05]

## diff degs
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
degs[, 'pct.diff' := pct.1 - pct.2]
degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'

degs_up = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]

tfs.deg = readRDS('MetaData/tfs_deg_state.rds')

useDAP = TRUE
show_tfs = NULL
for(cell_state0 in unique(degs$cell_state)){
  message(cell_state0, ':')
  edge_file = paste0('data/intermediate/TRN/TF_regulons_', cell_state0, '_useDAPs_useAllDEGs.tsv')
  node_file = paste0('data/intermediate/TRN/nodes_', cell_state0, '_useDAPs_useAllDEGs.tsv')
  if(!useDAP | cell_state0 == 'Interm_OxPhos_2'){
    ## because Interm_OxPhos_2 with DAP have very limited regulons
    edge_file = paste0('data/intermediate/TRN/TF_regulons_', cell_state0, '_peakFracCutoff0.1_useAllDEGs.tsv')
    node_file = paste0('data/intermediate/TRN/nodes_', cell_state0, '_peakFracCutoff0.1_useAllDEGs.tsv')
    
  }
  
  edges = fread(edge_file)
  nodes = fread(node_file)
  #tf by degree
  edges[, 'tf_degree' := .N, by = TF]
  tf_rk = subset(edges, select = c('TF', 'tf_degree')) %>% unique()
  tf_rk = tf_rk[order(-tf_degree)]
  tf_rk = tf_rk[1:min(10, nrow(tf_rk))]
  tf_top_degree  = tf_rk$TF
  
  
  sele_genes = degs_up[cell_state == cell_state0]$gene
  topg = sele_genes[1:min(100, length(sele_genes))]
  if(cell_state0 == 'MES_5') topg = sele_genes[1:min(200, length(sele_genes))]
  if(cell_state0 == 'ADRN_Calcium_0') topg = sele_genes[1:min(50, length(sele_genes))]
  edges = edges[gene_name %in% topg]
  

  diff.tf0 = diff.tf[cluster == cell_state0]
  tfs.deg0 =  tfs.deg[cell_state == cell_state0]
  
  
  #sele_tfs = tf_top_degree
  sig_tfs = intersect(diff.tf0$gene[1:min(50, nrow(diff.tf0))], 
                     tfs.deg0$gene[1:min(50, nrow(tfs.deg0))])
  #sele_tfs = sig_tfs
  message(paste0('sig_tfs: ', paste0(sig_tfs, collapse = ',')))
  sele_tfs = union(sig_tfs, tf_top_degree)
  
  
  tf_nodes = nodes[node_type == 'TF' & node_strength > 0 ]
  tf_nodes = tf_nodes[gene %in% sele_tfs]
  gene_nodes =  nodes[node_type != 'TF' & gene %in% topg]
  
  sub_nodes = rbind(tf_nodes, gene_nodes)
  sub_edges = edges[TF %in% tf_nodes$gene & gene_name %in% gene_nodes$gene]
  message(paste0(nrow(sub_edges), ' edges, ',
                 nrow(tf_nodes), ' TFs, '), nrow(gene_nodes), ' genes')
  message(paste0('TFs: ', paste0(tf_nodes$gene, collapse = ',')))
  
  sub_edge_file = paste0('data/intermediate/TRN/subnetworks/useDAP/sub_edges_', cell_state0, '_useDAPs_useAllDEGs.tsv')
  
  sub_node_file = paste0('data/intermediate/TRN/subnetworks/useDAP/sub_nodes_', cell_state0, '_useDAPs_useAllDEGs.tsv')
  if(!useDAP){
    sub_edge_file = paste0('data/intermediate/TRN/subnetworks/useFraction/sub_edges_', cell_state0, '_peakFracCutoff0.1_useAllDEGs.tsv')
    sub_node_file = paste0('data/intermediate/TRN/subnetworks/useFraction/sub_nodes_', cell_state0, '_peakFracCutoff0.1_useAllDEGs.tsv')
    
  }
  
  write.table(sub_edges, file = sub_edge_file, 
              sep = '\t', row.names = F, quote = F)
  write.table(sub_nodes, file = sub_node_file, 
              sep = '\t', row.names = F, quote = F)
  show_tfs = union(show_tfs, tf_nodes$gene)
}


