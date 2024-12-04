library(data.table)
library(magrittr)
library(stringr)
library(ggplot2)
`%notin%` = Negate(`%in%`)


# prepare input for liana dotplot ####
## rename adrenergic clusters
cl_state_map = data.table('cluster' = c(paste0('Malignant', 0:5), paste0('Macro', c(0:5,7, 8))),
                          'cell_state' = c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                                           'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5', 
                                           'THY1+Macro_0', 'HS3ST2+Macro_1',
                                           'F13A1+Macro_2', 'CCL4+Macro_3', 
                                           'IL18+Macro_4','VCAN+Macro_5',
                                           'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8'))

setkey(cl_state_map, cluster)

expr_cutoff = 0.1
cytotalk_dir = paste0('data/intermediate/cytotalk/secRound/exprCutoff', expr_cutoff, '/')


cell_pairs <- list.files(cytotalk_dir)
#cell_pairs <- cell_pairs[cell_pairs %notin% c('inputs', 'result_Adrenergic4_Macrophage0', 
#                                              'result_Adrenergic7_Macrophage0')]

df.cytotalk <- NULL
for(pair in cell_pairs){
  source_cl = str_split_i(pair, '_', 2 )
  target_cl = str_split_i(pair, '_', 3 )
  cfile = paste0(cytotalk_dir, pair, '/resultObj_', source_cl, '_', target_cl, '.rds')
  if(!file.exists(cfile)) next
  cytotalk_res = readRDS(cfile)
  df_LR = cytotalk_res$pathways$df_pval
  df_LR$ligand = toupper(df_LR$ligand)
  df_LR$receptor = toupper(df_LR$receptor)
  df_LR$ligand_type = cl_state_map[df_LR$ligand_type]$cell_state
  df_LR$receptor_type = cl_state_map[df_LR$receptor_type]$cell_state
  df_LR$source = df_LR$ligand_type
  df_LR$target = df_LR$receptor_type
  df_LR$int_specificity = -log10(df_LR$pval_potential + 1e-10)
  df_LR$L_R = paste0(df_LR$ligand, '_', df_LR$receptor)
  df_LR$L_R_source <- paste0(df_LR$ligand,'_',df_LR$ligand_type,'_',
                             df_LR$receptor,'_',df_LR$receptor_type)  
  
  # read cost from final network file
  pathway_file = paste0(cytotalk_dir, pair,  '/FinalNetwork.txt')
  if(!file.exists(pathway_file)) next
  
  df_pathway = fread(pathway_file, header = T)
  df_pathway <- df_pathway[is_ct_edge == TRUE, ]
  df_pathway$node1 <- toupper(df_pathway$node1)
  df_pathway$node2 <- toupper(df_pathway$node2)
  df_pathway$node1_type = cl_state_map[df_pathway$node1_type]$cell_state
  df_pathway$node2_type = cl_state_map[df_pathway$node2_type]$cell_state
  
  df_pathway$L_R <- paste0(df_pathway$node1,'_',df_pathway$node2)  #node1 - ligand; node2 - receptor
  df_pathway$L_R_source <- paste0(df_pathway$node1,'_',df_pathway$node1_type,'_',
                                  df_pathway$node2,'_',df_pathway$node2_type)  #node1 - ligand; node2 - receptor
  setkey(df_pathway, L_R_source)
  df_LR$cost = df_pathway[df_LR$L_R_source]$cost
  df_LR$cross_talk_score <- 1 - df_LR$cost
  
  
  
  #where is cross-talk score and 
  #what metric did they use to filter L-R/signaling pathway 
  #In general, cross-talk score/ node prize/edges / mutual connection between L_R and corresponding downstream nodes)
  df.cytotalk <- rbind(df.cytotalk, df_LR)
}

# plot manually ####
dt.cytotalk = data.table(df.cytotalk)
dt.cytotalk[, 'interaction' := paste(ligand, receptor, sep = '->')]

saveRDS(dt.cytotalk, file = paste0('data/intermediate/cytotalk/secRound/cci_summarized_res_exprCutoff', expr_cutoff, '.rds'))

dt.cytotalk0 = dt.cytotalk[source %in% c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                                         'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5')]
dt.cytotalk0 = dt.cytotalk0[target %in% c('THY1+Macro_0', 'HS3ST2+Macro_1',
                                          'F13A1+Macro_2', 'CCL4+Macro_3', 
                                          'IL18+Macro_4','VCAN+Macro_5',
                                          'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8')]
dt.cytotalk0$source = factor(dt.cytotalk0$source, levels = c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                                                             'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5'))
dt.cytotalk0$target = factor(dt.cytotalk0$target, levels = c('THY1+Macro_0', 'HS3ST2+Macro_1',
                                                             'F13A1+Macro_2', 'CCL4+Macro_3', 
                                                             'IL18+Macro_4','VCAN+Macro_5',
                                                             'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8'))

p0 <- ggplot(data = dt.cytotalk0, aes(x = source, y = interaction, col = cross_talk_score)) +
  geom_point(size = 4) + 
  scale_color_viridis(option = 'C')  + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(. ~ target,    space = "free",
                            scales ="free",
                            switch = "y")

dt.cytotalk1 = dt.cytotalk[target %in% c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                                         'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5')]
dt.cytotalk1 = dt.cytotalk1[source %in% c('THY1+Macro_0', 'HS3ST2+Macro_1',
                                          'F13A1+Macro_2', 'CCL4+Macro_3', 
                                          'IL18+Macro_4','VCAN+Macro_5',
                                          'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8')]
dt.cytotalk1$target = factor(dt.cytotalk1$target, levels = c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                                                             'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5'))
dt.cytotalk1$source = factor(dt.cytotalk1$source, levels = c('THY1+Macro_0', 'HS3ST2+Macro_1',
                                                             'F13A1+Macro_2', 'CCL4+Macro_3', 
                                                             'IL18+Macro_4','VCAN+Macro_5',
                                                             'C1QC+SPP1+Macro_7', 'ProliferatingMacro_8'))

p1 <- ggplot(data = dt.cytotalk1, aes(x = target, y = interaction, col = cross_talk_score)) +
  geom_point(size = 4) + 
  scale_color_viridis(option = 'C')  + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(. ~ source,    space = "free",
             scales ="free",
             switch = "y")

ggsave(p0, filename = paste0('Figures/RNA/dotplot_cci_malignant2macro_exprCutoff', expr_cutoff, '.pdf'), 
       device = 'pdf', width = 14, height = 8)
ggsave(p1, filename = paste0('Figures/RNA/dotplot_cci_macro2malignant_exprCutoff', expr_cutoff, '.pdf'),
       device = 'pdf',
       width = 14, height = 8)
