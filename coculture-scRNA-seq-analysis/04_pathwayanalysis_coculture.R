library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
library(harmony)
library(RColorBrewer)
library(enrichR)
library(viridis)

# CHLA15 ####

## enrichr up-regulated genes ####
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
sele.dbs <- c("GO_Biological_Process_2023", 'MSigDB_Hallmark_2020')


## should do it comparison by comparison ##
comparison0 = 'Coculture_vs_Monoculture'

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

degs = readRDS(paste0('data/intermediate/coculture/CHLA15_degs_', comparison0, 
                      '_cellStateAsCovariate.rds'))
if(comparison0 == 'Coculture_vs_Monoculture'){
  degs0 = degs[avg_log2FC > 0.25 & p_val_adj < 0.05]
}else{
  degs0 = degs[avg_log2FC > 0.1 & p_val < 0.05]
}


genes0 = unique(degs0$gene)

message(paste0(comparison0, ': ', length(genes0), ' genes'))
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
      message(paste0(comparison0, ': no sig terms!')) 
      next
    }
    tmp = tmp[order(P.value)]
    nmax = min(15, nrow(tmp))
    
    tmp = tmp[1:nmax, ]
    tmp$Term = factor(tmp$Term, levels = rev(tmp$Term))
     
    dir.create(paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/'),
               showWarnings = F, recursive = T)
    filekey = paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/CHLA15_', 
                     db0, '_', comparison0)
    fwrite(res.state0[Count >= 3 & P.value < 0.05], sep = '\t',
           file = paste0(filekey, '.tsv'))
  }
  
}



## enrichr down-regulated genes ####
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
sele.dbs <- c("GO_Biological_Process_2023", 'MSigDB_Hallmark_2020')


## should do it comparison by comparison ##
comparison0 = 'Coculture_vs_Monoculture'

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

degs = readRDS(paste0('data/intermediate/coculture/CHLA15_degs_', comparison0, 
                      '_cellStateAsCovariate.rds'))
if(comparison0 == 'Coculture_vs_Monoculture'){
  degs0 = degs[avg_log2FC < -0.25 & p_val_adj < 0.05]
}else{
  degs0 = degs[avg_log2FC < -0.1 & p_val < 0.05]
}


genes0 = unique(degs0$gene)

message(paste0(comparison0, ': ', length(genes0), ' genes'))
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
      message(paste0(comparison0, ': no sig terms!')) 
      next
    }
    tmp = tmp[order(P.value)]
    nmax = min(15, nrow(tmp))
    
    tmp = tmp[1:nmax, ]
    tmp$Term = factor(tmp$Term, levels = rev(tmp$Term))
    
    dir.create(paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/'),
               showWarnings = F, recursive = T)
    filekey = paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/CHLA15_down_', 
                     db0, '_', comparison0)
    
    fwrite(res.state0[Count >= 3 & P.value < 0.05], sep = '\t',
           file = paste0(filekey, '.tsv'))
  }
  
}



## plot hallmark up in coculture and down in treated, and vice versa ####
#### Hallmark
adj_thr = 0.05

enrichr_res_dir = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/')
comparisons = c('Coculture_vs_Monoculture', 'Afatinib_vs_Coculture', 'CRM_vs_Coculture')

### up and down
sele_terms = NULL

for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA15_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Afatinib_vs_Coculture', 'CRM_vs_Coculture')){
    efile = paste0(enrichr_res_dir, 'CHLA15_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  tmp$score = -log10(tmp$Adjusted.P.value)
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp = tmp[1:min(nrow(tmp), 10), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA15_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Afatinib_vs_Coculture', 'CRM_vs_Coculture')){
    efile = paste0(enrichr_res_dir, 'CHLA15_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  
  tmp$score = -log10(tmp$Adjusted.P.value)
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp$rk = order(-tmp$score)
  tmp[rk>10]$rk=10
  enrichr_res[[comparison0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'comparison'))

dt = dcast(dt, Term ~ comparison, value.var = 'score')

mat = dt[, -1]
mat = subset(mat, select = comparisons[comparisons %in% names(mat)]) ## order it
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap::pheatmap(mat, na_col = 'white', border_color = NA,
                        cluster_rows = T, cluster_cols = F, angle_col = 45)
#ggsave(aa, filename = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA15_heatmap_hallmark_UpDown.pdf'),
#       device = 'pdf', width = 6, height = 6)

### dotplot

enrichr_res$Term = factor(enrichr_res$Term, 
                          levels = rev(rownames(mat)[aa$tree_row$order]))


## rename for easy read
enrichr_res[comparison == 'Coculture_vs_Monoculture']$comparison = 'Up in Coculture'
enrichr_res[comparison == 'Afatinib_vs_Coculture']$comparison = 'Down in Afatinib'
enrichr_res[comparison == 'CRM_vs_Coculture']$comparison = 'Down in CRM'
enrichr_res$comparison = factor(enrichr_res$comparison, 
                                levels = c('Up in Coculture', 
                                           'Down in Afatinib', 
                                           'Down in CRM'))


enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$frac = enrichr_res$score
pu <- ggplot(data = enrichr_res, aes(x = comparison, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(option = "D") + theme_bw() +
  ggtitle(paste0('CHLA15-HALLMARK-p.adj<', adj_thr )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(p.adj)"), 
         size=guide_legend(title="-log10(p.adj)")) 

#ggsave(ph, filename = 'Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA15_dotplot_hallmark_UpDown.pdf',
#       device = 'pdf', width = 5, height = 6)

## down and up

#### Hallmark
enrichr_res_dir = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/')
comparisons = c('Coculture_vs_Monoculture', 'Afatinib_vs_Coculture', 'CRM_vs_Coculture')

### collect all top terms
sele_terms = NULL

for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA15_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Coculture_vs_Monoculture')){
    efile = paste0(enrichr_res_dir, 'CHLA15_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp = tmp[1:min(nrow(tmp), 10), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA15_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Coculture_vs_Monoculture')){
    efile = paste0(enrichr_res_dir, 'CHLA15_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp$rk = order(-tmp$score)
  tmp[rk>10]$rk=10
  enrichr_res[[comparison0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'comparison'))

dt = dcast(dt, Term ~ comparison, value.var = 'score')

mat = dt[, -1]
mat = subset(mat, select = comparisons[comparisons %in% names(mat)]) ## order it
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap::pheatmap(mat, na_col = 'white', border_color = NA,
                        cluster_rows = T, cluster_cols = F, angle_col = 45)
#ggsave(aa, filename = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA15_heatmap_hallmark_DownUp.pdf'),
#       device = 'pdf', width = 6, height = 6)

### dotplot

enrichr_res$Term = factor(enrichr_res$Term, levels = rev(rownames(mat)[aa$tree_row$order]))


## rename for easy read
enrichr_res[comparison == 'Coculture_vs_Monoculture']$comparison = 'Down in Coculture'
enrichr_res[comparison == 'Afatinib_vs_Coculture']$comparison = 'Up in Afatinib'
enrichr_res[comparison == 'CRM_vs_Coculture']$comparison = 'Up in CRM'
enrichr_res$comparison = factor(enrichr_res$comparison, 
                                levels = c('Down in Coculture', 
                                           'Up in Afatinib', 
                                           'Up in CRM'))


enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$frac = enrichr_res$score
pd <- ggplot(data = enrichr_res, aes(x = comparison, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(option = "D") + theme_bw() +
  ggtitle(paste0('CHLA15-HALLMARK-p.adj<', adj_thr )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(p.adj)"), 
         size=guide_legend(title="-log10(p.adj)")) 

#ggsave(ph, filename = 'Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA15_dotplot_hallmark_DownUp.pdf',
#       device = 'pdf', width = 5, height = 6)

ggsave(pu + pd, filename = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA15_dotplot_hallmark_padj',
                                  adj_thr, '.pdf'),
       device = 'pdf', width = 10, height = 6)








# CHLA20 ####

## enrichr up-regulated genes ####
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
sele.dbs <- c("GO_Biological_Process_2023", 'MSigDB_Hallmark_2020')


## should do it by comparison ##
comparison0 = 'Coculture_vs_Monoculture'

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

degs = readRDS(paste0('data/intermediate/coculture/CHLA20_degs_', comparison0, 
                      '_cellStateAsCovariate.rds'))
if(comparison0 == 'Coculture_vs_Monoculture'){
  degs0 = degs[avg_log2FC > 0.25 & p_val_adj < 0.05]
}else{
  degs0 = degs[avg_log2FC > 0.1 & p_val < 0.05]
}


genes0 = unique(degs0$gene)

message(paste0(comparison0, ': ', length(genes0), ' genes'))
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
      message(paste0(comparison0, ': no sig terms!')) 
      next
    }
    tmp = tmp[order(P.value)]
    nmax = min(15, nrow(tmp))
    
    tmp = tmp[1:nmax, ]
    tmp$Term = factor(tmp$Term, levels = rev(tmp$Term))
    
    dir.create(paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/'),
               showWarnings = F, recursive = T)
    filekey = paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/CHLA20_', 
                     db0, '_', comparison0)
    fwrite(res.state0[Count >= 3 & P.value < 0.05], sep = '\t',
           file = paste0(filekey, '.tsv'))
  }
  
}



## enrichr down-regulated genes ####
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
sele.dbs <- c("GO_Biological_Process_2023", 'MSigDB_Hallmark_2020')


## should do it by comparison ##
comparison0 = 'Coculture_vs_Monoculture'

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

degs = readRDS(paste0('data/intermediate/coculture/CHLA20_degs_', comparison0, 
                      '_cellStateAsCovariate.rds'))
if(comparison0 == 'Coculture_vs_Monoculture'){
  degs0 = degs[avg_log2FC < -0.25 & p_val_adj < 0.05]
}else{
  degs0 = degs[avg_log2FC < -0.1 & p_val < 0.05]
}


genes0 = unique(degs0$gene)

message(paste0(comparison0, ': ', length(genes0), ' genes'))
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
      message(paste0(comparison0, ': no sig terms!')) 
      next
    }
    tmp = tmp[order(P.value)]
    nmax = min(15, nrow(tmp))
    
    tmp = tmp[1:nmax, ]
    tmp$Term = factor(tmp$Term, levels = rev(tmp$Term))
   
    dir.create(paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/'),
               showWarnings = F, recursive = T)
    filekey = paste0('Figures/RNA/coculture/enrichR_', db0, '/adjustState/CHLA20_down_', 
                     db0, '_', comparison0)
    
    fwrite(res.state0[Count >= 3 & P.value < 0.05], sep = '\t',
           file = paste0(filekey, '.tsv'))
  }
  
}


## plot hallmark up in coculture and down in treated, and vice versa ####
#### Hallmark
adj_thr = 0.05

enrichr_res_dir = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/')
comparisons = c('Coculture_vs_Monoculture', 'Afatinib_vs_Coculture', 'CRM_vs_Coculture')

### up and down
sele_terms = NULL

for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA20_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Afatinib_vs_Coculture', 'CRM_vs_Coculture')){
    efile = paste0(enrichr_res_dir, 'CHLA20_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  tmp$score = -log10(tmp$Adjusted.P.value)
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp = tmp[1:min(nrow(tmp), 10), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA20_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Afatinib_vs_Coculture', 'CRM_vs_Coculture')){
    efile = paste0(enrichr_res_dir, 'CHLA20_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  
  tmp$score = -log10(tmp$Adjusted.P.value)
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp$rk = order(-tmp$score)
  tmp[rk>10]$rk=10
  enrichr_res[[comparison0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'comparison'))

dt = dcast(dt, Term ~ comparison, value.var = 'score')

mat = dt[, -1]
mat = subset(mat, select = comparisons[comparisons %in% names(mat)]) ## order it
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap::pheatmap(mat, na_col = 'white', border_color = NA,
                        cluster_rows = T, cluster_cols = F, angle_col = 45)
#ggsave(aa, filename = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA20_heatmap_hallmark_UpDown.pdf'),
#       device = 'pdf', width = 6, height = 6)

### dotplot

enrichr_res$Term = factor(enrichr_res$Term, 
                          levels = rev(rownames(mat)[aa$tree_row$order]))


## rename for easy read
enrichr_res[comparison == 'Coculture_vs_Monoculture']$comparison = 'Up in Coculture'
enrichr_res[comparison == 'Afatinib_vs_Coculture']$comparison = 'Down in Afatinib'
enrichr_res[comparison == 'CRM_vs_Coculture']$comparison = 'Down in CRM'
enrichr_res$comparison = factor(enrichr_res$comparison, 
                                levels = c('Up in Coculture', 
                                           'Down in Afatinib', 
                                           'Down in CRM'))


enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$frac = enrichr_res$score
pu <- ggplot(data = enrichr_res, aes(x = comparison, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(option = "D") + theme_bw() +
  ggtitle(paste0('CHLA20-HALLMARK-p.adj<', adj_thr )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(p.adj)"), 
         size=guide_legend(title="-log10(p.adj)")) 

#ggsave(ph, filename = 'Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA20_dotplot_hallmark_UpDown.pdf',
#       device = 'pdf', width = 5, height = 6)

## down and up

#### Hallmark
enrichr_res_dir = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/')
comparisons = c('Coculture_vs_Monoculture', 'Afatinib_vs_Coculture', 'CRM_vs_Coculture')

### collect all top terms
sele_terms = NULL

for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA20_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Coculture_vs_Monoculture')){
    efile = paste0(enrichr_res_dir, 'CHLA20_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp = tmp[1:min(nrow(tmp), 10), ]
  sele_terms = unique(c(sele_terms, tmp$Term))
}

enrichr_res = list()
for(comparison0 in comparisons){
  efile = paste0(enrichr_res_dir, 'CHLA20_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  if(comparison0 %in% c('Coculture_vs_Monoculture')){
    efile = paste0(enrichr_res_dir, 'CHLA20_down_MSigDB_Hallmark_2020_', comparison0,  '.tsv')
  }
  if(!file.exists(efile)) next
  tmp = fread(efile)
  tmp$comparison = comparison0
  tmp = tmp[tmp$Term %in% sele_terms, ]
  
  tmp = tmp[Adjusted.P.value < adj_thr]
  
  tmp$rk = order(-tmp$score)
  tmp[rk>10]$rk=10
  enrichr_res[[comparison0]] = tmp
}
enrichr_res = do.call('rbind', enrichr_res)
enrichr_res = enrichr_res[!is.na(score)]
dt = subset(enrichr_res, select = c('Term', 'score', 'comparison'))

dt = dcast(dt, Term ~ comparison, value.var = 'score')

mat = dt[, -1]
mat = subset(mat, select = comparisons[comparisons %in% names(mat)]) ## order it
rownames(mat) = dt$Term
mat[is.na(mat)] = 0
mat[mat > 10] = 10
aa = pheatmap::pheatmap(mat, na_col = 'white', border_color = NA,
                        cluster_rows = T, cluster_cols = F, angle_col = 45)
#ggsave(aa, filename = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA20_heatmap_hallmark_DownUp.pdf'),
#       device = 'pdf', width = 6, height = 6)

### dotplot

enrichr_res$Term = factor(enrichr_res$Term, levels = rev(rownames(mat)[aa$tree_row$order]))


## rename for easy read
enrichr_res[comparison == 'Coculture_vs_Monoculture']$comparison = 'Down in Coculture'
enrichr_res[comparison == 'Afatinib_vs_Coculture']$comparison = 'Up in Afatinib'
enrichr_res[comparison == 'CRM_vs_Coculture']$comparison = 'Up in CRM'
enrichr_res$comparison = factor(enrichr_res$comparison, 
                                levels = c('Down in Coculture', 
                                           'Up in Afatinib', 
                                           'Up in CRM'))


enrichr_res[enrichr_res$score > 10]$score = 10
enrichr_res$frac = enrichr_res$score
pd <- ggplot(data = enrichr_res, aes(x = comparison, y = Term, color = score)) +
  geom_point(aes(size = frac)) + 
  scale_color_viridis_c(option = "D") + theme_bw() +
  ggtitle(paste0('CHLA20-HALLMARK-p.adj<', adj_thr )) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = '', y = '') + 
  guides(color=guide_legend(title="-log10(p.adj)"), 
         size=guide_legend(title="-log10(p.adj)")) 

#ggsave(ph, filename = 'Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA20_dotplot_hallmark_DownUp.pdf',
#       device = 'pdf', width = 5, height = 6)

ggsave(pu + pd, filename = paste0('Figures/RNA/coculture/enrichR_MSigDB_Hallmark_2020/adjustState/CHLA20_dotplot_hallmark_padj',
                                  adj_thr, '.pdf'),
       device = 'pdf', width = 10, height = 6)
