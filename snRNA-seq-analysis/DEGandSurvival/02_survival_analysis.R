## ***do the survival analysis systematically using data from R2 database ***##

library(data.table)
library(magrittr)
library(survival)
library(ggplot2)
library(ggsurvfit)
library(ggfortify) ## enable autoplot
library(survminer)
library(RColorBrewer)
library(pheatmap)

osCoxPerGene <- function(interest_gene, R2_clinical, R2_rna, surv_cutoff = 0.75,
                         regrOutMycn = F, KMPlot = F){
  
  R2_clinical$StateScore = ifelse(R2_rna[interest_gene, ] > 
                                    quantile(R2_rna[interest_gene, ], surv_cutoff, na.rm = T), 
                                  'High', 'Low')
  if(var(R2_rna[interest_gene, ], na.rm = T) == 0) {
    coef = 0
    pv = 1
    return(c(pv, coef))
  }
  
  if(length(unique(R2_clinical$age)) > 3){
    R2_clinical$age = as.numeric(R2_clinical$age)
  }
  
  if(regrOutMycn){
    if(any(names(R2_clinical) == 'gender')){
      s3 <- coxph(Surv(os, os_status) ~ StateScore + age + mycn_status + gender, 
                  data = R2_clinical)
    }else{
      s3 <- coxph(Surv(os, os_status) ~ StateScore + age + mycn_status, 
                  data = R2_clinical)
    }
    
  }else{
    if(any(names(R2_clinical) == 'gender')){
      s3 <- coxph(Surv(os, os_status) ~ StateScore + age + gender,
                  data = R2_clinical) 
    }else{
      s3 <- coxph(Surv(os, os_status) ~ StateScore + age, data = R2_clinical)
    }
    
  }
  
  pv = coef(summary(s3))['StateScoreLow', 5]
  coef = round(coef(summary(s3))['StateScoreLow', 1], 2)
  
  if(KMPlot){
    s2 <- survfit(Surv(os, os_status) ~ StateScore, data = R2_clinical)
    
    gs_plot <- ggsurvplot(
      s2,
      data = R2_clinical,
      censor = F,
      size = 0.5,                 # change line size
      palette =
        c("#E7B800", "#2E9FDF"),# custom color palettes
      conf.int = FALSE,          # Add confidence interval
      pval = paste('pv =', round(pv, 3), '\ncox_coef =',  coef),              # Add p-value
      risk.table = TRUE,        # Add risk table
      risk.table.col = 'strata',# Risk table color by groups
      risk.table.height = 0.25, # Useful to change when you have multiple groups
      xlab = "Time in days",
      title = interest_gene,
      ggtheme = theme_classic()      # Change ggplot2 theme
    ) 
    gs_plot$plot <- gs_plot$plot + theme(legend.title=element_blank()) 
    return(list(pv, coef, gs_plot$plot))
  }else{
    return(c(pv, coef))
  }
  
}

efsCoxPerGene <- function(interest_gene, R2_clinical, R2_rna, surv_cutoff = 0.75,
                          regrOutMycn = F, KMPlot = F){
  
  R2_clinical$StateScore = ifelse(R2_rna[interest_gene, ] > 
                                    quantile(R2_rna[interest_gene, ], surv_cutoff, na.rm = T), 
                                  'High', 'Low')
  if(var(R2_rna[interest_gene, ], na.rm = T) == 0) {
    coef = 0
    pv = 1
    return(c(pv, coef))
  }
  if(length(unique(R2_clinical$age)) > 3){
    R2_clinical$age = as.numeric(R2_clinical$age)
  }
  if(regrOutMycn){
    if(any(names(R2_clinical) == 'gender')){
      s3 <- coxph(Surv(efs, efs_status) ~ StateScore + age + mycn_status + gender, 
                  data = R2_clinical)
    }else{
      s3 <- coxph(Surv(efs, efs_status) ~ StateScore + age + mycn_status, 
                  data = R2_clinical)
    }
    
  }else{
    if(any(names(R2_clinical) == 'gender')){
      s3 <- coxph(Surv(efs, efs_status) ~ StateScore + age + gender,
                  data = R2_clinical) 
    }else{
      s3 <- coxph(Surv(efs, efs_status) ~ StateScore + age, data = R2_clinical)
    }
    
  }
  
  pv = coef(summary(s3))['StateScoreLow', 5]
  coef = round(coef(summary(s3))['StateScoreLow', 1], 2)
  if(KMPlot){
    s2 <- survfit(Surv(efs, efs_status) ~ StateScore , 
                  data = R2_clinical)
    
    gs_plot <- ggsurvplot(
      s2,
      data = R2_clinical,
      censor = F,
      size = 0.5,                 # change line size
      palette =
        c("#E7B800", "#2E9FDF"),# custom color palettes
      conf.int = FALSE,          # Add confidence interval
      pval = paste('pv =', round(pv, 3), '\ncox_coef =',  coef),              # Add p-value
      risk.table = TRUE,        # Add risk table
      risk.table.col = 'strata',# Risk table color by groups
      risk.table.height = 0.25, # Useful to change when you have multiple groups
      xlab = "Time in days",
      ylab = 'EFS probability',
      title = interest_gene,
      ggtheme = theme_classic()      # Change ggplot2 theme
    ) 
    gs_plot$plot <- gs_plot$plot + theme(legend.title=element_blank()) 
    return(list(pv, coef, gs_plot$plot))
  }else{
    return(c(pv, coef))
  }
}


prepSignatures <- function(R2_rna_list, sig_genes, gene_weights = NULL, 
                           sig_name = 'MES-like'){
  for(data_name in names(R2_rna_list)){
    R2_rna = R2_rna_list[[data_name]]
    samples = colnames(R2_rna)
    if(!is.null(gene_weights)) gene_weights = gene_weights[sig_genes %in% rownames(R2_rna)]
    sig_genes = sig_genes[sig_genes %in% rownames(R2_rna)]
    sigs <- sapply(samples, function(x){
      sig0 = mean(R2_rna[sig_genes, x] * gene_weights, na.rm = T)
      sig0[is.na(sig0)] = 0
      return(sig0)
    })
    R2_rna = rbind(R2_rna, sigs)
    rownames(R2_rna)[nrow(R2_rna)] = sig_name
    R2_rna_list[[data_name]] = R2_rna
  }
  return(R2_rna_list)
}

# data preparation ####
R2_rna_list = list()
R2_clinical_list = list()
for(data_name in c('SEQC', 'Cangelosi')){
 
  if(data_name == 'Cangelosi'){
    R2_data = fread('data/raw/R2/ps_avgpres_dgc2102a786_dgc2102_box1695308502-datagrabber-.txt')
    R2_clinical = t(R2_data[1:10, -c(1,2)])
    R2_clinical = data.frame(R2_clinical)
    names(R2_clinical) = R2_data[1:10, ]$probeset
    R2_clinical = data.table(R2_clinical, keep.rownames = T)
    names(R2_clinical)[1] = 'sample'
    R2_clinical = R2_clinical[platform != 'rnaseq'] # the rna-seq overlaped with SEQC
    R2_rna = R2_data[-c(1:10), -(1:2)]
    R2_rna = apply(R2_rna, 2, as.numeric)
    rownames(R2_rna) = R2_data$`#H:hugo`[-(1:10)]
    R2_rna = R2_rna[, R2_clinical$sample]
    
    R2_clinical$os  = as.integer(as.numeric(R2_clinical$overall_survival) * 365)
    R2_clinical$efs = as.integer(as.numeric(R2_clinical$`event-free_survival`) * 365)
    R2_clinical$os_status = R2_clinical$event_overall
    R2_clinical$efs_status = R2_clinical$`event-free`
    R2_clinical[efs_status == 'na']$efs_status = NA
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'yes', 1, 0)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'yes', 1, 0)
    
    R2_clinical$inss = R2_clinical$inss_stage
    R2_clinical$mycn_status = R2_clinical$mycn_status
    R2_clinical[mycn_status == 'na']$efs_status = NA
    R2_clinical$age = R2_clinical$agr_group
    
  }
  
  if(data_name == 'SEQC'){
    R2_data = fread('data/raw/R2/ps_avgpres_gse62564geo498_seqcnb1_box1695408960-datagrabber-.txt')
    R2_clinical = t(R2_data[1:16, -c(1,2)])
    R2_clinical = data.frame(R2_clinical)
    names(R2_clinical) = R2_data[1:16, ]$probeset
    R2_clinical = data.table(R2_clinical, keep.rownames = T)
    names(R2_clinical)[1] = 'sample'
    R2_rna = R2_data[-c(1:16), -(1:2)]
    R2_rna = apply(R2_rna, 2, as.numeric)
    rownames(R2_rna) = R2_data$`#H:hugo`[-(1:16)]
    
    R2_clinical$os  = as.numeric(R2_clinical$os_day)
    R2_clinical$efs = as.numeric(R2_clinical$efs_day)
    R2_clinical$os_status = R2_clinical$os_bin
    R2_clinical$efs_status = R2_clinical$efs_bin
    R2_clinical[efs_status == 'na']$efs_status = NA
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'event', 1, 0)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'event', 1, 0)
    
    R2_clinical$inss = R2_clinical$inss_stage
    R2_clinical$mycn_status = R2_clinical$mycn_status
    R2_clinical[mycn_status == 'na']$mycn_status = NA
    R2_clinical$mycn_status = ifelse(R2_clinical$mycn_status == 'mycn_amp', 'yes', 'no')
    R2_clinical$age = R2_clinical$age_at_diagnosis
    R2_clinical$gender = R2_clinical$sex
  }
  
  R2_rna_list[[data_name]] = R2_rna
  R2_clinical_list[[data_name]] = R2_clinical
}
saveRDS(R2_rna_list, file = 'data/intermediate/R2/list_rna_bulk.rds')
saveRDS(R2_clinical_list, file = 'data/intermediate/R2/list_clinical.rds')


# Overall survival/EFS signature per cell state ####
degs = readRDS(file = 'data/intermediate/RNA/new_degs/degs_malignant_cluster.rds')
degs[, 'pct.diff' := pct.1 - pct.2]
degs = degs[avg_log2FC > 0.5 & p_val < 0.001 & pct.diff > 0]

degs$cell_state[degs$cluster == 4]  = 'ADRN-Proliferating'
degs$cell_state[degs$cluster == 0] = 'ADRN-Calcium'
degs$cell_state[degs$cluster == 1] = 'ADRN-Baseline'
degs$cell_state[degs$cluster == 2] = 'Interm-OxPhos'
degs$cell_state[degs$cluster == 3] = 'ADRN-Dopaminergic'
degs$cell_state[degs$cluster == 5] = 'MES'

R2_rna_list = readRDS('data/intermediate/R2/list_rna_bulk.rds')
R2_clinical_list = readRDS('data/intermediate/R2/list_clinical.rds')

## prepare signatures
for(cell_state0 in unique(degs$cell_state)){
  degs0 = degs[cell_state == cell_state0]
  degs0 = degs0[order(p_val)]
  genes0 = degs0$gene
  R2_rna_list <- prepSignatures(R2_rna_list, genes0, 
                                gene_weights = degs0$pct.diff * degs0$avg_log2FC,
                                sig_name = degs0$cell_state[1])
}

## < plot os survive curve for a cell state ####
regrOutMycn = T
for(data_name in c('SEQC', 'Cangelosi')){
  plot_lists = list()
  scutoff = 0.5
  for(state0 in c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                  'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')){
    plot_lists[[state0]] <- osCoxPerGene(state0, R2_clinical_list[[data_name]], R2_rna_list[[data_name]], surv_cutoff = scutoff,
                                         regrOutMycn = regrOutMycn, KMPlot = T)[[3]]
  }
  
  kmplot_file = paste0('Figures/survival/Revision/Tumor/KMplot_by_state_', data_name, '_cutoff', scutoff, '.pdf')
  if(regrOutMycn) kmplot_file = paste0('Figures/survival/Revision/Tumor/KMplot_by_state_mycnAdjusted_', data_name, 
                                       '_cutoff', scutoff, '.pdf')
  ggsave(gridExtra::grid.arrange(grobs = plot_lists, ncol = 3), 
         file = kmplot_file,
         device = 'pdf', width = 9, height = 6)
  
  
}

## < plot efs survival curve ####
regrOutMycn = T
for(data_name in c('SEQC', 'Cangelosi')){
  plot_lists = list()
  scutoff = 0.5
  for(state0 in c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                  'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')){
    plot_lists[[state0]] <- efsCoxPerGene(state0, R2_clinical_list[[data_name]], R2_rna_list[[data_name]], surv_cutoff = scutoff,
                                          regrOutMycn = regrOutMycn, KMPlot = T)[[3]]
  }
  
  kmplot_file = paste0('Figures/survival/Revision/Tumor/KMplot_efs_by_state_', data_name, '_cutoff', scutoff, '.pdf')
  if(regrOutMycn) kmplot_file = paste0('Figures/survival/Revision/Tumor/KMplot_efs_by_state_mycnAdjusted_', data_name, 
                                       '_cutoff', scutoff, '.pdf')
  ggsave(gridExtra::grid.arrange(grobs = plot_lists, ncol = 3), 
         file = kmplot_file,
         device = 'pdf', width = 9, heigh = 6)
  
}
   
## < assign state by max score ####
data_name = 'Cangelosi' ## or SEQ

tmp = R2_rna_list[[data_name]][c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                                 'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES'), ]
id_max = apply(tmp, 2, which.max)
assign_state = c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                 'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')[id_max]
tmp = rbind(tmp, 'assigned_state' = assign_state)
R2_clinical = R2_clinical_list[[data_name]]
R2_clinical$state = assign_state
if(length(unique(R2_clinical$age)) > 3) R2_clinical$age = as.numeric(R2_clinical$age)

s_os <- survfit(Surv(os, os_status) ~ state, data = R2_clinical)
s_efs <- survfit(Surv(efs, efs_status) ~ state, data = R2_clinical)

R2_clinical$state = factor(R2_clinical$state)
R2_clinical = dplyr::mutate(R2_clinical, state = relevel(state, ref="ADRN-Calcium"))


if(any(names(R2_clinical) == 'gender')){
  cox_os <- coxph(Surv(os, os_status) ~ state + age + gender + mycn_status, 
                  data = R2_clinical)
  cox_efs <- coxph(Surv(efs, efs_status) ~ state + age + gender + mycn_status, 
                   data = R2_clinical)
}else{
  cox_os <- coxph(Surv(os, os_status) ~ state + age + mycn_status, 
                  data = R2_clinical)
  cox_efs <- coxph(Surv(efs, efs_status) ~ state + age + mycn_status, 
                   data = R2_clinical)
}


state_colors = brewer.pal(6, 'Set1')
names(state_colors) = paste0('state=', c('ADRN-Calcium', 'ADRN-Baseline', 'Interm-OxPhs',
                                         'ADRN-Dopaminergic', 'ADRN-Proliferating', 'MES'))

os_plot <- ggsurvplot(
  s_os,
  data = R2_clinical,
  censor = F,
  size = 0.5,                 # change line size
  conf.int = FALSE,          # Add confidence interval
  risk.table = F,        # Add risk table
  risk.table.col = 'strata',# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  xlab = "Time in days",
  ggtheme = theme_classic()      # Change ggplot2 theme
) 
os_plot <- os_plot$plot + theme(legend.title=element_blank()) +
  scale_color_manual(values = state_colors)
os_plot
summary(cox_os)

efs_plot <- ggsurvplot(
  s_efs,
  data = R2_clinical,
  censor = F,
  size = 0.5,                 # change line size
  conf.int = FALSE,          # Add confidence interval
  risk.table = F,        # Add risk table
  risk.table.col = 'strata',# Risk table color by groups
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  xlab = "Time in days",
  ylab = 'EFS probability',
  ggtheme = theme_classic()      # Change ggplot2 theme
) 
efs_plot <- efs_plot$plot + theme(legend.title=element_blank()) +
  scale_color_manual(values = state_colors)
efs_plot
summary(cox_efs)
table(assign_state)

ggsave(os_plot, filename = paste0('Figures/survival/Revision/Tumor/KMplot_by_AssignedState_', data_name, '.pdf'),
       device = 'pdf', width = 5, height = 5)
ggsave(efs_plot, filename =paste0('Figures/survival/Revision/Tumor/KMplot_efs_by_AssignedState_', data_name, '.pdf'),
       device = 'pdf', width = 5, height = 5)

