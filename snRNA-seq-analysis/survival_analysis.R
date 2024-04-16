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

os_cox_gene <- function(R2_clinical_list, R2_rna_list, sele_genes, surv_cutoff = 0.75,
                        regrOutMycn = FALSE){
  coef_mat = pv_mat = matrix(0, length(sele_genes), length(R2_clinical_list))
  colnames(coef_mat) = colnames(pv_mat) = names(R2_clinical_list)
  rownames(coef_mat) = rownames(pv_mat) = sele_genes
  
  for(data_name in names(R2_clinical_list)){
    R2_clinical = R2_clinical_list[[data_name]]
    R2_rna = R2_rna_list[[data_name]]
    for(interest_gene in sele_genes){
      if(!any(rownames(R2_rna) == interest_gene) || all(is.na(R2_rna[interest_gene, ]))) {
        coef_mat[interest_gene, data_name] = pv_mat[interest_gene, data_name] = NA
        next
      }
      tmp_res = osCoxPerGene(interest_gene, R2_clinical, R2_rna, surv_cutoff, regrOutMycn)
      pv_mat[interest_gene, data_name] = tmp_res[1]
      coef_mat[interest_gene, data_name] = tmp_res[2]
    }
    
  }
  coef_mat[is.na(coef_mat)] = 0
  pv_mat[is.na(pv_mat)] = 1
  return(list('pv_mat' = pv_mat, 'coef_mat' = coef_mat))
}

efs_cox_gene <- function(R2_clinical_list, R2_rna_list, sele_genes, surv_cutoff = 0.75,
                         regrOutMycn = FALSE){
  coef_mat = pv_mat = matrix(0, length(sele_genes), length(R2_clinical_list))
  colnames(coef_mat) = colnames(pv_mat) = names(R2_clinical_list)
  rownames(coef_mat) = rownames(pv_mat) = sele_genes
  
  for(data_name in names(R2_clinical_list)){
    R2_clinical = R2_clinical_list[[data_name]]
    R2_rna = R2_rna_list[[data_name]]
    for(interest_gene in sele_genes){
      if(!any(rownames(R2_rna) == interest_gene) || all(is.na(R2_rna[interest_gene, ]))) {
        coef_mat[interest_gene, data_name] = pv_mat[interest_gene, data_name] = NA
        next
      }
      tmp_res = efsCoxPerGene(interest_gene, R2_clinical, R2_rna, surv_cutoff, regrOutMycn)
      pv_mat[interest_gene, data_name] = tmp_res[1]
      coef_mat[interest_gene, data_name] = tmp_res[2]
    }
    
  }
  coef_mat[is.na(coef_mat)] = 0
  pv_mat[is.na(pv_mat)] = 1
  return(list('pv_mat' = pv_mat, 'coef_mat' = coef_mat))
}


prepSignatures <- function(R2_rna_list, sig_genes, gene_weights = NULL, sig_name = 'MES-like_0'){
  for(data_name in names(R2_rna_list)){
    R2_rna = R2_rna_list[[data_name]]
    samples = colnames(R2_rna)
    sig_genes = sig_genes[sig_genes %in% rownames(R2_rna)]
    if(!is.null(gene_weights)) gene_weights = gene_weights[sig_genes %in% rownames(R2_rna)]
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
for(data_name in c('SEQC', 'Versteeg', 'Target', 'Oberthuer', 'Cangelosi', 'NRC')){
  if(data_name == 'Versteeg'){
    R2_data = fread('data/raw/R2/ps_avgpres_nbadam88_u133p2_box1695232170-datagrabber-.txt')
    R2_clinical = t(R2_data[1:15, -c(1,2)])
    R2_clinical = data.frame(R2_clinical)
    names(R2_clinical) = R2_data[1:15, ]$probeset
    R2_clinical = data.table(R2_clinical, keep.rownames = T)
    names(R2_clinical)[1] = 'sample'
    
    R2_clinical$os = as.numeric(R2_clinical$nti_surv_overall)
    R2_clinical$efs = as.numeric(R2_clinical$nti_surv_progrfree)
    R2_clinical$os_status = R2_clinical$nti_event_overall
    R2_clinical$efs_status = R2_clinical$nti_event_progrfree
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'yes', 1, 0)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'yes', 1, 0)
    
    R2_clinical$inss = R2_clinical$inss
    R2_clinical$mycn_status = R2_clinical$mycn_amp
    R2_clinical$age = R2_clinical$age_year
    R2_clinical$gender = R2_clinical$gender
    
    R2_rna = R2_data[-c(1:15), -(1:2)]
    R2_rna = apply(R2_rna, 2, as.numeric)
    rownames(R2_rna) = R2_data$`#H:hugo`[-(1:15)]
  }
  
  if(data_name == 'Cangelosi'){
    R2_data = fread('data/raw/R2/ps_avgpres_dgc2102a786_dgc2102_box1695308502-datagrabber-.txt')
    R2_clinical = t(R2_data[1:10, -c(1,2)])
    R2_clinical = data.frame(R2_clinical)
    names(R2_clinical) = R2_data[1:10, ]$probeset
    R2_clinical = data.table(R2_clinical, keep.rownames = T)
    names(R2_clinical)[1] = 'sample'
    R2_clinical = R2_clinical[platform != 'rnaseq']
    R2_rna = R2_data[-c(1:10), -(1:2)]
    R2_rna = apply(R2_rna, 2, as.numeric)
    rownames(R2_rna) = R2_data$`#H:hugo`[-(1:10)]
    R2_rna = R2_rna[, R2_clinical$sample]
    
    R2_clinical$os  = as.integer(as.numeric(R2_clinical$overall_survival) * 365)
    R2_clinical$efs = as.integer(as.numeric(R2_clinical$`event-free_survival`) * 365)
    R2_clinical$os_status = R2_clinical$event_overall
    R2_clinical$efs_status = R2_clinical$`event-free`
    R2_clinical[efs_status == 'na'] = NA
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'yes', 1, 0)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'yes', 1, 0)
    
    R2_clinical$inss = R2_clinical$inss_stage
    R2_clinical$mycn_status = R2_clinical$mycn_status
    R2_clinical[mycn_status == 'na'] = NA
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
    R2_clinical[efs_status == 'na'] = NA
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'event', 1, 0)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'event', 1, 0)
    
    R2_clinical$inss = R2_clinical$inss_stage
    R2_clinical$mycn_status = R2_clinical$mycn_status
    R2_clinical[mycn_status == 'na'] = NA
    R2_clinical$mycn_status = ifelse(R2_clinical$mycn_status == 'mycn_amp', 'yes', 'no')
    R2_clinical$age = R2_clinical$age_at_diagnosis
    R2_clinical$gender = R2_clinical$sex
  }
  
  if(data_name == 'Oberthuer'){
    R2_data = fread('data/raw/R2/ps_avgpres_nb251_amexp255_box1695308338-datagrabber-.txt')
    R2_clinical = t(R2_data[1:13, -c(1,2)])
    R2_clinical = data.frame(R2_clinical)
    names(R2_clinical) = R2_data[1:13, ]$probeset
    R2_clinical = data.table(R2_clinical, keep.rownames = T)
    names(R2_clinical)[1] = 'sample'
    
    R2_rna = R2_data[-c(1:13), -(1:2)]
    R2_rna = apply(R2_rna, 2, as.numeric)
    rownames(R2_rna) = R2_data$`#H:hugo`[-(1:13)]
    
    R2_clinical$os  = as.numeric(R2_clinical$overall_surv_days)
    R2_clinical$efs = as.numeric(R2_clinical$eventfree_surv_days)
    
    R2_clinical$os_status = R2_clinical$overall_surv_event
    R2_clinical$efs_status = R2_clinical$eventfree_surv_event
    
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'no event', 0, 1)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'no event', 0, 1)
    
    R2_clinical$inss = R2_clinical$inss
    R2_clinical$mycn_status = ifelse(R2_clinical$mycn_amp == 'amplified', 'yes', 'no')
    R2_clinical$age = R2_clinical$age_diagnosis_days
    R2_clinical$gender = R2_clinical$sex
  }
  
  if(data_name == 'NRC'){
    R2_data = fread('data/raw/R2/ps_avgpres_gse85047ageot283_huex10t_box1695413974-datagrabber-.txt')
    R2_clinical = t(R2_data[1:9, -c(1,2)])
    R2_clinical = data.frame(R2_clinical)
    names(R2_clinical) = R2_data[1:9, ]$probeset
    R2_clinical = data.table(R2_clinical, keep.rownames = T)
    names(R2_clinical)[1] = 'sample'
    
    R2_rna = R2_data[-c(1:9), -(1:2)]
    R2_rna = apply(R2_rna, 2, as.numeric)
    rownames(R2_rna) = R2_data$`#H:hugo`[-(1:9)]
    
    R2_clinical$os  = as.numeric(R2_clinical$`overall_survival_time(days)`)
    R2_clinical$efs = as.numeric(R2_clinical$`profression_free_survival_time(days)`)
    
    R2_clinical$os_status = R2_clinical$event_overall
    R2_clinical$efs_status = R2_clinical$event_progression_free
    
    R2_clinical$os_status = ifelse(R2_clinical$os_status == 'no', 0, 1)
    R2_clinical$efs_status = ifelse(R2_clinical$efs_status == 'no', 0, 1)
    
    R2_clinical$inss = R2_clinical$inss
    R2_clinical$mycn_status = R2_clinical$mycn_amplification
    R2_clinical$age = as.numeric(R2_clinical$`age_at_diagnosis(days)`)
    
  }
  
  if(data_name == 'Target'){
    #zscore or rpkm
    target_rna_sample = fread('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/bulk_rna/public/target/NB/cbioportal/nbl_target_2018_pub/data_RNA_Seq_mRNA_median_Zscores.txt')
    sample_inf = fread('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/bulk_rna/public/target/NB/cbioportal/nbl_target_2018_pub/data_clinical_sample.txt',
                       skip = 4)
    target_clinical = readRDS('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/bulk_rna/public/target/NB/nci/clinical/clinical.rds')
    target_clinical = data.table(target_clinical)
    #target_clinical$perc_tumor = stringr::str_split_i(target_clinical$PERCENT_TUMOR_VS_STROMA, '/', 1)
    
    sample_inf = sample_inf[SAMPLE_ID %in% names(target_rna_sample)]
    
    # create avg patient level expression
    target_rna = target_rna_sample[, 1:2]
    for(patientID in unique(sample_inf$PATIENT_ID)){
      sampleIDs = sample_inf[PATIENT_ID == patientID]$SAMPLE_ID
      target_rna[[patientID]] = rowMeans(subset(target_rna_sample, select = sampleIDs))
    }
    
    genes = target_rna$Hugo_Symbol
    target_rna = as.matrix(target_rna[, -(1:2)])
    rownames(target_rna) = genes
    
    target_rna = target_rna[, colnames(target_rna) %in% target_clinical$PATIENT_ID]
    names(target_clinical)[1] = 'sample'
    target_clinical = target_clinical[sample %in% colnames(target_rna), ]
    R2_clinical = target_clinical
    R2_rna = target_rna[, R2_clinical$sample]
    
    R2_clinical$os  = as.numeric(R2_clinical$OS_DAYS)
    R2_clinical$efs = as.numeric(R2_clinical$EFS_TIME)
    
    R2_clinical$os_status = as.integer(stringr::str_split_i(R2_clinical$OS_STATUS, ':', 1))
    R2_clinical$efs_status = R2_clinical$EFSCENS
    
    R2_clinical$inss = R2_clinical$INSS_STAGE
    
    R2_clinical$age = R2_clinical$AGE
    R2_clinical$gender = R2_clinical$SEX
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

degs$cell_state[degs$cluster == 4]  = 'ADRN_Proliferating_4'
degs$cell_state[degs$cluster == 0] = 'ADRN_Calcium_0'
degs$cell_state[degs$cluster == 1] = 'ADRN_Ribosome_1'
degs$cell_state[degs$cluster == 2] = 'Interm_OxPhos_2'
degs$cell_state[degs$cluster == 3] = 'ADRN_Dopaminergic_3'
degs$cell_state[degs$cluster == 5] = 'MES_5'

R2_rna_list = readRDS('data/intermediate/R2/list_rna_bulk.rds')
R2_clinical_list = readRDS('data/intermediate/R2/list_clinical.rds')

## prepare signatures
for(cell_state0 in unique(degs$cluster)){
  degs0 = degs[cluster == cell_state0]
  #if(cell_state0 != 1) degs0 = degs0[pct.diff > 0.1]
  degs0 = degs0[order(p_val)]
  #degs0 = degs0[1:min(50, nrow(degs0)), ]
  genes0 = degs0$gene
  R2_rna_list <- prepSignatures(R2_rna_list, genes0, gene_weights = degs0$pct.diff * degs0$avg_log2FC,
                                sig_name = degs0$cell_state[1])
}

## < plot survive curve for a cell state ####
regrOutMycn = T
data_name = 'Cangelosi' # SEQC or Oberthuer or Cangelosi
plot_lists = list()
scutoff = 0.75
for(state0 in c('ADRN_Ribosome_1', 'ADRN_Proliferating_4', 'ADRN_Calcium_0',
                'ADRN_Dopaminergic_3', 'Interm_OxPhos_2', 'MES_5')){
  plot_lists[[state0]] <- osCoxPerGene(state0, R2_clinical_list[[data_name]], R2_rna_list[[data_name]], surv_cutoff = scutoff,
                                      regrOutMycn = regrOutMycn, KMPlot = T)[[3]]
}

kmplot_file = paste0('Figures/survival/v3/states/KMplot_by_state_', data_name, '.pdf')
if(regrOutMycn) kmplot_file = paste0('Figures/survival/v3/states/KMplot_by_state_mycnAdjusted', data_name, '.pdf')
ggsave(gridExtra::grid.arrange(grobs = plot_lists, ncol = 3), 
       file = kmplot_file,
       device = 'pdf', width = 9, height = 6)
  
## < plot ef survival curve ####
regrOutMycn = T
data_name = 'Cangelosi' # SEQC or Cangelosi
plot_lists = list()
scutoff = 0.75
for(state0 in c('ADRN_Ribosome_1', 'ADRN_Proliferating_4', 'ADRN_Calcium_0',
                'ADRN_Dopaminergic_3', 'Interm_OxPhos_2', 'MES_5')){
  plot_lists[[state0]] <- efsCoxPerGene(state0, R2_clinical_list[[data_name]], R2_rna_list[[data_name]], surv_cutoff = scutoff,
                                       regrOutMycn = regrOutMycn, KMPlot = T)[[3]]
}

kmplot_file = paste0('Figures/survival/v3/states/KMplot_efs_by_state_', data_name, '.pdf')
if(regrOutMycn) kmplot_file = paste0('Figures/survival/v3/states/KMplot_efs_by_state_mycnAdjusted', data_name, '.pdf')
ggsave(gridExtra::grid.arrange(grobs = plot_lists, ncol = 3), 
       file = kmplot_file,
       device = 'pdf', width = 9, heigh = 6)
       
## < assign state by max score ####
data_name = 'Cangelosi' # SEQC or Oberthuer or NRC or Target or Versteeg or Cangelosi

tmp = R2_rna_list[[data_name]][c('ADRN_Ribosome_1', 'ADRN_Proliferating_4', 'ADRN_Calcium_0',
                         'ADRN_Dopaminergic_3', 'Interm_OxPhos_2', 'MES_5'), ]
id_max = apply(tmp, 2, which.max)
assign_state = c('ADRN_Ribosome_1', 'ADRN_Proliferating_4', 'ADRN_Calcium_0',
                 'ADRN_Dopaminergic_3', 'Interm_OxPhos_2', 'MES_5')[id_max]
tmp = rbind(tmp, 'assigned_state' = assign_state)
R2_clinical = R2_clinical_list[[data_name]]
R2_clinical$state = assign_state
if(length(unique(R2_clinical$age)) > 3) R2_clinical$age = as.numeric(R2_clinical$age)
s_os <- survfit(Surv(os, os_status) ~ state, data = R2_clinical)
s_efs <- survfit(Surv(efs, efs_status) ~ state, data = R2_clinical)

if(data_name == 'Target'){
  cox_os <- coxph(Surv(os, os_status) ~ state + age + gender, data = R2_clinical)
  cox_efs <- coxph(Surv(efs, efs_status) ~ state + age + gender, data = R2_clinical)
}else{
  if(any(names(R2_clinical) == 'gender')){
    cox_os <- coxph(Surv(os, os_status) ~ state + age + gender + mycn_status, data = R2_clinical)
    cox_efs <- coxph(Surv(efs, efs_status) ~ state + age + gender + mycn_status, data = R2_clinical)
  }else{
    cox_os <- coxph(Surv(os, os_status) ~ state + age + mycn_status, data = R2_clinical)
    cox_efs <- coxph(Surv(efs, efs_status) ~ state + age + mycn_status, data = R2_clinical)
  }
  
}


summary(cox_os)
state_colors = brewer.pal(6, 'Set1')
names(state_colors) = c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                        'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5')
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
os_plot <- os_plot$plot + theme(legend.title=element_blank()) 
os_plot

summary(cox_efs)
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
efs_plot <- efs_plot$plot + theme(legend.title=element_blank()) 
efs_plot

ggsave(os_plot, filename = paste0('Figures/survival/v3/states/KMplot_by_AssignedState_', data_name, '.pdf'),
       device = 'pdf', width = 5, height = 5)
ggsave(efs_plot, filename =paste0('Figures/survival/v3/states/KMplot_efs_by_AssignedState_', data_name, '.pdf'),
       device = 'pdf', width = 5, height = 5)
table(assign_state)

