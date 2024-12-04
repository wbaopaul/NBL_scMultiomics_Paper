library(Seurat)
library(data.table)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(survival)
library(ggsurvfit)
library(survminer)
library(ggpubr)

osCoxPerGene <- function(interest_gene, R2_clinical, R2_rna, surv_cutoff = 0.75,
                         regrOutMycn = F, KMPlot = F){
  
  R2_clinical$StateScore = ifelse(R2_rna[interest_gene, ] > 
                                    quantile(R2_rna[interest_gene, ], surv_cutoff, na.rm = T), 
                                  'High', 'Low')
  nhigh = nrow(R2_clinical[StateScore == 'High'])
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
    return(list(pv, coef, gs_plot$plot, nhigh))
  }else{
    return(c(pv, coef, nhigh))
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


## -- post civersortx analysis --##

## SEQC data ####
decov_result = fread('data/intermediate/RNA/CIBERSORTx/eachGroupDownTo10kCells_Data2_seqc_nbl/CIBERSORTx_Adjusted.txt')

### correlate with the clinical inf ####
seqc_data = fread('data/raw/R2/ps_avgpres_gse62564geo498_seqcnb1_box1701183946-datagrabber-rpm.txt')
seqc_clinical = t(seqc_data[1:16, -c(1,2)])
seqc_clinical = data.frame(seqc_clinical)
names(seqc_clinical) = seqc_data[1:16, ]$probeset
seqc_clinical = data.table(seqc_clinical, keep.rownames = T)
names(seqc_clinical)[1] = 'sample'
setkey(seqc_clinical, sample)
decov_result$inss_stage = seqc_clinical[decov_result$Mixture]$inss_stage

states = names(decov_result)[2:21]
malignant_frac = subset(decov_result, select = c('ADRN-Dopaminergic', 'MES',
                                                 'Interm-OxPhos', 'ADRN-Calcium',
                                                 'ADRN-Proliferating', 'ADRN-Baseline'))
decov_result$malignant = rowSums(malignant_frac)
malignant_states = c('ADRN-Dopaminergic', 'MES',
                     'Interm-OxPhos', 'ADRN-Calcium',
                     'ADRN-Proliferating', 'ADRN-Baseline')

decov_result[inss_stage == 'st1']$inss_stage = 'stage_1'
decov_result[inss_stage == 'st2']$inss_stage = 'stage_2'
decov_result[inss_stage == 'st3']$inss_stage = 'stage_3'
decov_result[inss_stage == 'st4']$inss_stage = 'stage_4'
decov_result[inss_stage == 'st4s']$inss_stage = 'stage_4s'

decov_result$inss_stage = factor(decov_result$inss_stage,
                                 levels = c('stage_1', 'stage_2', 'stage_3',
                                            'stage_4s', 'stage_4'))
plist = list()
for(state0 in states){
  decov_result$fraction = decov_result[[state0]]
 
  my_comparisons <- list( c("stage_1", "stage_4"), c("stage_2", "stage_4"), 
                          c("stage_3", "stage_4"),
                          c('stage_4s', 'stage_4'))
  sig_plot <- ggplot(data = decov_result, aes(x = inss_stage, y = fraction, 
                                              fill = inss_stage)) +
    geom_boxplot() + theme_classic() + 
    scale_fill_manual(values = brewer.pal(5, 'Reds')) +
    theme(legend.position = 'none') + ggtitle(state0) +
    stat_compare_means(comparisons = my_comparisons)
  plist[[state0]] = sig_plot
  
}

ggsave(gridExtra::grid.arrange(grobs = plist[c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                                               'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')], ncol = 6), 
       file = 'Figures/RNA/CIBERSORTx/ByStage/seqc_fraction_malignantStates_byStage.pdf',
       device = 'pdf', width = 20, height = 4)


ggsave(gridExtra::grid.arrange(grobs = plist[c('THY1+Macro', 'HS3ST2+Macro',
                                               'F13A1+Macro', 'CCL4+Macro', 
                                               'IL18+Macro','VCAN+Macro',
                                               'C1QC+SPP1+Macro', 'ProliferatingMacro')], 
                               ncol = 8), 
       file = 'Figures/RNA/CIBERSORTx/ByStage/seqc_fraction_macroStates_byStage.pdf',
       device = 'pdf', width = 20, height = 4)

decov_result$mycn_status = seqc_clinical[decov_result$Mixture]$mycn_status
decov_result = decov_result[mycn_status != 'na']
decov_result[mycn_status == 'mycn_amp']$mycn_status = 'yes'
decov_result[mycn_status == 'mycn_nonamp']$mycn_status = 'no'

plist = list()
for(state0 in states){
  decov_result$fraction = decov_result[[state0]]
  
  pv = wilcox.test(decov_result[mycn_status=='yes']$fraction, 
                   decov_result[mycn_status=='no']$fraction)$p.val/2
  
  sig_plot <- ggplot(data = decov_result, aes(x = mycn_status, y = fraction, 
                                              color = mycn_status)) +
    geom_boxplot() + theme_classic() + 
    scale_color_brewer(palette = 'Dark2') + 
    theme(legend.position = 'none') + ggtitle(state0) +
    annotate('text', x = 1.5, y = max(decov_result$fraction), 
             label = paste0('p = ', format(pv, digits = 2)))
  
  plist[[state0]] = sig_plot
}

ggsave(gridExtra::grid.arrange(grobs = plist[c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                                               'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')], ncol = 6), 
       file = 'Figures/RNA/CIBERSORTx/ByMycn/seqc_fraction_malignantStates_byMYCN.pdf',
       device = 'pdf', width = 12, height = 2)

ggsave(gridExtra::grid.arrange(grobs = plist[c('THY1+Macro', 'HS3ST2+Macro',
                                               'F13A1+Macro', 'CCL4+Macro', 
                                               'IL18+Macro','VCAN+Macro',
                                               'C1QC+SPP1+Macro', 'ProliferatingMacro')], ncol = 8), 
       file = 'Figures/RNA/CIBERSORTx/ByMycn/seqc_fraction_macroStates_byMYCN.pdf',
       device = 'pdf', width = 12, height = 2)


# Cangelosi data ####
### correlate with the clinical inf ####
decov_result = fread('data/intermediate/RNA/CIBERSORTx/eachGroupDownTo10kCells_Data2_cangelosi/CIBERSORTx_Adjusted.txt')

R2_rna_list = readRDS(file = 'data/intermediate/R2/list_rna_bulk.rds') ## prepared in survival analysis
R2_clinical_list = readRDS(file = 'data/intermediate/R2/list_clinical.rds')
cangelosi_clinical = R2_clinical_list[['Cangelosi']]
cangelosi_rna = R2_rna_list[['Cangelosi']]
setkey(cangelosi_clinical, sample)
decov_result$inss_stage = cangelosi_clinical[decov_result$Mixture]$inss_stage

states = names(decov_result)[2:21]
malignant_frac = subset(decov_result, select = c('ADRN-Dopaminergic', 'MES',
                                                 'Interm-OxPhos', 'ADRN-Calcium',
                                                 'ADRN-Proliferating', 'ADRN-Baseline'))
decov_result$malignant = rowSums(malignant_frac)
malignant_states = c('ADRN-Dopaminergic', 'MES',
                     'Interm-OxPhos', 'ADRN-Calcium',
                     'ADRN-Proliferating', 'ADRN-Baseline')
decov_result = decov_result[inss_stage != 'stage_na']

decov_result$inss_stage = factor(decov_result$inss_stage,
                                 levels = c('stage_1', 'stage_2', 'stage_3',
                                            'stage_4s', 'stage_4'))
plist = list()
for(state0 in states){
  decov_result$fraction = decov_result[[state0]]
  my_comparisons <- list( c("stage_1", "stage_4"), c("stage_2", "stage_4"), 
                          c("stage_3", "stage_4"),
                          c('stage_4s', 'stage_4'))
  
  sig_plot <- ggplot(data = decov_result, aes(x = inss_stage, y = fraction, 
                                              fill = inss_stage)) +
    geom_boxplot() + theme_classic() + 
    scale_fill_manual(values = brewer.pal(5, 'Reds')) +
    theme(legend.position = 'none') + ggtitle(state0) +
    stat_compare_means(comparisons = my_comparisons)
  
  plist[[state0]] = sig_plot
  state1 = gsub('+', '_', state0, fixed = T)
  ggsave(sig_plot, filename = paste0('Figures/RNA/CIBERSORTx/ByStage/cangelosi_fraction_', 
                                     state1, '_ByStage.pdf'), width = 5, height = 4)
  
}
ggsave(gridExtra::grid.arrange(grobs = plist[c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                                               'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')], ncol = 6), 
       file = 'Figures/RNA/CIBERSORTx/ByStage/cangelosi_fraction_malignantStates_byStage.pdf',
       device = 'pdf', width = 20, height = 4)

ggsave(gridExtra::grid.arrange(grobs = plist[c('THY1+Macro', 'HS3ST2+Macro',
                                               'F13A1+Macro', 'CCL4+Macro', 
                                               'IL18+Macro','VCAN+Macro',
                                               'C1QC+SPP1+Macro', 'ProliferatingMacro')], ncol = 8), 
       file = 'Figures/RNA/CIBERSORTx/ByStage/cangelosi_fraction_macroStates_byStage.pdf',
       device = 'pdf', width = 20, height = 4)

decov_result$mycn_status = cangelosi_clinical[decov_result$Mixture]$mycn_status
decov_result = decov_result[mycn_status != 'na']
plist = list()
for(state0 in states){
  decov_result$fraction = decov_result[[state0]]
  
  pv = wilcox.test(decov_result[mycn_status=='yes']$fraction, 
                   decov_result[mycn_status=='no']$fraction)$p.val/2
  
  sig_plot <- ggplot(data = decov_result, aes(x = mycn_status, y = fraction, 
                                                          color = mycn_status)) +
    geom_boxplot() + theme_classic() + 
    scale_color_brewer(palette = 'Dark2') + 
    theme(legend.position = 'none') + ggtitle(state0) +
    annotate('text', x = 1.5, y = max(decov_result$fraction), 
             label = paste0('p = ', format(pv, digits = 2)))
  plist[[state0]] = sig_plot
}

ggsave(gridExtra::grid.arrange(grobs = plist[c('ADRN-Baseline', 'ADRN-Proliferating', 'ADRN-Calcium',
                                               'ADRN-Dopaminergic', 'Interm-OxPhos', 'MES')], ncol = 6), 
       file = 'Figures/RNA/CIBERSORTx/ByMycn/cangelosi_fraction_malignantStates_byMYCN.pdf',
       device = 'pdf', width = 12, height = 2)
ggsave(gridExtra::grid.arrange(grobs = plist[c('THY1+Macro', 'HS3ST2+Macro',
                                               'F13A1+Macro', 'CCL4+Macro', 
                                               'IL18+Macro','VCAN+Macro',
                                               'C1QC+SPP1+Macro', 'ProliferatingMacro')], ncol = 8), 
       file = 'Figures/RNA/CIBERSORTx/ByMycn/cangelosi_fraction_macroStates_byMYCN.pdf',
       device = 'pdf', width = 12, height = 2)
