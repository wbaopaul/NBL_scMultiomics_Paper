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
library(ggalluvial)
`%notin%` = Negate(`%in%`)

myColors = brewer.pal(6, 'Set1')
names(myColors) = c('ADRN_Calcium_0', 'ADRN_Ribosome_1', 'Interm_OxPhos_2',
                    'ADRN_Dopaminergic_3', 'ADRN_Proliferating_4', 'MES_5')


## clean data ####
seurat.atac <- readRDS("Seurat_Objects/signac_obj_adrenergic_signac_modified_integration_v2.rds")
seurat.atac <- FindClusters(seurat.atac, resolution = 0.2)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.5)

seurat.atac = subset(seurat.atac, ATAC_snn_res.0.5 %in% c(0:11)) # rm tiny clusters and redo clustering
seurat.atac <- RunUMAP(seurat.atac, reduction = "integrated_lsi", 
                             dims = 2:30, n.neighbors = 50, min.dist = 0.1)
seurat.atac = FindNeighbors(seurat.atac, dims = 2:30, 
                                  reduction = 'integrated_lsi')
seurat.atac <- FindClusters(seurat.atac, resolution = 0.25)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.3)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.35)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.4)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.45)
seurat.atac <- FindClusters(seurat.atac, resolution = 0.5)

seurat.atac <- FindClusters(seurat.atac, resolution = 0.4)

pp <- DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.25', label = T) +
  DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.3', label = T) +
  DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.35', label = T) +
  DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.4', label = T) +
  DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.45', label = T) +
  DimPlot(seurat.atac, group.by = 'ATAC_snn_res.0.5', label = T) 

seurat.atac <- subset(seurat.atac, ATAC_snn_res.0.4 %in% c(0:10)) #cluster11 only two cells
saveRDS(seurat.atac, file = 'Seurat_Objects/signac_obj_adrenergic_signac_modified_integration_v3.rds')
ggsave(pp, filename = 'Figures/ATAC/State/umap_different_resls.pdf', device = 'pdf',
       width = 15, height = 6)

## seurat label transfer  ####
seurat.atac <- readRDS("Seurat_Objects/signac_obj_adrenergic_signac_modified_integration_v3.rds")
seurat.atac$ATAC_snn_res.0.4 = as.character(seurat.atac$ATAC_snn_res.0.4)
seurat.atac$seurat_clusters = seurat.atac$ATAC_snn_res.0.4

labelTransRes = readRDS('MetaData/labelTransferResult_cca_signac_obj_adrenergic_signac_modified_integration_v2.rds')
seurat.atac = AddMetaData(seurat.atac, metadata = labelTransRes)
seurat.atac$cell_state = seurat.atac$cell_state0 = seurat.atac$cell_state1 = NULL
saveRDS(seurat.atac, file = "Seurat_Objects/seurat_atac_malignant_withTransferLabels.rds")

pc <- DimPlot(seurat.atac, group.by = 'predicted.id') +
  scale_color_manual(values = myColors) 
ps <- FeaturePlot(seurat.atac, features = 'prediction.score.max',
                 max.cutoff = 'q99', raster = T) 

ggsave(pc + ps, filename = 'Figures/ATAC/State/umap_predicted_cell_state_byLabelTransfer.pdf',
       device = 'pdf', width = 14, height = 5)


tmp = subset(seurat.atac, prediction.score.max > 0.6)
freq <- table(tmp$predicted.id, tmp$ATAC_snn_res.0.4) %*% diag(1/table(tmp$ATAC_snn_res.0.4))
colnames(freq) = names(table(tmp$ATAC_snn_res.0.4))


## plot fraction
p0 <- pheatmap(freq, angle_col = 45, border_color = NA, 
               color = hcl.colors(12, "YlOrRd", rev = TRUE))
## plot enrichment
mat = table(tmp$predicted.id, tmp$ATAC_snn_res.0.4)
norm_mat = pearson_residual(mat)
pvs = pnorm(norm_mat, lower.tail = F)
plot_mat <- -log10(pvs)
plot_mat[is.na(plot_mat)] = 0
plot_mat[is.infinite(plot_mat)] = max(plot_mat[!is.infinite(plot_mat)]) + 1
plot_mat[plot_mat > quantile(plot_mat, 0.975)] = quantile(plot_mat, 0.975)

plot_mat_norm <- plot_mat %*% diag(1/colSums(plot_mat + 0.01)) 
colnames(plot_mat_norm) = colnames(plot_mat)

p1 <- pheatmap(plot_mat, angle_col = 45, border_color = NA, 
               color = hcl.colors(12, "YlOrRd", rev = TRUE))
p2 <- pheatmap(plot_mat_norm, angle_col = 45, border_color = NA, 
               color = hcl.colors(12, "YlOrRd", rev = TRUE))


mdata = data.table(subset(tmp@meta.data, select = c('seurat_clusters', 'predicted.id')))
mdata = mdata[, 'N':= .N, by = seurat_clusters]
mdata = mdata[, 'freq' := .N, by = list(seurat_clusters, predicted.id)] %>% unique()
mdata = mdata[, 'frac' := freq/N]

mdata$seurat_clusters = factor(mdata$seurat_clusters, levels = 0:10)
pal <- ggplot(data = mdata,
              aes(axis1 = seurat_clusters, axis2 = predicted.id, y = freq)) +
  geom_alluvium(aes(fill = predicted.id)) +
  geom_stratum(aes(fill = predicted.id)) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c('seurat_clusters', 'predicted.id'),
                   expand = c(0.15, 0.05), 
                   labels = c( "ATAC clusters", "RNA States")) +
  scale_fill_manual(values = myColors, na.value = 'grey90') + 
  theme_minimal() + theme(legend.position = 'none',
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())

#clean version
mdata1 = mdata[frac > 0.2]
#normalize to 1 for sum of frac per cluster
mdata1 = mdata1[, 'frac0' := sum(frac), by = seurat_clusters]
mdata1[, 'frac' := frac/frac0]
pal1 <- ggplot(data = mdata1,
              aes(axis1 = seurat_clusters, axis2 = predicted.id, y = frac)) +
  geom_alluvium(aes(fill = predicted.id)) +
  geom_stratum(aes(fill = predicted.id)) +
  geom_text(stat = "stratum", 
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c('seurat_clusters', 'predicted.id'),
                   expand = c(0.15, 0.05), 
                   labels = c( "ATAC clusters", "RNA States")) +
  scale_fill_manual(values = myColors, na.value = 'grey90') + 
  theme_minimal() + theme(legend.position = 'none',
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())


## differential tf activity ####
if(F){
  chromvar.obj = readRDS('chromVAR_Objects/seurat_atac_non_binary_vap20000.rds') # old version include all cells
  zscore = chromvar.obj@assays@data$z
  dscore = chromvar.obj@assays@data$deviations
  zscore = zscore[, colnames(seurat.atac)]
  dscore = dscore[, colnames(seurat.atac)]
  rownames(dscore) = rownames(zscore) = chromvar.obj@elementMetadata$name
}

#by signac
if(T){
  chromvar.obj <- readRDS('Seurat_Objects/seurat_allCells_chromvar_result_bySignac.rds')
  chromvar.obj <- readRDS('Seurat_Objects/seurat_malignant_chromvar_result_bySignac.rds') #with new peaks
  zscore = chromvar.obj@assays$chromvar@data
  tfs = readRDS('MetaData/snATAC/TF/motif_id_motif_name.rds')
  names(tfs) = gsub('_', '-', names(tfs))
  all(rownames(zscore) == names(tfs))
  names(tfs) = NULL
  rownames(zscore) = tfs
  zscore = zscore[, colnames(seurat.atac)]
  dscore = zscore
}

### < use seurat diff tf activity ####
seurat.atac[["TF"]] <- CreateAssayObject(counts = dscore)
seurat.atac@assays$TF@data = dscore
seurat.atac@assays$TF@scale.data = zscore
DefaultAssay(seurat.atac) = 'TF'
diff.tf = FindAllMarkers(seurat.atac, assay = 'TF', logfc.threshold = 0.05,
                         min.pct = 0.025, only.pos = T, max.cells.per.ident = 500,
                         mean.fxn = rowMeans,
                         fc.name = "avg_diff")
diff.tf = data.table(diff.tf)
saveRDS(diff.tf, file = 'data/intermediate/ATAC/seurat_diff_chromVAR_dscore_by_atac_clusters.rds')

## differential TF activity of predicted cell states ####
seurat.atac$predicted_cell_state = seurat.atac$predicted.id
seurat.atac$predicted_cell_state[seurat.atac$prediction.score.max < 0.6] = 'Unknown'
DefaultAssay(seurat.atac) = 'TF'
Idents(seurat.atac) = seurat.atac$predicted_cell_state
diff.tf0 = FindAllMarkers(seurat.atac, assay = 'TF', logfc.threshold = -Inf,
                         min.pct = 0, only.pos = F, max.cells.per.ident = 500,
                         mean.fxn = rowMeans, return.thresh = 1,
                         fc.name = "avg_diff") ##keep all TFs for calculating TF avg_diff
diff.tf0 = data.table(diff.tf0)
diff.tf0$cluster = as.character(diff.tf0$cluster)
saveRDS(diff.tf0, file = 'data/intermediate/ATAC/seurat_diff_chromVAR_dscore_by_atac_cell_states_allTFsKept.rds')

## filter redundant ones
diff.tf = diff.tf0[avg_diff > 0.5 & pct.1 > 0.1 & p_val < 0.01] # may need further filtering in downstream analysis
saveRDS(diff.tf, file = 'data/intermediate/ATAC/seurat_diff_chromVAR_dscore_by_atac_cell_states.rds')

## change btw timepoint: cell state level ####
seurat.atac$predicted_cell_state = seurat.atac$predicted.id
seurat.atac$predicted_cell_state[seurat.atac$prediction.score.max < 0.6] = 'Unknown'

## global composition per timepoint
mdata = data.table(subset(seurat.atac@meta.data, select = c('predicted_cell_state', 'Stage_Code')))
mdata = mdata[Stage_Code %in% c('DX', 'PTX') & predicted_cell_state != 'Unknown']
mdata[, 'n' := .N, by = list(predicted_cell_state, Stage_Code)]
mdata[, 'N' := .N, by = Stage_Code]
mdata[, 'percentage' := n/N]
mdata = subset(mdata, select = c( 'predicted_cell_state', 'Stage_Code', 'percentage')) %>% unique()

p2 <- ggplot(data = mdata, aes(x = Stage_Code, y = percentage, fill = predicted_cell_state)) +
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  scale_fill_manual(values = myColors) +
  xlab('') + ylab('Fraction') + guides(fill=guide_legend(title=""))
ggsave(p2, filename = 'Figures/ATAC/State/atac_malignant_cell_state_compositin_ByTime.pdf', device = 'pdf',
       width = 5, height = 5)



## change overtime per state
mdata0 = data.table(seurat.atac@meta.data)
mdata = subset(mdata0, select = c('predicted_cell_state', 'Stage_Code'))
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]

mdata[, 'n' := .N, by = list(predicted_cell_state, Stage_Code)]
mdata[, 'N' := .N, by = predicted_cell_state]
mdata[, 'percentage' := n/N]
mdata = subset(mdata, select = c( 'predicted_cell_state', 'Stage_Code', 'percentage')) %>% unique()
mdata = mdata[predicted_cell_state != 'Unknown']
mdata$predicted_cell_state = factor(mdata$predicted_cell_state, 
                                    levels = c('ADRN_Calcium_0', 'ADRN_Ribosome_1',
                                               'Interm_OxPhos_2', 'ADRN_Dopaminergic_3',
                                               'ADRN_Proliferating_4', 'MES_5') )
p3 <- ggplot(data = mdata, aes(x = predicted_cell_state, y = percentage, 
                               fill = Stage_Code)) +
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  scale_fill_brewer(palette = 'Set1', direction = -1) +
  xlab('') + ylab('Fraction') + guides(fill=guide_legend(title="")) +
  geom_hline(yintercept = nrow(mdata0[Stage_Code == 'PTX'])/nrow(mdata0),
             linetype = 2)
ggsave(p3, filename = 'Figures/ATAC/State/atac_malignant_cell_cluster_changeByTime.pdf', device = 'pdf',
       width = 6, height = 5)

## change btw timepoint: patient level comparison ####
mdata = seurat.atac@meta.data
mdata = data.table(mdata)

## remove patients that with no DX sample sequenced
DX_pts = unique(mdata[Stage_Code == 'DX']$Patient_No)
mdata = mdata[Patient_No %in% DX_pts]
mdata = mdata %>% subset(select = c(Patient_No, Stage_Code, predicted_cell_state))
mdata = mdata[Stage_Code %in% c('DX', 'PTX')]

mdata = mdata[predicted_cell_state != 'Unknown']
mdata[, 'N' := .N, by = list(Patient_No, Stage_Code)]

mdata1 = mdata[Stage_Code == 'DX']
mdata2 = mdata[Stage_Code == 'PTX']
mdata1[, 'n' := .N, by = list(Patient_No, predicted_cell_state)]
mdata2[, 'n' := .N, by = list(Patient_No, predicted_cell_state)]

mdata = rbind(mdata1, mdata2)
mdata[, 'frac' := n/N]
mdata = subset(mdata, select = c(Stage_Code, predicted_cell_state, frac, n, N, Patient_No)) %>% unique()
mdata = mdata[N > 100]

plist = list()
for(cl0 in c('ADRN_Calcium_0', 'ADRN_Ribosome_1',
             'Interm_OxPhos_2', 'ADRN_Dopaminergic_3',
             'ADRN_Proliferating_4', 'MES_5')){
  mdata0 = mdata[predicted_cell_state == cl0]
  setkey(mdata0, Patient_No)
  mdata0[, nSample := .N, by = Patient_No]
  mdata0 = subset(mdata0, select = c(Stage_Code, frac, Patient_No))
  mdata01 = mdata0[ Stage_Code == 'DX']
  mdata02 = mdata0[ Stage_Code == 'PTX']
  
  pts = union(mdata01$Patient_No, mdata02$Patient_No)
  y = mdata02[pts]$frac
  x = mdata01[pts]$frac
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  pv_cl = wilcox.test(y, x, paired = T)$p.value/2
  
  pdata = data.table('Patient_No' = c(pts, pts), 'frac' = c(x, y), 
                     'Stage_Code' = rep(c('DX', 'PTX'), each = length(pts)))
  pp <- ggplot(data = pdata, aes(x = Stage_Code, y = frac, colour = Stage_Code)) +
    geom_boxplot() + ggtitle(paste0(cl0)) + scale_color_brewer(palette = 'Set1', direction = -1) +
    geom_point(aes(fill = Stage_Code), size = 2) + theme_classic() + NoLegend() +
    xlab('') + ylab('Fraction') + scale_fill_brewer(palette = 'Set1', direction = -1) +
    geom_line(aes(group = Patient_No), color='gray', alpha=0.5, size = 1) +
    annotate('text', x = 1.5, y = max(pdata$frac), label = paste0('p = ', format(pv_cl, digits = 2)))
  cl1 = gsub(' ', '_', cl0)
  cl1 = gsub('+', '_', cl1, fixed = T)
  ggsave(pp, filename = paste0('Figures/ATAC/State/atac_cluster', cl1, '_frac_byTimepoint.pdf'), device = 'pdf',
         width = 3, height = 4)
  plist[[cl0]] = pp
}

## diff tf btw predicted cell state ####
diff.tf = readRDS('data/intermediate/ATAC/seurat_diff_chromVAR_dscore_by_atac_cell_states.rds')
diff.tf = diff.tf[!grepl(gene, pattern = '^ENSG')]
diff.tf = diff.tf[order(-avg_diff)]

## add TFs in top degs
tfs.deg = readRDS('MetaData/tfs_deg_state.rds')
tfs_sub_trn = readRDS('MetaData/tfs_in_subnetworks.rds')
sele.tfs = sig.tfs = NULL
for(cl0 in c('ADRN_Calcium_0', 'MES_5', 'ADRN_Proliferating_4',
             'Interm_OxPhos_2', 'ADRN_Dopaminergic_3', 'ADRN_Ribosome_1')){
  diff.tf0 = diff.tf[cluster == cl0]
  tfs.deg0 = tfs.deg[cell_state == cl0]
  sele.tfs = unique(c(sele.tfs, diff.tf0[1:min(20, nrow(diff.tf0))]$gene))
  sig.tfs = union(sig.tfs, intersect(diff.tf0$gene[1:min(50, nrow(diff.tf0))], 
                      tfs.deg0$gene[1:min(50, nrow(tfs.deg0))]))
  
}
sele.tfs = union(sele.tfs, show_tfs)


tmp = subset(seurat.atac, downsample = 200)

tmp_mat = seurat.atac@assays$TF@data[sele.tfs,]
tmp_mat[tmp_mat > 6] = 6
tmp_mat[tmp_mat < -2] = -2
pmat = sapply(c('ADRN_Calcium_0', 'MES_5', 'ADRN_Proliferating_4',
                'Interm_OxPhos_2', 'ADRN_Dopaminergic_3', 'ADRN_Ribosome_1'), 
              function(x) return(rowMeans(tmp_mat[, seurat.atac$predicted_cell_state == x])))
colnames(pmat) = c('ADRN_Calcium_0', 'MES_5', 'ADRN_Proliferating_4',
                   'Interm_OxPhos_2', 'ADRN_Dopaminergic_3', 'ADRN_Ribosome_1')
pmat[pmat > 2] = 2
pmat[pmat < -2] = -2

getPalette = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))
tmpColors = getPalette(39)

pb <- pheatmap::pheatmap(pmat, cluster_rows = T, border_color = NA,
                   cluster_cols = T, show_colnames = T, angle_col = 45,
                   breaks = sort(unique(c(seq(-2, 0, length.out = 20), 
                                          seq(0, 2, length.out = 20)))),
                   color = tmpColors)

tmp$predicted_cell_state = factor(tmp$predicted_cell_state, levels = c('ADRN_Calcium_0', 'ADRN_Ribosome_1',
                                                                       'Interm_OxPhos_2', 'ADRN_Dopaminergic_3',
                                                                       'ADRN_Proliferating_4', 'MES_5'))
ph <- DoHeatmap(tmp, 
                features = sele.tfs, group.colors = tmpColors,
                disp.max = 2, disp.min = -2) +  
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[-c(1, 11)])

ggsave(ph, filename = 'Figures/ATAC/Signatures/TF/heatmap_diff_tf_byCellState.pdf',
       width = 12, height = 12, device = 'pdf')
ggsave(pb, filename = 'Figures/ATAC/Signatures/TF/heatmap_diff_tf_byCellState_bulk.pdf',
       width = 7, height = 10, device = 'pdf')

