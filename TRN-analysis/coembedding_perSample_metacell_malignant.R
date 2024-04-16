
library(gridExtra)
library(patchwork)
`%notin%` = Negate(`%in%`)
# do normalization using log, tf-idf, or none, regress out confounds on pca or not 
doBasicSeurat_atac_updated <- function(mtx, npc = 30, top.variable = 5000, 
                                       norm_by = 'tf-idf',
                                       doScale = T, doCenter = T, assay = 'ATAC',
                                       reg.var = 'nCount_ATAC', regressOnPca = T,
                                       project = 'scATAC', 
                                       meta.data = NULL, excludePks.fromVAP = NULL){
  
  # top.variabl -- use top most variable features
  seurat.obj = CreateSeuratObject(mtx, project = project, assay = assay,
                                  names.delim = '-',
                                  meta.data = meta.data)
  
  if(norm_by == 'log') seurat.obj@assays[[assay]]@data <- log1p(mtx) / log(2)
  if(norm_by == 'tf-idf') seurat.obj@assays[[assay]]@data <- TF_IDF(mtx, verbose = F)
  if(norm_by == 'logNormalize') seurat.obj <- NormalizeData(
    object = seurat.obj,
    assay = assay,
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.obj[[paste0('nCount_', assay)]][, 1])
  ) ## treat the data as non-binary
  
  nvap = ifelse(top.variable > 1, top.variable, floor(top.variable * ncol(mtx)))
  seurat.obj <- FindVariableFeatures(object = seurat.obj,
                                     selection.method = 'vst',
                                     nfeatures = nvap)
  vaps = VariableFeatures(seurat.obj)
  
  #peak.frac = rowMeans(mtx > 0)
  #excludePks.fromVAP = names(which(peak.frac < vap.min.frac))
  
  if(!is.null(excludePks.fromVAP)) vaps = setdiff(vaps, excludePks.fromVAP)
  
  if(length(vaps) < 10) stop('Top few VAPs left!')
  ## redo normalization using vap
  if(norm_by == 'tf-idf'){
    mtx.norm = TF_IDF(mtx[vaps, ])
    tmp <- mtx[setdiff(rownames(mtx), vaps), ]
    data0 <- rbind(mtx.norm, tmp)
    seurat.obj[[assay]]@data = data0[rownames(mtx), ]
    rm(data0, tmp, mtx.norm)
  }
  
  
  if(regressOnPca){
    reg.var0 = NULL
  }else{
    reg.var0 = reg.var
  }
  VariableFeatures(seurat.obj) <- vaps
  seurat.obj <- ScaleData(object = seurat.obj,
                          features = vaps,
                          vars.to.regress = reg.var0, do.scale = doScale,
                          do.center = doCenter)
  
  seurat.obj <- RunPCA(object = seurat.obj,
                       features = vaps,
                       verbose = FALSE, npc = npc)
  if(length(reg.var) > 0 & regressOnPca) seurat.obj = regress_on_pca(seurat.obj, reg.var)
  
  return(seurat.obj)
}


args = commandArgs(T)
sampleID0 = args[1]

seuratAtacPath = 'Seurat_Objects/signac_obj_adrenergic_signac_modified_integration_v3.rds'
reduction2use = 'pca' # pca or harmony

## load scRNA data and cell type annotation ####
seurat.rna <- readRDS('Seurat_Objects/seurat_rna_malignant_harmony_clean.rds')
mdata <- readRDS('MetaData/metadata_seurat_rna_malignant_harmony_cellStateAnnotated.rds')
seurat.rna <- AddMetaData(seurat.rna, metadata = mdata)

seurat.rna = subset(seurat.rna, biospecimen_id == sampleID0 )

seurat.atac = readRDS(seuratAtacPath)
seurat.atac = subset(seurat.atac, sampleID == sampleID0)

## recreate seurat object using the sample only
mtx = seurat.atac@assays$ATAC@counts
npc = 30
rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.005))

nvap = 10000
inherited.mdata <- seurat.atac@meta.data
seurat.atac = doBasicSeurat_atac_updated(mtx, npc = npc, 
                                         norm_by = 'logNormalize',
                                         top.variable = nvap,
                                         regressOnPca = TRUE,
                                         reg.var = 'nFeature_ATAC', 
                                         excludePks.fromVAP = filtered.pks, 
                                         meta.data =inherited.mdata)
seurat.atac = RunUMAP(seurat.atac, reduction = 'pca', dim = 1:npc) %>% 
  FindNeighbors(reduction = 'pca', dim = 1:npc) %>% FindClusters(resolution = 0.2)
p0 = DimPlot(seurat.atac, raster = F) + ggtitle('ATAC')

seurat.atac$cell_state = NULL
## use GAS = promote + gene body accessibility for label transfer ####
if(all(names(seurat.atac@assays) != 'ACTIVITY')){
  activity.matrix = readRDS('Seurat_Objects/mtx_gene_activity_signac_new.rds')
  activity.matrix = activity.matrix[, colnames(seurat.atac)]
  seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
  
  seurat.atac <- NormalizeData(
    object = seurat.atac,
    assay = 'ACTIVITY',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.atac$nCount_ACTIVITY)
  )
  
}
DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

seurat.atac <- FindVariableFeatures(seurat.atac, nfeatures = 5000)
vags = VariableFeatures(seurat.atac)

refAssay = 'RNA'

## restrict RNA to the atac region if atac has 1 region but rna has 2
if(length(unique(seurat.rna$sample_name)) > length(unique(seurat.atac$region))){
  reg.atac = unique(seurat.atac$region)
  seurat.rna = subset(seurat.rna, sample_name == paste('NB_7767', sampleID0, 
                                                       reg.atac, sep = '_'))
}

## transfer label 
DefaultAssay(seurat.rna) = refAssay
mtx = seurat.rna[['RNA']]$counts
options(Seurat.object.assay.version = "v3")
seurat.rna <- CreateSeuratObject(mtx, meta.data = seurat.rna@meta.data)
seurat.rna = NormalizeData(seurat.rna)
seurat.rna = FindVariableFeatures(seurat.rna, nfeatures = 2200)
vegs = VariableFeatures(seurat.rna)
gfreq = rowMeans(seurat.rna[['RNA']]$counts > 0)
rgenes = names(which(gfreq < 0.0025))
vegs = setdiff(vegs, rgenes)
VariableFeatures(seurat.rna) <- vegs
seurat.rna = ScaleData(seurat.rna) %>% RunPCA(verbose = F, npcs = 30) %>% RunUMAP(reduction = 'pca', dim = 1:30)
seurat.rna = FindNeighbors(seurat.rna, reduction = 'pca', dims = 1:30) %>% FindClusters(resolution = 0.5)
p1 <- DimPlot(seurat.rna, group.by = 'seurat_clusters') + ggtitle('RNA')


genes4anchors = vegs

if(length(unique(seurat.rna$sample_name)) > length(unique(seurat.atac$region))) genes4anchors = intersect(vegs, vags)


transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = refAssay,
                                        normalization.method = 'LogNormalize',
                                        query.assay = "ACTIVITY",
                                        reduction = 'cca',
                                        k.anchor = 30, k.filter = 200)


## co-embedding ####
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <-seurat.rna[['RNA']]$data[genes4anchors %in% rownames(seurat.rna), ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
kweight = min(50, nrow(transfer.anchors@anchors) - 5)
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.atac[["pca"]],
                           dims = 1:ncol(seurat.atac[["pca"]]), k.weight = kweight)

# this line adds the imputed data matrix to the seurat.atac object
seurat.atac[[refAssay]] <- imputation
DefaultAssay(seurat.atac) = refAssay
#seurat.atac[['ATAC']] <- NULL

coembed <- merge(x = seurat.rna, y = seurat.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE, npcs = 30)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- FindNeighbors(coembed, dims = 1:30, reduction = 'pca')
coembed <- FindClusters(coembed, resolution = 0.1)

p2 <- DimPlot(coembed, group.by = 'tech', label = T) + ggtitle('Coembedded')
p3 <- DimPlot(coembed, group.by = 'seurat_clusters', label = T) + ggtitle('Coembedded')
p4 <- DimPlot(coembed, group.by = 'cell_state') + ggtitle('Coembedded')

message('Co-embedding done, working on metacell calling ... ')


## call metacell on coembedded obj ####
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

coembed <- RenameCells(coembed,  new.names = paste0(coembed$tech, '_', colnames(coembed)))

coembed <- SetupForWGCNA(
  coembed,
  gene_select = "fraction", # the gene selection approach
  #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "test", # the name of the hdWGCNA experiment
  features = genes4anchors
)

# construct metacells  in each group
coembed <- MetacellsByGroups(
  seurat_obj = coembed,
  group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 3, # maximum number of shared cells between two metacells
  min_cells = 100,
  reduction = 'pca',
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
coembed <- NormalizeMetacells(coembed)
metacell_obj <- GetMetacellObject(coembed)

if(ncol(metacell_obj) < 15){
  coembed <- MetacellsByGroups(
    seurat_obj = coembed,
    group.by = c("seurat_clusters"), # specify the columns in seurat_obj@meta.data to group by
    k = 20, # nearest-neighbors parameter
    max_shared = 3, # maximum number of shared cells between two metacells
    min_cells = 100, # ignore cluster with less than min_cells
    reduction = 'pca',
    ident.group = 'seurat_clusters'
  )
  # normalize metacell expression matrix:
  coembed <- NormalizeMetacells(coembed)
  metacell_obj <- GetMetacellObject(coembed)
  
}

coembed <- ScaleMetacells(coembed, features=VariableFeatures(seurat.rna))
message(paste('#of metacells:', ncol(metacell_obj)))

if(ncol(metacell_obj) > 50){
  coembed <- RunPCAMetacells(coembed, features=VariableFeatures(seurat.rna), npcs = 10)
  #coembed <- RunHarmonyMetacells(coembed, group.by.vars='Sample')
  coembed <- RunUMAPMetacells(coembed, reduction='pca', dims=1:10)
  
  p5 <- DimPlotMetacells(coembed, group.by = 'seurat_clusters') + ggtitle("cluster_metacell")
  p6 <- DimPlot(coembed, group.by = 'cell_state') + ggtitle('Coembedded')
  ggsave(p0 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3), 
         width = 15, height = 12, 
         filename = paste0('data/intermediate/coembedBySample/figures/new_coembed_', sampleID0, '.pdf'), 
         device = 'pdf')
  
}else{
  ggsave(p0 + p2 + p3 + p4 + plot_layout(ncol = 2, heights = c(2, 2)), 
         width = 15, height = 12, 
         filename = paste0('data/intermediate/coembedBySample/figures/new_coembed_', sampleID0, '.pdf'), 
         device = 'pdf')
}




message('Call metacells Done! Saving metadata matrices...')

## construct and save metacell rna and atac matrices 
## construct and save metacell rna and atac matrices 
merged_cells = metacell_obj$cells_merged
nrna_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'RNA_'))
natac_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'ATAC_'))

write(merged_cells, file = paste0('data/intermediate/coembedBySample/new_metacells/metacells_names_', sampleID0, '.txt'))


sele.metacells = merged_cells[nrna_cells >= 5 & natac_cells >= 5] ## ignore metacells which are mostly comprised of single modality
nrna_cells = nrna_cells[names(sele.metacells)]
natac_cells = natac_cells[names(sele.metacells)]


## manually
mask <- sapply(names(nrna_cells), function(x) colnames(coembed) %in%
                 unlist(strsplit(merged_cells[x], ',', fixed = T)))
rownames(mask) = colnames(coembed)
rna.mask <- sapply(names(nrna_cells), function(x) mask[, x]/nrna_cells[x])
rna.mtx <- coembed[['RNA']]$data %*% rna.mask # remember atac cells in RNA assay have 0 account

atac.mask <- sapply(names(natac_cells), function(x) mask[, x]/natac_cells[x])
atac.mtx <- coembed[['ATAC']]$data %*% atac.mask[colnames(coembed[['ATAC']]), ]



colnames(rna.mtx) = colnames(atac.mtx) = paste0(sampleID0, '_', names(sele.metacells))

saveRDS(rna.mtx, file = paste0('data/intermediate/coembedBySample/new_metacells/rna_metacell_mtx_', sampleID0, '.rds'))
saveRDS(atac.mtx, file = paste0('data/intermediate/coembedBySample/new_metacells/atac_metacell_mtx_', sampleID0, '.rds'))
message(paste('#of final metacells:', ncol(rna.mtx)))

message('All Done!')
