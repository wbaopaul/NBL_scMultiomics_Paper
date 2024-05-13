library(Signac)
library(Seurat)
library(data.table)
library(magrittr)
## pool all cells ####
dir0 = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/NB/integrated/reConstructed_matrix/'
dir1 = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/NB/output_chop/'
dir2 = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/NB/output_ucsd/'
samples = dir(dir0)
annotation = readRDS('/mnt/isilon/tan_lab/yuw1/R_work_dir/scGEL/data/annotation_hg38.rds')
seu_list = list()
for(sample0 in samples){
  mfile = paste0(dir0, sample0, '/matrix.rds')
  if(!file.exists(mfile)) next
  tmp = readRDS(mfile)
  #cs = colSums(tmp > 0)
  #tmp = tmp[, cs > 2000]
  rs = rowSums(tmp > 0)
  tmp = tmp[rs > 10, ]
  
  # add qc metric
  qc_stat_file = paste0(dir1, 'NB_7767_', sample0, '/summary/NB_7767_', sample0, '.MACS2.qc_per_barcode.txt')
  if(!file.exists(qc_stat_file)) qc_stat_file = paste0(dir2, 'NB_7767_', sample0, '/summary/NB_7767_', sample0, '.MACS2.qc_per_barcode.txt')
  
  if(file.exists(qc_stat_file)){
    qc_singlecell = fread(qc_stat_file)
    qc_singlecell = qc_singlecell[bc %in% colnames(tmp)]
    qc_singlecell = data.frame(qc_singlecell)
    rownames(qc_singlecell) = qc_singlecell$bc
    qc_singlecell$bc = NULL
    names(qc_singlecell) =  c("total.unique.frags", "frac.mito",  "frac.peak",
                              "frac.promoter", "frac.tss", "frac.enhancer", "tss_enrich_score")
    #seurat.obj <- AddMetaData(seurat.obj, metadata = qc_singlecell)
  }
  
  frag_file = paste0(dir1, 'NB_7767_', sample0, '/summary/NB_7767_', sample0, '.fragments.tsv.gz')
  if(!file.exists(frag_file)) frag_file = paste0(dir2, 'NB_7767_', sample0, '/summary/NB_7767_', sample0, '.fragments.tsv.gz')
  
  chrom_assay <- CreateChromatinAssay(
    counts = tmp,
    fragments = frag_file,
    annotation = annotation
  )
  # create a seurat obj for each sample
  seurat.obj = CreateSeuratObject(counts = chrom_assay, assay = 'ATAC', 
                                  meta.data = qc_singlecell)
  
  seurat.obj$biospecimenName = sample0
  seu_list[[sample0]] = seurat.obj
  message(paste(sample0, 'Done!'))
}

seurat.atac0 = merge(seu_list[[1]], seu_list[-1],
                     add.cell.ids = samples)

## add more metadata
sample.sheet = subset(seurat.atac0@meta.data, select = c(biospecimenName)) %>% unique()
sample.sheet = data.table(sample.sheet)
setkey(sample.sheet, biospecimenName)
sample.sheet[, 'sampleID' := unlist(strsplit(biospecimenName, '_'))[1], by = biospecimenName]
sample.sheet[, 'Region' := unlist(strsplit(biospecimenName, '_'))[2], by = biospecimenName]

seurat.atac0$sampleID = sample.sheet[seurat.atac0$biospecimenName]$sampleID
seurat.atac0$region = sample.sheet[seurat.atac0$biospecimenName]$Region

## add metadata from RNA-seq
mdata = readRDS('MetaData/snRNA/metadata_seurat_from_Yasin.rds')
mdata = data.table(subset(mdata, select = c(sample_name, Stage_Code, Patient_ID, Patient_No))) %>% unique()
mdata[, 'sampleID' := unlist(strsplit(sample_name, '_'))[3], by = sample_name]
mdata = subset(mdata, select = c(sampleID, Stage_Code, Patient_ID, Patient_No)) %>% unique()
setkey(mdata, sampleID)

dd = subset(sample.sheet, select = c('sampleID')) %>% unique()
setkey(dd, sampleID)
dd[, 'Stage_Code' := mdata[J(dd$sampleID)]$Stage_Code]
dd[, 'Patient_ID' := mdata[J(dd$sampleID)]$Patient_ID]
dd[, 'Patient_No' := mdata[J(dd$sampleID)]$Patient_No]

seurat.atac0$Stage_Code = dd[seurat.atac0$sampleID]$Stage_Code
seurat.atac0$Patient_No = dd[seurat.atac0$sampleID]$Patient_No


saveRDS(seurat.atac0, file = 'Seurat_Objects/seurat_atac_pooled.rds')

## prepare ACTIVITY assay ####
seurat.atac0 = readRDS(file = 'Seurat_Objects/seurat_atac_pooled.rds')
mtx.ga <- GeneActivity(seurat.atac0)
saveRDS(mtx.ga, file = 'Seurat_Objects/mtx_gene_activity_signac.rds')



