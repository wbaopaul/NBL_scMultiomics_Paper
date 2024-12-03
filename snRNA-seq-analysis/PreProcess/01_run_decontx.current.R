library('openxlsx')

args <- commandArgs(T)
cohort = args[1]
sample_name = args[2]
min_expressed_genes = as.numeric(args[3])
min_UMI = as.numeric(args[4])

#num_loaded_cells = 0
#num_loaded_cells = as.numeric(args[4])



#sample_name = "NB_7767_985_REG2"
#cohort = 'NB'


print(100)

seurat_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/scrna/seurat/',cohort, '/')

max_mit_pct = 10
#min_expressed_genes = 1000
max_expressed_genes = 10000
#min_UMI = 2000
max_UMI = 50000

filtering_str = paste0('filtered.' 
                       , 'max_mit_pct_', max_mit_pct
                       ,'.min_expressed_genes_', min_expressed_genes
                       ,'.max_expressed_genes_', max_expressed_genes
                       ,'.min_UMI_', min_UMI
                       ,'.max_UMI_', max_UMI)



plot_dir =  paste0('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/plots/scrna/',cohort, '/individual.sc_transformed.',
                   filtering_str,'/filtered/', sample_name, '/')
dir.create(plot_dir, recursive = T, showWarnings = F)

aln_summary_file = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/summary/scrna/htan_std/Merged_Alignment_Summary.2021_10_03.xlsx'
df_aln_summary = read.xlsx(aln_summary_file, 'Numeric_Tab')[1:98,]
frac_in_cells = df_aln_summary$Fraction.Reads.in.Cells
names(frac_in_cells) = paste0('NB_7767_', df_aln_summary$Sample)


singlet_matrices_dir = paste0(seurat_dir, '/filtered/singlet_matrices.', filtering_str, '/')
count_matrix_file = paste0(singlet_matrices_dir, '/count_matrix.',sample_name,'.rds')


#doublet_dir = paste0(seurat_dir, '/filtered/doublets_marked.', filtering_str, '/')
#decont_matrix_dir = paste0(seurat_dir, '/filtered/decont_matrices.', filtering_str, '/')
decont_matrix_dir = paste0(seurat_dir, '/filtered/decont_matrices.dynamic.', filtering_str, '/')
dir.create(decont_matrix_dir, recursive = T)
decont_matrix_file = paste0(decont_matrix_dir, '/decont_matrix.',sample_name,'.rds')


print('***********************************')
if(file.exists(decont_matrix_file))
{

  print('Decont already computed')
  print(decont_matrix_file)
  
}else{
  
  
  library(SingleCellExperiment)
  library(celda)
  
  print('Starting')
  
  frac_in_cells.sample = frac_in_cells[sample_name]
  print(count_matrix_file)
  count_matrix = readRDS(count_matrix_file)
  
  mean_count = mean(count_matrix[count_matrix>0])
  
  dim(count_matrix)
  count_matrix[1:5, 1:5]
  
  print('1')
  
  
  #sce_object = SingleCellExperiment(    assays = list(counts = count_matrix) )
  print('4')
  
  temp = as.matrix(count_matrix)
  temp[1:10, 1:10]
  storage.mode(temp) <- "integer"
  native_c = ceiling(20 * frac_in_cells.sample)
  contam_c = 20 - native_c
  if(native_c == 0 ){ native_c = 1 }
  delta_vec = as.numeric(c(native_c, contam_c))
  decontx_result = decontX(temp, delta = delta_vec, estimateDelta = F)
  
  #decontx_result = decontX(as.matrix(count_matrix))
  #sce_object = SingleCellExperiment(    assays = list(counts = temp) )
  
  #decontx_result = decontX(sce_object)
  print('5')
  
  #decont_matrix = decontx_result@assays@data$decontXcounts
  
  #decont_matrix = decontx_result$resList$estNativeCounts
  #decont_matrix.1 = decont_matrix
  decont_matrix = decontx_result$decontXcounts
  contam = decontx_result$contamination
  

  length(contam)
  
  print('10')
  
  head(decont_matrix[1:5, 1:5])

  #count_matrix['ALK', 1:10]
  
  #decont_matrix['ALK', 1:10]
  
  sum(decont_matrix - count_matrix)
  
  saveRDS(decont_matrix, decont_matrix_file)
  print(decont_matrix_file)
  
  cell_contamination_file = paste0(decont_matrix_dir, '/cell_contamination.',sample_name,'.rds')
  saveRDS(contam, cell_contamination_file)
}
# 
# 
# prominent_marker_genes = c('NRG1', 'ALK', 'PHOX2B', 'MYCN', 'THY1' , 
#                            'COL1A1', 'COL3A1', 'COL6A1', 'PDGFRB', 'DCN', 'TAGLN',
#                            'PECAM1', 'EPAS1', 'KDR', 'CDH5', 'ENG', 'CLDN5',
#                            'CD4', 'CD14', 'CD63', 'MRC1', 
#                            'MS4A1', 'CD22', 'PAX5', 
#                            'CD2', 'CD96', 'CD247', 'IL7R', 'LEF1', 'TCF7',
#                            'CYP11B1', 'CYP17A1',
#                            'SOX10', 'S100B')
# 
# prominent_marker_genes_short = c( 'ALK', 'PHOX2B',  
#                                   'COL1A1', 'COL3A1', 
#                                   'PECAM1',  'KDR', 'CDH5', 
#                                   'CD4', 'CD14', 
#                                   'MS4A1', 'CD22', 'PAX5', 
#                                   'CD2', 'CD96', 'CD247', 
#                                   'CYP11B1', 'CYP17A1',
#                                   'SOX10', 'S100B', 
#                                   'PKHD1', 'GLIS3')
# 
# decont_matrix.1 = decont_matrix
# decont_matrix.2 = decont_matrix
# 
# 
# seurat_object = CreateSeuratObject(decont_matrix)
# seurat_object = NormalizeData(seurat_object)
# seurat_object = FindVariableFeatures(seurat_object)
# seurat_object = ScaleData(seurat_object)
# seurat_object = RunPCA(seurat_object)
# seurat_object = RunUMAP(seurat_object, dims = 1:10)
# seurat_object = FindNeighbors(seurat_object, reduction = 'umap', dims = 1:2)
# seurat_object = FindClusters(seurat_object, resolution = 0.1)
# 
# 
# #seurat_object.1 = seurat_object
# seurat_object.2 = seurat_object
# #
# #DimPlot(seurat_object)
# DimPlot(seurat_object.2)
# 
# gg = CellChat::StackedVlnPlot(obj = seurat_object.1, features = prominent_marker_genes_short, #y.size = 3,
#                     plot.margin = unit(c(-5, -2, -5, -5), "cm"),
#                     group.by = 'seurat_clusters', assay = 'RNA')
# 
# print(gg)
#   
# 
# gg = CellChat::StackedVlnPlot(obj = seurat_object.1, features = prominent_marker_genes, #y.size = 3,
#                               plot.margin = unit(c(-5, -2, -5, -5), "cm"),
#                               group.by = 'seurat_clusters', assay = 'RNA')
# 
# print(gg)
