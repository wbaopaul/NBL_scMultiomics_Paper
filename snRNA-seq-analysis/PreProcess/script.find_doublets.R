args <- commandArgs(T)
cohort = args[1]
sample_name = args[2]
min_expressed_genes = as.numeric(args[3])
min_UMI = as.numeric(args[4])

num_loaded_cells = 0

library(ggplot2)
library(Seurat)

#sample_name = "NB_7767_3072_REG1"
#cohort = 'NB'
#cohort = 'HGG'
#cohort = 'Broad_NB'

seurat_dir = paste0('/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/scrna/seurat/',cohort, '/')

max_mit_pct = 10
#min_expressed_genes = 500
#min_expressed_genes = 1000
max_expressed_genes = 10000
#min_UMI = 1500
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

sc_dir = paste0(seurat_dir, '/filtered/individual.sc_transformed.', filtering_str, '/')
sc_file_path = paste0( sc_dir, "/seurat_object.",sample_name,".rds")

doublet_dir = paste0(seurat_dir, '/filtered/doublets_marked.', filtering_str, '/')
dir.create(doublet_dir)

singlet_dir = paste0(seurat_dir, '/filtered/singlets.', filtering_str, '/')
singlet_matrices_dir = paste0(seurat_dir, '/filtered/singlet_matrices.', filtering_str, '/')

dir.create(singlet_dir)
dir.create(singlet_matrices_dir)

singlet_file = paste0(singlet_dir, '/seurat_object.',sample_name,'.rds')

require(openxlsx)
num_loaded_cells_file = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/metadata/scrna/NB/Number_of_Cells_Loaded.xlsx'
df_num = read.xlsx(num_loaded_cells_file)
vec_num = df_num$Loaded_Cell_Count
names(vec_num) = df_num$Sample

num_loaded_cells = as.numeric(vec_num[sample_name])

if(is.null(num_loaded_cells))
{
    num_loaded_cells = 10000
}
if(is.na(num_loaded_cells))
{
    num_loaded_cells = 10000
}

print("----------------------------------------------------------------")
print(paste("Loaded cell count:", num_loaded_cells ))
print("----------------------------------------------------------------")

#sc_transformed_doublet_list <- foreach (i=1:num_objects) %dopar% {

  source('/mnt/isilon/tan_lab/uzuny/scripts/single_cell_matrix_processing/scDataAnalysis_Utilities.v10.R')
  source('/mnt/isilon/tan_lab/uzuny/scripts/10x/gene_expression_matrix/read_10x_matrix_functions.v04.R')
  #source('/mnt/isilon/tan_lab/uzuny/scripts/package_scripts/casper_functions.R')

  print('******************************************************')


  doublet_file = paste0(doublet_dir, '/seurat_object.',sample_name,'.rds')
  print(doublet_file)
    
  
  
  #if(!exists(doublet_file))
  #{
  
    
    seurat_object = readRDS( sc_file_path )
    
    
    
    library(Seurat)
    library(DoubletFinder)
    
    #DimPlot(seurat_object, reduction = 'umap')
    
    exp_doublet_rate = get_multiplet_rate(cells_loaded = num_loaded_cells)
    print('**********************************************************')
    print(num_loaded_cells)
    print(exp_doublet_rate)
    print('**********************************************************')
    
    seurat_object = FindDoublets(seurat_object, PCs = 1:30,
                             exp_rate = exp_doublet_rate, sct = TRUE)
    
    print(head(seurat_object@meta.data))
    print(paste('Finished sample ', sample_name) )
    print('******************************************************')
    
    #DimPlot(seurat_object, reduction = 'umap', group.by = 'Doublet_Singlet')
    

    saveRDS(seurat_object, file = doublet_file)
    #  seurat_object
    


seurat_object_singlet = subset(seurat_object, Doublet_Singlet == 'Singlet')
saveRDS(seurat_object, file = singlet_file)

print(singlet_file)

count_matrix_file = paste0(singlet_matrices_dir, '/count_matrix.',sample_name,'.rds')

count_matrix = seurat_object_singlet@assays$RNA@counts
saveRDS(count_matrix, count_matrix_file)


print(count_matrix_file)


  #}
 
