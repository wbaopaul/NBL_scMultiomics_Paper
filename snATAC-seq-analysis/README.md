# snATAC-seq analysis
Codes for pooling, integrating, annotating snATAC-seq data, and DAP, and TRN analyses

- combine_atac.R -- combine snATAC-seq matrices from all samples
- process_atac_non_binary.R -- initial process using standard seurat pipeline with minor midification
- atac_integration_allCells.R  -- integrate all cells by signac 
- transfer_label_allCells.R -- transfer rna cell type to cell in atac-seq data
- cell_type_ann_atac.R -- cell type annotation (by label transfer + cluster) 
- atac_int_malignant.R -- integrate malignant cells 
- cell_state_analysis_atac.R -- analysis of neoplastic cell states in atac-seq
- dap_analysis.R -- DAP analysis for different malignant cell states
