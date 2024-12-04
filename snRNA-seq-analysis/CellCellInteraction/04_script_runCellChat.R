library(Seurat)
library(data.table)
library(CellChat)
library(magrittr)

## prepare input: just need to run once
seurat.rna = readRDS(file = paste0('data/intermediate/liana/seurat_tumor_macro4liana.rds')) # same input as liana
labels <- Idents(seurat.rna)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

## Run cellchat ####
cellChat <- createCellChat(object = seurat.rna, group.by = "cell_state", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database

future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

saveRDS(cellchat, file = 'data/intermediate/cellchat/cellchat_res.rds')

## post analysis --plot ####
cellchat = readRDS(file = 'data/intermediate/cellchat/cellchat_res.rds')

df.net <- subsetCommunication(cellchat, sources.use = 'VCAN+',
                              targets.use = c("ADRN-Baseline", "ADRN-Calcium", 'ADRN-Dopaminergic',
                                              'ADRN-Proliferating', 'Interm-OxPhs', 'MES'))
df.net = data.table(df.net)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netAnalysis_contribution(cellchat, signaling = 'EGF')

## visualize HBEGF-ERBB4
pathways.show = 'EGF'
pairLR.HBEGF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- 'HBEGF_ERBB4' # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector

# Circle plot
pdf(paste0('Figures/RNA/', LR.show, '_signaling_byCellChat.pdf'), width = 6, height = 6)
netVisual_individual(cellchat, signaling = pathways.show, 
                     pairLR.use = LR.show, layout = "circle")
dev.off()

