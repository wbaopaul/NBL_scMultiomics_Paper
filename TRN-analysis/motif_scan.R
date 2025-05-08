library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(GenomicRanges)


## check motif overlapping with any peak -- new peaks ####
load('data/intermediate/coembedBySample/new_metacells/rna_atac_metacell_matrices.RData')
peaks = rownames(metacell.atac.mtx)

pks = data.table(tidyr::separate(data = data.table('peak_name' = peaks),
                                 col = 'peak_name', into = c('chr', 'start', 'end'),
                                 remove = F))
pks$start = as.integer(pks$start)
pks$end = as.integer(pks$end)
setkey(pks, chr, start)
# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = pks$chr,
                                ranges = IRanges::IRanges(start = pks$start,
                                                          end = pks$end))
pv.cutoff = 5e-05
motif_ix <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = pv.cutoff)

motif_name = motif_ix@colData$name
pk_inf = motif_ix@rowRanges
motif_pk_match <- motif_ix@assays@data$motifMatches
rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
names(motif_name) = NULL
colnames(motif_pk_match) = motif_name      

saveRDS(motif_pk_match, file = paste0('data/intermediate/TRN/peak_motif_match_pv', pv.cutoff, '.rds'))
