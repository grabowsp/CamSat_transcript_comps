# Script to make General DE gene lists from Camelina seed time course experiment
#   Get overlap of DESeq2-General, maSigPro, and splineTimeR-general

# Load 'R_GO_analysis' conda environment so all packages are available'

# LOAD FILES #
## Information about comparisons for loading and filtering data
comp_info_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/timecourse_info_for_combo_files.txt'
comp_info <- read.table(comp_info_file, header = T, sep = ' ', 
  stringsAsFactors = F)

## File with A.thaliana homologs
cs_to_at_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cam_to_AT_gene_table_v2.0.txt'
cs_to_at <- read.table(cs_to_at_file, header = T, stringsAsFactors = F, 
  sep = '\t')
# GO info
cs_go_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cam_GO_table.txt'
# Count data
count_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Camelina_transc_tot_counts_noBadLibs.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# Metadata
metadata_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_rnaseq_meta_for_DESeq2.txt'
meta <- read.table(metadata_file, header = T, stringsAsFactors = F, sep = '\t')

# SET VARIABLES #

deseq_dir <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_DESeq2/'

deseq_suf <- '_DESeq2_TC_genes.txt'

masig_dir <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_maSigPro/'
masig_suf <- '_TC_DE_genes.txt'

#control_treatment <- 'MT5'
#comp_treatment <- '8171'
#time_comps <- c('8DAF','12DAF')


# SET CONSTANTS #
min_floor_count_cut <- 5

# SET OUTPUT INFO #
out_dir <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_DE_combofiles/'
out_suf <- '_GeneList2_TimeDependent_DE_genes.txt'
#gen_DE_gene_out_file <- 'tmp_out'
#'/home/t4c1/WORK/grabowsk/data/Camelina_suneson_seed_TC/Csa_suneson_TC_GeneList1_DE_genes_v2.0.txt'


# LOAD PACKAGES #
library('topGO')
library('edgeR')

####################3
for(COMP_NUM in seq(nrow(comp_info))){
  comp_name <- comp_info$comp_name[COMP_NUM]
  control_treatment <- comp_info$control_treatment[COMP_NUM]
  comp_treatment <- comp_info$comp_treatment[COMP_NUM]
  time_comps <- unlist(strsplit(comp_info$time_comps[COMP_NUM], split = ','))

  deseq_tc_file <- paste(deseq_dir, comp_name, deseq_suf, sep = '')
  # General differences between the lines
  deseq_tc <- read.table(deseq_tc_file, header = T, stringsAsFactors = F, 
    sep = '\t')

  masig_tc_file <- paste(masig_dir, comp_name, masig_suf, sep = '')
  # General differences between the lines
  masig_tc <- read.table(masig_tc_file, header = T, stringsAsFactors = F, 
    sep = '\t')

  # add rank columns to result files
  deseq_tc$rank <- order(deseq_tc$adj_pval)
  masig_tc$rank <- order(masig_tc$p_val)

  # Generate "General" Gene list
  tc_names <- intersect(deseq_tc$gene, masig_tc$genes)

  # Make file with "general" genes
  tc_names_2 <- sort(tc_names)
  deseq_tc_ord <- deseq_tc[order(deseq_tc$gene),]
  masig_tc_ord <- masig_tc[order(masig_tc$genes),]

  # p-values are all multiple-testing corrected
  tc_DE_df <- data.frame(gene = tc_names_2, stringsAsFactors = F)

  # ASSIGN ARABIDOPSIS NAMES
  tc_inds_with_at <- which(tc_DE_df$gene %in% cs_to_at$cs_gene)
  at_to_incl_inds <- which(cs_to_at$cs_gene %in%tc_DE_df$gene)

  tc_DE_df$At_gene <- NA
  tc_DE_df$At_symbol <- NA
  tc_DE_df$At_description <- NA

  tc_DE_df$At_gene[tc_inds_with_at] <- cs_to_at$at_short[at_to_incl_inds]
  tc_DE_df$At_symbol[tc_inds_with_at] <- cs_to_at$at_symbol[at_to_incl_inds]
  tc_DE_df$At_description[tc_inds_with_at] <- cs_to_at$at_descr[at_to_incl_inds]

  # ASSIGN GENE NAME, P-VALUES AND RANKS FOR EACH APPROACH
  tc_DE_df$deseq_pval <- deseq_tc_ord$adj_pval[which(deseq_tc_ord$gene %in% 
    tc_DE_df$gene)]
  tc_DE_df$deseq_rank <- deseq_tc_ord$rank[which(deseq_tc_ord$gene %in% 
    tc_DE_df$gene)]
  tc_DE_df$masig_pval <- masig_tc_ord$p_val[which(masig_tc_ord$genes %in% 
    tc_DE_df$gene)]
  tc_DE_df$masig_rank <- masig_tc_ord$rank[which(masig_tc_ord$genes %in% 
    tc_DE_df$gene)]

  # ASSIGN GO TERMS
  geneID2GO <- readMappings(file = cs_go_file)

  tc_DE_df$gene_short <- unlist(lapply(
    strsplit(tc_DE_df$gene, split = '.', fixed = T), function(x) x[1]))

  tc_DE_df$go_terms <- NA
  for(i in seq(nrow(tc_DE_df))){
    tmp_go_terms <- unlist(geneID2GO[tc_DE_df$gene_short[i]])
    if(length(tmp_go_terms) > 0){tc_DE_df$go_terms[i] <- paste(tmp_go_terms, collapse = ';')}
  }

  # ASSIGN EXPRESSION LEVELS

  # adjust counts matrix
  counts_full <- counts[, -1]
  rownames(counts_full) <- counts[,1]

  samp_meta <- meta
  # Re-order metadata so samples are in same order as count data
  ph_order <- c()
  for(i in seq(ncol(counts_full))){
    temp_ind <- which(samp_meta$SampleName == colnames(counts_full)[i])
    ph_order <- c(ph_order, temp_ind)
  }
  samp_meta_full <- samp_meta[ph_order, ]

  # adjust count and metadata to only include required libraries
  control_meta_inds <- which(samp_meta_full$Treatment == control_treatment)
  comp_meta_inds <- which(samp_meta_full$Treatment == comp_treatment)

  meta_time_inds <- which(samp_meta_full$Time %in% time_comps)

  samp_meta_2 <- samp_meta_full[
    intersect(c(control_meta_inds, comp_meta_inds), meta_time_inds), ]

  counts_1 <- counts_full[, samp_meta_2$SampleName]

  norm_factors <- calcNormFactors(counts_1, method = 'TMM')
  counts_2 <- mapply(`/`, counts_1, norm_factors)
  rownames(counts_2) <- rownames(counts_1)

  for(i in seq(ncol(counts_2))){
    meta_ind <- which(samp_meta_2$SampleName == colnames(counts_2)[i])
    colnames(counts_2)[i] <- paste(samp_meta_2$SampleName[meta_ind], 
      samp_meta_2$Treatment[meta_ind], samp_meta_2$Time[meta_ind], sep = '_')
  }

  count_inds <- which(rownames(counts_2) %in% tc_DE_df$gene)
  tc_DE_df[, colnames(counts_2)] <- counts_2[count_inds, ]

  DE_out_file <- paste(out_dir, comp_name, out_suf, sep = '')

  write.table(tc_DE_df, file = DE_out_file, quote = F, sep = '\t', 
    row.names = F, col.names = T)
}

quit(save = 'no')

