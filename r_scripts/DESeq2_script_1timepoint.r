# Script for identifying differentially expressed genes in a time course
#  experiment of a Treatment vs Control using the DESeq2 R package
# This script has been adjusted for instances when only one time point is
#  compared 
#
# Arguments:
# [1]: (character) filename of read-count expression matrix; rows = genes, 
#        columns = read counts per library/sample; first column = gene names;
#        column names must have sample/library name
# [2]: (character) filename of metadata for sequencing libraries;
#        Must contain:
#        i) column named <SampleName> with names of technical samples/libraries
#             that corresponds to colunn names of read-count matrix
#        ii) column named <Treatment> which indicates which overall biological
#              sample, line, or treatment the library comes from 
#              ex: Control, Salt, 167 
# [3]: (character) prefix for the file names that will be returned/saved
# [4]: (character) the name of the Control line or treatment as written
#        in the <Treatment> column of the sample metadata
# [5]: (character) the name of the Comparison line or treatment as written
#        in the <Treatment> column of the sample metadata
# [6]: (character) the Time Points to be compared, as written in the <Time> 
#        column of the dample metadata; 
#        time pointes separated by ',' and no spaces
#################

args = commandArgs(trailingOnly = TRUE)

# INPUT FILES #
count_file <- args[1] 
#count_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Camelina_transc_tot_counts_noBadLibs.txt'
counts <- read.table(count_file, header = T, stringsAsFactors = F, sep = '\t')
# NOTE: first column of matrix should be the names of the genes 

samp_meta_file <- args[2]
#samp_meta_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_rnaseq_meta_for_DESeq2.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F,
               sep = '\t')

# SET VARIABLES #
out_pre <- args[3]
#out_pre <- 'test'

control_treatment <- args[4]
#control_treatment <- 'MT5'
comp_treatment <- args[5]
#comp_treatment <- '8171'

time_comps_in <- args[6]
#time_comps_in <- '8DAF,12DAF'
#time_comps <- unlist(strsplit(time_comps_in, split = ','))
time_comps <- time_comps_in

# SET CONSTANTS #
min_floor_count_cut <- 5
# any genes for which no libraries have a greater read count than this are
#   removed

p_cut <- 0.05
# the significance cutoff for the multiple testing-corrected p-values

# LOAD LIBRARYS #
library('DESeq2')

#####################
# SCRIPT
# adjust counts matrix
counts_full <- counts[, -1]
rownames(counts_full) <- counts[,1]

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

# Don't pre-normalize with DESeq2 because calculations require non-normalized
#   data

# Remove genes with low counts across samples
max_count <- apply(counts_1, 1, max)
low_rows <- which(max_count <= min_floor_count_cut)
counts_2 <- counts_1[-low_rows, ]
# REMOVE FOLLOWING LINE WHEN DONE TESTING!!!!
#counts_2 <- counts_2[c(1:500), ]

# format metadata
samp_meta_2$SampleName <- as.factor(samp_meta_2$SampleName)
alt_treats <- setdiff(unique(samp_meta_2$Treatment), control_treatment)
treat_order <- c(control_treatment, alt_treats)
samp_meta_2$Treatment <- factor(samp_meta_2$Treatment, levels = treat_order)
#samp_meta_2$Replicate <- as.factor(samp_meta_2$Replicate)
#samp_meta_2$Time <- factor(samp_meta_2$Time, levels = time_comps)

# Analysis looking for overall differences between lines, without modeling
#  different responses across time
dds <- DESeqDataSetFromMatrix(countData = counts_2, colData = samp_meta_2, 
         design = ~ Treatment)
dds <- DESeq(dds)
res <- results(dds)
resSig <- res[which(res$padj < p_cut), ]
resOrdered <- resSig[order(resSig$pvalue), ]
resOrderedDF <- as.data.frame(resOrdered)
gen_res_df <- data.frame(gene = rownames(resOrderedDF), resOrderedDF, 
                stringsAsFactors = F)
# Only outputting values for comparison of treatments -  the values for the 
#   other variables are not informative about the response of the genes 
#   across time points in the different lines
# Write file with gene names and p-values for Control vs Treatment
out_gen_gene_name <- paste(out_pre, 'DESeq2_general_genes.txt', sep = '_')
write.table(gen_res_df[ , c('gene', 'pvalue', 'padj')], 
  file = out_gen_gene_name, quote = F, sep = '\t', row.names = F, col.names = T)
# Write file with additional info
out_gene_full_file <- paste(out_pre, 'DESeq2_general_full_mat.txt', sep = '_')
write.table(gen_res_df, file = out_gene_full_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

quit(save = 'no')


############ This part is NOT run
# Analysis with modeling response to time; keep Time as a factor because
#   we do NOT expect fold-change to increase with time 
ddsTC <- DESeqDataSetFromMatrix(countData = counts_2, colData = samp_meta_2, 
           design = ~ Treatment + Time + Treatment:Time)
ddsTC <- DESeq(ddsTC, test = 'LRT', reduced = ~ Treatment + Time)
resTC <- results(ddsTC)
keep_inds <- which(resTC$padj < p_cut)
resTCSig <- resTC[keep_inds, ]
resTC_DF <- as.data.frame(resTCSig)

# data.frame with gene and p-value from LRT
TC_gene_df <- data.frame(gene = rownames(resTC_DF), pval = resTC_DF$pvalue, 
  adj_pval = resTC_DF$padj, stringsAsFactors = F) 

# Extract info from the different variables in the model
var_df <- data.frame(gene = rownames(resTC_DF), stringsAsFactors = F)
tot_df <- data.frame(gene = rownames(resTC_DF), baseMean = resTC_DF$baseMean, 
  LRT_pval = resTC_DF$pvalue, LRT_adj_pval = resTC_DF$padj, 
  stringsAsFactors = F)

var_names <- resultsNames(ddsTC)
var_names <- var_names[-grep('^Time', var_names)]
for(i in var_names){
  tmp_res <- results(ddsTC, name = i, test = 'Wald')
  tmp_resSig <- as.data.frame(tmp_res[keep_inds, -1])
  tmp_colnames <- paste(i, colnames(tmp_resSig), sep = '_')
  var_df[, tmp_colnames[1]] <- tmp_resSig[,1]
  tot_df[,tmp_colnames] <- tmp_resSig
}

# Write file with gene names and LRT p-values
out_gene_name <- paste(out_pre, 'DESeq2_TC_genes.txt', sep = '_')
write.table(TC_gene_df, file = out_gene_name, quote = F, sep = '\t', 
  row.names = F, col.names = T)

# Write file with log2foldchange for each of the non-TimevsTime variables
out_var_name <- paste(out_pre, 'DESeq2_TC_var_lfc.txt', sep = '_')
write.table(var_df, file = out_var_name, quote = F, sep = '\t', 
  row.names = F, col.names = T)

# Write file with all the info from the LRT and non-TimevsTime comparisons
out_tot_name <- paste(out_pre, 'DESeq2_TC_full_mat.txt', sep = '_')
write.table(tot_df, file = out_tot_name, quote = F, sep = '\t',     
  row.names = F, col.names = T)

quit(save = 'no')
