# Script for identifying differentially expressed genes between conditions
#  across a time course experiment using the maSigPro R package

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
#samp_meta_file <- samp_meta_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_rnaseq_meta_for_DESeq2.txt'
samp_meta <- read.table(samp_meta_file, header = T, stringsAsFactors = F,
               sep = '\t')

# SET VARIABLES #
out_pre <- args[3]
###out_pre <- 'test'

control_treatment <- args[4]
#control_treatment <- 'MT5'
comp_treatment <- args[5]
#comp_treatment <- '8171'

time_comps_in <- args[6]
#time_comps_in <- '8DAF,12DAF'
time_comps <- unlist(strsplit(time_comps_in, split = ','))

# SET CONSTANTS #
min_floor_count_cut <- 5
# any genes for which no libraries have a greater read count than this are
#   removed

# LOAD LIBRARYS #
library('edgeR')
library('maSigPro')

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
samp_meta_2$Treatment <- gsub('-', '_', samp_meta_2$Treatment)

counts_1 <- counts_full[, samp_meta_2$SampleName]

# generate edisign object
samp_meta_2$Condition <- as.numeric( factor( paste( samp_meta_2$Treatment, 
                           samp_meta_2$Time, sep = '_')))
edes <- data.frame(Time = samp_meta_2$Time, Replicate = samp_meta_2$Condition, 
          stringsAsFactors = F)
rownames(edes) <- samp_meta_2$SampleName

edes$Control <- as.numeric(samp_meta_2$Treatment == control_treatment)
alt_treat_vec <- setdiff(unique(samp_meta_2$Treatment), control_treatment)

for(t in alt_treat_vec){
  treat_name <- paste('treament', t, sep = '_')
  edes[,treat_name] <- as.numeric(samp_meta_2$Treatment == t)
}

# adjust time to be numeric
edes$Time <- as.numeric(gsub('DAF|H', '', edes$Time))

### Maybe add a number of time points and an input ###
# degree should be n_time_points - 1
design <- make.design.matrix(edes, degree = 1)

# Normalize read counts using edgeR function
norm_factors <- calcNormFactors(counts_1, method = 'TMM')
counts_2 <- mapply(`/`, counts_1, norm_factors)
rownames(counts_2) <- rownames(counts_1)

# Remove genes with low counts across samples
max_count <- apply(counts_2, 1, max)
low_rows <- which(max_count <= min_floor_count_cut)
counts_3 <- counts_2[-low_rows, ]

# estimate the GLM theta using edgeR function
est_theta <- 1/(estimateGLMCommonDisp(counts_3))

# run functions to find significant genes
fit <- p.vector(counts_3, design, Q = 0.05, MT.adjust = 'BH', counts = T, 
         min.obs = 1, theta = est_theta)
tstep <- T.fit(fit, step.method = 'backward', alfa = 0.05)
sigs_groups <- get.siggenes(tstep, rsq = 0.7, vars = 'groups')
sigs_each <- get.siggenes(tstep, rsq = 0.7, vars = 'each')

# Generate output files for the General (control vs alternate) DE genes

#file_types_names <- c('DE_genes', 'DE_betas', 'DE_full_mat')

out_list <- list()
#beta_tab <- sigs_groups$sig.genes[[2]]$coefficients
beta_tab <- sigs_each$sig.genes[[2]]$coefficients
#pval_tab <- sigs_groups$sig.genes[[2]]$sig.pvalues
pval_tab <- sigs_each$sig.genes[[2]]$sig.pvalues
# data.frame with gene names and overall p-value
out_list[['DE_genes']] <- data.frame(genes = rownames(pval_tab), 
                            p_val = pval_tab[,4], stringsAsFactors = F)
# data.frame with betas for all variables; can try to use for clustering
#   by behaviour across conditions
out_list[['DE_betas']] <- data.frame(genes = rownames(beta_tab), beta_tab, 
                            stringsAsFactors = F)
# data.frame with overall p-value, R^2, betas for each variable, and
#   p-values for each beta
out_list[['DE_full_mat']] <- data.frame(genes = rownames(pval_tab), 
                               pval_tab, beta_tab, stringsAsFactors = F)
for(df in seq(length(out_list))){
  out_name <- paste( paste(out_pre, 'General', 
                names(out_list)[df], sep = '_'), '.txt', sep = '')
  write.table(out_list[[df]], file = out_name, quote = F, sep = '\t', 
    row.names = F, col.names = T)
}

tc_out_list <- list()
tc_beta_tab <- sigs_each$sig.genes[[4]]$coefficients
tc_pval_tab <- sigs_each$sig.genes[[4]]$sig.pvalues
# data.frame with gene names and overall p-value
tc_out_list[['DE_genes']] <- data.frame(genes = rownames(tc_pval_tab),  
                            p_val = tc_pval_tab[,6], stringsAsFactors = F)
# data.frame with betas for all variables; can try to use for clustering
#   by behaviour across conditions
tc_out_list[['DE_betas']] <- data.frame(genes = rownames(tc_beta_tab), 
                               tc_beta_tab, stringsAsFactors = F)
# data.frame with overall p-value, R^2, betas for each variable, and
#   p-values for each beta
tc_out_list[['DE_full_mat']] <- data.frame(genes = rownames(tc_pval_tab),  
                               tc_pval_tab, tc_beta_tab, stringsAsFactors = F)
for(df in seq(length(tc_out_list))){
  tc_out_name <- paste( paste(out_pre, 'TC',  
                names(tc_out_list)[df], sep = '_'), '.txt', sep = '')
  write.table(tc_out_list[[df]], file = tc_out_name, quote = F, sep = '\t',  
    row.names = F, col.names = T)
}

quit(save = 'no')
