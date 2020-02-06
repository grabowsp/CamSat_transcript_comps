# Commands used to generate the count matrix

## Important Files
* Count matrix for all libraries
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Camelina_transc_tot_counts.txt`
* Count matrix for libraries used in DE analysis
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Camelina_transc_tot_counts_noBadLibs.txt`
* File with information about each library
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/transc_comp_samp_and_lib_info.tsv`
* List of libraries to be excluded from DE analysis because of poor
correlation with or lack of replicates
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/`
* Correlation heatmap using raw read counts for each library
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Cs_transc_correlation_heatmap.pdf`

## Notes
* Generated reads-per-gene counts for each library using STAR
* Need to combine the results into a single matrix
* The format of the `...ReadsPerGene.out.tab` file generated by STAR is:
  * Column 1: gene ID
  * Column 2: counts for unstranded RNA-seq
  * Column 3: counts for the 1st read strand aligned with RNA (forward strand)
  * Column 4: counts for the 2nd read strand aligned with RNA (reverse strand)
* JGI only includes the reverse strand counts in their tally
  * I will stick with that convention for now 
    * The vast majority of hits are reverse strand, anyhow

## Generate Initial Count Matrix from individual files
* Final matrix
  `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Camelina_transc_tot_counts.txt`
### Script for generating the matrix
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts

for LIB in `cat /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transcriptome_libs.txt`;
  do
  echo $LIB > tmp_head.txt;
  cut -f 4 ../Cs_transc_bams/$LIB'.star.ReadsPerGene.out.tab' | tail -n +5 | \
    cat tmp_head.txt - > tmp_$LIB'_counts.txt';
  done

echo GeneID > tmp_head.txt
cut -f 1 ../Cs_transc_bams/GCHPB.star.ReadsPerGene.out.tab | tail -n +5 | \
cat tmp_head.txt - > tmp_geneID.txt

paste tmp_geneID.txt tmp*counts.txt > Camelina_transc_tot_counts.txt
``` 

## Generate correlation matrix
* `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/Cs_transc_correlation_heatmap.pdf`
### Problem Libraries Identified with Correlation Heatmap
* HMT5_10DAF - won't be able to use
  * GCAOB: HMT5_10DAF - not correlated with replicate; very highly correlated with 15DAF libraries
  * GCANZ: HMT5_10DAF - not correlated with replicate; correlated with some 10DAF and 12 DAF libraries 
* HMT102_10DAF - won't be able to use
  * GCAPB: HMT102_10DAF - only one replicate; somewhat correlated with some 10DAF and 12 DAF libraries; 
* 8171_10DAF - won't be able to use because no reliable replication
  * GCAPU: 8171_10DAF - low correlation with replicates
  * GCAPW: 8171_10DAF - low correlation with replicates
  * GCAPX: 8171_10DAF - low correlation with replicates
* 8171_12DAF - only include other 2 replicates
  * GCAPY: 8171_12DAF - not well correlated with replicates
### File with problem libraries
* `/global/cscratch1/sd/grabowsp/CamSat_transcript/remove_libs.txt`
### Script to generate correlation heatmap
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

R

library(reshape2)
library(ggplot2)

count_file <- paste('/global/cscratch1/sd/grabowsp/CamSat_transcript/', 
  'Cs_transc_counts/Camelina_transc_tot_counts.txt', sep = '')
counts_0 <- read.table(count_file, header = T, stringsAsFactors = F, 
  sep = '\t')

info_file <- paste('/global/cscratch1/sd/grabowsp/CamSat_transcript/', 
  'Cs_transc_counts/transc_comp_samp_and_lib_info.tsv', sep = '')
info_0 <- read.table(info_file, header = T, sep = '\t', stringsAsFactors = F)

count_1 <- counts_0

info_1 <- info_0[order(info_0$Time), 1:3]

info_2 <- info_1[order(info_1$Sample), ]

counts_1 <- counts_0[, info_2$LibName]

descr_names <- paste(info_2$Sample, info_2$Time, info_2$LibName, sep = '_')

colnames(counts_1) <- descr_names

count_cor <- cor(x = counts_1, method = "pearson")

count_cor_round <- round(count_cor, digits = 2)

count_cor_melt <- melt(count_cor_round)

gg_heat <- ggplot(data = count_cor_melt, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = 'red', high = 'darkgreen', mid = 'white', 
    midpoint = 0.5, limit = c(0,1), space = 'Lab', 
    name = 'Pearson\nCorrelation') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 14, 
    hjust = 1), axis.text.y = element_text(size = 14)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = value), color = 'black', size = 4)

heat_file_1 <- paste('/global/cscratch1/sd/grabowsp/CamSat_transcript/',
  'Cs_transc_counts/Cs_transc_correlation_heatmap.pdf', sep = '')

pdf(heat_file_1, width = 12*4, height = 9*4)
gg_heat
dev.off()
```

## Generate Filtered Counts Matrix
* Libraries not used in DE analysis removed
```
module load python/3.7-anaconda-2019.07
source activate R_analysis

R

count_file <- paste('/global/cscratch1/sd/grabowsp/CamSat_transcript/',
  'Cs_transc_counts/Camelina_transc_tot_counts.txt', sep = '')
counts_0 <- read.table(count_file, header = T, stringsAsFactors = F,
  sep = '\t')

bad_lib_file <- paste('/global/cscratch1/sd/grabowsp/CamSat_transcript/', 
  'remove_libs.txt', sep = '')
bad_libs <- read.table(bad_lib_file, header = F, stringsAsFactors = F)

bad_col_inds <- c()
for(i in seq(nrow(bad_libs))){
  tmp_ind <- which(colnames(counts_0) == bad_libs[i, 1])
  bad_col_inds <- c(bad_col_inds, tmp_ind)
}

counts_filt <- counts_0[, -bad_col_inds]

filt_file_out <- paste('/global/cscratch1/sd/grabowsp/CamSat_transcript/',
  'Cs_transc_counts/Camelina_transc_tot_counts_noBadLibs.txt', sep = '')
write.table(counts_filt, file = filt_file_out, quote = F, sep = '\t', 
  row.names = F, col.names = T)
```

