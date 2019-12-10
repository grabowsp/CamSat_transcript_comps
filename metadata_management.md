# Management of Metadata for Project

## Sample and JGI sequencing info
* Tables are saved in this google Sheet
  * `https://docs.google.com/spreadsheets/d/1ARXuObWDpaJwKm2KNPXpso9qMpAfQFvml2cQ_G0fkIk/edit?usp=sharing`
* Combined file with sample and library info:
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/lib_and_samp_info.txt`
* List of library names
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transcriptome_libs.txt`
### Steps
* Manually input Experiment info from Chaofu's info sheet as "Sample Info"
* Process the JGI Library Sequencing data using R
* Combine the Sample and Sequencing info into single table
  * need to upload this to NERSC
### Process files
* In R
* I did all this on my PC because Cori was down for maintenance
  * all generated tables are saved in the Google Sheet noted above
```
lib_file <- 'Camelina_JGI_Lib_info.txt'
lib_info <- read.table(lib_file, header = F, stringsAsFactors = F, sep = '\t')

lib_info <- gsub('\u2028', 'SPLIT_HERE', lib_info)

lib_2 <- unlist(strsplit(lib_info, split = 'SPLIT_HERE'))

lib_list <- strsplit(lib_2, split = ' ')

# these files have a weird blank space that needs to be removed
mystery_char <- unlist(strsplit(lib_list[[1]][1], split = ''))[6]

lib_list_2 <- lapply(lib_list, function(x) gsub(mystery_char, '', x))

lib_df <- data.frame(lib_name = unlist(lapply(lib_list_2, function(x) x[1])),
  tot_frags = as.numeric(unlist(lapply(lib_list_2, function(x) x[2]))),
  mapped_frags = as.numeric(unlist(lapply(lib_list_2, function(x) x[3]))),
  assigned_frags = as.numeric(unlist(lapply(lib_list_2, function(x) x[4]))),
  unassignedAmbig = as.numeric(unlist(lapply(lib_list_2, function(x) x[5]))),
  unassignedNoFeatures = as.numeric(
    unlist(lapply(lib_list_2, function(x) x[6]))),
  unassignedSecondaryHits = as.numeric(
    unlist(lapply(lib_list_2, function(x) x[7]))),
  ratioStrandedness = as.numeric(unlist(lapply(lib_list_2, function(x) x[8]))),
  stringsAsFactors = F
)

lib_info_out <- 'CamSat_transcript_info.txt'

write.table(lib_df, file = lib_info_out, quote = F, sep = '\t', row.names = F,
  col.names = T)

samp_file <- 'transcriptome_samp_info.txt'

samp_info <- read.table(samp_file, header = T, stringsAsFactors = F, sep = '\t')

setdiff(lib_df$lib_name, samp_info$LibName)
 [1] "GCAPN" "GCAPO" "GCHOT" "GCHOU" "GCHOW" "GBCAS" "GBCAT" "GBCAU" "GBCAW"
[10] "GBCAX"
# These are dropped because their comparison libraries failed

# add sequencing info to sample info df

lib_df_2 <- lib_df[which(lib_df$lib_name %in% samp_info$LibName), ]

lib_df_2_sort <- lib_df_2[order(lib_df_2$lib_name), ]
samp_info_sort <- samp_info[order(samp_info$LibName), ]

#sum(samp_info_sort$LibName == lib_df_2_sort$lib_name)

combo_info <- cbind(samp_info_sort, lib_df_2_sort[, c(2:ncol(lib_df_2_sort))])

combo_out <- 'lib_and_samp_info.txt'

write.table(combo_info, file = combo_out, quote = F, sep = '\t', row.names = F,
  col.names = T)
```

## Generate list of libraries
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript

cut -f 1 lib_and_samp_info.txt | tail -n +2  > Cs_transcriptome_libs.txt
```
