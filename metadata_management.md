# Management of Metadata for Project

## Get JGI sequencing info

```
lib_file <- 'Camelina_JGI_Lib_info.txt'
lib_info <- read.table(lib_file, header = F, stringsAsFactors = F, sep = '\t')

lib_info <- gsub('\u2028', 'SPLIT_HERE', lib_info)

lib_2 <- unlist(strsplit(lib_info, split = 'SPLIT_HERE'))

lib_list <- strsplit(lib_2, split = ' ')


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
```
