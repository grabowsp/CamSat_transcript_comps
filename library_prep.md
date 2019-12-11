# CamSat Transcriptome Library Prep

## Steps
1. Make file with library names
2. Prepare directory
3. Generate lib.config file
4. Make sure included libraries did not FAIL the QC steps
3. Run prepper
* use CamSat_snp_calling as template for this

## Prepare Directory for for libarary prepping
### Overview
* Directory on NERSC:
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps`
* Use Chris's prepper program to prepare the libraries for SNP calling
  * Info about prepper found in shared Google Doc
    * `https://docs.google.com/document/d/1ISRXR1s6WHnNDBQSN0eTu1-11gajfdHRouUMX03Ev3Y/edit`
* Requires setting up many subdirectories
  * Ex: CONFIG, linkerSeq, contaminant, directories with each library code...
  * Info found in GoogleDoc
### Generate Subdirectories
* Generate necessary subdirectories and copy files to `linkerSeq` directory
```
bash
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
mkdir CONFIG
mkdir contaminant
mkdir linkerSeq
cp /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/PREP_TESTING/linkerSeq/* ./linkerSeq
for i in `cat ../Cs_transcriptome_libs.txt`; do mkdir ./$i; done
```

## Generate lib.config file
### Overview
* Requires using the jamo module/function to get information about the \
status and location of the libraries at JGI
* Used `jamo info LIB_CODE` to get information about the libraries
  * Any 'PURGED' libraries had to be 'fetched' using \
`jamo fetch library LIB_CODE`
* Use `jamo report short_form library LIB_CODE` to get info about the \
libraries that could be put into the `lib.config` file
  * Need to include a `grep UNPROCESSED` term so that info about all \
libraries tied to a `LIB_CODE` are included in the prepping
### Get status of libraries
```
bash
module load jamo
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
for i in `cat ../Cs_transcriptome_libs.txt`;
do jamo info library $i >> transcriptome_jamo_lib_info.txt; done
```
### Fetch any purged libraries
* Use `grep PURGED transcriptome_jamo_lib_info.txt` to get `LIB_CODE` of any \
libraries that were not in the main storage of JGI
* Use `jamo fetch library LIB_CODE` for those libraries
* ex:
```
for i in `grep PURGED transcriptome_jamo_lib_info.txt`;
do jamo fetch library $i; done
```
* Need to wait until `jamo info library LIB_CODE` shows a state of \
`BACKUP_COMPLETE` or `RESTORED` before going on to next step
```
for i in `cat ../Cs_transcriptome_libs.txt`;
do jamo info library $i; done
```
### Generate lib.config file
```
bash
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
for i in `cat transcriptome_jamo_lib_info.txt`;
do jamo report short_form library $i | grep UNPROCESSED >> ./CONFIG/lib.config;
done
```
### Look for and adjust LIBs sequenced 2+ times
* There are no duplicated libraries
* in R
```
config_file <- '/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps/CONFIG/lib.config'

config <- read.table(config_file, header = F, sep = '\t', stringsAsFactors = F)

dup_inds <- which(duplicated(config[,1]))
# none
```

## Check that libraries pass QC
### Overview
* use jamo to see if any `LIB_CODE`s are linked to fastq files that failed \
the JGI QC test
* compare the Failed QC test files to those in lib.config
### Looked for Failed runs
* need to load the jamo module
```
module load jamo
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
for i in `cat ../Cs_transcriptome_libs.txt`;
do jamo report short_form library $i | grep 'passedQC=F' >> qc_F_report.txt;
echo $i; done
```
* the file is empty, so none of the libraries failed QC
### Check if failed runs are included in lib.config
* Did not have to run this
```
#qc_report_file <- '/global/projectb/scratch/grabowsp/Csativa_reseq/snp_call_1/qc_F_report.txt'
qc_report <- read.table(qc_report_file, header = F, sep = ' ',
  stringsAsFactors = F)
#
#lib_config_file <- '/global/projectb/scratch/grabowsp/Csativa_reseq/snp_call_1/CONFIG/lib.config'
lib_config <- read.table(lib_config_file, header = F, sep = '\t',
  stringsAsFactors = F)
#
f_info <- strsplit(qc_report[,4], split = '/')
f_files <- unlist(lapply(f_info, function(x) x[[length(x)]]))
#
f_in_lc_list <- sapply(f_files, function(x) grep(x, lib_config[,6], fixed = T))
f_in_lc_length <- unlist(lapply(f_in_lc_list, length))
f_in_lc <- which(f_in_lc_length > 0)
length(f_in_lc) > 0
# [1] FALSE
```
* no overlap between QC failure and files in lib.config

## Submit prepping jobs
### Overview
* Generated sub-lists of 45 `LIB_CODE`s because launching the prepper submits \
many jobs for each `LIB_CODE`, and Chris recommended prepping 45 libraries \
at a time if use '-q 100' flag to stay withing the rules of NERSC
### Generate sub-lists
* Generate sub-lists of 45 `LIB_CODE` entries for submitting jobs
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript
split -l 45 -d Cs_transcriptome_libs.txt Cs_transc_small_list_
```
### Submit jobs
* have submitted: `Cs_transc_small_list_00`; `Cs_transc_small_list_00`
```
bash
module load python
source activate /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/ANACONDA_ENVS/PREP_ENV/
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps

for i in `cat ../Cs_transc_small_list_01`;
  do
  cd ./$i;
  python3 /global/dna/projectdirs/plant/geneAtlas/HAGSC_TOOLS/PREP_TESTING/splittingOPP.py \
/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps $i -q 100;
  cd ..;
  sleep 5s;
  done
```

## Checking that preps have finished
### Overview
* check that there is a `prepComplete` file for each library
* for those that didn't finish, check what happened
### Check for `prepComplete` file
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
for i in `cat ../Cs_transc_small_list_01`;
do ls -oh ./$i/prepComplete; done
```
* all finished for: `Cs_transc_small_list_00`;, `Cs_transc_small_list_01`







