# Need to split libraries into R1 and R2 Reads to Align using STAR

## Example of command used to split the R1 and R2 reads
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps/GCAHO

bzip2 -dc GCAHO.prepped.fastq.bz2 | head -200 | \
awk '{if(NR%8 >= 1 && NR%8 < 5) print}' | bzip2 -c > test_R1.fastq.bz2

bzip2 -dc GCAHO.prepped.fastq.bz2 | head -200 | \
awk '{if(NR%8 == 0 || NR%8 >= 5) print}' | bzip2 -c > test_R2.fastq.bz2
```

## Full splitting job
### Command
* On NERSC
	* `/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps/split_prepped_fastq_files.sh`
### Submit job
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
sbatch split_prepped_fastq_files.sh
```
### Need to break it down into smaller sample sets
* Make smaller jobs for the final 27 libraries
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps
sbatch split_prep_fastq_sub01.sh
sbatch split_prep_fastq_sub02.sh
sbatch split_prep_fastq_sub03.sh
```



