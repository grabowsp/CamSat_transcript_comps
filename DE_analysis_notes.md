# Notes for DE analsysis of the different comparisons

## MT5 vs 8171
* start by copying r script from rnaseq analysis before
```
cp /home/grabowsky/tools/workflows/rnaseq_TC_contrasts/r_scripts/DESeq2_TC_script.r /home/grabowsky/tools/workflows/CamSat_transcript_comps/r_scripts/DESeq2_script.r
```

```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_DESeq2

sbatch run_DESeq2_01.sh

```


```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_DESeq2

module load python/3.7-anaconda-2019.07
source activate R_RNAseq

Rscript \
/global/homes/g/grabowsp/tools/CamSat_transcript_comps/r_scripts/\
DESeq2_script_1timepoint.r \
/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_counts/\
Camelina_transc_tot_counts_noBadLibs.txt \
/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_rnaseq_meta_for_DESeq2.txt \
MT5_vs_HMT102 MT5 HMT102 15DAF



``` 
