# Process for aligning libraries and generating count matrix

## Testing
* Make a subfile
* Try aligning the subfile

## Get Genome and Annotation
* Downloaded for JGI
* Genome
  * `/global/homes/g/grabowsp/data/CamSat/Cs_genome_v2.fasta`
* Annotation
  * `/global/homes/g/grabowsp/data/CamSat/Cs_genes_v2_annot_corr.gff3` 

## Index Genome for/with STAR using gff3-specific flags
* Need to include a few extra flags at the end to accomodate format 
differences between gtf and gff3
  * alternatively, I could convert the annotation file to gtf using gffread
    * `https://github.com/gpertea/gffread`
* Index directory
  * `/global/cscratch1/sd/grabowsp/CamSat_transcript/starIndex_v3`
### Generate index
* I ran this interactively on Cori
```
module load python/3.7-anaconda-2019.07
source activate STAR_aligner

STARINDEX=/global/cscratch1/sd/grabowsp/CamSat_transcript/starIndex_v3
REFFASTA=/global/homes/g/grabowsp/data/CamSat/Cs_genome_v2.fasta
GTFFILE=/global/homes/g/grabowsp/data/CamSat/Cs_genes_v2_annot_corr.gff3

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $STARINDEX \
--genomeFastaFiles $REFFASTA --sjdbGTFfile $GTFFILE --sjdbOverhang 100 \
--sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript ID  \
--sjdbGTFtagExonParentGene Parent
```

## Map/Align Libraries with STAR and count gene hits
### Example
```
module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/STAR_aligner

cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_bams/

STARINDEX=/global/cscratch1/sd/grabowsp/CamSat_transcript/starIndex_v3
REFFASTA=/global/homes/g/grabowsp/data/CamSat/Cs_genome_v2.fasta
GTFFILE=/global/homes/g/grabowsp/data/CamSat/Cs_genes_v2_annot_corr.gff3

PREP_R1=/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps/GCAHO/GCAHO.prepped_R1.fastq.bz2
PREP_R2=/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_preps/GCAHO/GCAHO.prepped_R2.fastq.bz2
OUT_COUNT_PREFIX=/global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_bams/GCAHO.star.

STAR --runThreadN 16 --genomeDir $STARINDEX --readFilesCommand bzcat \
--readFilesIn $PREP_R1 $PREP_R2 --outFilterMultimapNmax 7 \
--outFilterMismatchNmax 4 --outFileNamePrefix $OUT_COUNT_PREFIX \
--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif \
--quantMode GeneCounts

```
### Generate submit scripts
* I generated `TESTLIB_star_v3.sh` during the testing phase
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_bams

for LIB in `cat /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transcriptome_libs.txt`;
  do
  sed 's/TESTLIB/'"$LIB"'/g' TESTLIB_star_v3.sh > $LIB'_star.sh';
  done
```
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/CamSat_transcript/Cs_transc_bams

for SHJOB in `ls *star.sh`;
  do 
  sbatch $SHJOB;
  done
```


