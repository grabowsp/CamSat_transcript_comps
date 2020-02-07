# Samples used for each comparison and notes about omited samples
* note: The libraries were manually transcribed so there may be some \
transcribing errors - double check other documents if there is a question

## MT5 vs 8171
* 8DAF and 12DAF
### 8DAF
* MT5: GCASC, GCASG
* 8171: GCAPP, GCAPS, GCAPT
### 10DAF - Not included
* Not included because there is low correlation between all 3 8171_10DAF \
replicates
  * DESeq2 requires replication
  * Without any correlated replicates, is impossible to tell which library(s) \
are reliable
### 12DAF
* MT5: GCASP, GCASS, GCAST
* 8171: GCAPZ, GCASA
  * GCAPY was removed because of low correlation with other replicates

## HMT5 vs HMT102 - Not Run
* Did not do this comparison because there is low correlation between the \
two HMT5_10DAF libraries and because the HMT5_15DAF libraries are actually \
from MT5, so would be repeating that comparison

## MT5 vs HMT102
* 15DAF
### 10DAF - Not Run
* Did not use this time point because HMT102 has only 1 10DAF replicate
  * DESeq2 requires replication
  * Without any correlated replicates, is impossible to tell which library(s) \
are reliable
### 15DAF
* MT5: GCAOC, GCAOG, GCAOH
* HMT102: GCAPC, GCAPG

## PR33-Sh vs PS69-Sh
* 0H and 72H
### 0H
* PR33-Sh: GCASU, GCASW, GCASX
* PS69-Sh: GCATB, GCATC, CGATG
### 72H
* PR33-Sh: GCAUH, GCAUN, GCAUO
* PS69-Sh: GCAHS, GCAHT, GCAHU

## PR33-R vs PS69-R
* 0H and 72H
### 0H
* PR33-R: GCASY, GCASZ, GCATA
* PS69-R: GCATH, GCATN, GCATO
### 72H
* PR33-R: GCAUP, GCAHO, GCAHP
* PS69-R: GCAHW, GCAHX, GCAHY

## NR130-Sh vs NS233-Sh
* 0H and 72H
### 0H
* NR130-Sh: GCATP, GCATS, GCATT
* NS233-Sh: GCATY, GCATZ, GCAUA
### 72H
* NR130-Sh: GCAHZ, GCANA, GCANB
* NS233-Sh: GCANN, GCANO, GCANP

## NR130-R vs NS233-R
* 0H and 72H
### 0H
* NR130-R: GCATU, GCATW, GCATX
* NS233-R: GCAUB, GCAUC, GCAUG
### 72H
* NR130-R: GCANC, GCANG, GCANH 
* NS233-R: GCANS, GCANT, GCANU

