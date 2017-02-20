#!/bin/bash

#TSS pickle contains locations of TSS with most H3K4me3 reads for each gene
tags=2000
TSS=H3K4me3_t${tags}_mrc3_all.pickle

#min tags per profile
min_tags=40
#binning for profile
bins=100
#moving window averaging
avg=200
#downstream TSS distance for profiles and counts
pdist=3000
#format either "BED" or "SAM"
form=SAM

#make profile directory
PROFDIR=../results/tss_profiles/
mkdir $PROFDIR

#make counts directory
COUNTDIR=../results/tss_counts/
mkdir $COUNTDIR

#loop through profiling and counts for all conditions
for x in HD109n1 HD109n4 HD 60n5 HD60n8 HD33n1 HD28n6 HD21n1
do
    SAM="../data/${x}_sorted_fixseq.sam"
    prof_dir=$PROFDIR/${x}
    mkdir $prof_dir
    prof_pre=${prof_dir}/TSStags${tags}_tags${min_tags}_bin${bins}_avg${avg}_pdist${pdist}_

    CMD="python tss_profiles_predeftss_multiread.py --pdist $pdist --min_tags $min_tags --format $form --norm density --mv_avg $avg --bin $bins --TSS $TSS $SAM $prof_pre"
    echo $CMD
    $CMD

    OUT=$COUNTDIR/$x
    CMD="python tss_count_predeftss_multiread.py --format $form --TSS $TSS --pdist $pdist $SAM $OUT"
    echo $CMD
    $CMD
done
