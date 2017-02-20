#!/bin/bash

#TSS pickle
tags=2000
TSS=/nfs/latdata/cwng/ChIP/iPSC_MSN/H3K4me3/results/primary_tss/H3K4me3_t${tags}_mrc3_all.pickle

#min tags per profile
min_tags=40
#binning for profile
bins=100
#moving window averaging
avg=100
#format either "BED" or "SAM"
form=SAM

#make profile directory
mkdir ../results/tss_profiles_tse/

#make counts directory
COUNTDIR=../results/tss_counts_t${tags}_mrc3_all/
mkdir $COUNTDIR

#loop through profiling and counts for all conditions
jobids=
for x in HD109n1  HD109n4 HD21n1  HD28n6  HD33n1  HD60n5  HD60n8
do
    SAM="../data/${x}_sorted_fixseq.sam"
    prof_dir=../results/tss_profiles_tse/${x}
    mkdir $prof_dir
    prof_pre=${prof_dir}/TSStags${tags}_tags${min_tags}_bin${bins}_avg${avg}_pdist${pdist}_

    CMD="python tss_profiles_tss_to_tse.py --min_tags $min_tags --format $form --norm density --mv_avg $avg --bin $bins --TSS $TSS $SAM $prof_pre"
    echo $CMD
    $CMD

    OUT=$COUNTDIR/$x
    CMD="wqsub.py python tss_count_tss_to_tse.py --format $form --TSS $TSS $SAM $OUT"
    echo $CMD
    $CMD
done
