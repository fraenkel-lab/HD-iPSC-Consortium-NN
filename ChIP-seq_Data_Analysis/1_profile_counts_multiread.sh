#!/bin/bash

#TSS pickle
tags=2000
TSS=/nfs/latdata/cwng/ChIP/iPSC_MSN/H3K4me3/results/primary_tss/H3K4me3_t${tags}_mrc3_all.pickle

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
#mkdir ../results/tss_profiles_t${tags}_mrc3_all/
mkdir ../results/tss_profiles/

#make counts directory
COUNTDIR=../results/tss_counts/
mkdir $COUNTDIR

#loop through profiling and counts for all conditions
jobids=
for x in HD109n1 HD109n4 HD 60n5 HD60n8 HD33n1 HD28n6 HD21n1
do
    rep1="../data/${x}_sorted_fixseq.sam"
    SAM="$rep1"
    prof_dir=../results/tss_profiles/${x}
    mkdir $prof_dir
    prof_pre=${prof_dir}/TSStags${tags}_tags${min_tags}_bin${bins}_avg${avg}_pdist${pdist}_

    CMD="wqsub.py python tss_profiles_predeftss_multiread.py --pdist $pdist --min_tags $min_tags --format $form --norm density --mv_avg $avg --bin $bins --TSS $TSS $SAM $prof_pre"
    echo $CMD
    id=`$CMD`
    jobids="$id $jobids"

    OUT=$COUNTDIR/$x
    CMD="wqsub.py python tss_count_predeftss_multiread.py --format $form --TSS $TSS --pdist $pdist $SAM $OUT"
    echo $CMD
    id=`$CMD`
    jobids="$id $jobids"
done
########################################################
wait_for_jobid.py $jobids
########################################################
