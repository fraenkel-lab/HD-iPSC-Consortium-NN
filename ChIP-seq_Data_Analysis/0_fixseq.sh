#!/bin/bash

for x in HD180n3 HD109n1  HD109n4  HD109n5  HD180n1 HD21n1  HD28n6  HD33n1  HD60n5  HD60n8
do
    echo $x
    # The first step is to create the read count histogram for Fixseq
    # input.  These flags exclude unmapped reads, secondary alignments,
    # and QC failures.  You may also consider excluding high-copy
    # sequences like mitochondrial reads.
    BAM="../data/${x}_sorted.bam"
    echo 'Getting counts histogram from alignment file...'
    samtools view -F 772 $BAM | cut -f 3,4 | uniq -c | awk '{print $1}' | sort | uniq -c | sort -n | awk '{print $1","$2}' > input.csv

    # Next, run methods.r and go through the resulting output.csv file.
    # Choose which rounding scheme you'd like to employ, and create a
    # remapping count file.

    echo 'Running fixseq methods.r...'
    Rscript /Software/fixseq/thashim-fixseq-b7f0e45600be/methods.r
    cut -f 2,5 -d "," output.csv | tr "," "\t" | tail -n +2 | sort -n > counts_map.txt

    # Now use BEDTools and this small Python script to remap read counts
    # by expanding/contracting a bed file.  IMPORTANT: The original bam
    # must be coordinate-sorted, since the remapper expects duplicate
    # reads to be adjacent to each other in the bam (then bed).

    echo 'Re-mapping reads according to fixseq output...'
    pyth=/Software/fixseq/thashim-fixseq-b7f0e45600be/bed_count_remapper.py
    bamToBed -i $BAM | python $pyth counts_map.txt > output.bed
    echo 'Converting output.bed to BAM...'
    bedToBam -i output.bed -g /nfs/genomes/human_gp_feb_09/hg19.chrom.sizes > output.bam
    echo 'Sort and index BAM file...'
    newname=$(basename $BAM .bam)_fixseq
    samtools sort output.bam $newname
    samtools index ${newname}.bam
    samtools view -h ${newname}.bam > ${newname}.sam

    rm output.bed
    rm output.bam
    rm input.csv
    rm output.csv
done
# Finally, run the peak caller or other processing algorithm of your
# choice on output.bed.  You could also pipe the above string directly
# into the peak caller in place of its bed input.

