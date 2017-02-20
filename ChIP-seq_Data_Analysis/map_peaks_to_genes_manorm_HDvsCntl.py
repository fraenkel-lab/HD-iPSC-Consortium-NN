#/usr/bin/python
#Amanda Daigle

import os

treat_peak_file = '../results/MAnorm/HDvsCtrl/output_tables/biased_peaks_of_HD_all_peaks.bed'
treat_out_file = '../results/MAnorm/HDvsCtrl/output_tables/biased_peaks_of_HD_all_peaks_genes.txt'

control_peak_file = '../results/MAnorm/HDvsCtrl/output_tables/biased_peaks_of_Ctrl_all_peaks.bed'
control_out_file = '../results/MAnorm/HDvsCtrl/output_tables/biased_peaks_of_Ctrl_all_peaks_genes.txt'

print'Mapping peaks from HD biased peaks...\n'
os.system('python map_peaks_to_known_genes_AJK.py --tss --upstream-window=10000 --downstream-window=10000 --map-output=%s --detail /nfs/genomes/human_gp_feb_09/anno/knownGene-hg19-2011-01-04.txt /nfs/genomes/human_gp_feb_09/anno/kgXref-hg19-2011-01-04.txt %s'%(treat_out_file, treat_peak_file))

print'Mapping peaks from control biased peaks...\n'
os.system('python map_peaks_to_known_genes_AJK.py --tss --upstream-window=10000 --downstream-window=10000 --map-output=%s --detail /nfs/genomes/human_gp_feb_09/anno/knownGene-hg19-2011-01-04.txt /nfs/genomes/human_gp_feb_09/anno/kgXref-hg19-2011-01-04.txt %s'%(control_out_file, control_peak_file))
