#Amanda Daigle

cd ../results/MAnorm/

python ../../src/MAnorm/MAnorm.py -v --p1 ../MACS/HD/HD_no180s_mfold10,30_pval1e-5_peaks.bed --r1 ../../data/allHD.bed --p2 ../MACS/Ctrl/Ctrl_no180s_mfold10,30_pval1e-5_peaks.bed --r2 ../../data/allcontrols.bed --l1 40 --l2 40 -o HDvsCtrl

cd ../../src/
