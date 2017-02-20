#Amanda Daigle
#8/4/14

import matplotlib
matplotlib.use('pdf')
from matplotlib import pyplot as plt
matplotlib.rcParams['pdf.fonttype']=42

mark='H3K4me3'
treatfile = open('../results/MAnorm/HDvsControl/output_tables/%s_HD_genes_GO_cluster.txt'%mark, 'r')
controlfile = open('../results/MAnorm/HDvsControl/output_tables/%s_Ctrl_genes_GO_cluster.txt'%mark, 'r')

#list of lists for each replicate
top10t = []
while len(top10t) < 10:
    #List for each row of the table
    tline = treatfile.readline()
    if tline.startswith('Category'):
	tline = treatfile.readline()
        twords = tline.strip().split('\t')
        if len(twords) >1:
            try:
                tterm = twords[1].split('~')[1]
            except IndexError:
                tterm = twords[1]
            if tterm.startswith('IPR'):
                tterm = tterm.split(':')[1]
            #concatenate long terms to fit in column
            if len(tterm) > 40:
                tterm = tterm[:40]
            tbon = '%.3g' %float(twords[12])
            #Only include in table if significant
            if float(tbon) > 0.05:
                tterm = ''
                tbon = ''
        top10t.append([tterm,tbon])
treatfile.close()

top10c = []
while len(top10c) < 10:
    cline = controlfile.readline()
    if cline.startswith('Category'):
        cline = controlfile.readline()
        cwords = cline.strip().split('\t')
        try:
            cterm = cwords[1].split('~')[1]
        except IndexError:
            cterm = cwords[1]
        if cterm.startswith('IPR'):
            cterm = cterm.split(':')[1]
        if len(cterm) > 40:
            cterm = cterm[:40]
        cbon = '%.3g' %float(cwords[12])
        if float(cbon) > 0.05:
           cterm = ''
           cbon = ''
        top10c.append([cterm,cbon])
controlfile.close()
    
top10 = [top10t[i]+top10c[i] for i in range(0,10)]

#make table
fig=plt.figure()
ax = fig.add_subplot(111)
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
ax.axis('off')
colLabels=['GO terms for genes with higher %s in HD'%mark, 'FDR', 'GO terms for genes with lower %s in HD'%mark, 'FDR']
the_table = ax.table(cellText=top10, 
                     colLabels=colLabels, 
                     colWidths=[0.3,0.09,0.3,0.09],
                     loc='center')
the_table.auto_set_font_size(False)
the_table.set_fontsize(6.5)
the_table.scale(1.3,1.3)
plt.title('Gene Ontology results for differential %s peaks'%(mark))
plt.savefig('../results/MAnorm/HDvsControl/output_tables/%s_GO_table.pdf'%mark)
