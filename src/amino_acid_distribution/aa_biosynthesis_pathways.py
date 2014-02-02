import matplotlib.pyplot as plt
import numpy as np
import csv
from download_aa_sequences import calculate_aa_dist_per_proteome as csp
from download_aa_sequences import normalize_aa_dist as nd

AA_LETTERS = sorted("ACEDGFIHKMLNQPSRTWVY")
aa_biosythesis_reader = csv.reader(open('aa_biosynthesis.csv'), delimiter='\t')

for row in aa_biosythesis_reader:
    print row[0,:] 


#bio_path_UPID = []
#aa_dist_csv_reader = csv.reader(open('aa_dist_by_UP_ID.csv'), delimiter='\t')
## skip the first header row
#aa_dist_csv_reader.next()
#
#load_UPID_to_aa_dist
#
#
#UPID_aa_dist = {}
#aa_dist_path_genomic = np.zeros((1, len(AA_LETTERS)))
#for row in aa_dist_csv_reader:
#    if row[0] in bio_path_UPID:
#        UPID_aa_dist[row[0]] = np.array([float(x) for x in row[1:]])
#        distribution = [float(x) for x in row[1:len(AA_LETTERS)+1]]
#        aa_dist_path_genomic += distribution
#aa_dist_path_genomic_normed = aa_dist_path_genomic/sum(aa_dist_path_genomic)
#a = csp('Ecoli_19_Conditions_Proteomics.csv', UPID_aa_dist)
##print a
#a_normed = nd(a,20)
#
#AA_PATH = sorted(['P_'+ x for x in AA_LETTERS], reverse=True)
#plt.clf()
#plt.imshow(a_normed, interpolation="nearest")
#plt.colorbar()
#plt.show()
#N = len(AA_LETTERS)
#ind = np.arange(N)  # the x locations for the groups
#width = 0.1        # the width of the bars
#ax = plt.axes()
#ax.set_ylabel('Relative abundance')
#ax.set_title('Lysine (K) biosynthesis pathway')
#ax.set_xticks(ind+width)
#ax.set_xticklabels((AA_LETTERS))
#ax.set_yticks(ind+width)
#ax.set_yticklabels((AA_PATH))

#fig = plt.figure()

#rects1 = ax.bar(ind, aa_dist_path_genomic_normed.T, width, color='b')
##rects2 = ax.bar(ind+width,aa_dist_path_proteomic_normed.T , width, color='g')
## add some
##ax.legend( (rects1[0]),('Proteome'))