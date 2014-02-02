import matplotlib.pyplot as plt
import numpy as np
import csv
from download_aa_sequences import calculate_aa_dist_per_proteome as csp

biosynthetic_path_reader = csv.reader(open('K_biosynthesis.csv'), delimiter='\t')
biosynthetic_path_reader.next()
biosynthetic_path_reader.next()

bio_path_bnumber = []
bio_path_UPID = []

for row in biosynthetic_path_reader:
    bio_path_bnumber.append(row[1]) 
for row in open('all_ecoli_genes.txt'):
    if row[0:5] in bio_path_bnumber:
        bio_path_UPID.append(row[48:54]) 

aa_dist_csv_reader = csv.reader(open('aa_dist_by_UPID.csv'), delimiter='\t')
# skip the first header row
AA_LETTERS = aa_dist_csv_reader.next()[1:21]

UPID_aa_dist = {}
aa_dist_path_genomic = np.zeros((1, len(AA_LETTERS)))
for row in aa_dist_csv_reader:
    if row[0] in bio_path_UPID:
        UPID_aa_dist[row[0]] = np.array([float(x) for x in row[1:]])
        distribution = [float(x) for x in row[1:len(AA_LETTERS)+1]]
        aa_dist_path_genomic += distribution
aa_dist_path_genomic_normed = aa_dist_path_genomic/sum(aa_dist_path_genomic)
a = csp('Ecoli_19_Conditions_Proteomics.csv', UPID_aa_dist)

N = len(AA_LETTERS)
ind = np.arange(N)  # the x locations for the groups
width = 0.35        # the width of the bars
fig = plt.figure()
ax = plt.axes()
rects1 = ax.bar(ind, aa_dist_path_genomic_normed.T, width, color='b')
#rects2 = ax.bar(ind+width,aa_dist_path_proteomic_normed.T , width, color='g')
# add some
ax.set_ylabel('Relative abundance')
ax.set_title('Lysine (K) biosynthesis pathway')
ax.set_xticks(ind+width)
ax.set_xticklabels((AA_LETTERS))
#ax.legend( (rects1[0]),('Proteome'))