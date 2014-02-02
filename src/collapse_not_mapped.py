from Bio import ExPASy
from Bio import SeqIO
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.backends.backend_pdf import PdfPages


AA_LETTERS = sorted("ACEDGFIHKMLNQPSRTWVY")

# list all proteins in ecoli by uniprot_ID - aprse from .txt file
UPID_list = []
for row in open('all_ecoli_genes.txt'):
    if row: 
        UPID_list.append(row[48:54])
    
# write output file                
output_csv = csv.writer(open('aa_dist_by_UPID.csv', 'w'), 
                        delimiter='\t')

# header for output file
output_csv.writerow(['UPID'] + AA_LETTERS)

# a dictionary containing the aa_dist for each uniprot ID
UPID_to_aa_dist = {}

# the distribution of amino acid for the entire genome
total_genomic_aa_dist = np.zeros((1, len(AA_LETTERS)))

for i, UPID in enumerate(UPID_list[0:50]):  
    print i, "\t", UPID

    # initialize a dictionary for amino acids frequency in each protein
    aa_dist = dict([(aa, 0) for aa in AA_LETTERS])
    
    # call for aa sequence for each uniprot from swiss prot - biopython tool
    handle = ExPASy.get_sprot_raw(UPID)
    seq_record = SeqIO.read(handle, "swiss")
    
    # count frequency for each aa in each UPID
    # update aa_frequency in aa_dict - to avoid bugs where for example an aa seq from
    # swiss prot may contain weired letters such as 'X'
    for aa, count in Counter(seq_record).items():
        if aa in AA_LETTERS:
            aa_dist[aa] = int(count)
    
    UPID_to_aa_dist[UPID] = np.array([aa_dist[aa] for aa in AA_LETTERS])
    output_csv.writerow([UPID] + UPID_to_aa_dist[UPID].tolist())
    
    # sum UPID_to_aa_dist to total_genomic_aa_dist
    total_genomic_aa_dist += UPID_to_aa_dist[UPID]

# read protein abundances ant create a dictionary with uniprot_ID keys
# and abundances array as values
proteomics_csv_reader = csv.reader(open('Ecoli_19_Conditions_Proteomics.csv'), delimiter='\t')
# skip the first empty row
proteomics_csv_reader.next()

conditions = proteomics_csv_reader.next()[10:29]

UPID_to_abundance_vectors = {}
for row in proteomics_csv_reader:
    if row[0]:
        # 19 different growth conditions
        UPID = row[1]        
        UPID_to_abundance_vectors[UPID] = [float(x) for x in row[10:29]]

# calculate the total amino acid distribution of the proteome (i.e. weighted by
# protein abundances), for each of the 19 conditions
proteomic_aa_dist_matrix = np.zeros((len(conditions), len(AA_LETTERS)))
proteomic_aa_dist_normed = np.zeros((len(conditions), len(AA_LETTERS)))
for i, condition in enumerate(conditions):
    total_proteomic_aa_dist = np.zeros((1, len(AA_LETTERS)))
    for UPID in UPID_to_aa_dist.keys():
        aa_dist = UPID_to_aa_dist[UPID]
        if UPID in UPID_to_abundance_vectors:
            abundance = UPID_to_abundance_vectors[UPID]
            total_proteomic_aa_dist += abundance[i] * aa_dist
    
    proteomic_aa_dist_matrix[i, :] = total_proteomic_aa_dist
    proteomic_aa_dist_normed[i, :] = total_proteomic_aa_dist / total_proteomic_aa_dist.sum()

pp = PdfPages('aa_dist.pdf')
for i, condition in enumerate(conditions):
    fig = plt.figure(figsize=(6, 4))
    print condition
    plt.bar(range(len(AA_LETTERS)), proteomic_aa_dist_normed[i, :].T, figure=fig)
    plt.title(unicode(condition), figure=fig)
    pp.savefig(fig)
pp.close()