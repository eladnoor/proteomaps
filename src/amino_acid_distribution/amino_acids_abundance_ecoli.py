from Bio import ExPASy
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import csv
from matplotlib.backends.backend_pdf import PdfPages


def main(UPID_list_fname, abundance_fname, cutoff, output_fname):
    AA_LETTERS = sorted("ACEDGFIHKMLNQPSRTWVY")
    
    # list all proteins in ecoli by uniprot_ID - aprse from .txt file
    UPID_list = []
    for row in open(UPID_list_fname, 'r'):
        if row:
            UPID_list.append(row[48:54])
    
    if cutoff > 0:
        UPID_list = UPID_list[0:min(cutoff, len(UPID_list))]
    
    # write output file                
    output_csv = csv.writer(open(output_fname, 'w'), 
                            delimiter='\t')
    
    # header for output file
    output_csv.writerow(['UPID'] + AA_LETTERS)
    
    # a dictionary containing the aa_dist for each uniprot ID
    UPID_to_aa_dist = {}
    
    # the distribution of amino acid for the entire genome
    total_genomic_aa_dist = np.matrix(np.zeros((1, len(AA_LETTERS))))
    for i, UPID in enumerate(UPID_list):  
        print i, "\t", UPID
    
        # initialize a dictionary for amino acids frequency in each protein
        aa_dist = dict([(aa, 0) for aa in AA_LETTERS])
        
        # call for aa sequence for each uniprot from swiss prot - biopython tool
        handle = ExPASy.get_sprot_raw(UPID)
        seq_record = SeqIO.read(handle, "swiss")
        
        # count frequency for each aa in each UPID
        # update aa_frequency in aa_dict - to avoid bugs where for example an aa seq from
        # swiss prot may contain weired letters such as 'X'
        for aa in list(seq_record):
            if aa in AA_LETTERS:
                aa_dist[aa] += 1
        
        UPID_to_aa_dist[UPID] = np.array([aa_dist[aa] for aa in AA_LETTERS])
        output_csv.writerow([UPID] + list(UPID_to_aa_dist[UPID]))
        
        # sum UPID_to_aa_dist to total_genomic_aa_dist
        total_genomic_aa_dist += np.array(UPID_to_aa_dist[UPID])
    
    # read protein abundances ant create a dictionary with uniprot_ID keys
    # and abundances array as values
    proteomics_csv_reader = csv.reader(open(abundance_fname, 'r'), delimiter='\t')
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
    total_proteomic_aa_dist = np.zeros((len(conditions), len(AA_LETTERS)))
    total_proteomic_aa_dist_normed = np.zeros((len(conditions), len(AA_LETTERS)))
    for i, condition in enumerate(conditions):
        for UPID in UPID_to_aa_dist.keys():
            aa_dist = UPID_to_aa_dist[UPID]
            if UPID in UPID_to_abundance_vectors:
                abundance = UPID_to_abundance_vectors[UPID]
                total_proteomic_aa_dist[i, :] += abundance[i] * aa_dist
        total_proteomic_aa_dist_normed[i, :] = total_proteomic_aa_dist[i, :] / sum(total_proteomic_aa_dist[i, :])
        
    
    total_genomic_aa_dist_normed = total_genomic_aa_dist / total_genomic_aa_dist.sum(1)
    N = len(AA_LETTERS)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35        # the width of the bars
    
    for i, condition in enumerate(conditions):
        fig = plt.figure()
        ax = plt.axes()
        rects1 = ax.bar(ind, total_genomic_aa_dist_normed.T, width, color='r')
        rects2 = ax.bar(ind+width,total_proteomic_aa_dist_normed[i, :].T , width, color='y')
        # add some
        ax.set_ylabel('Relative abundance')
        ax.set_title(condition)
        ax.set_xticks(ind+width)
        ax.set_xticklabels((AA_LETTERS))
        ax.legend( (rects1[0], rects2[0]),('Genome','Proteome'))
    
        fig.savefig('aa_dist_%s.svg' % condition)
        #sys.exit(0)
        
if __name__ == "__main__":
    main('all_ecoli_genes.txt', 'Ecoli_19_Conditions_Proteomics.csv', 0, 'aa_dist_by_UP_ID.csv')
