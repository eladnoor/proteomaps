from Bio import ExPASy
from Bio import SeqIO
import numpy as np
import csv
AA_LETTERS = sorted("ACEDGFIHKMLNQPSRTWVY")

# list all proteins in ecoli by uniprot_ID - aprse from .txt file
def download_aa_dist_per_gene(UPID_list_fname, cutoff):
    UPID_list = []
    for row in open(UPID_list_fname, 'r'):
        if row:
            UPID_list.append(row[48:54])
    
    if cutoff > 0:
        UPID_list = UPID_list[0:min(cutoff, len(UPID_list))]    
    
    # a dictionary containing the aa_dist for each uniprot ID
    UPID_to_aa_dist = {}
    
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
    return UPID_to_aa_dist

def load_UPID_to_aa_dist(aa_dist_by_gene_fname):
    input_csv = csv.reader(open(aa_dist_by_gene_fname, 'r'), delimiter='\t')
    input_csv.next()
    UPID_to_aa_dist = {}
    for row in input_csv:
        UPID = row[0]
        UPID_to_aa_dist[UPID] = np.array([float(x) for x in row[1:]])
    return UPID_to_aa_dist

def calculate_aa_dist_per_genome(UPID_to_aa_dist):
    genomic_aa_dist = np.zeros((1, len(AA_LETTERS)))
    for aa_dist in UPID_to_aa_dist.values():
        genomic_aa_dist += aa_dist
    return genomic_aa_dist
    
def write_to_tsv(header, dictionary, output_fname):
    # write output file                
    output_csv = csv.writer(open(output_fname, 'w'), delimiter='\t')
    # header for output file
    output_csv.writerow(header)
    for key in dictionary.keys():
        output_csv.writerow([key] + list(dictionary[key]))
            
def calculate_aa_dist_per_proteome(proteome_fname, UPID_to_aa_dist):
    proteomics_csv_reader = csv.reader(open(proteome_fname, 'r'), delimiter='\t')
    # skip the first empty row
    proteomics_csv_reader.next()
    conditions = proteomics_csv_reader.next()[10:29]
    UPID_to_abundance_vectors = {}
    total_proteomic_aa_dist = np.zeros((len(conditions), len(AA_LETTERS)))    
    for row in proteomics_csv_reader:
        if row[0]:
            # 19 different growth conditions
            UPID = row[1]        
            UPID_to_abundance_vectors[UPID] = [float(x) for x in row[10:29]]
    for i, condition in enumerate(conditions):
        for UPID in UPID_to_aa_dist.keys():
            aa_dist = UPID_to_aa_dist[UPID]
            if UPID in UPID_to_abundance_vectors:
                abundance = UPID_to_abundance_vectors[UPID]
                total_proteomic_aa_dist[i, :] += abundance[i] * aa_dist
    return total_proteomic_aa_dist
                
def normalize_aa_dist(total_proteomic_aa_dist, conditions):
    total_proteomic_aa_dist_normed = a_normed = np.zeros((conditions, len(AA_LETTERS)))
    for i, row in enumerate(total_proteomic_aa_dist): 
        total_proteomic_aa_dist_normed[i] = total_proteomic_aa_dist[i] / sum(total_proteomic_aa_dist[i])
    return total_proteomic_aa_dist_normed
        
if __name__ == "__main__":
    UPID_to_aa_dist = download_aa_dist_per_gene('all_ecoli_genes.txt',20)
    write_to_tsv(['UPID'] + AA_LETTERS, UPID_to_aa_dist, 'aa_dist_by_UP_ID.csv')
    aa_dist_genome = calculate_aa_dist_per_genome(UPID_to_aa_dist)
    total_proteomic_aa_dist = calculate_aa_dist_per_proteome('Ecoli_19_Conditions_Proteomics.csv', UPID_to_aa_dist)        
    a = normalize_aa_dist(total_proteomic_aa_dist)
    print a
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    