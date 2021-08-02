#!/usr/bin/env python

'''
Given a prediction of the effect of each single mutation and the sequence of the reference in amino acids and nucleotides, average the effect of the given measurement for each position
'''

import argparse
import copy
from Bio import AlignIO
import re
import numpy as np
from statistics import mean 
import sys

genetic_code = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

nucleotides = ['A', 'T', 'G', 'C']

def main(aa_file, nt_file, prediction_file, column, output_file, mode, averaged):
    
#     # Load sequences
#     with open(aa_file) as aa_handler:
#         aa_seq = list(AlignIO.parse(aa_handler, "fasta"))[0]
#         aa_record = aa_seq[0]
#         aa_str = str(aa_record.seq)
# #         aa_str = list(aa_record.seq)
#         
#     with open(nt_file) as nt_handler:
#         nt_seq = list(AlignIO.parse(nt_handler, "fasta"))[0]
#         nt_record = nt_seq[0]
# #         nt_str = list(nt_record.seq)
#         nt_str = str(nt_record.seq)
        
    # TODO: Check nt seq translated is identical to aa seq
#     nt_translated = nt_record.translate()
#     preds_ind = load_predictions(prediction_file, 4)
#     preds_dca = load_predictions(prediction_file, 6)
    
    if mode == 'all':
        avg = average_all_variants(prediction_file, column)
    else:
        avg = average_codons(aa_file, nt_file, prediction_file, column, mode)
#     avg_dca = average_codons(aa_str, nt_str, prediction_file, 6, mode)
    
    with open(output_file, 'w') as fh:
#         print('position\tdelta_IND\tdelta_DCA', file = fh)
        for pos in avg:
            print(str(pos+1) + "\t" + str(avg[pos]), file = fh)
    


'''
Make DCA delta energy average for each position taking the average of all 19 possible variants for each position
'''
def average_all_variants(prediction_file, column, averaged):
    
#     delta_pos = []
    delta_pos = {}
    current_pos = 1
    delta_list = []
    delta_mean = 0
    with open(prediction_file, 'r') as fh:
        for line in fh:
            match = re.search("^position", line)
            if match:
                continue
            
            fields = line.split("\t")
            pos = int(fields[0])
            aa_mut = fields[2]
            aa_wt = fields[1]
            if aa_mut == '-':
                continue
            
#             if(pos != current_pos):
#                 delta_pos.append(mean(delta_list))
#                 current_pos = pos
#                 delta_list = []
            
            delta = float(fields[column])
            if aa_wt != aa_mut:
                delta_list.append(delta)
            
            if(pos != current_pos):
                delta_mean = mean(delta_list)
#                 if len(delta_list) > 0:
#                     delta_mean = mean(delta_list)
#                 delta_pos.append(delta_mean)
                delta_pos[current_pos-1] = np.round(delta_mean,4)
                current_pos = pos
                delta_list = []
    
    # Last position
#     delta_pos.append(delta_mean)
    delta_pos[current_pos-1] = np.round(delta_mean,4)
    
    if averaged:
        average = abs(np.mean(list(delta_pos.values())))
        for i in range(len(delta_pos)):
            delta_pos[i] = np.round(delta_pos[i]/average, 4)
    
#     return np.round(delta_pos,4)
    return delta_pos

'''
Make DCA delta energy average for each position considering the codons and looking at the reachable amino acids. It requires the nucleotide sequence
Available modes:
- 'reachable': Uniform distribution on the reachable amino acids
- 'weighted': Biased to consider the multiplicity of codons, e.g. if a amino acid have several reachable codons then is more likely
'''
def average_codons(aa_file, nt_file, prediction_file, column, mode = 'reachable', averaged = False):
    preds = load_predictions(prediction_file, column)
    
    # Load sequences
    with open(aa_file) as aa_handler:
        aa_seq = list(AlignIO.parse(aa_handler, "fasta"))[0]
        aa_record = aa_seq[0]
        aa_str = list(aa_record.seq)
        
    with open(nt_file) as nt_handler:
        nt_seq = list(AlignIO.parse(nt_handler, "fasta"))[0]
        nt_record = nt_seq[0]
        nt_str = list(nt_record.seq)
    
    return compute_average_codons(aa_str, nt_str, preds, mode, averaged)


def compute_average_codons(aa_str, nt_str, preds, mode = 'reachable', averaged = False):
    
    if mode == 'reachable': 
        pass 
    elif mode == 'weighted':
        pass
    else:
        sys.exit("ERROR: Unrecognized mode: " + mode + ". Please use reachable (uniform distribution on the reachable amino acids) or weighted (biases to consider the multiplicity of codons, e.g. if a amino acid have several reachable codons then is more likely)")
    
    delta_pos = {}
    for pos in preds:
        pos_py = pos-1
        aa_wt = aa_str[pos_py]
        codon = nt_str[(pos_py*3):(pos_py*3+3)]
        
        delta_list = []
        aas_reachable = []
        aas_num_codons_reachable = {}
        for i in range(3):
            nt_wt = codon[i]
            for nt_mut in nucleotides:
                if nt_wt == nt_mut:
                    continue
            
                codon_mut = list(codon)
                codon_mut[i] = nt_mut
                
                codon_mut_str = ''.join(codon_mut)
                aa_mut = genetic_code[codon_mut_str]
                
                if aa_mut == aa_wt or aa_mut == '*':
                    continue
                
                # Check the mutation is in the predictions (it may come from
                # experimental data where not all mutations are measured
                if aa_mut not in preds[pos]:
                    continue
                
                if aa_mut not in aas_reachable:
                    aas_reachable.append(aa_mut)
                    delta_list.append(preds[pos][aa_mut])
                elif mode == 'weighted':
                    delta_list.append(preds[pos][aa_mut])
                    
                if aa_mut not in aas_num_codons_reachable:
                    aas_num_codons_reachable[aa_mut] = 1
                else:
                    aas_num_codons_reachable[aa_mut] += 1
                    
#         print(str(pos) + "\t" + str(len(aas_reachable)))
#         for aa_mut in aas_num_codons_reachable:
#             print(str(pos) + "\t" + aa_wt + "\t" + aa_mut + "\t" + str(aas_num_codons_reachable[aa_mut]))
        
        delta_pos[pos_py] = np.round(mean(delta_list), 4)
        
    if averaged:
        average = abs(np.mean(list(delta_pos.values())))
        for i in range(len(delta_pos)):
            delta_pos[i] = np.round(delta_pos[i]/average, 4)
        
    return delta_pos
    
    
    

def load_predictions(prediction_file, column):
    
    delta_pos = {}
    current_pos = 1
    delta_list = []
    delta_mean = 0
    with open(prediction_file, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
            match = re.search("^position", line)
            if match:
                continue
            
            fields = line.split("\t")
            pos = int(fields[0])
            aa_mut = fields[2]
            aa_wt = fields[1]
            if aa_mut == '-' or aa_mut == '*':
                continue
            
            pred = fields[column]
            if pred == 'NA':
                continue
            
            if pos in delta_pos:
                delta_pos[pos][aa_mut] = float(pred)
            else:
                delta_pos[pos] = {}
                delta_pos[pos][aa_mut] = float(pred)
                
            pass
            
# #             if(pos != current_pos):
# #                 delta_pos.append(mean(delta_list))
# #                 current_pos = pos
# #                 delta_list = []
#             
#             delta = float(fields[column])
#             if aa_wt != aa_mut:
#                 delta_list.append(delta)
#             
#             if(pos != current_pos):
#                 delta_mean = mean(delta_list)
# #                 if len(delta_list) > 0:
# #                     delta_mean = mean(delta_list)
#                 delta_pos.append(delta_mean)
#                 current_pos = pos
#                 delta_list = []
#     
#     # Last position
#     delta_pos.append(delta_mean)
    
    return delta_pos
    

if __name__ == '__main__':
    
    parser=argparse.ArgumentParser(description='Given a prediction of the effect of each single mutation and the sequence of the reference in amino acids and nucleotides, average the effect of the given measurement for each position')
         
    parser = argparse.ArgumentParser()
    parser.add_argument("-aa", "--aa_file", type=str, help="Input file with the amino acid sequence in FASTA format", required=True)
    parser.add_argument("-nt", "--nt_file", type=str, help="Input file with the nucleotide sequence in FASTA format")
    parser.add_argument("-m", "--mode", type=str, help="Optional. Select the mode to make the position average for each position from single mutants predictions. Options: all (average of 19 possible variants), reachable (average from reachable amino acids with one mutation), weighted (average from reachable amino acids with one mutation weigthed by the number of codons that codified each amino acid). Reachable and weighted require the nucleotide sequence. Default: all.", default = "all")
       
    parser.add_argument("-i", "--prediction_file", type=str, help="File with the predictions of the effect of the of every single mutations", required=True)
    parser.add_argument("-c", "--column", type=str, help="Column of the prediction file containing the prediction of interest", required=True)
    parser.add_argument('--averaged', dest='averaged', help="Optional. If True, prediction are corrected by the average per domain (divided by the average per domain). False by default", default=False, action='store_true')
   
    parser.add_argument("-o","--output_file",  help="Output file", type=str, required=True)
    
    args = parser.parse_args()
      
    aa_file = args.aa_file
    nt_file = args.nt_file
    prediction_file = args.prediction_file
    mode = args.mode
    column = int(args.column)-1
    averaged = args.averaged
    output_file = args.output_file
    
    # Test
#     aa_file = "/home/fenix/projects/sars_cov2/temp/Spike.bCoV_S1_RBD.fasta"
#     nt_file = "/home/fenix/projects/sars_cov2/temp/Spike.bCoV_S1_RBD.nucleotides.fasta"
#     prediction_file = "/home/fenix/projects/sars_cov2/temp/Spike.bCoV_S1_RBD.single_muts"
# #     column = args.column
#     mode = "reachable"
#     output_file = "/home/fenix/projects/sars_cov2/temp/preds_averaged"

    # Test
#     aa_file = "/media/sakura_6TB/projects/sars_cov2/data/pierre/models/Results_052021/binding/lambda_0.856/Spike.bCoV_S1_RBD.fasta"
#     mode = "all"
#     prediction_file = "/media/sakura_6TB/projects/sars_cov2/data/pierre/models/Results_052021/binding/lambda_0.856/Spike.bCoV_S1_RBD.not_NA.single_muts"
#     column = 4
#     output_file = "/media/sakura_6TB/projects/sars_cov2/data/pierre/models/Results_052021/binding/lambda_0.856/Spike.bCoV_S1_RBD.not_NA.single_muts.averaged_all"
#     nt_file = ""
    
    main(aa_file, nt_file, prediction_file, column, output_file, mode, averaged)
#     main(aa_file, nt_file, prediction_file, mode, output_file)

    
    
    