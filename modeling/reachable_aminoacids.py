#!/usr/bin/env python

'''
Given a sequence of nucleotides, it returns the reachable amino acids 
'''
import argparse
from Bio import AlignIO
import os
import random

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

'''
Given the amino acid sequences of a protein domain and the nucleotide sequence of a gene or genome, it find the corresponding nucleotide sequence
It returns True if the matching is perfect and the nucleotide sequence corresponding to the domain is written in the output_file 
'''
def reachable_aminoacids(nt_file):
    
    
    with open(nt_file) as nt_handler:
        nt_seq = list(AlignIO.parse(nt_handler, "fasta"))[0]
        nt_record = nt_seq[0]
        seq = nt_record.seq
    
    if len(seq) % 3 != 0:
        print("Warning: the nucleotides sequences is not multiple of 3")
    
    reachable = {}
    
    for i in range(int(len(seq)/3)):
        codon = seq[(i*3):(i*3+3)]
        aa_wt = genetic_code[codon]
        
        print(str(i+1) + "\t" + codon + "\t" + aa_wt + "\t", end = "")
        
        # Make all the possible mutations and record the reachable aminoacids
        # including how many times they are reached
#         aas_reachable = []
        aas_num_codons_reachable = {}
        for j in range(3):
            nt_wt = codon[j]
            
            for nt_mut in nucleotides:
                if nt_wt == nt_mut:
                    continue
            
                codon_mut = list(codon)
                codon_mut[j] = nt_mut
                
                codon_mut_str = ''.join(codon_mut)
                aa_mut = genetic_code[codon_mut_str]
                
                if aa_mut == aa_wt or aa_mut == '*':
                    continue
                
#                 if aa_mut not in aas_reachable:
#                     aas_reachable.append(aa_mut)
                    
                if aa_mut not in aas_num_codons_reachable:
                    aas_num_codons_reachable[aa_mut] = 1
                else:
                    aas_num_codons_reachable[aa_mut] += 1
        
        count = 0
        for aa in aas_num_codons_reachable:
            if count == 0:
                print(aa + ":" + str(aas_num_codons_reachable[aa]), end = "")
                count += 1
            else:
                print("|" + aa + ":" + str(aas_num_codons_reachable[aa]), end = "")
        
        print()


if __name__ == '__main__':
    
#     parser=argparse.ArgumentParser(description='Given a sequence of amino acids of a protein domain and sequence of nucleotides (of a gene or genome), it extract the sequence of nucleotides that codify for the sequence of amino acids')
#       
#     parser = argparse.ArgumentParser()
#     parser.add_argument("-aa", "--aa_file", type=str, help="Input file with the amino acid sequence in FASTA format", required=True)
#     parser.add_argument("-nt", "--nt_file", type=str, help="Input file with the nucleotide sequence in FASTA format", required=True)
#     parser.add_argument("-f", "--frame", type=int, help="Frame (offset) for reading the nucleotide sequence", default = 0)
#     
#     parser.add_argument("-o","--output_file",  help="Output file", type=str, required=True)
#  
#     args = parser.parse_args()
#     
#     aa_file = args.aa_file
#     nt_file = args.nt_file
#     output_file = args.output_file
#     frame = args.frame
    
    # Test
    nt_file = "/media/sakura_6TB/projects/sars_cov2/data/pipeline/03_2021/sequences/domains/Spike.bCoV_S1_RBD.nt.fasta"
    
    reachable_aminoacids(nt_file)






