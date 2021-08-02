#!/usr/bin/env python

'''
Given a sequence of amino acids of a protein domain and sequence of nucleotides (of a gene or genome), it extracts the sequence of nucleotides that codify for the sequence of amino acids 
'''
import argparse
from Bio import AlignIO
import os
import random


'''
Given the amino acid sequences of a protein domain and the nucleotide sequence of a gene or genome, it find the corresponding nucleotide sequence
It returns True if the matching is perfect and the nucleotide sequence corresponding to the domain is written in the output_file 
'''
def extract_seq_nucleotides(aa_file, nt_file, output_file, frame = 0, keep = False):
    
    
    with open(aa_file) as aa_handler:
        aa_seq = list(AlignIO.parse(aa_handler, "fasta"))[0]
        aa_record = aa_seq[0]
        
    with open(nt_file) as nt_handler:
        nt_seq = list(AlignIO.parse(nt_handler, "fasta"))[0]
        nt_record = nt_seq[0]
        if frame > 0:
            nt_record.seq = nt_record.seq[frame:]
                
    nt_translated = nt_record.translate()
    
    ali_file_in = "temp_ali_in" + str(random.random())
    ali_file_out = "temp_ali_out" + str(random.random())
    with open(ali_file_in, 'w') as out_handler:
        print("> Protein domain amino acid sequence\n" + aa_record.seq, file=out_handler)
        # We have to replace stop codons or will lose the numeric correspondence (mafft ignores them)
        # They are replace by Tryptophane as it is the less common amino acid and more unlike to create any problem
        # It should not be problems as long as the nucleotide sequence codifies exactly the amino acid sequence,
        # as we expect
        print("> Translated nucleotides sequence\n" + str(nt_translated.seq).replace('*', 'X'), file=out_handler)    
    
    os.system("mafft --auto " + ali_file_in + " > " + ali_file_out)
    
    # Load resulting alignment
    with open(ali_file_out) as aa_handler:
        ali_seq = list(AlignIO.parse(aa_handler, "fasta"))[0]
        ali_seq_aa = ali_seq[0].seq
        ali_seq_nt = ali_seq[1].seq
        if not keep:
            os.system("rm " + ali_file_in + " " + ali_file_out)
        
    # The nucleotide sequence should codify for amino acid sequence,
    # so we expect all gaps except for the codifying region, check this
    
#     length = len(ali_record_aa.seq)
    chars_aa = list(ali_seq_aa)
    chars_nt = list(ali_seq_nt)
    cod_region = False
    start = 0
    for i in range(len(chars_aa)):
        if chars_aa[i] == '-':
            if cod_region:
                end = i
                break
            continue
        elif chars_aa[i] != chars_nt[i]:
#             print("WARNING: The translated nucleotide sequence is not identical to the amino acid sequence. "
#             + "Check alignment file: " + ali_file_out)
            return False
        elif chars_aa[i] == chars_nt[i]:
            if not cod_region:
                start = i
                cod_region = True
            continue
    
    target_len = len(aa_record.seq)
    if end - start == target_len:
        with open(output_file, 'w') as out_hand:
            print(">", aa_record.id, file=out_hand)
            print(nt_record.seq[start*3:end*3], file=out_hand)
        
        return True
    else:
        
        return False
        

    

if __name__ == '__main__':
    
    parser=argparse.ArgumentParser(description='Given a sequence of amino acids of a protein domain and sequence of nucleotides (of a gene or genome), it extract the sequence of nucleotides that codify for the sequence of amino acids')
       
    parser = argparse.ArgumentParser()
    parser.add_argument("-aa", "--aa_file", type=str, help="Input file with the amino acid sequence in FASTA format", required=True)
    parser.add_argument("-nt", "--nt_file", type=str, help="Input file with the nucleotide sequence in FASTA format", required=True)
    parser.add_argument("-f", "--frame", type=int, help="Frame (offset) for reading the nucleotide sequence", default = 0)
    parser.add_argument("-k", "--keep", dest='keep', default=False, action='store_true', help="Do not remove the alignment between the translated sequence and amino acid sequence")
     
    parser.add_argument("-o","--output_file",  help="Output file", type=str, required=True)
  
    args = parser.parse_args()
     
    aa_file = args.aa_file
    nt_file = args.nt_file
    output_file = args.output_file
    keep = args.keep
    frame = args.frame
    
    # Test
#     aa_file = "/media/sakura_6TB/projects/sars_cov2/data/pipeline/03_2021/sequences/domains/Spike.bCoV_S1_RBD.fasta"
#     nt_file = "/home/fenix/projects/sars_cov2/temp/complete_genome_nts.frame_2.fasta"
#     frame = 0
#     output_file = "/home/fenix/projects/sars_cov2/temp/test"
    
    success = extract_seq_nucleotides(aa_file, nt_file, output_file, frame, keep)
    if success:
        print("Nucleotide sequence successfully recovered")
    else: 
        print("ERROR: The translated nucleotide sequence was not completely recovered. "
              + "Check alignment file")






