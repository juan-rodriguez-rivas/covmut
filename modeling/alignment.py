#!/usr/bin/env python

import subprocess
from Bio import SeqIO
from Bio import AlignIO
from glob import glob
import os
import argparse
from numpy import zeros

def compare_seqs(str1, str2, discard_gaps=True):
    for i in range(len(str1)):
        if(discard_gaps):
            if(str1[i] == '-' and str2[i] != '-'):
                continue
            elif(str2[i] == '-' and str1[i] != '-'):
                continue
        if(str1[i] != str2[i]):
            return False
    return True
        

# Remove identical sequences (keeping one copy)
# Remove sequences than have been observed fewer than min_observation times
def remove_identical_sequences(input_fasta, output_fasta, min_observations = 1):
    with open(input_fasta, "r") as input_handle:
        ali = list(SeqIO.parse(input_handle, "fasta"))
        tam = len(ali)
        strs = []
        for i in range(tam):
            strs.append(str(ali[i].seq))
        
        index_to_keep = []
        index_to_check = list(range(0,len(ali))) 
        
        while(index_to_check):
            
            for i in index_to_check:
                count_observations = 1
                index_to_check.remove(i)
                
                index_to_remove = []
                for j in index_to_check:
                    if(strs[i] == strs[j]):
                        count_observations = count_observations + 1
                        index_to_remove.append(j)
                
                for j in index_to_remove:
                    index_to_check.remove(j)
                
                if count_observations >= min_observations:
                    index_to_keep.append(i)
            
        
        out_ali = []
        for i in index_to_keep:
            out_ali.append(ali[i])
        
        with open(output_fasta, 'w') as output_handle:
            SeqIO.write(out_ali, output_handle, "fasta-2line")
              

def generate_alignment(input_fasta, iterations, database_path, output_root, num_cores="1"):
    
    filename = os.path.basename(input_fasta)
    
#     if os.path.isfile(output_root + ".a2m.xz"):
    if os.path.isfile(output_root + ".a2m"):
        print("Alignment for " + filename + " already computed, skipping")
    else:
        print("Genereting alignment for " + filename)
        subprocess.run(["jackhmmer", "-o", "/dev/null", "--cpu", num_cores, "-N", iterations, "--domtblout", output_root + ".domtblout",
                        "--chkali", output_root, "--chkhmm", output_root, "-A", output_root + ".sto", 
                        input_fasta, database_path])
        
        sto_to_fasta(output_root + ".sto", output_root + ".a2m")
            
        for i in range(1,int(iterations)+1):
            input_sto = output_root + "-" + str(i) + ".sto"
            if os.path.isfile(input_sto):
                sto_to_fasta(input_sto, output_root + "-" + str(i) + ".a2m")
        
#         subprocess.run(["xz"] + glob(output_root + "*"))
        print("Alignment for " + filename + " generated")
        
    
    
def sto_to_fasta(input_file, output_file):
        with open(input_file, "r") as input_handle:
            ali = AlignIO.parse(input_handle, "stockholm")
            AlignIO.write(ali, output_file, "fasta-2line")    
    
    

if __name__ == '__main__':
    
    parser=argparse.ArgumentParser(description='Generate an alignment using jackhmmer starting with the given sequence. It stores all intermediate iterations. All output files are compressed using xz program')
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fasta", type=str, help="Input sequence in FASTA format", required=True)
    parser.add_argument("-it", "--iterations", type=str, help="The number of jackhmmer iterations. By default, 5 iterations", default = 5)
    parser.add_argument("-db", "--database_path", help="The sequence database to be used in the search", type=str, required=True)
    parser.add_argument("-o","--output_root",  help="Output root path for the output files", type=str, required=True)
    parser.add_argument("-c", "--num_cores", help="Number of cores to used in the jackhmmer search. By default, 1 core is used", type=str, default = "1")
    args = parser.parse_args()
    
    
    generate_alignment(args.input_fasta, str(args.iterations), args.database_path, args.output_root, str(args.num_cores))
    
    
    
