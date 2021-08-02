#!/usr/bin/env python

import modeling.utils as utils #@UnresolvedImport
import modeling.alignment as alignment #@UnresolvedImport
import modeling.position_average as average #@UnresolvedImport

# import re
import sys
# from statistics import mean 
import argparse
import os.path
import subprocess
from shutil import copyfile
from Bio import SeqIO
# import numpy as np
# import pandas as pd
 


def make_dir(path):
    try:
        os.makedirs(path, exist_ok=True)
    except OSError:
        print("Creation of the directory %s failed" % path)
        exit()



def generate_alignment(id_protein, input_fasta, output_dir, iterations, database_path, num_cores=1, include_refseq=False):
    
    # Generate jackhmer alignment 
    out_root_dir = output_dir + "/alignment/raw" 
    make_dir(out_root_dir)
    out_root_ali = out_root_dir + "/" + id_protein
    alignment.generate_alignment(input_fasta, iterations, database_path, out_root_ali, num_cores="1")
     
    # Uncompress the target alignment
#     ali_file_compress = out_root_ali + ".a2m.xz"
    ali_file_uncompress = out_root_ali + ".a2m"
#     subprocess.run(["unxz", ali_file_compress])
    
    # Remove insertion respect to the reference sequence
    temp_dir = output_dir + "/alignment/clean"
    make_dir(temp_dir)
    ali_file_no_gaps = temp_dir + "/" + id_protein + ".no_gaps.a2m" 
    subprocess.run(["reformat.pl", "a2m", "fas", ali_file_uncompress, ali_file_no_gaps, "-M", "first", "-r", "-l", "30000"])
    
    # Clean alignment
    filtered_ali = temp_dir + "/" + id_protein + ".filtered.a2m"
    utils.filterMSA(input_fasta, ali_file_no_gaps, filtered_ali, include_refseq)
    no_identical_ali = temp_dir + "/" + id_protein + ".no_identical.a2m"
    alignment.remove_identical_sequences(filtered_ali, no_identical_ali)
    
    # Copy final alignment of distant sequences
    out_ali_file = output_dir + "/alignment/" + id_protein + ".a2m"
    copyfile(no_identical_ali, out_ali_file)
    

def main(id_protein, input_fasta, input_ali, output_dir, database_path, cpus = 1, iterations = 3, mode = 'all', input_nt_fasta = ''):
    # Test
#     input_fasta = "/home/fenix/projects/sars_cov2/data/pipeline/02_2021/alignments/sequences/domains/Spike_S1.bCoV_S1_RBD.fasta"
#     id_protein = "spike.RBD"
#     output_dir = "test"
#     input_ali = "/home/fenix/projects/sars_cov2/src/covmut/test/alignment/spike.RBD.a2m"
    
    plmDCA_inference_path = os.path.dirname(__file__) + "/modeling/plmDCA_inference.jl"
    single_muts_path = os.path.dirname(__file__) + "/modeling/single_muts.jl"
        
    # Config
#     iterations = "3"
#     num_cores = "1"
#     database_path = "/data2/teams/weigt/juan/sars_cov2/data/input/databases/uniref90.vipr.ncbi_virus.vv_mers.fasta"
#     database_path = "/media/sakura_6TB/projects/sars_cov2/data/input/sequences/databases/viprbrc.09_2020.fasta"
#     database_path = "/data2/teams/weigt/juan/databases/uniprot/uniref100.fasta"
#     gisaid_path = ""
#     input_fasta = "/data2/teams/weigt/juan/sars_cov2/src/covmut/modeling/spike/Spike.fasta"
#     output_root = "."
#     id_protein = "spike"

    if not os.path.isdir(output_dir):
        print("Output dir: " + output_dir + " does not exist. Creating output directory")
        make_dir(output_dir)

    # Generate the alignment if it is not given, skips if it's already present
    out_ali_file = output_dir + "/alignment/" + id_protein + ".a2m"
    if not os.path.isfile(out_ali_file):
        if not input_ali:
            generate_alignment(id_protein, input_fasta, output_dir, iterations, database_path)
        else:
            make_dir(output_dir + "/alignment/")
            copyfile(input_ali, out_ali_file)
    
    # Check the alignment exists and there are sequences in the alignment
    if not os.path.exists(out_ali_file):
        sys.exit("The alignment is not found, the predictions cannot be computed")
    elif os.stat(out_ali_file).st_size == 0:
        sys.exit("The alignment is empty, not homologs were found. The predictions cannot be computed")  
    
    
    # Compute the model
    print("Computing model")
    model_dir = output_dir + "/model/"
    make_dir(model_dir)
    stdout_file = model_dir + id_protein + ".out"
    outroot_file = model_dir + id_protein
#     if not os.path.exists(outroot_file + ".h.gz") or not os.path.exists(outroot_file + ".J.gz"):
    if not os.path.exists(outroot_file + ".h") or not os.path.exists(outroot_file + ".J"):
        with open(stdout_file, 'w') as fh:
            subprocess.run([plmDCA_inference_path, "-i", out_ali_file, "-o", outroot_file], stdout=fh)

    # TODO: Internally in python, avoiding calling the julia script
    # Generate delta energy all single mutations
    # Uncompress parameters files
    print("Computing predictions")
    hfile_comp = outroot_file + ".h.gz"
    Jfile_comp = outroot_file + ".J.gz"
#     subprocess.run(["gzip", "-d", hfile_comp, Jfile_comp])
    outfile_delta = output_dir + "/" + id_protein + ".single_muts"
    # Compute single mutations
    if not os.path.exists(outfile_delta):
        subprocess.run([single_muts_path, "-r", input_fasta, "-i", out_ali_file, "-p", outroot_file, "-o", outfile_delta])
    # Delete temporal uncompressed files
    hfile_uncomp = outroot_file + ".h"
    Jfile_uncomp = outroot_file + ".J"
#     subprocess.run(["rm", hfile_uncomp, Jfile_uncomp])
#     subprocess.run(["gzip", "-f", hfile_uncomp, Jfile_uncomp])
    
    # Read reference
    records_ref = list(SeqIO.parse(open(input_fasta,'r'), "fasta"))
    ref_seq = str(records_ref[0].seq)
    ref_name = str(records_ref[0].id)
    
    # Read DCA parameters
#     npz_h = outroot_file + ".h.npz.tar.gz"
#     npz_J = outroot_file + ".J.npz.tar.gz"
#     h,J = utils.read_dca_par(npz_h, npz_J)
    
    # Out final prediction
    # Site entropy, only conservation
#     site_entropy = np.round(utils.compute_entropy_context_ind(out_ali_file),4)
    # Local entropy, conservation + epistasis
#     local_entropy = np.round(utils.compute_entropy_context_dep(ref_seq, h,J),4)
    
    
    averaged = True
#     averaged = False
    # Average single mutations effect over positions, giving a delta energy (effect mutations) for each position
    if mode == 'all':
        delta_position_dca = average.average_all_variants(outfile_delta, 6, averaged)
        delta_position_ind = average.average_all_variants(outfile_delta, 4, averaged)
    elif mode == 'reachable':
        delta_position_dca = average.average_codons(input_fasta, input_nt_fasta, outfile_delta, 6, mode, averaged)
        delta_position_ind = average.average_codons(input_fasta, input_nt_fasta, outfile_delta, 4, mode, averaged)
    elif mode == 'weighted':
        delta_position_dca = average.average_codons(input_fasta, input_nt_fasta, outfile_delta, 6, mode, averaged)
        delta_position_ind = average.average_codons(input_fasta, input_nt_fasta, outfile_delta, 4, mode, averaged)
    else:
        sys.exit("ERROR: The mode " + mode + " is not recognized. The available options are: all, reachable or weighted")
    
    index = [i for i in range(1,len(ref_seq)+1)]
    print("Writing predictions")
    output_file = output_dir + "/" + id_protein + ".covmut"
    with open(output_file, 'w') as fh:
#         print("position\tsite_entropy\tlocal_entropy\tdelta_IND\tdelta_DCA",file=fh)
#         for i in range(len(ref_seq)):
#             print(str(index[i]) + "\t" + str(site_entropy[i]) + "\t" + str(local_entropy[i]) + "\t" + str(delta_position_ind[i]) + "\t" + str(delta_position_dca[i]), file=fh)
        
        print("position\tdelta_IND\tdelta_DCA",file=fh)
        for i in range(len(ref_seq)):
            print(str(index[i]) + "\t" + str(delta_position_ind[i]) + "\t" + str(delta_position_dca[i]), file=fh)

    

     


if __name__ == '__main__':
    
    parser=argparse.ArgumentParser(description='Given a sequence in FASTA format generate all covmut predictions, mutability by position and effect of single mutations')
      
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_fasta", type=str, help="Input amino acid sequence in FASTA format", required=True)
    parser.add_argument("-id", "--id_protein", type=str, help="Identifier of the protein. It will be used to name the output files", required=True)
    parser.add_argument("-db", "--database", type=str, help="Path to the sequence database in FASTA format", required = True)
    parser.add_argument("-it", "--iterations", type=str, help="Options. Number of jackhmmer iterations to perform", default = 3)
    parser.add_argument("-ali", "--input_ali", type=str, help="Optional. Multiple sequence alignment in FASTA to build the model. If it is not given, it is generated", default = "")
    parser.add_argument("-m", "--mode", type=str, help="Optional. Select the mode to make the position average for each position from single mutants predictions. Options: all (average of 19 possible variants), reachable (average from reachable amino acids with one mutation), weighted (average from reachable amino acids with one mutation weigthed by the number of codons that codified each amino acid). Reachable and weighted require the nucleotide sequence. Default: all.", default = "all")
    parser.add_argument("-nt", "--input_nt_fasta", type=str, help="Optional. Input nucleotide sequence in FASTA format", default = "None")
    parser.add_argument("-cpus", "--num_cpus", type=str, help="Optional. Number of CPUs to use during jackhmmer searches", default = 1)    
    parser.add_argument("-o","--output_dir",  help="Output directory for the output files", type=str, required=True)
 
    args = parser.parse_args()
    
    id_protein = args.id_protein
    input_fasta = args.input_fasta
    database_path = args.database
    input_ali = args.input_ali
    iterations= str(args.iterations)
    mode = args.mode
    input_nt_fasta = args.input_nt_fasta
    num_cpus = str(args.num_cpus)
    
    # Use absolute paths for output_cd dir, reformat.pl seems to have a bug which is problematic with relative paths
    output_dir = os.path.abspath(args.output_dir)
    
    # Check input
    if not os.path.isfile(input_fasta):
        sys.exit("ERROR: The input FASTA file does not exist")
        
    if mode == 'all' or mode == 'All':
        pass
    elif mode == 'reachable' or mode == 'Reachable':
        if not os.path.isfile(input_nt_fasta):
            sys.exit("ERROR: With the strategy " + mode + " a FASTA file with the nucleotide sequence has to be provided")
    elif mode == 'weighted' or mode == 'Weighted':
        if not os.path.isfile(input_nt_fasta):
            sys.exit("ERROR: With the strategy " + mode + " a FASTA file with the nucleotide sequence has to be provided")
    else:
        sys.exit("ERROR: The mode " + mode + " is not recognized. The available options are: all, reachable or weighted")
     
       
    main(id_protein, input_fasta, input_ali, output_dir, database_path, num_cpus, iterations, mode, input_nt_fasta)
    

