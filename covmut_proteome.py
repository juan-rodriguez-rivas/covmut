#!/usr/bin/env python


import modeling.alignment as alignment  # @UnresolvedImport
import modeling.proteome as prot  # @UnresolvedImport
import modeling.utils as utils # @UnresolvedImport
import covmut_sequence as covmut_sequence

# from shutil import copyfile
import sys
from Bio import AlignIO
import subprocess
# import os
import os.path
import re
import argparse





def compute_distant_alignment(protein_id, domain_id, iterations, database_path, num_cores, output_dir):
    
    both_id = protein_id + "." + domain_id
    input_fasta_dir = output_dir + "/sequences"
    input_domain_fasta = input_fasta_dir + "/domains/" + both_id + ".fasta"
    
    # Final alignment file, skip if it is already computed
    ali_dir = output_dir + "/alignments/distant"
    covmut_sequence.make_dir(ali_dir)
    output_ali_file = output_dir + "/alignments/" + both_id + ".only_distant.a2m"
    if os.path.exists(output_ali_file):
        return
    
    ali_proteins_dir = ali_dir + "/proteins"
    ali_domains_dir = ali_dir + "/domains"
    covmut_sequence.make_dir(ali_proteins_dir)
    covmut_sequence.make_dir(ali_domains_dir)
    
    # Compute the alignment for the protein
    ali_prot_root = ali_proteins_dir + "/" + protein_id
    ali_prot = ali_prot_root + ".a2m"
    ali_prot_ungapped = ali_prot_root + ".no_gaps.a2m"
    if protein.length < 2000 and not os.path.exists(ali_prot_ungapped):
        input_fasta = input_fasta_dir + "/proteins/" + protein_id + ".fasta"
        alignment.generate_alignment(input_fasta, iterations, database_path, ali_prot_root, num_cores)
        # Remove gaps (insertions) compare to the reference sequence
        os.system("reformat.pl a2m fas " + ali_prot + " " + ali_prot_ungapped + " -M first -r -l 30000")
    
    
    # Compute the alignment for the domain
    root_domain = ali_domains_dir + "/" + both_id
    ali_domain = root_domain + ".a2m"
    
    if not os.path.exists(ali_domain):
        alignment.generate_alignment(input_domain_fasta, iterations, database_path, root_domain, num_cores)
    
    if os.stat(ali_domain).st_size == 0:
        sys.exit("The alignment is empty, not homologs were found. The predictions cannot be computed")
   
    # Clean alignment and run hhfilter to find out the number of effective sequences
    clean_domain_ali = root_domain + ".clean.a2m"
    hhfilter_domain_ali = root_domain + ".clean.id80.cov80.a2m"
    if not os.path.exists(hhfilter_domain_ali):
        clean_alignment(input_domain_fasta, root_domain, clean_domain_ali)
        subprocess.run(["hhfilter", "-i", clean_domain_ali, "-o", hhfilter_domain_ali, "-id", "80", "-cov", "80" ,"-M" ,"first", "-maxseq", "10000000"])


    domain_ali_object = list(AlignIO.parse(hhfilter_domain_ali, "fasta-2line"))[0]
   
    # Trim protein alignment to the domain
    ali_proteins_domains_dir = ali_proteins_dir + "/domains"
    covmut_sequence.make_dir(ali_proteins_domains_dir)
    root_prot_domain = ali_proteins_domains_dir + "/" + both_id
    clean_prot_domain_ali = root_prot_domain + ".clean.a2m"
    hhfilter_prot_domain_ali = root_prot_domain + ".clean.id80.cov80.a2m"
    if os.path.exists(ali_prot_ungapped):
        with open(ali_prot_ungapped, "r") as input_handle:
            ali_prot = list(AlignIO.parse(input_handle, "fasta-2line"))[0]
            ali_prot_domain = ali_prot[:,(domain.start-1):(domain.end)]
           
            ali_prot_domain_file = root_prot_domain + ".a2m" 
            with open(ali_prot_domain_file, 'w') as output_handle:
                # Substitute the id of the protein for the domain of the reference sequence
                # It includes which domain we are working with
                ali_prot_domain[0].id = domain_ali_object[0].id
                AlignIO.write(ali_prot_domain, output_handle, "fasta-2line")
           
            clean_alignment(input_domain_fasta, root_prot_domain, clean_prot_domain_ali)
            subprocess.run(["hhfilter", "-i", clean_prot_domain_ali, "-o", hhfilter_prot_domain_ali, "-id", "80", "-cov", "80" ,"-M" ,"first", "-maxseq", "10000000"]) 
           
        # Load alignment and compare number of sequences
        num_seqs_domain_ali = len(list(AlignIO.parse(hhfilter_domain_ali, "fasta-2line"))[0])
        num_seqs_prot_domain_ali = len(list(AlignIO.parse(hhfilter_prot_domain_ali, "fasta-2line"))[0])
       
        # Select the one with more effective sequences, remove reference sequence
        if num_seqs_prot_domain_ali > num_seqs_domain_ali:
            with open(clean_prot_domain_ali, "r") as input_handle:
                ali = list(AlignIO.parse(input_handle, "fasta-2line"))[0]
                ali_out = ali[1:,:]
                with open(output_ali_file, "w") as output_handle:
                    AlignIO.write(ali_out, output_handle, "fasta-2line")
        else:
            with open(clean_domain_ali, "r") as input_handle:
                ali = list(AlignIO.parse(input_handle, "fasta-2line"))[0]
                ali_out = ali[1:,:]
                with open(output_ali_file, "w") as output_handle:
                    AlignIO.write(ali_out, output_handle, "fasta-2line")
           
    else: # If there is not protein alignment, use the domain alignment
        with open(clean_domain_ali, "r") as input_handle:
            ali = list(AlignIO.parse(input_handle, "fasta-2line"))[0]
            ali_out = ali[1:,:]
            with open(output_ali_file, "w") as output_handle:
                AlignIO.write(ali_out, output_handle, "fasta-2line")



def compute_alignment_including_close(protein_id, domain_id, gisaid_path, num_cores, coverage, output_dir):

    both_id = protein_id + "." + domain_id

    # Generate GISAID alignments
    output_close_dir = output_dir + "/alignments/close"
    covmut_sequence.make_dir(output_close_dir)
            
    final_output_ali = output_dir + "/alignments/" + both_id + ".including_close.a2m"
    if os.path.exists(final_output_ali):
        return
            
    # Close sequences, 1 iteration is enough
    input_fasta = output_dir + "/sequences/domains/" + both_id + ".fasta"
    
    out_root =  output_close_dir + "/" + both_id
    alignment.generate_alignment(input_fasta, "1", gisaid_path, out_root, num_cores)
    out_file = out_root + ".a2m"
    
    clean_file = out_root + ".no_gaps.cov" + str(coverage) + ".clean.a2m"
    if not os.path.exists(clean_file):
        no_gaps_file = out_root + ".no_gaps.a2m"
        os.system("reformat.pl a2m fas " + out_file + " " + no_gaps_file + " -M first -r -l 30000")
        hhfilter_file = out_root + ".no_gaps.cov" + str(coverage) + ".a2m"
        subprocess.run(["hhfilter", "-i", no_gaps_file, "-o", hhfilter_file, "-id", "100", "-cov", str(coverage),"-M" ,"first", "-maxseq", "10000000"])
        utils.filterMSA_gisaid(input_fasta, hhfilter_file, clean_file)
    
    no_dups_file = out_root + ".no_gaps.cov" + str(coverage) + ".clean.no_dups.a2m"
    alignment.remove_identical_sequences(clean_file, no_dups_file, 1)
    
    # Load and concatenate the distant alingment and the gisaid alignment, remove duplicates adn write
    ali_distant_file = output_dir + "/alignments/" + both_id + ".only_distant.a2m"
    combined_ali_file = out_root + ".no_gaps.cov" + str(coverage) + ".clean.no_dups.combined.a2m"
    
    with open(combined_ali_file, 'w') as fh:
        subprocess.run(["cat", no_dups_file, ali_distant_file], stdout=fh)
    
    combined_ali = list(AlignIO.parse(combined_ali_file, "fasta-2line"))[0]
    with open(combined_ali_file, 'w') as output_fh:
        AlignIO.write(combined_ali, output_fh, "fasta-2line")
    
    alignment.remove_identical_sequences(combined_ali_file, final_output_ali)
            
            


def run_prediction(protein_id, domain_id, input_fasta, ali_file, output_preds_dir, database_path, cpus, iterations, mode, input_nt_fasta):
    both_id = protein_id + "." + domain_id
    output_file = output_preds_dir + "/" + both_id + ".covmut"
    
    if os.path.exists(output_file):
        return
        
    if not os.path.exists(ali_file):
        print("Alignment " + ali_file + " not found. Skipping running predictions")
        return
            
    # Run prediction
    covmut_sequence.main(both_id, input_fasta, ali_file, output_preds_dir, database_path, cpus, iterations, mode, input_nt_fasta)


def clean_alignment(input_fasta, input_root_ali, output_file):
    
    if os.path.exists(output_file):
        return
     
    # Uncompress the target alignment
    ali_file_uncompress = input_root_ali + ".a2m"
    
    # Remove insertion respect to the reference sequence
    ali_file_no_gaps = input_root_ali + ".no_gaps.a2m"
    os.system("reformat.pl a2m fas " + ali_file_uncompress + " " + ali_file_no_gaps + " -M first -r -l 30000") 
#     subprocess.run(["reformat.pl", "a2m", "fas", ali_file_uncompress, ali_file_no_gaps, "-M", "first", "-r", "-l", "30000"])
        
    # Clean alignment
    filtered_ali = input_root_ali + ".filtered.a2m"
    utils.filterMSA(input_fasta, ali_file_no_gaps, filtered_ali)
    alignment.remove_identical_sequences(filtered_ali, output_file)


def load_domain_predictions(input_file, protein, domain, type_ali):
           
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            match = re.search("^position", line)
            if match:
                continue
            
            line = line.replace('\n', '')
            fields = line.split("\t")
            
            pos = int(fields[0])-1
#             s_ent = fields[1]
#             l_ent = fields[2]
#             d_ind = fields[3]
#             d_dca = fields[4]
            d_ind = fields[1]
            d_dca = fields[2]
            
            
            protein_pos = pos + domain.start-1
            protein.pos_domain[protein_pos] = domain.id
            
            if(type_ali == 'distant'):
                if protein.distant_site_entropy[protein_pos] != 'NA':
                    print("Warning, position " + str(protein_pos) + " in domain " + domain.id + " of the protein " + protein.id + " already predicted. There are overlaps between domains")
                    continue
#                 protein.distant_site_entropy[protein_pos] = s_ent
#                 protein.distant_local_entropy[protein_pos] = l_ent
                protein.distant_score_IND[protein_pos] = d_ind
                protein.distant_score_DCA[protein_pos] = d_dca
            elif(type_ali == 'close'):
                if protein.close_site_entropy[protein_pos] != 'NA':
                    print("Warning, position " + str(protein_pos) + " in domain " + domain.id + " of the protein " + protein.id + " already predicted. There are overlaps between domains")
                    continue
#                 protein.close_site_entropy[protein_pos] = s_ent
#                 protein.close_local_entropy[protein_pos] = l_ent
                protein.close_score_IND[protein_pos] = d_ind
                protein.close_score_DCA[protein_pos] = d_dca
            else:
                print("Warning: unrecognized type of alignment. Received: " + type_ali + ". Expected: distant or close")
            

def load_mutability(input_file, protein, domain):
    
    mutability = []
    
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            match = re.search("^position", line)
            if match:
                continue
            
            line = line.replace('\n', '')
            fields = line.split("\t")
            
#             pos = int(fields[0])-1
            num_muts = fields[3]
            mutability.append(num_muts)
            
    return mutability


def collapse_predictions(proteome):
    
    # Load predictions, collapse domains prediction into proteins
    for protein_id in proteome.proteins_id:
        protein = proteome.proteins_id[protein_id]
        
        for domain_id in protein.domains_id:
#             domain = protein.domains_id[domain_id]
            both_id = protein_id + "." + domain_id
            domain = protein.domains_id[domain_id]
            
            domain_preds_file = output_dir + "/predictions/distant/" + both_id + ".covmut"
            load_domain_predictions(domain_preds_file, protein, domain, 'distant')
            domain_preds_file = output_dir + "/predictions/close/" + both_id + ".covmut"
            load_domain_predictions(domain_preds_file, protein, domain, 'close')
            

if __name__ == '__main__':
    
    parser=argparse.ArgumentParser(description='Given a proteome in FASTA format and file with a set of HMM profiles (e.g. pfam), generate predictions of the mutability of each position.')
      
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--fasta_proteome", type=str, help="Input proteome in FASTA format", required=True)
    parser.add_argument("-hmm", "--hmm_file", type=str, help="The number of jackhmmer iterations. By default, 5 iterations", required=True)
    parser.add_argument("-db", "--database", type=str, help="Path to the sequence database in FASTA format", required = True)
    parser.add_argument("-g", "--database_gisaid", type=str, help="Path to the GISAID sequence database in FASTA format", required = True)
    parser.add_argument("-it", "--iterations", type=str, help="Options. Number of jackhmmer iterations to perform", default = 3)
    parser.add_argument("-m", "--mode", type=str, help="Optional. Select the mode to make the position average for each position from single mutants predictions. Options: all (average of 19 possible variants), reachable (average from reachable amino acids with one mutation), weighted (average from reachable amino acids with one mutation weigthed by the number of codons that codified each amino acid). Reachable and weighted require the nucleotide sequence. Default: all", default = "all")
    parser.add_argument("-cpus", "--num_cpus", type=str, help="Optional. Number of CPUs to use during jackhmmer searches", default = 1)
    parser.add_argument("-nt", "--nucleotide_genome", type=str, help="Optional. Input nucleotide genome sequence in FASTA format. Required if modes reachable or weighted are used", default = "None")
    parser.add_argument("-o","--output_dir",  help="Output directory for the output files", type=str, required=True)
    args = parser.parse_args()
  
    fasta_proteome = args.fasta_proteome
    hmm_file = args.hmm_file
    database_path = args.database
    iterations= str(args.iterations)
    mode = args.mode
    nucleotide_genome = args.nucleotide_genome
    num_cpus = str(args.num_cpus)
    gisaid_path = args.database_gisaid
    coverage_gisaid = 80
    
    # Use absolute paths for output_dir, reformat.pl seems to have a bug which is problematic with relative paths
    output_dir = os.path.abspath(args.output_dir)
    
    # Check input
    if not os.path.exists(fasta_proteome):
        sys.exit("Input proteome file not found. Please, provide a proteome file in FASTA format")
        
    if not os.path.exists(hmm_file):
        sys.exit("HMM file not found. Please, provide a HMM file")
        
    if not os.path.exists(database_path):
        sys.exit("Sequence database file not found. Please, provide a database file in FASTA format")
    
    if not os.path.exists(gisaid_path):
        sys.exit("GISAID sequence database file not found. Please, provide a database file in FASTA format")
        
    if mode == 'reachable' or mode == 'weighted':
        if nucleotide_genome == 'None':
            sys.exit("Nucleatide genome not provided. If modes reachable or weighted are used, a file containing the genome sequence in nucleotides have to be provide")
        elif not os.path.exists(nucleotide_genome):
            sys.exit("Nucleatide genome not found. If modes reachable or weighted are used, a file containing the genome sequence in nucleotides have to be provide")
        
    
#     # Generate alignments for proteins and domains, chose the one with highest number of effective sequences
#     fasta_proteome = "/home/fenix/projects/sars_cov2/data/genome/proteome.Wuhan-Hu-1.MN908947.fasta"
# #     fasta_proteome = "/home/fenix/projects/sars_cov2/data/pipeline/test/input/for_testing.proteome.fasta"
#     hmm_file = "/home/fenix/projects/sars_cov2/data/pipeline/test/input/all.hmm"
#     output_dir = "/media/sakura_6TB/projects/sars_cov2/data/pipeline/06_2021_min_obs1"
#     database_path = "/media/sakura_6TB/projects/sars_cov2/data/input/sequences/databases/viprbrc.09_2020.fasta"
#     gisaid_path = "/media/sakura_6TB/projects/sars_cov2/data/input/sequences/gisaid/01_2021/allprot0119.fasta/allprot0119.fasta"
#     iterations = "1"
#     num_cores = "1"
#     mutability = True
#     nucleotide_genome = "/home/fenix/projects/sars_cov2/data/genome/complete_genome_nts.fasta"
#     mode = "reachable"
#     coverage_gisaid = 80
    
    # Make output dir if it dos not exist
    if not os.path.isdir(output_dir):
        print("Output dir: " + output_dir + " does not exist. Creating output directory")
        covmut_sequence.make_dir(output_dir)

    # Obtain protein and domain sequences from the proteome and HMM profiles
    seqs_dir = output_dir + "/sequences"
    covmut_sequence.make_dir(seqs_dir)
    
    filename = os.path.basename(hmm_file)
    hmm_outfile = seqs_dir +  "/" + filename + ".domtblout" 
    proteome = prot.load_proteome(fasta_proteome, hmm_file, hmm_outfile)
    
    seqs_domains_dir = seqs_dir + "/domains"
    seqs_proteins_dir = seqs_dir + "/proteins"
    covmut_sequence.make_dir(seqs_domains_dir)
    covmut_sequence.make_dir(seqs_proteins_dir)
    proteome.write_domain_fastas(seqs_domains_dir)
    proteome.write_protein_fastas(seqs_proteins_dir)
    
    if not nucleotide_genome == "None":
        prot.add_nucleotide_sequence(proteome, seqs_domains_dir, nucleotide_genome)
    
    # Compute alignment (for sequences shorter than 2000 aas) for proteins and domains
    # Remove redundancy and select between domains and protein based alignments the one with more effective sequences
    ali_dir = output_dir + "/alignments"
    preds_dir = output_dir + "/predictions"
    covmut_sequence.make_dir(preds_dir)
    

    for protein_id in proteome.proteins_id:
        protein = proteome.proteins_id[protein_id]
        
        for domain_id in protein.domains_id:
            domain = protein.domains_id[domain_id]
            both_id = protein_id + "." + domain_id
            
            # Compute distant alignment
            
            compute_distant_alignment(protein_id, domain_id, iterations, database_path, num_cpus, output_dir)
            # Compute prediction for the distant alignment
            input_fasta = seqs_domains_dir + "/" + both_id + ".fasta"
            input_nt_fasta = seqs_domains_dir + "/" + both_id + ".nt.fasta"
            ali_file = ali_dir + "/" + both_id + ".only_distant.a2m"
            output_preds_dir = preds_dir + "/distant"
            covmut_sequence.make_dir(output_preds_dir)
#             run_prediction(protein_id, domain_id, input_fasta, ali_file, output_preds_dir, mode, input_nt_fasta)
            run_prediction(protein_id, domain_id, input_fasta, ali_file, output_preds_dir, database_path, num_cpus, iterations, mode, input_nt_fasta)
            
            # Compute alignment including close sequences
            compute_alignment_including_close(protein_id, domain_id, gisaid_path, num_cpus, coverage_gisaid, output_dir)
            # Compute prediction for ali including close seqs
            ali_file = ali_dir + "/" + both_id + ".including_close.a2m"
            output_preds_dir = preds_dir + "/close"
            covmut_sequence.make_dir(output_preds_dir)
#             run_prediction(protein_id, domain_id, input_fasta, ali_file, output_preds_dir, mode, input_nt_fasta)
            run_prediction(protein_id, domain_id, input_fasta, ali_file, output_preds_dir, database_path, num_cpus, iterations, mode, input_nt_fasta)
            
    # Collapse domain predictions into protein
    collapse_predictions(proteome)
                    
    # Write prediction files
    output_file_distant = output_dir + "/only_distant.covmut"
    
    
    with open(output_file_distant, 'w') as output_handle:
        print("protein", "domain", "position_protein", "position_domain", "aa_reference", "score_IND", "score_DCA", "mutability_may21", "mutability_december20", "mutability_july20", sep="\t", file=output_handle)
        for protein in proteome.proteins_pos:
            protein_id = protein.id
            for i in range(protein.length):
                domain_id = protein.pos_domain[i]
                pos_domain = 'NA'
                if domain_id != 'NA':
                    pos_domain = str(i + 2 - protein.domains_id[domain_id].start)
                    
                print(protein_id, protein.pos_domain[i], str(i+1), pos_domain, protein.seq[i], str(protein.distant_score_IND[i]), str(protein.distant_score_DCA[i]), sep="\t", file=output_handle)
                
    output_file_close = output_dir + "/including_close.covmut"        
    with open(output_file_close, 'w') as output_handle:
        print("protein", "domain", "position_protein", "position_domain", "aa_reference", "score_IND", "score_DCA", "mutability_may21", "mutability_december20", "mutability_july20", sep="\t", file=output_handle)
        for protein in proteome.proteins_pos:
            protein_id = protein.id
            for i in range(protein.length):
                domain_id = protein.pos_domain[i]
                pos_domain = 'NA'
                if domain_id != 'NA':
                    pos_domain = str(i + 2 - protein.domains_id[domain_id].start)

                print(protein_id, protein.pos_domain[i], str(i+1), pos_domain, protein.seq[i], str(protein.close_score_IND[i]), str(protein.close_score_DCA[i]), sep="\t", file=output_handle)
 
    
    