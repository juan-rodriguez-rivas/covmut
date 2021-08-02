#!/usr/bin/env python

import modeling.extract_ntseq as extract_ntseq #@UnresolvedImport

from Bio import SeqIO
import subprocess
import os.path
import math
from _ast import Or


''' From a list of domain, it removes overlapping domains (keeping the one with the best evalue)
 and sort them according to their positions in the protein
'''
def sort_domains(domains_list: list) -> list:
    
    overlapping = []
    domain_dict = {}
    for i in range(0, len(domains_list)):
        domain = domains_list[i]
        domain_acc = domain.id + "_" + str(domain.start) + "_" + str(domain.end)
        domain.acc = domain_acc
        domain.size = domain.end - domain.start
        
        domain_dict[domain_acc] = domain
        
    
    for i in range(0, len(domains_list)-1):
        domain1 = domains_list[i]
        start1 = domain1.start
        end1 = domain1.end
        
        if not domain1.acc in domain_dict:
            continue
        
        for j in range(i+1, len(domains_list)):
            domain2 = domains_list[j]
            start2 = domain2.start
            end2 = domain2.end
            
            if not domain2.acc in domain_dict:
                continue
            
            # Check overlap, if start1 or end1 is between start2 and end2, 
            # then there is overlap. In that case, discard the one with the
            # worse evalue
            if((start1 >= start2 and start1 <= end2) or
               (end1 >= start2 and end1 <= end2) or 
               (start1 >= start2 and end1 <= end2) or 
               (start2 >= start1 and end2 <= end1)):
#                 if domain2.evalue < domain1.evalue:
#                     overlapping.append(i)
#                 else:
#                     overlapping.append(j)
                # If the domain have a similar evalue (log10 diff < 3), use the size as criteria
                diff_log = math.log10(domain2.evalue) - math.log10(domain1.evalue)
                 
                if diff_log < -3:
                    # domain2 better, remove domain1
                    del domain_dict[domain1.acc]
                elif diff_log > 3:
                    # domain1 better, remove domain2
                    del domain_dict[domain2.acc]
                else:
                    if(domain1.size > domain2.size):
                        # domain1 larger, remove domain2
                        del domain_dict[domain2.acc]
                    else:
                        # domain2 larger, remove domain1
                        del domain_dict[domain1.acc]
    
#     # New list keeping only the non overlapping domains
#     unsorted = []
#     for i in range(0, len(domains_list)):
#         if not i in overlapping:
#             unsorted.append(domains_list[i])
    
    # Iterate over the list sorting elements till there is not element in unsorted
    unsorted = list(domain_dict.values())
    sorted_domains = []
    while len(unsorted) > 0:
        domain_to_sort = unsorted.pop(0)
        if not sorted_domains: # first element
            sorted_domains.append(domain_to_sort)
        else:
            for i in range(0, len(sorted_domains)):
                if domain1.start < domain2.start:
                    sorted_domains.insert(i,domain_to_sort)
                    break
                elif (i == len(sorted_domains)):
                    sorted_domains.append(domain_to_sort)
                    break
    
    return sorted_domains

def add_nucleotide_sequence(proteome, seq_dir, nucleotide_genome):
    
    for protein_id in proteome.proteins_id:
        protein = proteome.proteins_id[protein_id]
        for domain_id in protein.domains_id:
            both_id = protein_id + "." + domain_id
            aa_file = seq_dir + "/" + both_id + ".fasta"
            output_file = seq_dir + "/" + both_id + ".nt.fasta"
            
            if not os.path.isfile(output_file):
                if extract_ntseq.extract_seq_nucleotides(aa_file, nucleotide_genome, output_file, 0):
                    pass
                elif extract_ntseq.extract_seq_nucleotides(aa_file, nucleotide_genome, output_file, 1):
                    pass
                elif extract_ntseq.extract_seq_nucleotides(aa_file, nucleotide_genome, output_file, 2):
                    pass
                else:
                    print("The nucleotide sequence was not found for the domain " + domain_id + " of the protein " + protein_id)
            
    

def add_domains(hmm_outfile: str, proteins_id: str, evalue_threshold: float = 0.00001):
    
    with open(hmm_outfile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.split()
                identifier = fields[0]
                accession = fields[1]
                length = int(fields[5])
                evalue = float(fields[12])
                score = float(fields[13])
                bias = float(fields[14])
                start = int(fields[17])
                end = int(fields[18])
                
                protein_id = fields[3]
                
                if evalue < evalue_threshold:
                    seq = proteins_id[protein_id].seq[start-1:end]
                    domain = Domain(identifier, accession, evalue, score, bias, start, end, seq, length)
                    proteins_id[protein_id].domains_id[identifier] = domain
                    
    for protein_id in proteins_id:
        protein = proteins_id[protein_id]
        
        domains_sorted = sort_domains(list(protein.domains_id.values()))
        proteins_id[protein_id].domains_list = domains_sorted
        # Remove previous  dict of domains, it may contain overlapping domains
        domains_id = {}
        
        for i in range(0, len(domains_sorted)):
            domain = domains_sorted[i]
            domains_id[domain.id] = domain
        
        proteins_id[protein_id].domains_id = domains_id
                    


class Domain:
    """
    It stores one domain

    Attributes
    ----------
    id: str
        Identifier of the domain
    acc: str
        Accession of the domain
    record: SeqRecord
        SeqRecord object of the domain
    seq: str
        Sequence of the domain
    length: int
        Length of the domain
    
    position: int
        Position of the protein in the proteome
    domain_id: dict
        Dict containing the domains of the protein
    domain_list:
        List containing the domains of the protein, sorte by the positions of the domain in the protein    
    """
    def __init__(self, identifier, accession, evalue, score, bias, start, end, seq, length):
        self.id = identifier
        self.acc = accession
        self.length = length
        self.evalue = evalue
        self.score = score
        self.bias = bias
        self.start = start
        self.end = end
        self.seq = seq
        

class Protein:
    """
    It stores one protein

    Attributes
    ----------
    id: str
        Identifier of the protein
    record: SeqRecord
        SeqRecord object of the protein
    seq: str
        Sequence of the protein
    length: int
        Length of the sequence
    position: int
        Position of the protein in the proteome
    domain_id: dict
        Dict containing the domains of the protein
    domain_list:
        List containing the domains of the protein, sorte by the positions of the domain in the protein    
    """
    def __init__(self, record, position, identifier, seq):
        self.record = record
        self.seq = seq
        self.length = len(seq)
        self.position = position
        self.id = identifier
        
        self.domains_id = {}
        self.domains_list = []
        
        # Store the domain in each positions
        self.pos_domain = ["NA"]*self.length
        
        # Store prediction for each position using the distant alingment
        self.distant_site_entropy = ["NA"]*self.length
        self.distant_local_entropy = ["NA"]*self.length
        self.distant_score_IND = ["NA"]*self.length
        self.distant_score_DCA = ["NA"]*self.length
        
        # Store prediction for each position using the alignmkent including close sequences
        self.close_site_entropy = ["NA"]*self.length
        self.close_local_entropy = ["NA"]*self.length
        self.close_score_IND = ["NA"]*self.length
        self.close_score_DCA = ["NA"]*self.length



class Proteome:    
    """
    Load from a fasta file and HMM profiels the set of proteins and domains in that proteome

    Attributes
    ----------
    proteins_pos : list
        List of proteins
    proteins_id : Dict
        Dict of proteins

    """
    def __init__(self, fasta_file, hmm_file, hmm_outfile, evalue_threshold = 0.00001):
        
        self.proteins_pos = []
        self.proteins_id = {}
        position = 1
        
        with open(fasta_file, "r") as input_handle:
            
            for record in SeqIO.parse(input_handle, "fasta"):
              
                protein = Protein(record, position, record.id, str(record.seq))
                self.proteins_pos.append(protein)
                self.proteins_id[record.id] = protein
                                
                position += 1
        
        # If hmm_outfile (domtblout type of HMMER output files) exists, load domains from that file
        # Otherwise, generate the output file using hmmscan and load domains
        if not os.path.isfile(hmm_outfile):
#             hmm_outfile = hmm_file + ".domtblout"
            hmm_pressed = hmm_file + ".h3m"
            if not os.path.isfile(hmm_pressed):
                subprocess.run(["hmmpress", hmm_file])
            subprocess.run(["hmmscan", "--domtblout", hmm_outfile, hmm_file, fasta_file])
    
        add_domains(hmm_outfile, self.proteins_id, evalue_threshold)
                
                
    def get_protein(self, identifier):
        if identifier in self.proteins_id:
            return self.proteins_id[identifier]
        else:
            print("Protein with ID " + identifier + " not found")
    
    def write_protein_fastas(self, output_dir):
        for protein in self.proteins_id:
            output_file = output_dir + "/" + protein + ".fasta"
            with open(output_file, 'w') as f:
                f.write("> " + protein + "\n" + self.proteins_id[protein].seq)

            
    def write_domain_fastas(self, output_dir):
        for protein in self.proteins_id:
            for domain in self.proteins_id[protein].domains_id:
                output_file = output_dir + "/" + protein + "." + domain + ".fasta"
                with open(output_file, 'w') as f:
                    f.write("> " + protein + "|" + domain + "\n" + self.proteins_id[protein].domains_id[domain].seq)



'''Load a proteome
Input:
    - FASTA of the target proteome
    - HMM of the domains in the proteome
    - Path to the output hmmscan file
    - Optional. Threshold for the evalue. Domains with higher evalues are discarded. Default value, 0.00001
Output
    - An object Proteome containing a set of proteins, indexed by id or position in the proteome.\n Each object protein contain the corresponding domain objects
'''
def load_proteome(fasta_file: str, hmm_file: str, hmm_outfile: str) -> Proteome:
    proteome = Proteome(fasta_file, hmm_file, hmm_outfile)
    
    return proteome



# Test
# fasta_file = "/home/fenix/projects/sars_cov2/data/genome/SARS_CoV2.USCS.fasta.ol"
# hmm_file = "/home/fenix/projects/sars_cov2/data/pfam/hmms/all/all.hmm"
# hmm_outfile = "/home/fenix/projects/sars_cov2/data/pfam/hmms/all/all.hmm.domtblout.sorted.selected"
# 
# proteome = load_proteome(fasta_file, hmm_file, hmm_outfile)
# output_dir_domain = "/home/fenix/sars_cov2/data/pipeline/02_2021/alignments/sequences/domains"
# output_dir_protein = "/home/fenix/sars_cov2/data/pipeline/02_2021/alignments/sequences/proteins"
# proteome.write_domain_fastas(output_dir_domain)
# proteome.write_protein_fastas(output_dir_protein)

# count = 1
# for protein in proteome.proteins_id:
#     for domain in proteome.proteins_id[protein].domains_id:
#         print(count)
#         count += 1



    

