import numpy as np
import gzip
from Bio import SeqIO
from pathlib import Path
import os
import subprocess
import tarfile
from io import BytesIO
#for parallel computing
from joblib import Parallel, delayed
import multiprocessing
num_cores_energy = multiprocessing.cpu_count()
from tqdm import tqdm
import pandas as pd
import sys


valid_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-']
aa_3= ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR','-']

d_aa_num= {a:i for i,a in enumerate(valid_aa)}
d_3to1 = {a3:a1 for a3,a1 in zip(aa_3,valid_aa)}



def read_dca_par(path_h_DCA, path_J_DCA):
    ''' read compressed DCA file '''
    tar = tarfile.open(path_h_DCA, "r:gz")
    for member in tar.getmembers():
        f = tar.extractfile(member)
        if f is not None:
            content = f.read()
            load_bytes = BytesIO(content)
            h = np.load(load_bytes)
    tar = tarfile.open(path_J_DCA, "r:gz")
    for member in tar.getmembers():
        f = tar.extractfile(member)
        if f is not None:
            content = f.read()
            load_bytes = BytesIO(content)
            J = np.load(load_bytes)
    return h,J

def compute_sm_energy_dict(seq, h ,J):
    ''' for SINGLE MUTANTS, return a dictionary d['idx', 'mutated_aa'] = energy - energy_wild_type '''
    ''' it can be VERY SLOW and d_sm BIG(all possible sm ~ 21*L) '''
    ''' see below to speed it up '''
    E0 = compute_energy(seq,h,J)
    d_sm = {}
    for i in range(0, len(seq)):
        print(i, len(seq))
        #add also the gap
        for aa in valid_aa:
            new_seq = seq[:i] + aa + seq[(i+1):]
            E = compute_energy(new_seq,h,J)
            print(E)
            d_sm[i,aa] = np.round(E-E0,4)
    return d_sm

def compute_sm_energy(seq, h ,J, idx, aa ):
    ''' for SINGLE MUTANTS, given the ref_seq,h,J and idx(pos_mutations) aa(mutated_aa)
    return energy_sum_single_mutants - energy_wild_type '''
    E0 = compute_energy(seq,h,J)
    E_sum_sm = 0
    for i,a_i in zip(idx, aa):
        new_seq = seq[:i] + a_i + seq[(i+1):]
        E = compute_energy(new_seq,h,J)
        E_sum_sm += E
    return np.round(E_sum_sm,4)

def compute_energy(seq, h, J, parallel = False):
    if all_standard_aa(seq):
        if(parallel == True):
            #DO NOT USE FOR NOW!!!
            #something weird... E_parallel != E_non_parallel
            # parallel actually slower than non parallel (execution time limited by memory access and not processor time??)
            E = 0
            all_ei = Parallel(n_jobs=num_cores_energy)(delayed(compute_energy_given_ai)(seq, h, J, idx_ai) for idx_ai in range(0,len(seq)))
            E = np.sum(all_ei)
            return E
        if(parallel == False):
            E = 0
            for idx_aa1 in range(0, len(seq)):
                aa1 = seq[idx_aa1]
                E -= h[d_aa_num[aa1], idx_aa1]
                for idx_aa2 in range(idx_aa1+1, len(seq)):
                    aa2 = seq[idx_aa2]
                    E -= J[d_aa_num[aa1], d_aa_num[aa2], idx_aa1, idx_aa2]
            return E

def compute_energy_given_ai(seq,h,J, idx_ai):
    '''e.g. idx_ai=1; computing E_1 = h_1 + J_12 + J_13 etc. (good for parallelization)'''
    ai = seq[idx_ai]
    #print("**", idx_ai, ai)
    ei = h[d_aa_num[ai], idx_ai]
    for idx_aj in range(idx_ai+1, len(seq)):
        aj = seq[idx_aj]
        #print(idx_aj, aj)
        ei -= J[d_aa_num[ai], d_aa_num[aj], idx_ai, idx_aj]
    return ei

def compute_entropy_context_ind(path_msa):
    ''' compute context-independent entropy (from msa)'''
    fi = compute_freq(path_msa)
    S = compute_entropy_from_freq(fi)
    return S

def compute_entropy_from_freq(fi, remove_gaps = True, base2 = True):
    if remove_gaps:
        fi = (fi[:20,:])/np.sum(fi[:20,:], axis = 0)
    qq, N = fi.shape
    S = []
    for i in range(0,N):
        si = 0
        for q in range(0,qq):
            si -= fi[q,i]*np.log(fi[q,i])
        if base2:
            si /= np.log(2)
        S.append(si)
    return S

def compute_entropy_context_dep(ref_seq, h,J ):
    ''' compute context-DEPENDENT entropy (from hhblits ref_seq, h, J)'''
    q, N = h.shape
    fi_plm = np.zeros(h.shape)
    #same conventions than in Eq.5.8 (PhD thesis)
    for i in range(0,N):
        #compute denominator
        denom = 0
        for b in range(0,q):
            arg_denom = h[b,i]
            for j in range(0,N):
                if(j!=i):
                    aj = d_aa_num[ref_seq[j]]
                    arg_denom += J[b, aj ,i, j]
            denom += np.exp(arg_denom)
        # compute numerator
        for ai in range(0,q):
            arg_num = h[ai,i]
            for j in range(0,N):
                if(j!=i):
                    aj = d_aa_num[ref_seq[j]]
                    arg_num += J[ai, aj ,i, j]
            num = np.exp(arg_num)
            fi_plm[ai,i] = num/denom
    #return the entropy
    S = compute_entropy_from_freq(fi_plm)
    return S

def compute_num_gap(seq):
    '''return the number of gaps in a sequence '''
    num_gap = 0
    for _,char in enumerate(seq):
        if(char == '-'):
            num_gap += 1
    return num_gap

def compute_gap_fraction(seq):
    num_gap = compute_num_gap(seq)
    frac_gap = (num_gap+0.0)/len(seq)
    return frac_gap

def compute_diff(ref_seq, seq):
    ''' compute the mutations between two strings, return idx_mut, aa_first_seq(wt), aa_second_seq(mutant)'''
    vec_idx = []
    vec_aa1 = []
    vec_aa2 = []
    for idx, aa in enumerate(zip(ref_seq,seq)):
        aa1 = aa[0]
        aa2 = aa[1]
        if (aa1.lower() != aa2.lower()):
            vec_idx.append(idx)
            vec_aa1.append(aa1)
            vec_aa2.append(aa2)
    return vec_idx, vec_aa1, vec_aa2

def compute_dist(ref_seq, seq):
    distance = sum([1 for x, y in zip(ref_seq, seq) if x.lower() != y.lower()])
    return distance

def compute_dist_excluding_gaps(ref_seq, seq):
#     distance = sum([1 for x, y in zip(ref_seq, seq) if ( x.lower() != y.lower() or x == '-' or y == '-' )])
    distance = 0
    for x, y in zip(ref_seq, seq):
        if x == '-':
            continue
        elif y == '-':
            continue      
        elif x.lower() != y.lower():
            distance += 1
            
    return distance

def compute_seqid(ref_seq, seq):
    '''return the sequence identity (seqid) '''
    distance = compute_dist_excluding_gaps(ref_seq,seq)
    distance /= len(seq)
    seqid = 1 - distance
    return seqid

def compute_freq(path_msa):
    ''' compute single point frequencies of an MSA '''
    records_msa = list(SeqIO.parse(open(path_msa,'r'), "fasta"))
    fi = np.zeros(( len(d_aa_num), len(records_msa[0].seq) ))
    for idx_rec, rec in enumerate(records_msa):
        seq = rec.seq
        for idx_aa, amino_a in enumerate(seq):
            fi[d_aa_num[amino_a], idx_aa] += 1
    #add (small) pseudocount to take into account 0 frequencies (0*log(0))
    alpha = 0.0001
    fi = (1-alpha)*fi + alpha/2
    #normalize
    fi /= fi.sum(axis = 0)
    return fi

def all_standard_aa(seq):
    '''return True if sequence contains only standard-aa'''
    for char in seq:
        if((char not in valid_aa) and char !='-'):
            #print("seq containing non standard aa: "+char)
            return False
            break
    return True

def split_proteome(path_ref_proteome, name_ref_proteome, tmp_path):
    ''' simple function to split the reference proteome in reference proteins'''
    with open(os.path.join(path_ref_proteome, name_ref_proteome), "r") as input_handle:
        for record_ref in SeqIO.parse(input_handle, "fasta"):
            name = record_ref.id
            seq_ref = str(record_ref.seq)
            #save tmp file with the seq of the reference
            name_tmp_file = "ref_"+name
            f_tmp = open(os.path.join(tmp_path,name_tmp_file),"w")
            f_tmp.write(">"+name+"\n")
            f_tmp.write(seq_ref)
            f_tmp.close()
        return 0

def run_hhblits(path_hhblits, path_ref_prot, path_db, path_msa_out, num_cores):
        ''' run hhblits, get the distant homologs MSA, return the number of sequences '''
        #1) run hhblits
        FNULL = open(os.devnull, 'w')
        subprocess.run([path_hhblits, '-i', path_ref_prot, '-d', path_db, '-oa3m', path_ref_prot+".a3m", '-cpu' , str(num_cores)], stdout=FNULL, stderr=subprocess.STDOUT)
        #num of sequences
        file_out = open(path_msa_out, 'w')
        #2) parse and filter the hhblits msa
        with open(path_ref_prot+".a3m", "r") as input_handle:
            for idx_record, record in enumerate(SeqIO.parse(input_handle, "fasta")):
                seq = str(record.seq)
                #hhblits ouput is a3m format, to make it a fasta remove dot and lower
                seq = ''.join(char for char in seq if (char.isupper() or char =='-'))
                # 2.1) do the  filtering
                records_ref = list(SeqIO.parse(open(path_ref_prot,'r'), "fasta"))
                ref_seq = str(records_ref[0].seq)
                # - remove sequences which are to gapped (i.e. gap_fraction mst be less than 10% gap)
                # - remove sequence which are CLOSE to the reference sequence (i.e. sequence_identity must be LESS than 90%)
                # - remove sequences containing non standard aa
                if( (compute_gap_fraction(seq) < 0.1) and (compute_seqid(ref_seq, seq) < 0.9) and all_standard_aa(seq)):
                    file_out.write(str(">"+record.id)+'\n')
                    file_out.write(str(seq)+'\n')
            file_out.close()
        return 0

def filterMSA(path_ref_prot, path_msa_in, path_msa_out, include_refseq=True, max_grap_fraction = 0.2, max_seq_id = 0.9):
        file_out = open(path_msa_out, 'w')
        # parse and filter the msa
        records_ref = list(SeqIO.parse(open(path_ref_prot,'r'), "fasta"))
        ref_seq = str(records_ref[0].seq)
        with open(path_msa_in, "r") as input_handle:
            count = 1
            for idx_record, record in enumerate(SeqIO.parse(input_handle, "fasta")):
                seq = str(record.seq)
                #remove dot and lower
                seq = ''.join(char for char in seq if (char.isupper() or char =='-'))
                #  do the filtering
                # - remove the sequences which are to gapped (i.e. sequence must contain less than 10 gap)
                # - remove sequence which are close the reference sequence (i.e. sequence_identity must be less than 90%)
                # - remove sequences containing non standard aa
                if include_refseq and count == 1: # Keep the first seq, i.e., the reference sequence
                    file_out.write(str(">"+record.id)+'\n')
                    file_out.write(str(seq)+'\n')
                    count += 1 
                elif( (compute_gap_fraction(seq) < max_grap_fraction) and (compute_seqid(ref_seq, seq) < max_seq_id) and all_standard_aa(seq)):
                    file_out.write(str(">"+record.id)+'\n')
                    file_out.write(str(seq)+'\n')
        file_out.close()
        return 0
    
def filterMSA_gisaid(path_ref_prot, path_msa_in, path_msa_out, max_grap_fraction = 0.2, min_seq_id = 0.9):
        file_out = open(path_msa_out, 'w')
        # parse and filter the msa
        records_ref = list(SeqIO.parse(open(path_ref_prot,'r'), "fasta"))
        ref_seq = str(records_ref[0].seq)
        with open(path_msa_in, "r") as input_handle:
            count = 1
            for idx_record, record in enumerate(SeqIO.parse(input_handle, "fasta")):
                seq = str(record.seq)
                #remove dot and lower
                seq = ''.join(char for char in seq if (char.isupper() or char =='-'))
                #  do the filtering
                # - remove the sequences which are to gapped (i.e. sequence must contain less than 10 gap)
                # - remove sequence which are far the reference sequence (i.e. sequence_identity must greater than 90%)
                # - remove sequences containing non standard aa
                if count == 1: # Keep the first seq, i.e., the reference sequence
                    file_out.write(str(">"+record.id)+'\n')
                    file_out.write(str(seq)+'\n')
                    count += 1 
                elif( (compute_gap_fraction(seq) < max_grap_fraction) and (compute_seqid(ref_seq, seq) > min_seq_id) and all_standard_aa(seq)):
                    file_out.write(str(">"+record.id)+'\n')
                    file_out.write(str(seq)+'\n')
        file_out.close()
        return 0


def do_DCA_inference(path_msa, path_dca_par, min_num_seq, num_cores):
        #1) number of lines (sequences). N.b. in Uniclust30 Meff ~ M
        M = len(open(path_msa).readlines())/2
        #2) do the inference with DCA
        #only for msa with more than min_num_seq sequences
        if( M > min_num_seq):
            #import julia (Try to include julia variables into python => TO DO, see 'import julia') ---> doesn't work, I gave up...
#             filename = 
            out_file = path_msa + '.out'
#             f_julia= open(os.path.join(path_dca_par,out_file), 'a')
            f_julia= open(out_file, 'w')
            f_julia.write(path_msa.split("/")[-1]+".fa"+'\n')
            f_julia.close() #close and re-open. Otherwise it writes the msa only at the end (after plmDCA info) (function waits subprocess to finish)
#             f_julia= open(os.path.join(path_dca_par,out_file), 'a')
            f_julia= open(out_file, 'w')
            subprocess.run(["julia",'-p', str(num_cores), './src/plmDCA_inference.jl',path_msa, path_dca_par], stdout=f_julia, stderr=subprocess.STDOUT)
            f_julia.close()
        else:
            print('... ERROR! too few seqs (M={0})!'.format(str(M)))
        return 0

def run_phmmer(path_phmmer, path_ref_prot, path_db, path_msa_out, path_tmp_stockholm, path_tmp_msa, num_cores):
        ''' run phmmer, get the local homologs MSA (form E coli strains) '''
        file_out = open(path_msa_out, 'w')
        #1) run phmmer
        FNULL = open(os.devnull, 'w')
        subprocess.run([path_phmmer, '-A', path_tmp_stockholm, '--cpu', str(num_cores), path_ref_prot, path_db], stdout=FNULL, stderr=subprocess.STDOUT)
        #2) convert stockholm to fasta
        subprocess.run(['./src/stockholm2fasta.pl', '-g',  path_tmp_stockholm ], stdout=open(path_tmp_msa,'w'),stderr=subprocess.STDOUT)
        #3) parse and filter the hhblits msa
        with open(path_tmp_msa, "r") as input_handle:
            for idx_record, record in enumerate(SeqIO.parse(input_handle, "fasta")):
                seq = str(record.seq)
                #remove dot and lower
                seq = ''.join(char for char in seq if (char.isupper() or char =='-'))
                # 2.1) do the filtering
                records_ref = list(SeqIO.parse(open(path_ref_prot,'r'), "fasta"))
                ref_seq = str(records_ref[0].seq)
                # - remove the sequences which are to gapped (i.e. sequence must contain less than 10 gap)
                # - remove sequence which are FAR the reference sequence (i.e. sequence_identity must be MORE than 90%)
                # - remove sequences containing non standard aa
                if( (compute_num_gap(seq) < 10) and (compute_seqid(ref_seq, seq) > 0.9) and all_standard_aa(seq)):
                    file_out.write(str(">"+record.id)+'\n')
                    file_out.write(str(seq)+'\n')
        file_out.close()
        #rm prefiltering msa
        subprocess.run(['rm' , path_tmp_msa])
        #rm stockholm (too big!)
        subprocess.run(['rm' , path_tmp_stockholm])
        return 0

def compute_energy_local_msa(ref_prot_file, output_file, ali_file,h,J, verbose ):
    ''' load DCA model, compute energy of strains (also e_sum_sm and postions mutated '''
    records_ref = list(SeqIO.parse(open(ref_prot_file,'r'), "fasta"))
    ref_seq = str(records_ref[0].seq)
    ref_name = str(records_ref[0].id)
    E0 = 0
    #0. compute energy single mutants -> NO need for it ( and very slow, seq function compute_sm_energy(seq, h, J, idx, aa)
    #d_sm = compute_sm_energy(ref_seq, h, J)
    path_msa_local = os.path.join(ali_file)
    all_seq_name = []
    all_seq_num_occurences= []
    all_seq_e = []
    all_seq_dist = []
    all_seq_ref_prot = []
    all_seq_sum_sm = []
    all_seq_mut_idx= []
    all_seq_mut_aa = []
    E0 = compute_energy(ref_seq, h, J)
    with open(path_msa_local,"r") as f:
        for record in tqdm(SeqIO.parse(f,"fasta")):
            seq = str(record.seq)
            E = compute_energy(seq, h, J)
            idx, aa1, aa2 = compute_diff(ref_seq,seq)
            # sum of energies single mutants
            E_sum_sm = compute_sm_energy(ref_seq, h, J, idx, aa2)
            #num mutations
            dist = len(idx)
            name_seq_list = (str(record.id).split('/')[0]).split('-')
            if(verbose == True):
                for name_seq in name_seq_list:
                        all_seq_ref_prot.append(ref_name)
                        all_seq_name.append(name_seq)
                        all_seq_e.append(np.round(E,4))
                        all_seq_dist.append(int(dist))
                        all_seq_sum_sm.append(np.round(E_sum_sm,4))
                        all_seq_mut_idx.append(idx)
                        all_seq_mut_aa.append(aa2)
                all_seq_e_e0 = np.round(all_seq_e - E0,4)
            if(verbose == False):
                all_seq_ref_prot.append(ref_name)
                all_seq_num_occurences.append(len(name_seq_list))
                all_seq_e.append(np.round(E,4))
                all_seq_dist.append(int(dist))
                all_seq_sum_sm.append(np.round(E_sum_sm,4))
                all_seq_e_e0 = np.round(all_seq_e - E0 ,4)
                all_seq_mut_idx.append(idx)
                all_seq_mut_aa.append(aa2)
    if(verbose == True):
        df = pd.DataFrame({'ref':all_seq_ref_prot, 'seq_name':all_seq_name, 'e':all_seq_e, 'e-e0':all_seq_e_e0, 'e_sum_sm': all_seq_sum_sm, 'dist':all_seq_dist, 'idx_mut':all_seq_mut_idx, 'aa_mut': all_seq_mut_aa})
    if(verbose == False):
        df = pd.DataFrame({'ref':all_seq_ref_prot, 'num_occurences':all_seq_num_occurences, 'e':all_seq_e, 'e-e0':all_seq_e_e0, 'e_sum_sm': all_seq_sum_sm,'dist':all_seq_dist, 'idx_mut':all_seq_mut_idx, 'aa_mut': all_seq_mut_aa})
    df.to_csv(os.path.join(output_file), index = False)
    return 0

def compute_energy_ind_msa(ref_prot_file, ali_file, output_file, h,J ):
    ''' load DCA model, compute energy of mutation sampled from the profile model'''
    records_ref = list(SeqIO.parse(open(ref_prot_file,'r'), "fasta"))
    ref_seq = str(records_ref[0].seq)
    ref_name = str(records_ref[0].id)
    path_msa_local_ind = os.path.join(ali_file)
    all_seq_ref_prot = []
    all_seq_e = []
    all_seq_e_e0 = []
    all_seq_dist = []
    E0 = compute_energy(ref_seq, h, J)
    with open(path_msa_local_ind,"r") as f:
        for record in tqdm(SeqIO.parse(f,"fasta")):
            seq = str(record.seq)
            dist = compute_dist(ref_seq,seq)
            E = compute_energy(seq, h, J)
            all_seq_ref_prot.append(ref_name)
            all_seq_e.append(np.round(E,4))
            all_seq_dist.append(int(dist))
        all_seq_e_e0 = np.round(all_seq_e - E0 ,4)
        df_ind = pd.DataFrame({'ref':all_seq_ref_prot, 'e':all_seq_e, 'e-e0':all_seq_e_e0, 'dist':all_seq_dist})
    df_ind.to_csv(os.path.join(output_file), index = False)
    return 0

def compute_energy_rand_msa(ref_prot_file, ali_file, output_file, h,J):
    ''' load DCA model, compute energy of mutation sampled from the random model'''
    records_ref = list(SeqIO.parse(open(ref_prot_file,'r'), "fasta"))
    ref_seq = str(records_ref[0].seq)
    ref_name = str(records_ref[0].id)
    path_msa_local_rand= os.path.join(ali_file)
    all_seq_ref_prot = []
    all_seq_e = []
    all_seq_e_e0 = []
    all_seq_dist = []
    E0 = compute_energy(ref_seq, h, J)
    with open(path_msa_local_rand,"r") as f:
        for record in tqdm(SeqIO.parse(f,"fasta")):
            seq = str(record.seq)
            dist = compute_dist(ref_seq,seq)
            E = compute_energy(seq, h, J)
            all_seq_ref_prot.append(ref_name)
            all_seq_e.append(np.round(E,4))
            all_seq_dist.append(int(dist))
        all_seq_e_e0 = np.round(all_seq_e - E0 ,4)
        df_rand = pd.DataFrame({'ref':all_seq_ref_prot, 'e':all_seq_e, 'e-e0':all_seq_e_e0, 'dist':all_seq_dist})
    df_rand.to_csv(os.path.join(output_file, 'e_'+ref_name+'_rand.csv'), index = False)
    return 0

def compute_all_entropies(ref_prot_file, ali_file, ali_file_local, output_file, h, J ):
    ''' compute s_ind, s_dep , s_strains '''
    records_ref = list(SeqIO.parse(open(ref_prot_file,'r'), "fasta"))
    ref_seq = str(records_ref[0].seq)
    ref_name = str(records_ref[0].id)
    ####################################################################################################
    #context IND entropy (from msa_hhblits)
    path_msa_hhblits = os.path.join(ali_file)
    S_ind = np.round(compute_entropy_context_ind(path_msa_hhblits),4)
    #context DEP entropy (from ref_seq, h, J)
    S_dep = np.round(compute_entropy_context_dep(ref_seq, h,J),4)
    #compute entropy in MSA_local (hhblits) (i.e. observed polymorphism?)
    path_msa_local = os.path.join(ali_file_local)
    S_local_obs = np.round(compute_entropy_context_ind(path_msa_local),4)
    all_seq_ref_prot = [ref_name for i in range(0,len(ref_seq))]
    all_seq_idx= [i for i in range(0,len(ref_seq))]
    df_s = pd.DataFrame({'ref':all_seq_ref_prot, 'idx': all_seq_idx, 's_ind':S_ind, 's_dep':S_dep, 's_local_obs':S_local_obs})
    df_s.to_csv(os.path.join(output_file), index = False)
    return 0

