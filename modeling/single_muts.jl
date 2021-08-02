#!/usr/bin/env julia

using GaussDCA
using FastaIO
using DelimitedFiles
using ArgParse

# Some needed scripts
# PATH_SCRIPT = "/home/jrodrig5/landscapes_coevolution/src/giancarlo"
# include(joinpath(PATH_SCRIPT,"read_script.jl"))
# include(joinpath(PATH_SCRIPT,"compute_fi_fij.jl"))

# PATH_SCRIPT = "/home/fenix/projects/sars_cov2/src/energies"
include(joinpath(@__DIR__,"tools_energy_analysis.jl"))
include(joinpath(@__DIR__,"compute_fi_fij.jl"))
# @__DIR__

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--refseq", "-r"
            help = "Path to the reference sequence in FASTA format"
            required = true
        "--ali", "-i"
            help = "Path to the alignment in FASTA format"
            required = true
        "--root_parameters", "-p"
            help = "Root path to the gzip compressed files with the h and J parameters. They should be at root.h.gz and root.J.gz"
            required = true
        "--outfile", "-o"
            help = "Path for the output file"
            required = true
    end

    return parse_args(s)
end

let alphabet = [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
    # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt8(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end


aas = Dict(1 => "A",
    2 => "C",
    3 => "D",
    4 => "E",
    5 => "F",
    6 => "G",
    7 => "H",
    8 => "I",
    9 => "K",
    10 => "L",
    11 => "M",
    12 => "N",
    13 => "P",
    14 => "Q",
    15 => "R",
    16 => "S",
    17 => "T",
    18 => "V",
    19 => "W",
    20 => "Y",
    21 => "-")
    
function read_par(path_par::AbstractString)

   # read fields
   h_name = path_par*".h"
   h = readdlm(h_name)

   q, N = size(h)

   J_name = path_par*".J"
   J_m = readdlm(J_name)
   J_m += J_m'

   J = zeros(q,q,N,N);

   for i=1:N
      for j = (i+1):N
           J[:,:,i,j] = J_m[(1+q*(i-1)):(q+q*(i-1)),(1+q*(j-1)):(q+q*(j-1))]
      J[:,:,j,i] = (J_m[(1+q*(i-1)):(q+q*(i-1)),(1+q*(j-1)):(q+q*(j-1))])'
       end
   end

   return h,J

end


function fasta2matrix(filename::AbstractString, max_gap_fraction::Real)

    f = FastaReader(filename)

    max_gap_fraction = Float64(max_gap_fraction)

    # pass 1

    seqs = Int[]
    inds = Int[]
    fseqlen = 0

    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end

    length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2
    # Z = Array{Int8}(fseqlen, length(seqs))
    Z = zeros(Int8, (fseqlen, length(seqs)))

    seqid = 1
    for (name, seq) in f
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        seqid += 1
    end
    @assert seqid == length(seqs) + 1

    close(f)

    return Z'
end


function main()
    parsed_args = parse_commandline()

    reference_sequence_file = parsed_args["refseq"]
    input_ali = parsed_args["ali"]
    rootname_parameters = parsed_args["root_parameters"]
    output_file = parsed_args["outfile"]
    
    # Load info for IND model
    Z = GaussDCA.read_fasta_alignment(input_ali,0.9)
    theta = 0.2
    N, M = size(Z)
    q = Int(maximum(Z))
    fi_tens, _, _, _, Meff, _ = compute_new_fi_fij(Z,q,theta)
    fi_tens_pseudo = add_pseudocounts(fi_tens, 0.01)
    
    # Load info for DCA model
    h,J = read_par(rootname_parameters)
    S = fasta2matrix(reference_sequence_file,0.9)
    
    # To vector
    ref_seq = S[1,:]
    
    # Compute energies
    energy_ref_dca = compute_energy_single_sequence(h,J,ref_seq)
    energy_ref_ind = compute_energy_independent_model(fi_tens_pseudo, ref_seq)
        
    file_out = open(output_file,"w")
    tam = length(ref_seq)
    # effects = Array{Float64,2}(tam,19)
    effects = zeros(Float64,tam,19)
    println(file_out, "position\taa_wildtype\taa_mutation\tenergy_IND\tscore_IND\tenergy_DCA\tscore_DCA")
    for pos in 1:tam
        original = ref_seq[pos]
        for mut in 1:21
            mut_seq = copy(ref_seq)
            mut_seq[pos] = mut
            energy_mut_dca = compute_energy_single_sequence(h,J,mut_seq)
            energy_mut_ind = compute_energy_independent_model(fi_tens_pseudo, mut_seq)

            diff_energy_dca = -(energy_mut_dca - energy_ref_dca)
            diff_energy_ind = -(energy_mut_ind - energy_ref_ind)
            println(file_out, pos,"\t", aas[original],"\t",aas[mut],"\t",round(energy_mut_ind, digits=4),"\t",round(diff_energy_ind, digits=4),"\t",round(energy_mut_dca, digits=4),"\t",
            round(diff_energy_dca, digits=4))
        end
#         println(file_out, pos,"\t", aas[original],"\t","*","\t","NA","\t","NA")
    end
end

main()


