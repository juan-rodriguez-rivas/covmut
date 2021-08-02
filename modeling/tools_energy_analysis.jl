#!/usr/bin/env julia

#compute energy of a single sequence
#h=(q,N) and J=(q,N)
function compute_energy_single_sequence(h::Array{Float64,2},
                                        J::Array{Float64,4},
                                        S::Vector)

    N = size(h)[2]
    q = size(h)[1]
    E = 0.0
    for i = 1:N
        E -= h[S[i],i]
        for j = (i+1):N
			E -= J[S[i],S[j],i,j]
		end
	end
return E
end

function make_random_mutations_compute_energy(h, J, num_muts, ref_seq, length_ref_seq, freqs_pseudo)

		unsorted_columns = randperm(length_ref_seq)
		j = 1
		muts_done = 0
		mutated_sequence = copy(ref_seq)
		while muts_done <= num_muts
			column = unsorted_columns[j]
			if ref_seq[column] != 21 # If it is gap, we not consider the position
				unsorted_residues = randperm(20)
				if ref_seq[column] == unsorted_residues[1] # If the first random residue is the same as in the wild type, take the second. Take the first otherwise
					mut = unsorted_residues[2]
				else
					mut = unsorted_residues[1]
				end
				mutated_sequence[column] = mut
				muts_done = muts_done + 1
				# println("muts done "*string(muts_done))
			end
			j = j + 1
			# println(j)
		end
	return compute_energy_single_sequence(h,J,mutated_sequence)
end


function make_random_mutations_single_sequence(num_muts, ref_seq)

	unsorted_columns = randperm(length(ref_seq))
	j = 1
	muts_done = 0
	mutated_sequence = copy(ref_seq)
	while muts_done <= num_muts
		column = unsorted_columns[j]
		if ref_seq[column] != 21 # If it is gap, we not consider the position
			unsorted_residues = randperm(20)
			if ref_seq[column] == unsorted_residues[1] # If the first random residue is the same as in the wild type, take the second. Take the first otherwise
				mut = unsorted_residues[2]
			else
				mut = unsorted_residues[1]
			end
			mutated_sequence[column] = mut
			muts_done = muts_done + 1
		end
			j = j + 1
	end
	return mutated_sequence
end



function make_random_mutations_alignment(mutations, ref_seq)

	num_seqs = length(mutations[:,1])
	tam = length(ref_seq)
	alignment = Array{Int8}(tam, num_seqs)

	for seq in 1:num_seqs
		num_muts = length(mutations[seq,1])

		mutated_sequence = make_random_mutations_single_sequence(num_muts, ref_seq)
		for j in 1:tam
			alignment[j, seq] = mutated_sequence[j]
		end
	end
	return alignment
end



function compute_energy_independent_model(freqs, seq)
	energy_ind = 0
	for col in 1:length(freqs[1,:])
		aa = seq[col]
		energy_ind += log(freqs[aa,col])
	end
	return -energy_ind
end


function compute_diff_independent_only_muts(freqs, ref_seq, seq)
	diff_ind = 0
	for col in 1:length(seq)
		aa = seq[col]
		aa_ref_seq = ref_seq[col]
		if aa != aa_ref_seq
			diff_ind += log(freqs[aa,col]) - log(freqs[aa_ref_seq,col])
		end
	end
	return -diff_ind
end


function compute_diff_independent_only_muts_sequences(freqs, ref_seq, sequences)
	num_seqs = length(sequences[1,:])
	diff_ind = Vector(num_seqs)
	for i in 1:num_seqs
		seq = sequences[:,i]
		diff_ind[i] = compute_diff_independent_only_muts(freqs, ref_seq, seq)
	end
	return diff_ind
end


function compute_energy_independent_model_sequences(freqs, sequences)
	num_seqs = length(sequences[1,:])
	energy_ind = Vector(num_seqs)
	for i in 1:num_seqs
		seq = sequences[:,i]
		energy_ind[i] = compute_energy_independent_model(freqs, seq)
	end
	return energy_ind
end


function add_pseudocounts(freqs, alpha)
	tam_i = length(freqs[1,:])
	tam_j = length(freqs[:,1])
# 	freqs_pseudo = Array{Float64,2}(tam_j,tam_i)
    freqs_pseudo = zeros(Float64, tam_j,tam_i)
	for i in 1:tam_i
		for j in 1:tam_j
			freqs_pseudo[j,i] = (1-alpha)*freqs[j,i] + alpha/21
		end
	end

	return freqs_pseudo
end


# function make_mutation_according_to_IND(ref_seq, pos, freqs)
# 	aa = ref_seq[pos]
# 	mut_aa = aa
#
# 	if aa != 21 # If it is a gap, we do not change anything, just return the gap
# 		random_num = rand()
# 		found = 0
# 		while(!found)
#
# 			found = 1
# 		end
# 	end
#
# 	return mut_aa
# end


function normalized_cummulative_without_gaps(freqs)
	num_cols = length(freqs[1,:])
	normalized_freqs = Array{Float64,2}(20, num_cols)
	cummulative_normalized_freqs = Array{Float64,2}(20, num_cols)
	for col in 1:num_cols
		sum_freq = 0
		for aa in 1:20 # to avoid gap = 21
			sum_freq = sum_freq + freqs[aa,col]
		end
		# Normalize without gap
		for aa in 1:20 # to avoid gap = 21
			normalized_freqs[aa,col] = freqs[aa,col]/sum_freq
		end
		cummulative_normalized_freqs[1,col] = normalized_freqs[1,col]
		for aa in 2:20 # to avoid gap = 21
			cummulative_normalized_freqs[aa,col] = cummulative_normalized_freqs[aa-1,col]+normalized_freqs[aa,col]
		end
	end

	return cummulative_normalized_freqs
end


function get_mutations(ali, ref_seq)
	num_seqs = length(ali[1,:])
	num_pos = length(ali[:,1])

	mutations = []
	for seq in 1:num_seqs
		mutations_this_seq = []
		for pos in 1:num_pos
			if ali[pos, seq] != ref_seq[pos]
				push!(mutations_this_seq, [pos, ali[pos, seq]])
			end
		end
		push!(mutations, mutations_this_seq)
	end

	return mutations
end


function compute_delta_energy_mut_separately(mutations, ref_seq)
	num_seqs = length(mutations[:,1])
	diff_energy = zeros(num_seqs)
	for seq in 1:num_seqs
		num_muts = length(mutations[seq,1])
		sum_delta = 0
		for mut in 1:num_muts
			mutation = mutations[seq,1][mut]
			mut_seq = copy(ref_seq)
			mut_pos = mutation[1]
			mut_aa = mutation[2]
			mut_seq[mut_pos] = mut_aa
			energy_mut = compute_energy_single_sequence(h,J,mut_seq)
			delta = energy_mut - energy_ref
			sum_delta = sum_delta + delta
		end
		diff_energy[seq] = sum_delta
	end
	return diff_energy
end


function compute_delta_energy_mut_return_separately(mutations, ref_seq)
	num_seqs = length(mutations[:,1])
	diff_energy = zeros(num_seqs)

	effects = []
	for seq in 1:num_seqs
		num_muts = length(mutations[seq,1])
		sum_delta = 0
		mutations_this_seq = []
		for mut in 1:num_muts
			mutation = mutations[seq,1][mut]
			mut_seq = copy(ref_seq)
			mut_pos = mutation[1]
			mut_aa = mutation[2]
			mut_seq[mut_pos] = mut_aa
			energy_mut = compute_energy_single_sequence(h,J,mut_seq)
			delta = energy_mut - energy_ref
			sum_delta = sum_delta + delta
			push!(mutations_this_seq, delta)
		end
		diff_energy[seq] = sum_delta
		push!(effects, mutations_this_seq)
	end
	return effects
end
