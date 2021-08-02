
# it take Z as Gauss DCA 
# and returns 
# - fi_tens (q,N), fij_tens(q,q,N,N)
# - fi (q*N),fij (q*N*q*N)
# - Meff (effective numb of seq)
# - W (the vector of M weights)

# [normal gaussDCA returns N*(q-1), frequencies]

q = 21

function compute_new_fi_fij(Z::Matrix{Int8}, q, theta)

	#G.C
	#q+1 to get sequences which are 21X21

    W, Meff = GaussDCA.compute_weights(Z, q+1, theta) 
    Pi_matr, Pij_matr = compute_all_freqs(Z, W, Meff,q+1)
    Pi_tens, Pij_tens = freq_tens_form(Pi_matr, Pij_matr, Z)

    return Pi_tens, Pij_tens, Pi_matr, Pij_matr, Meff, W
end


function compute_all_freqs(Z::Matrix{Int8}, W::Vector{Float64}, Meff::Float64,q::Int)
    N, M = size(Z)
    s = q - 1

    Ns = N * s

    Pij = zeros(Ns, Ns)
    Pi = zeros(Ns)

    ZZ = Vector{Int8}[vec(Z[i,:]) for i = 1:N]

    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        for k = 1:M
            a = Zi[k]
            a == q && continue
            Pi[i0 + a] += W[k]
        end
        i0 += s
    end
    Pi /= Meff

    i0 = 0
    for i = 1:N
        Zi = ZZ[i]
        j0 = i0
        for j = i:N
            Zj = ZZ[j]
            for k = 1:M
                a = Zi[k]
                b = Zj[k]
                (a == q || b == q) && continue
                Pij[i0+a, j0+b] += W[k]
            end
            j0 += s
        end
        i0 += s
    end
    for i = 1:Ns
        Pij[i,i] /= Meff
        for j = i+1:Ns
            Pij[i,j] /= Meff
            Pij[j,i] = Pij[i,j]
        end
    end

    return Pi, Pij
end


function freq_tens_form(fi, fij, Z::Matrix{Int8})

    N, M = size(Z)

	#adapt the frequencies
	fi_tens = zeros(q,N)

	for i = 1:N
		fi_tens[:,i] = fi[(1+(q)*(i-1)):((q)+(q)*(i-1))]
	end


	fij_tens= zeros(q,q,N,N);
	for i = 1:N
		for j = 1:N
			fij_tens[:,:,i,j] = fij[(1+(q)*(i-1)):((q)+(q)*(i-1)),(1+(q)*(j-1)):((q)+(q)*(j-1))]
		end
	end

	return fi_tens, fij_tens
end

