#!/usr/bin/env julia

#import PlmDCA
using PlmDCA
using LinearAlgebra
using DelimitedFiles
using ArgParse
using NPZ

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "Input alignment"
            required = true
        "--output", "-o"
            help = "Output root"
            required = true
        "--theta", "-t"
            help = "Reweighting value, default=0.2 (reweighting at 80% sequence identity)"
            default = 0.2
            arg_type = Float64
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    input = parsed_args["input"]
    output = parsed_args["output"]
    theta = parsed_args["theta"]

    a = plmdca(input, verbose = false, theta = theta)

    h = a.htensor
    J = a.Jtensor
    scores = a.score

    output_scores = output*".scores"
    output_stream = open(output_scores, "w")
    for pred in scores
        println(output_stream, string(pred[1]) * "\t" * string(pred[2]) * "\t" * string(pred[3]))
    end
    close(output_stream)


    # J is too big! Reduce is size
    # switch to matrix form of J
    q,N = size(h)
    J_m = zeros(q*N,q*N)
    for i=1:N
    	for j = 1:N
    		J_m[(1+q*(i-1)):(q+q*(i-1)),(1+q*(j-1)):(q+q*(j-1))] = J[:,:,i,j]
    	end
    end

    # Remove lower_triangular_part (use the fact that J_m = J_m')
    J_m -= LowerTriangular(J_m)

    # Reduce precision to 8 digits
    J_m = round.(J_m,digits=8)

    # Write a out files
    J_name = output*".J"
    writedlm(J_name, J_m)

    h_name = output*".h"
    writedlm(h_name,h)
    
    # Save also in NPZ format    
#     name_h = h_name*".npz"
#     name_J = J_name*".npz"
#     
#     npzwrite(name_h, h)
#     npzwrite(name_J, J)
#     
#     name_h_compressed = name_h*".tar.gz"
#     name_J_compressed = name_J*".tar.gz"
#     run(`tar -czvf $name_h_compressed $name_h`)
#     run(`tar -czvf $name_J_compressed $name_J`)
#     run(`rm $name_h`)
#     run(`rm $name_J`)
end

main()
