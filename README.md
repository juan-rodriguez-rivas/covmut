# CovMut

## Introduction
Covmuts is a software to build sequence-based models in order to predict the mutability of each position in protein domains. It computes both an indenpendent model and an epistitatic model. The former assumes independence between positions of the domains while the latter considerer pairwise dependecies between positions, providing, in the vast majority of cases, a better prediction of the mutability.

More information about the application of CovMut for SARS-CoV-2 can be found at https://github.com/GiancarloCroce/DCA_SARS-CoV-2 and https://giancarlocroce.github.io/DCA_SARS-CoV-2/ including the data generated and/or reproduce the results.

Two different interfaces are given, covmut_sequence.py makes a prediction for a single sequence while covmut_proteome.py performs the predicition for an entire proteome.

covmut_sequence.py requires an input FASTA file with the query sequence, an identifier (root for output files), a target sequence database, and an output directory. More information can be found by running `covmut_sequence.py -h`. The usage summary is:
> usage: covmut_sequence.py [-h] -i INPUT_FASTA -id ID_PROTEIN -db DATABASE [-it ITERATIONS] [-ali INPUT_ALI] [-m MODE] [-nt INPUT_NT_FASTA] [-cpus NUM_CPUS] -o OUTPUT_DIR

covmut_proteome.py requires an input FASTA file with the query proteome, a HMM file with the protein domain profiles of interest, two sequence databases (one with a set of diverse proteome, e.g. uniref90, and one for local sequences close to the query, e.g. GISAID for SARS-CoV-2), and an output directory. More information can be found by running `covmut_proteome.py -h`. The usage summary is:
> usage: covmut_proteome.py [-h] -p FASTA_PROTEOME -hmm HMM_FILE -db DATABASE -g DATABASE_GISAID [-it ITERATIONS] [-m MODE] [-cpus NUM_CPUS] [-nt NUCLEOTIDE_GENOME] -o OUTPUT_DIR

In the subdirectoty testing, a small dataset is provided to be able to check that it is running properly.


## Installation
CovMut is wrote in python (tested in 3.9.5 and 3.8.5) and julia (tested in 1.6.1). They can be installed using the package manager:
```
sudo apt install python3 julia
```
In case of problems with the default version of julia, the specific version (1.6.X) can be found at https://julialang.org/downloads/

MAFFT, HHSUITE and HMMER are required to run CovMut:
```
sudo apt install hmmer hhsuite mafft
```

reformat.pl from hhsuite is used, check if it's in the path or add to the path (e.g. in ~/.bashrc, the default path typically is /usr/share/hhsuite/scripts/):
```
export PATH=$PATH:/usr/share/hhsuite/scripts
```

Install the required python packages (e.g. using conda packge manager, https://docs.conda.io/en/latest/miniconda.html):
```
conda install biopython pandas numpy
```

Add covmut to python library path (e.g. in ~/.bashrc):
```
export PYTHONPATH=$PYTHONPATH:/path/covmut
```

Add required packages in julia:
```
julia
]
add FastaIO
add NPZ
add ArgParse
add https://github.com/pagnani/PlmDCA
add https://github.com/carlobaldassi/GaussDCA.jl
```

Any feedback is most welcome, do not hesitate to create an issue or contact me





