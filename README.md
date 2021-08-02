# covmut

Covmuts is a software to build sequence-based models in order to predict the mutability of each position in protein domains. It computes both an indenpendent model and an epistitatic model. The former assumes independence between positions of the domains while the latter considerer pairwise dependecies between positions, providing, in the vast majority of cases, a better prediction of the mutability.

Two different interfaces are given, covmut_sequence.py makes a prediction for a single sequence while covmut_proteome.py performs the predicition for an entire proteome.

covmut_sequence.py requires an input FASTA file with the query sequence, an identifier (root for output files), a target sequence database, and an output directory. More information can be found by running `covmut_sequence.py -h`. The usage summary is:
> usage: covmut_sequence.py [-h] -i INPUT_FASTA -id ID_PROTEIN -db DATABASE [-it ITERATIONS] [-ali INPUT_ALI] [-m MODE] [-nt INPUT_NT_FASTA] [-cpus NUM_CPUS] -o OUTPUT_DIR

covmut_proteome.py requires an input FASTA file with the query proteome, a HMM file with the protein domain profiles of interest, two sequence databases (one with a set of diverse proteome, e.g. uniref90, and one for local sequences close to the query, e.g. GISAID for SARS-CoV-2), and an output directory. More information can be found by running `covmut_proteome.py -h`. The usage summary is:
> usage: covmut_proteome.py [-h] -p FASTA_PROTEOME -hmm HMM_FILE -db DATABASE -g DATABASE_GISAID [-it ITERATIONS] [-m MODE] [-cpus NUM_CPUS] [-nt NUCLEOTIDE_GENOME] -o OUTPUT_DI

In the subdirectoty testing, a small dataset is provided to be able to check that it is running properly.


