# covmut

Covmuts is a software to build sequence-based model in order to predict the mutability of each position in protein domains. It computes both an indenpent model and a epistitatic model. The latter almost always improves the mutability prediction compared to the former.

Two different interfaces are given, covmut_sequence.py makes a prediction for a single sequence while covmut_proteome.py performs the predicition for an entire proteome.

> ./covmut_sequence.py -h
> usage: covmut_sequence.py [-h] -i INPUT_FASTA -id ID_PROTEIN -db DATABASE [-it ITERATIONS] [-ali INPUT_ALI] [-m MODE] [-nt INPUT_NT_FASTA] [-cpus NUM_CPUS] -o OUTPUT_DIR
>
>optional arguments:
>  -h, --help            show this help message and exit
>  -i INPUT_FASTA, --input_fasta INPUT_FASTA
>                        Input amino acid sequence in FASTA format
>  -id ID_PROTEIN, --id_protein ID_PROTEIN
>                        Identifier of the protein. It will be used to name the output files
>  -db DATABASE, --database DATABASE
>                        Path to the sequence database in FASTA format
>  -it ITERATIONS, --iterations ITERATIONS
>                        Options. Number of jackhmmer iterations to perform
>  -ali INPUT_ALI, --input_ali INPUT_ALI
>                        Optional. Multiple sequence alignment in FASTA to build the model. If it is not given, it is generated
>  -m MODE, --mode MODE  Optional. Select the mode to make the position average for each position from single mutants predictions. Options: all (average of 19 >possible variants), reachable (average from reachable amino acids with one
>                        mutation), weighted (average from reachable amino acids with one mutation weigthed by the number of codons that codified each amino acid). >Reachable and weighted require the nucleotide sequence. Default: all.
>  -nt INPUT_NT_FASTA, --input_nt_fasta INPUT_NT_FASTA
>                        Optional. Input nucleotide sequence in FASTA format
>  -cpus NUM_CPUS, --num_cpus NUM_CPUS
>                        Optional. Number of CPUs to use during jackhmmer searches
>  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
>                        Output directory for the output files
