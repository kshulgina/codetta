# Codetta v2.0


## Description

Codetta is a Python program for predicting the genetic code (codon table) of an 
organism from nucleotide sequence data.   

The analysis consists of three steps:

1. Align the input nucleotide sequence to a database of profile Hidden Markov 
models (HMMs) of proteins (such as the Pfam database)
2. Summarize the resulting alignments
3. Infer the genetic code from the alignment summary file

If you are looking to reproduce results from 
[Shulgina & Eddy (2021)](https://elifesciences.org/articles/71402), please follow 
the README for [Codetta v1.0](https://github.com/kshulgina/codetta/releases/tag/v1.0).

If you encounter any problems in the installation or usage of Codetta, please 
leave a Github issue or email me at shulgina@g.harvard.edu


## Download and setup

Codetta was developed for Python versions 3.5+ on Linux and MacOS. Required 
Python packages are `numpy` and `scipy`.

Once inside the Codetta directory, run the setup script
	
	./setup.sh

This will check Python package requirements and set up a local version of HMMER.

Codetta additionally requires:

- `wget` and `gzip`: on Mac, use install commands `brew install wget` and `brew 
install gzip`. For Linux, you'll have to use your system's package management 
tool.


## Setting up a profile HMM database

Codetta aligns the input sequence to a database of profile HMM models of proteins. 

By default, Codetta will use the Pfam database. You can download a version of 
Pfam 35.0 specially built for Codetta from our website. From the `codetta` directory, 
use this command to download and uncompress it (note: this will take about 3 Gb of 
disk space)

	cd resources/
	wget http://eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz
	tar xf Pfam-A_enone.tar.gz
	rm Pfam-A_enone.tar.gz
	cd ..

If you want to build your own custom profile HMM database, this will be described in a
later section. 


## Usage

The program `codetta.py` rolls the three main analysis steps into a single script. 
Usage for `codetta.py` is

	./codetta.py [input sequence file] [optional arguments]

The three steps can also be run separately using the programs:

- `codetta_align.py`: Align profile HMMs to the input nucleotide sequence.
- `codetta_summary.py`: Summarize profile HMM alignments into an alignment summary 
file.
- `codetta_infer.py`: Infer the genetic code from the alignment summary file.

General usage for these programs is

	./[program name] [input file or prefix] [optional arguments]

For any of these programs, type `./[program name] --help` for complete usage details.


## Tutorials
### Starting out...

The simplest way to run Codetta is by using the `codetta.py` program. 

This program performs the three analysis steps in order. All you have 
to do is specify an input nucleotide sequence file.

In the `examples/` directory, you will find the file `GCA_000442605.1.fna`
 which contains the genome of the bacterium _Nasuia deltocephalinicola_.

We can predict the genetic code of this bacterium simply by running

	./codetta.py examples/GCA_000442605.1.fna

You will see outputs written to the terminal indicating that each of the three 
Codetta steps is executing. 

At the end, the inferred genetic code is printed 

	Genetic code: FFLLSSSSYY??CCWWL?L?PPPPHHQQ????I?IMTTT?NNKKSSRRV?VVAAAADDEEGGGG

This corresponds to the inferred translation of each of the 64 codons, in order 
from 'UUU, UUC, UUA, UUG, UCU, UCC, ..., GGA, GGG' (iterating 3rd, 2nd, then 
1st base through UCAG). 

A detailed summary of the analysis can be found at 
`examples/GCA_000442605.1.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_excl-mtvuy.genetic_code.out`. 
The long file extension specifies the inference parameters.

This file contains a detailed summary of the genetic code inference results:

	# Analysis arguments
	alignment_prefix   examples/GCA_000442605.1.fna
	profile_database   Pfam-A_enone.hmm
	results_summary    None
	evalue_threshold   1e-10
	prob_threshold     0.9999
	max_fraction       0.01
	excluded_pfams     mtvuy
	#
	# Codon inferences                      Consensus columns
	# codon   inference  std code  diff?    N aligned  N used
	TTT       F          F         N        546        546       
	TTC       F          F         N        41         41        
	TTA       L          L         N        712        712       
	TTG       L          L         N        87         87        
	TCT       S          S         N        306        306   
	...
	GGA       G          G         N        267        267       
	GGG       G          G         N        66         66        
	#
	# Log decoding probabilities
	# codon      logP(A)      logP(C)      logP(D)      logP(E)      logP(F)      logP(G) ...   
	TTT         -1301.6631   -1465.3643   -2182.8130   -1914.5273       0.0000   -2091.3867 ...
	#
	# Final genetic code inference
	FFLLSSSSYY??CCWWL?L?PPPPHHQQ????I?IMTTT?NNKKSSRRV?VVAAAADDEEGGGG

You might choose to change some parameters of the analysis. Some commonly used options are

- Use the `-p` argument to specify a different profile HMM database. Note that the database 
must be located in the `resources/` directory.
- Use `-m -t -v -u -y` to change which groups of problematic Pfam domains are included. For 
instance, if you're analyzing a mitochondrial genome you may want to use `-m`  to include 
Pfams commonly found in mitochondrial genomes. Likewise, use the `-v` and `-t` flags for 
viral genomes, and the `-u` and `-y` flags if you want to include selenocysteine and 
pyrrolysine-containing domains, respectively. By default, all of these domains are excluded.
- You can specify your own name for the inference output file with the `--inference_output` 
argument.
- Use `-e` to change the profile HMM hit e-value threshold (default is 1e-10). Use `-r` 
to change the probability threshold to call an amino acid meaning for a codon (default 
is 0.9999). Use `-f` to change the maximum fraction of observations for a codon coming 
from a single profile HMM position (default is 0.01).
- If you're running several Codetta analyses, you can use `--results_summary` to specify 
a file to which a one-line summary of the results will be appended to.

See the full list of options with `./codetta.py --help`


### Now that you're comfortable...

The basic usage of Codetta is to use `codetta.py`, which rolls the three main steps of 
the analysis into a single program. See the "Starting out..." section above for an 
example and description of the output.

If, for any reason, you want to run the three Codetta analysis steps separately, you
can use the programs `codetta_align.py`, `codetta_summary.py`, and `codetta_infer.py`.

The example in the previous section

	./codetta.py examples/GCA_000442605.1.fna

can be also run by breaking the steps up as

	./codetta_align.py examples/GCA_000442605.1.fna
	./codetta_summary.py examples/GCA_000442605.1.fna
	./codetta_infer.py examples/GCA_000442605.1.fna
	
You might want to do this if:

- You want to infer the genetic code multiple times using different parameters without 
having to re-align the profile HMMs to your input sequence. In this case, you would 
run `codetta_align.py` and `codetta_summary.py` once and then `codetta_infer.py` several
times with different parameters.
- You want to parallelize the computationally intensive alignment step on a computing 
cluster. See the next section.


### Parallelizing the alignment step on a computing cluster 

Codetta is default set up to run locally. However, the `codetta_align` step can 
start to take a long time for longer genomes. For instance, a 4 Mb genome takes 
about 10 minutes on a MacBook Pro. If you're analyzing a large genome or many 
sequences, we recommended parallelizing the analysis on a computing cluster.

This is simple for clusters using a SLURM job scheduler. Install Codetta 
following the usual steps. Then, manually edit the 
`resources/template_jobarray.sh` file to have the right parameters for your 
computing cluster. The `template_jobarray.sh` file looks like this:

	#!/bin/bash
	
	#SBATCH -p [SPECIFY PARTITION NAME]      # partition
	#SBATCH --time=8:30:00                   # wall-clock time (mins:secs)
	#SBATCH -c 1                             # requesting 1 core
	#SBATCH --mem=4000M
	#SBATCH -o /dev/null                     # file to which STDOUT + STDERR will be written

You'll probably only need to specify a partition name.

Then, when you run the `codetta_align` step, add the argument 
`--parallelize_hmmscan 's'`. This will cause the `hmmscan` jobs to be sent as a 
SLURM job array to the specified partition instead of being executed locally. 
Make sure to wait for the SLURM jobs to finish before proceeding to 
`codetta_summary`. That's it!

If your cluster uses a different job scheduler, you will have to manually 
modify the job array template file and the code in `codetta.py` (`processing_genome` 
function) that writes and sends the job array file.


### Building a custom profile HMM database

Pfam domains are expected to align to about ~50% of coding sequence. If you're focused
on a specific clade, you could increase the proportion of aligned sequence by using a 
customized profile HMM database. Starting from a set of multiple sequence alignments 
of expected proteins in the organism of interest, you can build a custom profile HMM
database.

In this example, we will build a profile HMM database tailored for metazoan 
mitochondria which have only 13 protein-coding genes. In the 
`examples/mito-db/` directory, you can find multiple sequence alignment files 
in Stockholm format for each of the 13 genes. *  

The first step is to use `hmmbuild` to create profile HMMs from each of the 
alignment files.

	cd examples/mito-db/
	ls metazoan_mito*.msa | xargs -I {} ../../hmmer-3.3.2/bin/hmmbuild --enone {}.hmm {}

Then we concatenate all of these `.hmm` files into a single database and 
finally run `hmmpress` to finish

	cat metazoan_mito*.hmm > metazoan_mito_proteins.hmm
	../../hmmer-3.3.2/bin/hmmpress metazoan_mito_proteins.hmm
	mv metazoan_mito_proteins.hmm* ../../resources/.

At the end, the profile HMM database is moved into the resources directory, 
where Codetta expects to find all profile HMM databases. When running a Codetta 
analysis, you will have to specify the custom profile HMM database 
with `-p [name of profile HMM db]`.

*These alignments were created by searching the human mitochondrial protein 
sequences against SwissProt with `jackhmmer`


### Bonus: downloading nucleotide sequences from GenBank

We have also provided a simple program for a downloading FASTA file from 
GenBank by specifying either a genome assembly accession or a nucleotide 
accession. 

- `codetta_download`: Download a genome assembly or nucleotide sequence from 
GenBank

We can download the _Nasuia deltocephalinicola_ genome using the Genbank 
genome assembly accession by

	./codetta_download.py GCA_000442605.1 a --sequence_file examples/GCA_000442605.1.fna

This will download a FASTA file containing the GCA_000442605.1 sequence into 
`examples/GCA_000442605.1.fna`. The argument `a` specifies that this is an assembly 
database accession and not a nucleotide accession (which would be `c`).

We can download the mitochondrial genome of the green algae 
_Pycnococcus provasolii_, which is under NCBI nucleotide accession NC_013935.1 by

	./codetta_download.py NC_013935.1 c --sequence_file examples/NC_013935.1.fna

Notice the argument `c` specifying that this is a nucleotide 
database accession.


### Summary, with one more example

Now let's pull it all together by predicting the genetic code of the 
_P. provasolii_ mitochondrial genome:

	./codetta_download.py NC_013935.1 c --sequence_file examples/NC_013935.1.fna
	./codetta.py examples/NC_013935.1.fna -m

Notice how we specified the `-m` argument to indicate that we do not want to 
exclude Pfam domains associated with mitochondrial genomes. 

Alternatively, we could also run the same analysis as 

	./codetta_download.py NC_013935.1 c --sequence_file examples/NC_013935.1.fna
	./codetta_align.py examples/NC_013935.1.fna --align_output examples/Pprovasolii_mito
	./codetta_summary.py examples/Pprovasolii_mito
	./codetta_infer.py examples/Pprovasolii_mito -m --inference_output examples/Pprovasolii_mito_Pfam_genetic_code.out

Here we are using the `--align_output` argument to write the 
alignment output files with a more informative prefix `Pprovasolii_mito` 
and the `--inference_output` argument to write the inference output file 
`examples/Pprovasolii_mito_Pfam_genetic_code.out`.

The output genetic code is:

	FF??SSSSYY??CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	
Comparing to the standard genetic code (below), you can see that two codons 
have alternative meanings: the stop codon UGA is now tryptophan codon and the 
isoleucine codon AUA is now a methionine codon. Some codons are uninferred (?) 
due to few aligned Pfam consensus columns (look at the inference output file 
for more detail).

	P mt code : FF??SSSSYY??CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	std code  : FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	                          ^                   ^

This alternative genetic code in _Pycnococcus_ mitochondria is summarized in 
[Noutahi et al (2019)](https://pubmed.ncbi.nlm.nih.gov/30698742/).










