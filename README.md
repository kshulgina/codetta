# Codetta v2.0


## Description

Codetta is a Python program for predicting the genetic code (codon table) of 
an organism from nucleotide sequence data.   

The analysis consists of three steps:

1. Align the input nucleotide sequence to a database of profile Hidden Markov 
models (HMMs) of proteins (such as the Pfam database)
2. Collate the resulting alignments into a single output file
3. Infer the genetic code from the alignment output file

A detailed explanation of the underlying probability model can be found in 
[Shulgina & Eddy (2021)](https://elifesciences.org/articles/71402). If you 
are looking to reproduce results from Shulgina & Eddy (2021), please follow 
the README for 
[Codetta v1.0](https://github.com/kshulgina/codetta/releases/tag/v1.0).

If you encounter any problems in the installation or usage of Codetta, please 
leave a Github issue or email me at shulgina@g.harvard.edu


## Download and setup

Codetta was developed for Python versions 3.5+ on Linux and MacOS. Required 
Python packages are `numpy` and `scipy`.

Once inside the Codetta directory, run the setup script
	
	./setup.sh

This script will check Codetta requirements and set up a local version of 
HMMER. If you are missing `wget` or `gzip`, on MacOS, you can use 
`brew install`. (For Linux, you'll have to use your system's package 
management tool.) Note that on MacOS, you'll also need XCode installed so 
that HMMER can compile.


## Setting up a profile HMM database

Codetta aligns the input sequence to a database of protein profile HMMs. 

By default, Codetta will use the Pfam database. You can download a version of 
Pfam 35.0 specially built for Codetta from our website. From the `codetta`
directory, use these commands to download and uncompress it (note: this will
 take about 3 Gb of disk space)

	cd resources/
	wget http://eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz
	tar xf Pfam-A_enone.tar.gz
	rm Pfam-A_enone.tar.gz
	cd ..

If you want to build your own custom profile HMM database, this is be described in a later section. 


## Usage

The program `codetta.py` rolls the three main analysis steps into a single 
script. Usage for `codetta.py` is

	./codetta.py [input sequence file] [optional arguments]

The three steps can also be run separately using the programs:

- `codetta_align.py`: Align profile HMMs to the input nucleotide sequence.
- `codetta_summary.py`: Summarize profile HMM alignments into an alignment 
output file.
- `codetta_infer.py`: Infer the genetic code from the alignment output file.

General usage for these programs is

	./[program name] [input file or prefix] [optional arguments]

For any of these programs, type `./[program name] --help` for complete usage 
details.


## Tutorials
### Starting out...

The simplest way to run Codetta is by using the `codetta.py` program. 

This program performs the three analysis steps in order. All you have to do 
is specify a nucleotide sequence input file. This file should contain 
nucleotide sequences from a single organism in FASTA format. This can be a 
genome, transcriptome, collection of genes, etc.

Make sure you have first set up a profile HMM database (see above)!

In the `examples/` directory, you will find the file `GCA_014211875.1.fna`
 which contains the genome of the bacterium _Nasuia deltocephalinicola_.

We can predict the genetic code of this bacterium simply by running

	./codetta.py examples/GCA_014211875.1.fna

You will see outputs written to the terminal indicating that each of the 
three Codetta steps is executing. Codetta will create five files in the 
directory containing the input sequence file.

- Processed sequence file: `examples/GCA_014211875.1.fna.sequence_pieces.fna`
- Preliminary translation file + ssi index: 
`examples/GCA_014211875.1.fna.preliminary_translation.faa`
- Alignment output file: 
`examples/GCA_014211875.1.fna.Pfam-A_enone.hmm.alignment_output.txt`
- Inference output file: 
`examples/GCA_014211875.1.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_excl-mtvuy.genetic_code.out`

At the end, the inferred genetic code is printed to the terminal

	Genetic code: FFLLSSSSYY??CCWWL?LLPPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG

This corresponds to the inferred translation of each of the 64 codons, in 
order of 'UUU, UUC, UUA, UUG, UCU, UCC, ..., GGA, GGG' (iterating 3rd, 2nd, 
then 1st base through UCAG).

The ?s correspond to codons that had no inferred amino acid meaning-- this 
means that there was insufficient or ambiguous information about the codon to 
make a confident inference. This is also the expected inference for stop 
codons (since Codetta does not explicitly predict stop codons).

A detailed overview of the analysis can be found at 
`examples/GCA_014211875.1.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_excl-mtvuy.genetic_code.out`. 
The long file extension specifies the inference parameters. You can specify 
your own (simpler) name for the inference output file using the 
`--inference_output` argument.


This file contains a detailed overview of the genetic code inference results:

	# Analysis arguments
	alignment_prefix   examples/GCA_014211875.1.fna
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
	TTC       F          F         N        48         48        
	TTA       L          L         N        703        703       
	TTG       L          L         N        82         82        
	TCT       S          S         N        307        307      
	...
	GGA       G          G         N        248        248       
	GGG       G          G         N        73         72        
	#
	# Log decoding probabilities
	# codon      logP(A)      logP(C)      logP(D)      logP(E)      logP(F)      logP(G)     ...
	TTT         -1250.3530   -1453.8418   -2107.4121   -1850.9243       0.0000   -2034.8423   ...
	...
	#
	# Final genetic code inference
	FFLLSSSSYY??CCWWL?LLPPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG

You might choose to change some parameters of the analysis. Some commonly 
used options are

- Use the `-p` argument to specify a different profile HMM database. Note 
that the database must be located in the `resources/` directory.
- The `excluded_pfams` line in the above output file refers to which groups 
of problematic Pfam domains are excluded from the analysis (`m` 
mitochondrial, `t` transposon and other mobile genetic element, `v` viral, `u`
 selenocysteine-containing, and `y` pyrrolysine-containing). By default, all 
of these groups are excluded. Use `-m -t -v -u -y` to include these groups. 
For instance, if you're analyzing a mitochondrial genome you will want to use 
`-m`  to include Pfams commonly found in mitochondrial genomes. Likewise, use 
the `-v` and `-t` flags when analyzing a viral genome. 
- You can specify your own name for the inference output file with the 
`--inference_output` argument.
- Use `-e` to change the profile HMM hit e-value threshold (default is 
1e-10). Use `-r` to change the probability threshold to call an amino acid 
meaning for a codon (default is 0.9999). Use `-f` to change the maximum 
fraction of observations for a codon coming from a single profile HMM 
position (default is 0.01).
- If you're running several Codetta analyses, you can use `--results_summary` 
to specify a file to which a one-line summary of the results will be appended 
to.

See the full list of options with `./codetta.py --help`


### Now that you're comfortable...

The basic usage of Codetta is to use `codetta.py`, which rolls the three main 
steps of the analysis into a single program. See the "Starting out..." 
section above for an example and description of the output.

If, for any reason, you want to run the three Codetta analysis steps 
separately, you can use the programs `codetta_align.py`, `codetta_summary.py`
, and `codetta_infer.py`.

The example in the previous section

	./codetta.py examples/GCA_014211875.1.fna

can be also run by breaking the steps up as

	./codetta_align.py examples/GCA_014211875.1.fna
	./codetta_summary.py examples/GCA_014211875.1.fna
	./codetta_infer.py examples/GCA_014211875.1.fna
	
You might want to do this if:

- You want to infer the genetic code multiple times using different 
parameters without having to re-align the profile HMMs to your input 
sequence. In this case, you would run `codetta_align.py` and 
`codetta_summary.py` once and then `codetta_infer.py` several times with 
different parameters.
- You want to parallelize the computationally intensive alignment step on a 
computing cluster. See the next section.


### Parallelizing the alignment step on a computing cluster 

Codetta is default set up to run locally. However, the `codetta_align` step 
can start to take a long time for longer genomes. For instance, a 4 Mb 
bacterial genome takes about 10 minutes on a MacBook Pro. If you're analyzing 
a large genome or many sequences, we recommended parallelizing the 
`codetta_align` step on a computing cluster.

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
`--parallelize_hmmscan 's'`. This will cause the `hmmscan` jobs to be sent as 
a SLURM job array to the specified partition instead of being executed 
locally. Make sure to wait for the SLURM jobs to finish before proceeding to 
`codetta_summary`. That's it!

If your cluster uses a different job scheduler, you will have to manually 
modify the job array template file and the code in `codetta.py` 
(`processing_genome` function) that writes and sends the job array file.


### Building a custom profile HMM database

Pfam domains are expected to align to 
[about ~50% of coding sequence](https://pubmed.ncbi.nlm.nih.gov/30357350/). 
If you're focused on a specific clade, you could increase the proportion of 
aligned sequence by using a customized profile HMM database. Starting from a 
set of multiple sequence alignments of expected proteins in the organism of 
interest, you can build a custom profile HMM database.

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
where Codetta expects to find all profile HMM databases. When running a 
Codetta analysis, you will have to specify the custom profile HMM database 
with `-p [name of profile HMM db]`.

\*These alignments were created by searching the human mitochondrial protein 
sequences against SwissProt with `jackhmmer`


### Bonus: downloading nucleotide sequences from GenBank

We have also provided a simple program for a downloading FASTA file from 
GenBank by specifying either a genome assembly accession or a nucleotide 
accession. 

- `codetta_download`: Download a genome assembly or nucleotide sequence from 
GenBank

We can download the _Nasuia deltocephalinicola_ genome using the Genbank 
genome assembly accession by

	./codetta_download.py GCA_014211875.1 a --sequence_file examples/GCA_014211875.1.fna

This will download a FASTA file containing the GCA_014211875.1 sequence into 
`examples/GCA_014211875.1.fna`. The argument `a` specifies that this is an 
assembly database accession and not a nucleotide accession (which would be 
`c`).

We can download the mitochondrial genome of the green algae 
_Pycnococcus provasolii_, which is under NCBI nucleotide accession 
NC_013935.1 by

	./codetta_download.py NC_013935.1 c --sequence_file examples/NC_013935.1.fna

Notice the argument `c` specifying that this is a nucleotide database 
accession.


### Summary, with one more example

Now let's pull it all together by predicting the genetic code of the 
_P. provasolii_ mitochondrial genome:

	./codetta_download.py NC_013935.1 c --sequence_file examples/NC_013935.1.fna
	./codetta.py examples/NC_013935.1.fna -m

Notice how we specified the `-m` argument to indicate that we do not want to 
exclude Pfam domains associated with mitochondrial genomes, which is expected to have two codon reassignments. 

Alternatively, we could also run the same analysis as 

	./codetta_download.py NC_013935.1 c --sequence_file examples/NC_013935.1.fna
	./codetta_align.py examples/NC_013935.1.fna --align_output examples/Pprovasolii_mito
	./codetta_summary.py examples/Pprovasolii_mito
	./codetta_infer.py examples/Pprovasolii_mito -m --inference_output examples/Pprovasolii_mito_Pfam_genetic_code.out

Here we are also showing how to use the `--align_output` argument to write 
the alignment output files with the more informative prefix `Pprovasolii_mito`
 and the `--inference_output` argument to write the inference output file 
 `examples/Pprovasolii_mito_Pfam_genetic_code.out`.

The output genetic code is:

	FF??SSSSYY??CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	
Comparing to the standard genetic code (below), you can see that two codons 
have alternative meanings: the stop codon UGA is now tryptophan codon and the 
isoleucine codon AUA is now a methionine codon. Some codons are uninferred 
(?) due to few aligned Pfam consensus columns (look at the inference output 
file for more detail).

	P mt code : FF??SSSSYY??CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	std code  : FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	                          ^                   ^

This alternative genetic code in _Pycnococcus_ mitochondria is summarized in 
[Noutahi et al (2019)](https://pubmed.ncbi.nlm.nih.gov/30698742/).










