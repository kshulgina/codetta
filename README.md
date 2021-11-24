# Codetta

WARNING: this branch is under active development, use at your own risk.

## Description

Codetta is a Python program for predicting the genetic code (codon table) of an 
organism from nucleotide sequence data.   

The analysis consists of three steps:

1. Aligning the input nucleotide sequence a database of profile Hidden Markov 
models (HMMs) of proteins (such as the Pfam database)
2. Generating an alignment summary file
3. Inferring the genetic code from the alignment summary file

Step 1 (profile HMM alignment) is by far the more computationally intensive 
step of the analysis. By default, the alignment is set up to run locally, which 
for a 4 Mb genome would take about 10 minutes (on a MacBook Pro). If you plan 
to analyze many genomes (or large genomes), we recommend parallelizing across 
many machines on a computing cluster. We provide instructions on how to do this 
below.

If you are looking to reproduce results from [Shulgina & Eddy (2021)]
(https://elifesciences.org/articles/71402), please follow the README for 
[Codetta v1.0](https://github.com/kshulgina/codetta/releases/tag/v1.0). We can 
also provide the alignment summary files for any genomes analyzed in Shulgina & 
Eddy (2021) upon request.

If you encounter any problems in the installation or usage of Codetta, please 
leave a Github issue or email me at shulgina@g.harvard.edu!

## Download and setup

<!---
### Cloning the GitHub repo
Clone the Codetta repository from GitHub with the command 

	git clone https://github.com/kshulgina/codetta
	cd codetta
--->

### Python version and packages
Codetta was developed for Python versions 3.5+ on Linux and MacOS. Required 
Python packages are `numpy` and `scipy`.

Set up Codetta by running the setup script
	
	./setup.sh

This will check Python requirements and install a local version of HMMER.

### Additional requirements
Codetta additionally requires:

- `wget` and `gzip`: on Mac, use install commands `brew install wget` and `brew 
install gzip`. For Linux, you'll have to use your system's package management 
tool.
 
<!---
- HMMER v3 and Easel library: the commands shown below will install these 
programs into `codetta/hmmer-3.3.2`. For more detail on installation, see the 
[HMMER user's guide](http://eddylab.org/software/hmmer/Userguide.pdf).

		wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz
		tar xf hmmer-3.3.2.tar.gz
		rm hmmer-3.3.2.tar.gz 
		cd hmmer-3.3.2
		pwd | xargs -I {} ./configure --prefix={}
		make
		make install
		cd easel; make install
		cd ../..

	Codetta will expect to find the HMMER and Easel programs in the directory 
	`codetta/hmmer-3.3.2/bin`.
--->

### Building a local version of the Pfam database
By default, Codetta will assume that the Pfam database is the source of profile 
HMMs. Unless you intend to use a custom profile HMM database, you will need to 
download and build a local version of the Pfam database.

Download current Pfam database into the `resources` directory. This may take a few 
minutes because this a ~140 Mb file.

	cd resources
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz
	gunzip Pfam-A.seed.gz

Then, use HMMER to build a searchable database, using the `--enone` flag to 
turn off entropy weighting. This will also take a few minutes. This process 
creates 3 Gb worth of files, so make sure you have sufficient disk space.

	../hmmer-3.3.2/bin/hmmbuild --enone Pfam-A_enone.hmm Pfam-A.seed
	../hmmer-3.3.2/bin/hmmpress Pfam-A_enone.hmm
	rm Pfam-A.seed
	cd ..

Pfam 35.0 download from website.

### Building a custom profile HMM database
[TBD]

<!---

To make a custom profile HMM database from a set of multiple sequence 
alignments, you need to follow a similar series of steps.

[DRAFT]

Here is an example showing how to build a custom profile HMM database from a 
set of multiple sequence alignments. As an example, we will use metazoan 
mitochondria, which have only 12 protein coding genes. In the 
`codetta\examples` directory, you can find a multiple sequence alignment files 
in Stockholm format for each of the 12 genes. 

The first step is to use `hmmbuild` to create profile HMMs from each of the 
alignment files.

	cd examples
	ls metazoan_mito*.msa | xargs -I {} hmmbuild --enone {}

Then we concatenate all of these `.hmm` files into a single database and 
finally run `hmmpress` to finish

	cat metazoan_mito*.hmm > metazoan_mito_proteins.hmm
	hmmpress metazoan_mito_proteins.hmm

When running a Codetta analysis, you will have to specify the custom profile 
HMM database with the `-p` option.
-->

## Usage

Codetta consists of three main programs.

- `codetta_align`: Align profile HMMs to the input nucleotide sequence.
- `codetta_summary`: Summarize profile HMM alignments into an alignment summary 
file.
- `codetta_infer`: Infer the genetic code from the alignment summary file.

General usage for these programs is

	./[program name].py [input file or prefix] [optional arguments]

For any of these programs, type `./[program name].py --help` for complete usage details.

## Tutorial

### Super quick example

The following commands will predict the genetic code of bacteria _Nasuia 
deltocephalinicola_, whose genome can be found in 
`examples/GCA_000442605.1.fna`.

	./codetta_align.py examples/GCA_000442605.1.fna
	./codetta_summary.py examples/GCA_000442605.1.fna
	./codetta_infer.py examples/GCA_000442605.1.fna

The output genetic code (in a one-line representation) is:

	FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTT?NNKKS?RRVVVVAAAADDEEGGGG

This corresponds to the inferred translation of each of the 64 codons, in order 
from 'UUU, UUC, UUA, UUG, UCU, UCC, ..., GGA, GGG' (iterating 3rd, 2nd, then 
1st base through UCAG). 

An output file with a detailed summary of the analysis can be found at 
`examples/GCA_000442605.1.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_excl-mtvuy.genetic_code.out`. 
The long file extension specifies the inference parameters.

### Example with more explanations

In the `examples` directory, there is a FASTA file called `GCA_000442605.1.fna` 
with the _Nasuia deltocephalinicola_ genome sequence (assembly accession 
GCA_000442605.1). _N. deltocephalinicola_ is known to have reassigned the 
canonical stop codon UGA to tryptophan. Let's see if we can predict this 
reassignment. 

The first step is to create a six-frame standard genetic code translation of 
the genome and align it to the entire Pfam database. (To align against a custom 
profile HMM database, use the `-p` argument.) We can do this with

	./codetta_align.py examples/GCA_000442605.1.fna

The input nucleotide sequence must a valid FASTA file as a DNA sequence 
(T instead of U).

This step may take a while, depending on the size of the input nucleotide 
sequence. Rough estimate of about an hour on a single CPU core to analyze a 
typical 6 Mb bacterial genome. However, _N. deltocephalinicola_ has a small 112 
Kb genome, so this will take only a minute.

If you intend to analyze many sequences (or longer sequences), we recommend 
parallelizing the computationally-intensive `hmmscan` step of the analysis over 
many machines on a computing cluster. Instructions can be found in the next 
section.

This Python program creates several files which are used by the subsequent step 
to generate an alignment summary file. The default location of these files is 
the same as the input sequence file, with different file extension. However, an 
alternative location for the alignment output files can be specified with the 
`--align_output` argument.

The next step is to process these files into an alignment summary file. The 
only required argument is the location of the alignment output files (either 
the default or what was specified by the `--align_output` argument.

	./codetta_summary.py examples/GCA_000442605.1.fna

Then, we can infer the genetic code of _N. deltocephalinicola_ with default 
parameters using

	./codetta_infer.py examples/GCA_000442605.1.fna

The output is a one line representation of the genetic code

	Genetic code: FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTT?NNKKS?RRVVVVAAAADDEEGGGG

This corresponds to the inferred translation of each of the 64 codons, in order 
from 'UUU, UUC, UUA, UUG, UCU, UCC, ..., GGA, GGG' (iterating 3rd, 2nd, then 
1st base through UCAG). This same one line representation of the genetic code 
is used on the [NCBI Genetic Codes page]
(https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). 

Notice that the 14th codon (corresponding to UGA) is W instead of ?. This means 
that we have correctly predicted the UGA reassignment to tryptophan in this 
bacterial genome.

Additionally, a file with detailed information about the run is created, named 
`examples/GCA_000442605.1.fna.Pfam-A_enone.hmm.1e-10_0.9999_0.01_
excl-mtvuy.genetic_code.out`. The long file extension specifies the inference 
parameters. You can specify an alternative output file name using the 
`--inference_output` argument. 

This file contains a detailed summary of the genetic code inference results:

	# Analysis arguments
	alignment_prefix   examples/GCA_000442605.1.fna
	profile_database   Pfam-A_enone.hmm
	output_summary     None
	evalue_threshold   1e-10
	prob_threshold     0.9999
	max_fraction       0.01
	excluded_pfams     mtvuy
	#
	# Codon inferences                      Consensus columns
	# codon   inference  std code  diff?    N aligned  N used
	TTT       F          F         N        433        433       
	TTC       F          F         N        36         36        
	TTA       L          L         N        572        572       
	TTG       L          L         N        77         77        
	TCT       S          S         N        257        257    
	...
	GGA       G          G         N        236        236       
	GGG       G          G         N        60         60                               
	#
	# Log decoding probabilities
	# codon      logP(A)      logP(C)      logP(D)      logP(E)      logP(F)      logP(G)   ...
	TTT         -1015.4800   -1150.8584   -1756.1719   -1531.7520       0.0000   -1670.3149 ...  
	...
	#
	# Final genetic code inference
	FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG


If you would like to change the default parameters of the `codetta_infer` 
analysis, you can do so by specifying the arguments `--evalue` to change the 
profile HMM hit e-value threshold, `--probability_threshold` to change the 
probability threshold used to call amino acid inferences, `--max_fraction` to 
change the maximum fraction a single consensus column can contribute to a 
single codon, and `-m -t -v -u -y` to change which groups of problematic Pfam 
domains are excluded. You can specify your own name for the inference output 
file with the `--inference_output` argument.

Other options you might choose to use:

- If you want to use a custom profile HMM database, specify it using the 
`-p` argument in `codetta_align`, `codetta_summary`, and `codetta_infer`
- If you're analyzing a large genome or many genomes on a computing cluster, 
you can parallelize the `codetta_align` step with the `--parallelize_hmmscan` 
argument (see more detail below)
- If you plan on running a large number of analyses, you can use 
`--results_summary` for `codetta_infer` to specify a file to which a one-line 
summary of the results will be appended to.
- If you're analyzing a mitochondrial genome, remember to use the `-m` flag in 
`codetta_infer` to turn off the exclusion of common mitochondrial Pfam domains. 
Likewise, use the `-v` and `-t` flags for viral genomes, and the `-u` and `-y` 
flags if you want to include selenocysteine and pyrrolysine-containing domains, 
respectively.

### Parallelizing hmmscan jobs on a computing cluster 

Codetta is default set up to run locally. However, the `codetta_align` step can 
start to take a long time for longer genomes. For instance, a 4 Mb genome takes 
about 10 minutes on a MacBook Pro. If you're analyzing a large genome or many 
sequences, we recommended parallelizing the analysis on a computing cluster.

This is simple for clusters using a SLURM job scheduler. Install Codetta 
following the same directions as above. Then, manually edit the 
`resources/template_jobarray.sh` file to have the right parameters for your 
computing cluster. The `template_jobarray.sh` file looks like this:

	#!/bin/bash
	                                                                                                                                                                                              
	#SBATCH -p [SPECIFY PARTITION NAME]      # queue
	#SBATCH --time=8:30:00                   # wall-clock time (mins:secs) 
	#SBATCH -c 1                             # requesting 1 core      
	#SBATCH --mem=4000M
	#SBATCH -o /dev/null                     # File to which STDOUT + STDERR 
	will be written  

You'll probably only need to specify a partition name.

Then, when you run the `codetta_align` step, add the argument 
`--parallelize_hmmscan 's'`. This will cause the `hmmscan` jobs to be sent as a 
SLURM job array to the specified partition instead of being executed locally. 
Make sure to wait for the SLURM jobs to finish before proceeding to 
`codetta_summary`. That's it!

If your cluster uses a different job scheduler, you will have to manually 
modify the job array template file and the code in `codetta.py` (`hmmscan_jobs` 
function) that writes and sends the job array file.

### Bonus: downloading nucleotide sequences from GenBank

We have also provided a simple program for a downloading FASTA file from 
GenBank by specifying either a genome assembly accession or a nucleotide 
accession. 

- `codetta_download`: Download a genome assembly or nucleotide sequence from 
GenBank

Let's use this to download the mitochondrial genome of the green algae 
_Pycnococcus provasolii_, which is under NCBI nucleotide accession GQ497137.1

	./codetta_download.py GQ497137.1 c --sequence_file examples/GQ497137.1.fna

This will download a FASTA file containing the GQ497137.1 sequence into 
`examples/GQ497137.1.fna`. The argument `c` specifies that this is a nucleotide 
database accession and not an assembly accession (which would be `a`).

### Summary, with one more example

Now let's pull it all together by predicting the genetic code of the 
_P. provasolii_ mitochondrial genome:

	./codetta_align.py examples/GQ497137.1.fna --align_output examples/Pprovasolii_mito
	./codetta_summary.py examples/Pprovasolii_mito
	./codetta_infer.py examples/Pprovasolii_mito -m --inference_output examples/Pprovasolii_mito_Pfam_genetic_code.out

Notice how we specified that the alignment output files are written with a more 
informative prefix `Pprovasolii_mito` and the inference output file is written 
to `examples/Pprovasolii_mito_Pfam_genetic_code.out`.
The `-m` argument indicates that we do not want to exclude Pfam domains 
associated with mitochondrial genomes. The output genetic code is:

	FF??S?SSYY??CCWWLLLLP?PPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	
Comparing to the standard genetic code (below), you can see that two codons 
have alternative meanings: the stop codon UGA is now tryptophan codon and the 
isoleucine codon AUA is now a methionine codon. Some codons are uninferred (?) 
due to few aligned Pfam consensus columns (look at the inference output file 
for more detail).

	P mt code : FF??S?SSYY??CCWWLLLLP?PPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	std code  : FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	                          ^                   ^

This alternative genetic code in _Pycnococcus_ mitochondria is summarized in 
[Noutahi et al (2019)](https://pubmed.ncbi.nlm.nih.gov/30698742/).










