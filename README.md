# Codetta

Codetta is a Python program for predicting the genetic code (codon table) of an organism from nucleotide sequence data.   

The analysis consists of three steps:

1. Aligning the input nucleotide sequence a database of profile Hidden Markov models (HMMs) of proteins (such as the Pfam database)
2. Generating an alignment summary file
3. Inferring the genetic code from the alignment summary file

Step 1 (profile HMM alignment) is by far the more computationally intensive step of the analysis. If you plan to analyze many genomes (or large genomes), we recommend parallelizing across many machines on a computing cluster. We provide instructions on how to do this below.

If you are looking to reproduce results from [Shulgina & Eddy (2021)](https://www.biorxiv.org/content/10.1101/2021.06.18.448887v1), please follow the README at [this git tag](https://github.com/kshulgina/codetta/releases/tag/v1.0) which has more specific instructions. We can also provide the alignment summary files for any genomes analyzed in Shulgina & Eddy (2021) upon request. 

## Download and setup

### Cloning the GitHub repo
Clone the Codetta repository from GitHub with the command 

	git clone https://github.com/kshulgina/codetta
	cd codetta

### Python version and packages
Codetta was developed for Python version 3.7-3.9 on Linux and MacOS. 

Type `python --version` into your terminal to check which version of Python you have. If you don't want to update your version of Python, you try using `conda` to create a virtual Python 3.9 environment using the commands 
	
	conda create --name py39 python=3.9
	source activate py39

To ensure that the correct Python package versions are installed, use the command

	python -m pip install -r requirements.txt

Otherwise, you can manually install the packages listed in the `requirements.txt` file.


### Additional requirements
Codetta additionally requires:

- `wget` and `gzip`: on Mac, use install commands `brew install wget` and `brew install gzip`. For Linux, you'll have to use your system's package management tool. 
- HMMER (v3.1b2) and Easel library: the commands shown below will install these programs into `codetta/hmmer-3.1b2`. For more detail on installation, see the [HMMER user's guide](http://eddylab.org/software/hmmer/Userguide.pdf). This specific version of HMMER (v3.1b2) is required to reproduce the results in Shulgina & Eddy (2021). Codetta may work with other versions of HMMER3, as long as the same version of HMMER is used to build the Pfam database. 

		wget http://eddylab.org/software/hmmer/hmmer-3.1b2.tar.gz
		tar xf hmmer-3.1b2.tar.gz
		rm hmmer-3.1b2.tar.gz 
		cd hmmer-3.1b2
		pwd | xargs -I {} ./configure --prefix={}
		make
		make install
		cd easel; make install
		cd ../..

	Codetta will expect to find the HMMER and Easel programs in the directory `codetta/hmmer-3.1b2/bin` unless otherwise specified as an argument.


If you plan on analyzing your own nucleotide sequences (Step 1), then you will also need:

- `gtar`: on Mac, use install command `brew install gnutar`. `gtar` is the default version of tar on most Linux machines. You can check which version of `tar` you have by typing `man tar` and looking at the first line. If you have gnutar but the command `gtar` does not work, you can map it by adding an alias `alias gtar='tar'` in your `~/.bashrc` file (or equivalent for your shell) and restarting your terminal.

### Building a local version of the Pfam database
You also need to download and build a local version of the Pfam database.

Download Pfam database into the `resources` directory. This may take a few minutes because this a ~19 Mb file.

	cd resources
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz
	gunzip Pfam-A.seed.gz

Then, use HMMER to build a searchable database, using the `--enone` flag to turn off entropy weighting. This will also take a few minutes. This process creates 3 Gb worth of files, so make sure you have sufficient disk space. If you use either a different version of HMMER or a different version of the Pfam database, then you will not be able to use the alignment summary files provided by us without unexpected errors.

	../hmmer-3.1b2/bin/hmmbuild --enone Pfam-A_enone.hmm Pfam-A.seed
	rm Pfam-A.seed
	../hmmer-3.1b2/bin/hmmpress Pfam-A_enone.hmm
	cd ..

Now you're ready to predict some genetic codes!

## Usage

Codetta consists of three main programs.

- `codetta_align`: Align profile HMMs to the input nucleotide sequence.
- `codetta_summary`: Summarize profile HMM alignments into an alignment summary file.
- `codetta_infer`: Infer the genetic code from the alignment summary file.

General usage for these programs is

	python [program name].py [input file prefix] [optional arguments]

For any of these programs, type `python [program name].py --help` for complete usage details.

If you want to be able to run Codetta from anywhere on your machine (without having to invoke python):

	chmod +x codetta_*.py
	export PATH=$PATH:$(pwd)

## Tutorial

### Super quick example

The following commands will predict the genetic code of bacteria _Nasuia deltocephalinicola_, whose genome can be found in `examples/GCA_000442605.1.fna`.

	python codetta_align.py examples/GCA_000442605.1
	python codetta_summary.py examples/GCA_000442605.1
	python codetta_infer.py examples/GCA_000442605.1

Notice that we specify the prefix of the genome, without the `.fna` file extension. The output genetic code (in a one-line representation) is:

	FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG

An output file with a detailed summary of the analysis can be found at `examples/GCA_001661245.1.inference_output_1e-10_0.9999_0.01_excl-mtvuy.out`. The long file extension specifies the inference parameters.

### Example with more explanations

In the `examples` directory, there is a FASTA file called `GCA_000442605.1.fna` with the _Nasuia deltocephalinicola_ genome sequence (assembly accession GCA_000442605.1). _N. deltocephalinicola_ is known to have reassigned the canonical stop codon UGA to tryptophan. Let's see if we can predict this reassignment. 

The first step is to create a six-frame standard genetic code translation of the genome and align it to the entire Pfam database. We can do this with

	python codetta_align.py examples/GCA_000442605.1

Notice that the input argument is the prefix of the fasta file, without the `.fna`. The input nucleotide sequence must a valid FASTA file as a DNA sequence (T instead of U).

This step may take a while, depending on the size of the input nucleotide sequence. Rough estimate of about an hour on a single CPU core to analyze a typical 6 Mb bacterial genome. However, _N. deltocephalinicola_ has a small 112 Kb genome, so this will take only a minute.

If you intend to analyze many sequences (or longer sequences), we recommend parallelizing the computationally-intensive `hmmscan` step of the analysis over many machines on a computing cluster. To do this, you will need to modify the code in the `codetta.py` file. Instructions can be found in the comment at the bottom of the `hmmscan_jobs()` function. 

This Python program creates several files in the `examples` directory with the specified prefix, which are used by the subsequent step, which is to generate an alignment summary file. This can be simply done with

	python codetta_summary.py examples/GCA_000442605.1

This creates an alignment summary file in the `examples` directory with the specified prefix.

We can infer the genetic code of _N. deltocephalinicola_ with default parameters using

	python codetta_infer.py examples/GCA_000442605.1

The output is a one line representation of the genetic code

	FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG

This corresponds to the inferred translation of each of the 64 codons, in order from 'UUU, UUC, UUA, UUG, UCU, UCC, ..., GGA, GGG' (iterating 3rd, 2nd, then 1st base through UCAG). This same one line representation of the genetic code is used on the [NCBI Genetic Codes page](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). 

Notice that the 14th codon (corresponding to UGA) is W instead of ?. This means that we have correctly predicted the UGA reassignment to tryptophan in this bacterial genome.

Additionally, a file is created, named `examples/GCA_000442605.1.inference_output_1e-10_0.9999_0.01_excl-mtvuy.out`. The long file extension specifies the inference parameters. This file contains a detailed summary of the genetic code inference results:

	# Analysis arguments
	prefix            examples/GCA_000442605.1
	output_summary    None
	evalue_threshold  1e-10
	prob_threshold    0.9999
	max_fraction      0.01
	excluded_pfams    mtvuy
	#
	# Codon inferences
	# codon   inference   N consensus columns   N column types subsampled
	TTT       F           443                             
	TTC       F           35                              
	TTA       L           575                             
	TTG       L           73                              
	TCT       S           264 
	...
	GGA       G           244                             
	GGG       G           61                              
	#
	# Log decoding probabilities
	# codon      logP(A)      logP(C)      logP(D)      logP(E)      logP(F)      logP(G)   ...
	TTT         -1093.2048   -1227.7571   -1912.3596   -1650.0266       0.0000   -1822.9573 ...  
	...
	#
	# Final genetic code inference
	FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG


If you would like to change the default parameters of the `codetta_infer` analysis, you can do so by specifying the arguments `--evalue` to change the profile HMM hit e-value threshold, `--probability_threshold` to change the probability threshold used to call amino acid inferences, `--max_fraction` to change the maximum fraction a single consensus column can contribute to a single codon, and `-m -t -v -u -y` to change which groups of problematic Pfam domains are excluded.

If you plan on running a large number of analyses, you can use `--results_summary` specify a file to which a one-line summary of the results will be appended to.


### Bonus: downloading nucleotide sequences from GenBank

We have also provided a simple program for a downloading FASTA file from GenBank by specifying either a genome assembly accession or a nucleotide accession. 

- `codetta_download`: Download a genome assembly or nucleotide sequence from GenBank

Let's use this to download the mitochondrial genome of the green algae _Pycnococcus provasolii_, which is under NCBI nucleotide accession GQ497137.1

	python codetta_download.py GQ497137.1 c --prefix examples/GQ497137.1

This will download a FASTA file containing the GQ497137.1 sequence into `examples/GQ497137.1.fna`. Again, notice that the `.fna` file extension should not included in the `--prefix` argument. The argument `c` specifies that this is a nucleotide database accession and not an assembly accession (which would be `a`).

### Summary

Now let's pull it all together by predicting the genetic code of the _P. provasolii_ mitochondrial genome from scratch:

	python codetta_download.py GQ497137.1 c --prefix examples/GQ497137.1
	python codetta_align.py examples/GQ497137.1
	python codetta_summary.py examples/GQ497137.1
	python codetta_infer.py examples/GQ497137.1 -m

The `-m` argument indicates that we do not want to exclude Pfam domains associated with mitochondrial genomes. The output genetic code is:

	FF??S?SSYY??CCWWLLLLP?PPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG

Comparing to the standard genetic code (below), you can see that two codons have alternative meanings: the stop codon UGA is now tryptophan codon and the isoleucine codon AUA is now a methionine codon. Some codons are uninferred (?) due to few aligned Pfam consensus columns (look at the inference output file for more detail).

	P. provasolii mt code : FF??S?SSYY??C?WWLLLLP?PPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	standard genetic code : FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	                                      ^                   ^

This alternative genetic code in _Pycnococcus_ mitochondria has been previously reported by [Noutahi et al (2019)](https://pubmed.ncbi.nlm.nih.gov/30698742/).







