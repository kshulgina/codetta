# Codetta

Codetta is a Python program for cracking the genetic code of an organism from nucleotide sequence data.   

The analysis consists of two main steps:

1. Aligning the input nucleotide sequence to the entire Pfam database, generating an alignment summary file.
2. Inferring the genetic code from the alignment summary file.

Step 1 (Pfam alignment) is by far the more computationally intensive step of the analysis. We can provide the alignment summary files for any subset of the bacterial and archaeal genomes analyzed in Shulgina & Eddy (2021) upon request.

## Download and setup

### Cloning the GitHub repo
Clone the Codetta repository from GitHub with the command 

	git clone https://github.com/kshulgina/codetta
	cd codetta

### Python version and packages
Codetta was developed for Python version 3.7-3.9 on Linux and MacOS. 

Type `python --version` into your terminal to check which version of Python you have. If you don't want to update your version of Python, you try using `conda` to create a virtual Python 3.9 environment using the commands 
	
	conda create --name py39 python=3.9
	conda activate py39

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


If you plan on analyzing your own nucleotide sequences, then you will also need:

- `gtar`: on Mac, use install command `brew install gnutar`. `gtar` is the default version of tar on most Linux machines. You can check which version of `tar` you have by typing `man tar` and looking at the first line. If you have gnutar but the command `gtar` does not work, you can map it by adding an alias `alias gtar='tar'` in your `~/.bashrc` file (or equivalent for your shell) and restarting your terminal.

### Building a local version of the Pfam database
You also need to download and build a local version of the Pfam database. We used Pfam version 32.0.

Download Pfam database into the `resources` directory. This may take a few minutes because this a ~19 Mb file.

	cd resources
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.seed.gz
	gunzip Pfam-A.seed.gz

Then, use HMMER to build a searchable database, using the `--enone` flag to turn off entropy weighting. This will also take a few minutes. This process creates 3 Gb worth of files, so make sure you have sufficient disk space. If you use either a different version of HMMER or a different version of the Pfam database, then you will not be able to use the alignment summary files provided by us without unexpected errors.

	../hmmer-3.1b2/bin/hmmbuild --enone Pfam-A_enone.hmm Pfam-A.seed
	rm Pfam-A.seed
	../hmmer-3.1b2/bin/hmmpress Pfam-A_enone.hmm
	cd ..

Now you're ready to predict some genetic codes!

## Usage

Codetta consists of three main programs. To infer the genetic code from an alignment summary file, you will use:

- `codetta_infer`: Infer the genetic code from the alignment summary file.

To analyze your own nucleotide sequence and generate an alignment summary file, you will use:

- `codetta_align`: Align Pfam profile HMMs to the input nucleotide sequence.
- `codetta_summary`: Summarize Pfam alignments into an alignment summary file.

General usage for these programs is

	python [program name].py [input file prefix] [optional arguments]

For any of these programs, type `python [program name] --help` for complete usage details.

If you want to be able to run Codetta from anywhere on your machine (without having to invoke python):

	chmod +x codetta_*.py
	export PATH=$PATH:$(pwd)

## Tutorial

Let's see a few examples of how Codetta is used.

### Genetic code inference from an alignment summary file

In the `examples` directory, you will find `GCA_001661245.1.hmmscan_summary.txt.gz`, an alignment summary file for the yeast _Pachysolen tannophilus_ (GenBank assembly GCA_001661245.1). This file summarizes the result of aligning the entire Pfam database against a six-frame translation of the entire genome.

_P. tannophilus_ is known to have reassigned the canonical leucine codon CUG to alanine. Let's see if we can predict this reassignment.

We can infer the genetic code of _P. tannophilus_ with default parameters using

	python codetta_infer.py examples/GCA_001661245.1

Notice that the input argument is the prefix of the alignment summary file, without a file extension.

The output is a one line representation of the genetic code

	FFLLSSSSYY??CC?WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG

This corresponds to the inferred translation of each of the 64 codons, in order from 'UUU, UUC, UUA, UUG, UCU, UCC, ..., GGA, GGG' (iterating 3rd, 2nd, then 1st base through UCAG). This same one line representation of the genetic code is used on the [NCBI Genetic Codes page](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). 

Notice that the 19th codon (corresponding to CUG) is A instead of L. This means that we have correctly predicted the CUG reassignment to alanine in this yeast genome.

Additionally, a file is created, named `examples/GCA_001661245.1.inference_output_1e-10_0.9999_0.01_excl-mtvuy.out`. The long file extension specifies the inference parameters. This file contains a detailed summary of the genetic code inference results:

	# Analysis arguments
	# arguments
	prefix            examples/GCA_001661245.1
	output_summary    None
	evalue_threshold  1e-10
	prob_threshold    0.9999
	max_fraction      0.01
	excluded_pfams    mtvuy
	#
	# Codon inferences
	# codon   inference   N consensus columns   N column types subsampled
	TTT       F           21345                           
	TTC       F           6945                            
	TTA       L           34424                           
	TTG       L           16958                           
	TCT       S           13015    
	...
	GGA       G           9874                            
	GGG       G           2725                            
	#
	# Log decoding probabilities
	# codon      logP(A)      logP(C)      logP(D)      logP(E)      logP(F)      logP(G) ...
	TTT        -74503.3438  -76349.0312 -113527.5000 -105288.6562       0.0000 -101927.0625 ...
	...
	#
	# Final genetic code inference
	FFLLSSSSYY??CC?WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG


If you would like to change the default parameters of the analysis, you can do so by specifying the arguments `--evalue` to change the profile HMM hit e-value threshold, `--probability_threshold` to change the probability threshold used to call amino acid inferences, `--max_fraction` to change the maximum fraction a single consensus column can contribute to a single codon, and `-m -t -v -u -y` to change which groups of problematic Pfam domains are excluded.

If you plan on running a large number of analyses, you can use `--results_summary` specify a file to which a one-line summary of the results will be appended to.


### Genetic code inference from a nucleotide sequence

Generating an alignment summary file from a nucleotide sequence requires two additional steps.

In the `examples` directory, there is a FASTA file called `GCA_000442605.1.fna` with the _Nasuia deltocephalinicola_ genome sequence (assembly accession GCA_000442605.1). The input nucleotide sequence must a valid FASTA file as a DNA sequence (T instead of U).

The first step is to create a six-frame standard genetic code translation of the genome and align it to the entire Pfam database. We can do this with

	python codetta_align.py examples/GCA_000442605.1

Notice that the input argument is the prefix of the fasta file, without the `.fna`.

This step may take a while, depending on the size of the input nucleotide sequence. Rough estimate of about an hour on a single CPU core to analyze a typical 6 Mb bacterial genome. However, _N. deltocephalinicola_ has a small 112 Kb genome, so this will take only a minute.

If you intend to analyze many sequences (or longer sequences), we recommend parallelizing the computationally-intensive `hmmscan` step of the analysis over many machines on a computing cluster. To do this, you will need to modify the code in the `codetta.py` file. Instructions can be found in the comment at the bottom of the `hmmscan_jobs()` function. 

This python program creates several files, which are used by the subsequent step:

- `examples/GCA_000442605.1.sequence_pieces.fna`: input nucleotide sequence, broken into pieces <100,000 nt
- `examples/GCA_000442605.1.preliminary_translation.faa`: six-frame standard genetic code translation
- `examples/GCA_000442605.1.temp_files/` is a temporary directory which stores the `hmmscan` result files

Next step is to process the `hmmscan` alignment files into an alignment summary file. This can be simply done with

	python codetta_summary.py examples/GCA_000442605.1

This process creates an alignment summary file is created (called `examples/GCA_000442605.1.hmmscan_summary.txt.gz`) and deletes the temporary directory `examples/GCA_000442605.1.temp_files/`.

To infer the genetic code (as described in the first example), just use 

	python codetta_infer.py examples/GCA_000442605.1

And you will see the output inferred genetic code:

	FFLLSSSSYY??CCW?L?L?PPPPHHQQ????I?IMTTTTNNKKSSRRVVVVAAAADDEEGGGG

And the output file `GCA_000442605.1.inference_output_1e-10_0.9999_0.01_excl-mtvuy.out` is generated.

### Bonus: downloading nucleotide sequences from GenBank

We have also provided a simple program for a downloading FASTA file from GenBank by specifying either a genome assembly accession or a nucleotide accession. 

- `codetta_download`: Download a genome assembly or nucleotide sequence from GenBank

Let's use this to download the mitochondrial genome of the green algae _Pycnococcus provasolii_, which is under NCBI nucleotide accession GQ497137.1

	python codetta_download.py GQ497137.1 c --prefix examples/GQ497137.1

This will download a FASTA file containing the GQ497137.1 sequence into `examples/GQ497137.1.fna`. Again, notice that the `.fna` file extension should not included in the `--prefix` argument. The argument `c` specifies that this is a nucleotide database accession and not an assembly accession (which would be `a`).

### Summary

Now let's pull it all together by predicting the genetic code of the _Pycnococcus provasolii_ mitochondrial genome from scratch:

	python codetta_download.py GQ497137.1 c --prefix examples/GQ497137.1
	python codetta_align.py examples/GQ497137.1
	python codetta_summary.py examples/GQ497137.1
	python codetta_infer.py examples/GQ497137.1 -m

The `-m` argument indicates that we do not want to exclude Pfam domains associated with mitochondrial genomes. The output genetic code is:

	FF??S?SSYY??CCWWLLLLP?PPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG

Comparing to the standard genetic code (below), you can see that two codons have alternative meanings: the stop codon UGA is now tryptophan codon and the isoleucine codon AUA is now a methionine codon. Some codons are uninferred (?) due to few aligned Pfam consensus columns (look at inference output file for more detail).

	P. provasolii mt code : FF??S?SSYY??C?WWLLLLP?PPHHQQRRRRIIMMTTTTNNKKSSR?V?VVAAAADDEEGGGG
	standard genetic code : FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
	                                      ^                   ^

This alternative genetic code in _Pycnococcus_ mitochondria has been previously reported by [Noutahi et al (2019)](https://pubmed.ncbi.nlm.nih.gov/30698742/).

 









