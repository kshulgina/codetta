# Codetta

Codetta is a Python program for cracking the genetic code of an organism from nucleotide sequence data.   

The analysis can be performed either 

- directly from an alignment summary file (easily downloaded from link) or 
- from your own nucleotide sequence (more complicated and computationally intensive). 

We provide the alignment summary files for over 250,000 bacterial and archaeal genomes analyzed in Shulgina & Eddy (2021). See our preprint for more detail: [link to preprint]


## Download and dependencies

Codetta was developed for Python version 3.x-3.9 on Linux and MacOS. Type `python --version` into your terminal to check which version of Python you have.

We use the following Python packages, with the recommended versions:

	scipy >= 1.6.2
	argparse >= 1.1
	numpy >= 1.20.1
	subprocess
	os
	sys
	re
	random
	itertools
	pickle
	ftplib
	datetime



Clone the Codetta repository from GitHub with the command 

	git clone https://github.com/kshulgina/codetta

Codetta additionally requires:

- HMMER (v3.1) and Easel library: follow "Quick-ish" installation instructions in [HMMER user's guide](http://eddylab.org/software/hmmer/Userguide.pdf). 
- `gzip`: on Mac, use install command `brew install gzip`

If you plan on analyzing your own nucleotide sequences, then you will also need:

- `gtar`: on Mac, use install command `brew install gnutar`
- `wget`

Before we're done, we also need to download and build a local version of the Pfam database. We used Pfam version 32.0.

Download Pfam database into the `resources` directory. This may take a few minutes because this a ~19 Mb file.

	cd codetta/resources
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.seed.gz
	gunzip Pfam-A.seed.gz

Then, use HMMER to build a searchable database, using the `--enone` flag to turn off entropy weighting. This will also take a few minutes.

	hmmbuild --enone Pfam-A_enone.hmm Pfam-A.seed
	hmmpress Pfam-A_enone.hmm

Now you're ready to predict some genetic codes!

## Usage

Codetta consists of four separate programs. If you plan on inferring the genetic code from an alignment summary file, then you only need to use:

- `codetta_predict`: Infer the genetic code from the alignment summary file.

If you plan on analyzing your own nucleotide sequence, you will also use:

- `codetta_summary`: Summarize Pfam alignments into an alignment summary file.
- `codetta_align`: Align Pfam profile HMMs to the input nucleotide sequence.

General usage for these programs is

		python [program name].py [input file prefix] [optional arguments]

For any of these programs, type `python [program name] --help` for complete usage details.

## Tutorial

Let's see a few examples of how Codetta is used.

### Genetic code inference from an alignment summary file

In the `examples` directory, you will find `xxx`, an alignment summary file for species name (GenBank assembly GCA_xxx). This file summarizes the result of aligning the entire Pfam database against a six-frame translation of the entire genome.

We can infer the genetic code of xx using default parameters using

	python codetta_predict.py prefix

Notice that the input argument is the prefix of the alignment summary file, without a file extension.

The output is a string representation of the genetic code

`FFLLSSSS....`

Additionally, a file is created `examples/xxx.inference_xxx.txt`, which contains a detailed summary of the genetic code inference results.

[show preview of that file]

If you would like to change the default parameters of the analysis, you can do so by specifying the arguments `--evalue` to change the profile HMM hit e-value threshold, `--probability_threshold` to change the probability threshold used to call amino acid inferences, `--max_fraction` to change the maximum fraction a single consensus column can contribute to a single codon, and `-m -t -v -u -y` to change which groups of problematic Pfam domains are excluded.

If you plan on running a large number of analyses, you can use `--summary` specify a file to which a one-line summary of the results will be appended to.


### Genetic code inference from a nucleotide sequence

Predicting the genetic code from a nucleotide sequence requires two additional steps.

In the `examples` directory, there is a FASTA file called xxx with the _Mycoplasma pneumonia_ genome sequence. The input nucleotide sequence must a valid FASTA file as a DNA sequence (T instead of U).

The first step is to create a six-frame standard genetic code translation of the genome and align it to the entire Pfam database. We can do this with

	python codetta_align.py prefix

Notice that the input argument is the prefix of the fasta file, without the '.fna'.

This python program create several files, which are used by subsequent programs:

- `xxx.sequence_pieces.fna` and `.ssi` index: input nucleotide sequence, broken into pieces <100,000 nt
- `xxx.preliminary_translation.faa` and `.ssi` index: six-frame standard genetic code translation
- `xxx.temp_files/` is a temporary directory which stores the `hmmscan` result files

This step may take a while, depending on the size of the input nucleotide sequence. Rough estimate of about an hour on a single CPU core to analyze a ~12 Maa six-frame translation of a typical 6 Mb bacterial genome.   However, _Mycoplasma pneumonia_ has a small x Mb genome, so this will take only a few minutes.

If you intend to analyze a large number of genomes, we recommend parallelizing the computationally-intensize `hmmscan` step of the analysis over many machines on a computing cluster. To do this, you will need to modify the code in the `Codetta.py` file. Instructions can be found in the comment at the bottom of the `hmmscan_jobs()` function.

Next step is to process the `hmmscan` alignment files into an alignment summary file. This can be simply done with

	python codetta_summary.py prefix
	
And that's it! An alignment summary file is created (called `xxx.alignment_summary.txt.gz`) that can be used as input for the `codetta_predict` program as described above. The `xxx.temp_files/` has also been cleaned and deleted.


### Bonus: downloading nucleotide sequences from GenBank

We have provided a simple program for downloading FASTA files from GenBank for either assembly accessions or single nucleotide database accessions. 

- `codetta_download`: Download a genome assembly or nucleotide sequence from GenBank

Let's use this to download the mitochondrial genome of xxx

	python codetta_download.py xxxx c --prefix examples/xxxx

This will download a FASTA file containing the xxx sequence into `examples/xxxx`. The argument `c` specifies that this is a nucleotide database accession and not an assembly accession (which would be `a`).

### Summary

Now let's pull it all together and analyze the genetic code of this mitochondrial genome from scratch:

	python codetta_download.py xxxx c --prefix examples/xxxx
	python codetta_align.py prefix
	python codetta_summary.py prefix
	python codetta_predict prefix -m

The `-m` argument indicates that we do not want to exclude Pfam domains associated with mitochondrial genomes. The output genetic code is:


Which is consistent with the known alternative genetic code used in xx mitochondria.

 









