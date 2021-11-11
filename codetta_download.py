#!/usr/bin/env python

from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="download GenBank genome assembly or nucleotide sequences")
    parser.add_argument('identifier', help='GenBank genome assembly accession or GenBank nucleotide accession', type=str)
    parser.add_argument('download_type', help='specify whether download is for GenBank genome assembly accession (a) or GenBank nucleotide accession (c)', 
                        type=str, choices=['a', 'c'])
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('--sequence_file', help='specify where to download the FASTA file. This can include a path. (default: [IDENTIFIER].fna)')
    parser.add_argument('--resource_directory', help='directory where resource files can be found (default: [script dir]/resources)', type=str)
    parser.add_argument('--hmmer_directory', help='directory where HMMER and Easel executables can be found (default: [script dir]/hmmer-3.1b2/bin)', type=str)

    return parser.parse_args()

def main():
    
    args = argument_parsing()
    
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    args.resource_directory = os.path.normpath(args.resource_directory)
    
    if args.hmmer_directory == None:
        args.hmmer_directory = os.path.join(os.path.dirname(__file__), 'hmmer-3.1b2/bin')
    args.hmmer_directory = os.path.normpath(args.hmmer_directory)
    
    # initialize genetic code with command line args and download genome
    args.align_output = None
    args.inference_output = None
    args.profiles = 'Pfam-A_enone.hmm'
    args.results_summary = None
    args.evalue = None
    args.probability_threshold = None
    args.max_fraction = None
    args.mito_pfams = None
    args.transposon_pfams = None
    args.viral_pfams = None
    args.selenocysteine_pfams = None
    args.pyrrolysine_pfams = None
    initialize_globals()
    gc = GeneticCode(args)
    
    # infer genetic code
    gc.get_genome()


if __name__ == "__main__":
    sys.exit(main())
