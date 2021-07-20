#!/usr/bin/env python

from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="infer the genetic code of an organism from an intermediate file summarizing profile HMM alignments")
    parser.add_argument('prefix', help='specify prefix to hmmscan alignment summary input file (ie [PREFIX].hmmscan_summary.txt.gz). This can include \
                                          a path. Inference output will be written to [PREFIX].inference_[string of parameters].txt')
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('-p', '--profiles', help='profile HMM database file, must be in located in resource directory (default: Pfam-A_enone.hmm)')
    parser.add_argument('-e', '--evalue', help='Pfam e-value threshold (default: 1e-10)', type=float, default=1e-10)
    parser.add_argument('-r', '--probability_threshold', help='threshold for decoding probabilities (default: 0.9999)', type=float, default=0.9999)
    parser.add_argument('-f', '--max_fraction', help='maximum fraction of observations for a codon coming from a single Pfam position (default: 0.01)', type=float, default=0.01)
    parser.add_argument('-s', '--results_summary', help='file path to append one-line result summary', type=str, default=None)
    parser.add_argument('-m', '--mito_pfams', help='flag to include Pfam domains commonly found in mitochondria', action="store_true", default=False)
    parser.add_argument('-t', '--transposon_pfams', help='flag to include Pfam domains associated with transposons and other mobile genetic elements', action="store_true", default=False)
    parser.add_argument('-v', '--viral_pfams', help='flag to include Pfam domains associated with viruses', action="store_true", default=False)
    parser.add_argument('-u', '--selenocysteine_pfams', help='flag to include Pfam domains known to contain selenocysteine', action="store_true", default=False)
    parser.add_argument('-y', '--pyrrolysine_pfams', help='flag to include Pfam domains known to contain pyrrolysine', action="store_true", default=False)
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
    
    if args.profiles == None:
        args.profiles = 'Pfam-A_enone.hmm'

    args.download_type = None
    args.identifier = None
    
    # initialize genetic code with command line args and download genome
    initialize_globals(args.resource_directory, args.profiles)
    gc = GeneticCode(args)
    
    # infer genetic code
    gc.compute_decoding_probabilities()


if __name__ == "__main__":
    sys.exit(main())
