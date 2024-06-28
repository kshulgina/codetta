#!/usr/bin/env python3

from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="infer the genetic code of an organism from an intermediate file summarizing profile HMM alignments")
    parser.add_argument('align_output', help='specify prefix of files created by codetta_align and codetta_summary. This can include a path. Alignment output file \
                                              should be located at [ALIGN_OUTPUT].[PROFILES FILE].alignment_output.txt')
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('--inference_output', help='output file for codetta_infer step. Default is [align_output].[profiles].[inference parameters].genetic_code.out')
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
    parser.add_argument('--bad_profiles', help='list of profiles that should be excluded from the analysis', type=str)
        
    return parser.parse_args()

def main():
    
    args = argument_parsing()
    
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    args.resource_directory = os.path.normpath(args.resource_directory)
    
    if args.profiles == None:
        args.profiles = 'Pfam-A_enone.hmm'

    args.sequence_file = args.align_output  # it's not used here so doesn't matter
    args.download_type = None
    args.identifier = None
    args.parallelize_hmmscan = None
    
    # initialize genetic code with command line args and download genome
    initialize_globals()
    initialize_emissions_dict(args.resource_directory, args.profiles, args.bad_profiles)
    gc = GeneticCode(args)
    
    # infer genetic code
    gc.compute_decoding_probabilities()


if __name__ == "__main__":
    sys.exit(main())
