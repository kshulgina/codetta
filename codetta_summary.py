#!/usr/bin/env python

from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="process hmmscan outputs into an alignment output file that can be used for genetic code inference")
    
    parser.add_argument('align_output', help='specify prefix of files created by codetta_align. This can include a path. Hmmscan alignment output \
                                              file will be written to [ALIGN_OUTPUT].[PROFILES FILE].alignment_output.txt.gz')
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('-p', '--profiles', help='profile HMM database file, must be in located in resource directory (default: Pfam-A_enone.hmm)')
    parser.add_argument('-e', '--evalue', help='profile HMM hit e-value threshold (default: 1e-10)', type=float, default=1e-10)
    parser.add_argument('--resource_directory', help='directory where resource files can be found (default: [script dir]/resources)', type=str)
        
    return parser.parse_args()

def main():
    
    args = argument_parsing()
    
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    args.resource_directory = os.path.normpath(args.resource_directory)
    
    if args.profiles == None:
        args.profiles = 'Pfam-A_enone.hmm'
    
    # initialize genetic code with command line args and download genome
    args.sequence_file = args.align_output  # it's not used here so doesn't matter
    args.parallelize_hmmscan = None
    args.inference_output = None
    args.results_summary = None
    args.identifier = None
    args.download_type = None
    args.probability_threshold = None
    args.max_fraction = None
    args.mito_pfams = None
    args.transposon_pfams = None
    args.viral_pfams = None
    args.selenocysteine_pfams = None
    args.pyrrolysine_pfams = None
    initialize_globals()
    initialize_emissions_dict(args.resource_directory, args.profiles)
    gc = GeneticCode(args)
    
    gc.process_hmmscan_results()
    print('Done.')

if __name__ == "__main__":
    sys.exit(main())
