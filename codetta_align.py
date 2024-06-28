#!/usr/bin/env python3

from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="aligns profile HMM database to six-frame translation of input nucleotide sequence")
    parser.add_argument('sequence_file', help='specify the input nucleotide sequence file in FASTA format.')
    
    # remaining arguments all are set optionally, otherwise default values    
    parser.add_argument('--align_output', help='specify prefix of files created by codetta_align. This can include a path. Temporary files will be written \
                                                to [ALIGN_OUTPUT].[PROFILES FILE].temp_files/ and output files will be written to \
                                                [ALIGN_OUTPUT].[PROFILES FILE].[extensions]. (default: [SEQUENCE_FILE])')
    parser.add_argument('-p', '--profiles', help='profile HMM database file, must be in located in resource directory (default: Pfam-A_enone.hmm)')
    parser.add_argument('--parallelize_hmmscan', help='send hmmscan jobs to computing cluster. Specify SLURM (s), and modify template file in resources directory accordingly.', type=str, choices=['s'])
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
    args.inference_output = None
    args.results_summary = None
    args.identifier = None
    args.download_type = None
    args.evalue = None
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
    
    gc.processing_genome()
    print('Done.')

if __name__ == "__main__":
    sys.exit(main())
