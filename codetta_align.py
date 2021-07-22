#!/usr/bin/env python

from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="aligns profile HMM database to six-frame translation of input nucleotide sequence")
    parser.add_argument('prefix', help='specify the file prefix for the input nucleotide sequence in FASTA format (ie [PREFIX].fna). This can include \
                                          a path. Temporary files will be written to [PREFIX]_[PROFILE DB FILE].temp_files/ and output files will be written to [PREFIX]_[PROFILE DB FILE].[extensions]')
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('-p', '--profiles', help='profile HMM database file, must be in located in resource directory (default: Pfam-A_enone.hmm)')
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
    
    # initialize genetic code with command line args and download genome
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
    gc.create_preliminary_translation()
    gc.hmmscan_jobs()

if __name__ == "__main__":
    sys.exit(main())
