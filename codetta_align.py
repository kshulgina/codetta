from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="aligns profile HMM database to six-frame translation of input nucleotide sequence")
    parser.add_argument('prefix', help='specify the file prefix for the input nucleotide sequence in FASTA format (ie [PREFIX].fna). This can include \
                                          a path. Temporary files will be written to [PREFIX]_temp_files/ and output files will be written to [PREFIX].[extensions]')
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('-b', '--resource_directory', help='directory where resource files can be found (default: [script dir]/resources)', type=str)
    
    return parser.parse_args()

def main():
    
    args = argument_parsing()
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    
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
    initialize_globals(args.resource_directory)
    gc = GeneticCode(args)
    
    gc.processing_genome()
    gc.create_preliminary_translation()
    gc.hmmscan_jobs()

if __name__ == "__main__":
    sys.exit(main())
