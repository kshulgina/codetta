from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="process hmmscan outputs into a summary file that can be used for genetic code inference")
    parser.add_argument('prefix', help='specify prefix of output files from gc_alignment (ie [PREFIX].preliminary_translation.faa). This can include \
                                          a path. Hmmscan alignment summary file will be written to [PREFIX].hmmscan_summary.txt.gz')
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('-e', '--evalue', help='Pfam e-value threshold (default: 1e-10)', type=float, default=1e-10)
    parser.add_argument('-b', '--resource_directory', help='directory where resource files can be found (default: [script dir]/resources)', type=str)
    
    return parser.parse_args()

def main():
    
    args = argument_parsing()
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    args.resource_directory = os.path.normpath(args.resource_directory)
    
    # initialize genetic code with command line args and download genome
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
    initialize_globals(args.resource_directory)
    gc = GeneticCode(args)
    
    gc.process_hmmscan_results()

if __name__ == "__main__":
    sys.exit(main())
