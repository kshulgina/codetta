from codetta import *

def argument_parsing():
    # initialize parser
    parser = argparse.ArgumentParser(description="download GenBank genome assembly or nucleotide sequences")
    parser.add_argument('identifier', help='GenBank genome assembly accession or GenBank nucleotide accession', type=str)
    parser.add_argument('download_type', help='specify whether download is for GenBank genome assembly accession (a) or GenBank nucleotide accession (c)', 
                        type=str, choices=['a', 'c'])
    
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument('--prefix', help='specify prefix of where to download the FASTA file (ie [PREFIX].fna). This can include \
                                          a path. (default: ./results/[IDENTIFIER].fna)')
    parser.add_argument('-b', '--resource_directory', help='directory where resource files can be found (default: [script dir]/resources)', type=str)
    
    return parser.parse_args()

def main():
    
    args = argument_parsing()
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    args.resource_directory = os.path.normpath(args.resource_directory)
    
    # initialize genetic code with command line args and download genome
    initialize_globals(args.resource_directory)
    args.results_summary = None
    args.evalue = None
    args.probability_threshold = None
    args.max_fraction = None
    args.mito_pfams = None
    args.transposon_pfams = None
    args.viral_pfams = None
    args.selenocysteine_pfams = None
    args.pyrrolysine_pfams = None
    gc = GeneticCode(args)
    
    # infer genetic code
    gc.get_genome()


if __name__ == "__main__":
    sys.exit(main())
