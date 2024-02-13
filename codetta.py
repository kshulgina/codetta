#!/usr/bin/env python

import scipy.special
import argparse
import os
from subprocess import call, run, Popen, PIPE, CalledProcessError
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import sys
import random
import itertools
import pickle

# load functions and datasets for specific genetic code inference tasks
from helper_functions import *

# silence numpy overflow errors (in np.exp, especially)
np.seterr(over='ignore')


def argument_parsing():
    parser = argparse.ArgumentParser(description="infer genetic code used by an organism from nucleotide sequence")
    parser.add_argument(
        'sequence_file',
        help='specify the input nucleotide sequence file in FASTA format.')
    # remaining arguments all are set optionally, otherwise default values
    parser.add_argument(
        '-p', '--profiles',
        help='profile HMM database file, must be in located in resource directory (default: Pfam-A_enone.hmm)')
    parser.add_argument(
        '-s', '--results_summary', type=str, default=None,
        help='file path to append one-line result summary ')
    #parser.add_argument('-i', '--identifier', help='GenBank genome assembly accession or GenBank nucleotide accession', type=str)
    #parser.add_argument('-d', '--download_type', help='specify whether download is for GenBank genome assembly accession (a) or GenBank nucleotide accession (c)', type=str, choices=['a', 'c'])
    parser.add_argument(
        '-e', '--evalue', type=float, default=1e-10, 
        help='profile HMM hit e-value threshold (default: 1e-10)')
    parser.add_argument(
        '-r', '--probability_threshold', type=float, default=0.9999,
        help='threshold for decoding probabilities (default: 0.9999)')
    parser.add_argument(
        '-f', '--max_fraction', type=float, default=0.01,
        help='maximum fraction of observations for a codon coming from a single Pfam position (default: 0.01)')
    parser.add_argument(
        '-m', '--mito_pfams', action="store_true", default=False,
        help='flag to include Pfam domains commonly found in mitochondria')
    parser.add_argument(
        '-t', '--transposon_pfams', action="store_true", default=False,
        help='flag to include Pfam domains associated with transposons and other mobile genetic elements')
    parser.add_argument(
        '-v', '--viral_pfams', action="store_true", default=False,
        help='flag to include Pfam domains associated with viruses')
    parser.add_argument(
        '-u', '--selenocysteine_pfams', action="store_true", default=False,
        help='flag to include Pfam domains known to contain selenocysteine')
    parser.add_argument(
        '-y', '--pyrrolysine_pfams', action="store_true", default=False,
        help='flag to include Pfam domains known to contain pyrrolysine')
    parser.add_argument(
        '--align_output',
        help='prefix for files created by codetta_align and codetta_summary. This can include a path. (default: [SEQUENCE_FILE])')
    parser.add_argument(
        '--inference_output',
        help='output file for codetta_infer step. Default is [ALIGN_OUTPUT].[PROFILES FILE].[inference parameters].genetic_code.out')
    parser.add_argument(
        '--resource_directory', type=str,
        help='directory where resource files can be found (default: [script dir]/resources)')
    parser.add_argument(
        '--parallelize_hmmscan', nargs="?", choices=['s', 'l'], default=None,
        help='parallelize hmmscan jobs [l]ocally or by sending to [s]LURM computing cluster.  Remember to modify the template file in resources directory accordingly. Default: None')
    parser.add_argument(
        '--njobs', default=2, type=int,
        help='number of hmmscan jobs to run in parallel, if parallelizing locally. Note that each hmmscan job uses 2 CPUs.')

    return parser.parse_args()


def initialize_globals():
    # standard genetic code dictionary used for translating nucleotide triplets to amino acids
    global gencode 
    gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        'ata':'I', 'atc':'I', 'att':'I', 'atg':'M', 'aca':'T', 'acc':'T', 'acg':'T', 'act':'T',
        'aac':'N', 'aat':'N', 'aaa':'K', 'aag':'K', 'agc':'S', 'agt':'S', 'aga':'R', 'agg':'R',
        'cta':'L', 'ctc':'L', 'ctg':'L', 'ctt':'L', 'cca':'P', 'ccc':'P', 'ccg':'P', 'cct':'P',
        'cac':'H', 'cat':'H', 'caa':'Q', 'cag':'Q', 'cga':'R', 'cgc':'R', 'cgg':'R', 'cgt':'R',
        'gta':'V', 'gtc':'V', 'gtg':'V', 'gtt':'V', 'gca':'A', 'gcc':'A', 'gcg':'A', 'gct':'A',
        'gac':'D', 'gat':'D', 'gaa':'E', 'gag':'E', 'gga':'G', 'ggc':'G', 'ggg':'G', 'ggt':'G',
        'tca':'S', 'tcc':'S', 'tcg':'S', 'tct':'S', 'ttc':'F', 'ttt':'F', 'tta':'L', 'ttg':'L',
        'tac':'Y', 'tat':'Y', 'taa':'_', 'tag':'_', 'tgc':'C', 'tgt':'C', 'tga':'_', 'tgg':'W'}
    
    # numerical indices corresponding to amino acids in most applications, such as column indices in profile HMM models
    global aa_indices
    aa_indices = {
            -2:'?', -1:'*', 0:'A', 1:'C', 2:'D', 3:'E', 4:'F', 5:'G', 6:'H', 7:'I', 8:'K', 
            9:'L', 10:'M', 11:'N', 12:'P', 13:'Q', 14:'R', 15:'S', 16:'T', 17:'V', 18:'W', 19:'Y', 20:'?'}
    
    # dictionary where keys are codons and values are ints, corresponds to order of codons in standard 
    # codon permutation such as in NCBI, and also order of codons in decoding probability lists
    global codons
    codons = list(itertools.product('TCAG',repeat=3))
    #
    global codon_order
    codon_order = {''.join(codons[0]):0}
    for cod in range(1, 64):
        codon_order[''.join(codons[cod])] = cod
    
    global std_gen_code
    std_gen_code = ''.join([gencode[''.join(codon)] for codon in codons]).replace('_', '*')


def initialize_emissions_dict(resource_dir, profile_db):
    # dictionary where keys are profile HMM names and values are emission probabilities for every position
    # load dictionary if it exists
    global emissions
    hmm_dictionary_file = '%s/%s.emissions_dict.p' % (resource_dir, profile_db)

    # get time of file creation to update emissions file if profile has been changed
    try:
        dict_time = int(datetime.datetime.fromtimestamp(os.path.getmtime(hmm_dictionary_file)).strftime("%Y%m%d%H%M%S%f"))
        hmm_time = int(datetime.datetime.fromtimestamp(os.path.getmtime('%s/%s' % (resource_dir, profile_db))).strftime("%Y%m%d%H%M%S%f"))
    except FileNotFoundError:
        dict_time = 0
        hmm_time = 1

    # load emissions if pickled file exists already
    if os.path.isfile(hmm_dictionary_file) and dict_time >= hmm_time:
        with open(hmm_dictionary_file, 'rb') as fp:
            emissions = pickle.load(fp)
    # make sure profile HMM database file is in the expected location
    elif not os.path.isfile('%s/%s' % (resource_dir, profile_db)):
        sys.exit('ERROR: the profile HMM database does not exist')
    # make a look-up dictionary for all the emission probabilities for the profile HMMs
    else:
        print('Profile HMM database emissions dictionary has not been created yet; creating...')
        emissions = {'profile_hmm': np.zeros(shape = (100, 20))}
        
        # open file containing HMM information
        f = open('%s/%s' % (resource_dir, profile_db))
        line = f.readline()
        
        # reading line one at a time until the file is empty, store emission probabiltiies into the dictionary
        # order of the amino acids is: A C D E F G H I K L M N P Q R S T V W Y
        while line:
            info = line.split()
            if info[0] == 'NAME':
                name = info[1]
                while info[0]!='LENG':
                    line = f.readline()
                    info = line.split()
                length = int(info[1])
                em_matrix = np.zeros(shape = (length, 20))
                while info[0]!='HMM':
                    line = f.readline()
                    info = line.split()
                for k in range(0, 5):
                    line = f.readline()
                while line.split()[0] != "//":
                    site = int(line.split()[0])-1
                    em_matrix[site,:] = line.split()[1:21]
                    line = f.readline()
                    line = f.readline()
                    line = f.readline()
                emissions[name] = em_matrix
            line = f.readline()
        
        f.close()
        
        # save dictionary 
        with open(hmm_dictionary_file, 'wb') as fp:
            pickle.dump(emissions, fp, protocol=2)
    
    # if profile database is Pfam, then remove list of bad Pfam domains
    if profile_db == 'Pfam-A_enone.hmm':
        if not os.path.isfile("%s/bad_pfams.txt" % resource_dir):
            sys.exit('ERROR: bad Pfams file cannot be found in the resource directory')
        with open("%s/bad_pfams.txt" % resource_dir, "r") as rf:
            pfam_rem = rf.read().splitlines()
        for pfam in pfam_rem:
            try:
                del emissions[pfam]
            except KeyError:
                pass


def _exec_script(path):
    """Chmod and execute shell script, returning True if success"""
    try:
        run(['chmod', '777', path], check=True)
        run([path], check=True)
        return True
    except CalledProcessError as err:
        print('error executing hmmscan shell script: %s' % err)
        return False


def _exec_script_parallel(paths, njobs):
    """Execute list of shell scripts in parallel

    Parameters
    ----------
    paths : list
        Paths to shell scripts
    jobs : int
        Number of scripts to run in parallel
    """
    with ProcessPoolExecutor(max_workers=njobs) as executor:
        futures = {
            executor.submit(_exec_script, path) : path
            for path in paths
        }
        for future in as_completed(futures):
            path = futures[future]
            if future.result():
                print('hmmscan shell script %s executed successfully' % path)
            else:
                print('hmmscan shell script %s failed' % path)


class GeneticCode:
    """
    Class to store parameters and inference of the genetic code of a set of 
    nucleotide sequences and functions to do the inference.
    """
    
    def __init__(self, args):
        """
        Initializes the object with parameters from the arguments given to the Python script
        """
        #
        if args.results_summary != None:
            summ_dir = os.path.dirname(args.results_summary)
            if summ_dir != '' and not os.path.isdir(summ_dir):
                sys.exit('ERROR: the path leading up to output summary file has not been created')
        self.profiles = args.profiles
        self.summary_file = args.results_summary
        self.resource_dir = args.resource_directory
        self.hmmer_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), 'hmmer-3.3.2/bin'))
        self.identifier = args.identifier
        self.parallelize_hmmscan = args.parallelize_hmmscan
        self.njobs = args.njobs
        if args.evalue != None and args.evalue < 0:
            sys.exit('ERROR: e-value threshold must be positive')
        else:
            self.e_value_threshold = args.evalue
        if args.probability_threshold != None and (args.probability_threshold < 0 or args.probability_threshold > 1):
            sys.exit('ERROR: probability threshold must be in the range [0, 1]')
        else:
            self.probability_threshold = args.probability_threshold
        if args.max_fraction != None and (args.max_fraction < 0 or args.max_fraction > 1):
            sys.exit('ERROR: max fraction parameter must be in the range [0, 1]')
        else:
            self.max_fraction = args.max_fraction
        if args.download_type != None:
            self.download = args.download_type
        
        # string designating which Pfam domain groups are excluded
        self.excluded_string = ''
        
        # if profile database isn't even Pfam, then don't exclude sets of Pfam domains
        if self.profiles != 'Pfam-A_enone.hmm':
            if not args.mito_pfams or not args.transposon_pfams or not args.viral_pfams or not args.selenocysteine_pfams or not args.pyrrolysine_pfams:
                print("Warning: Specified profile HMM database is not Pfam, so will not exclude Pfam domain sets from analysis")
        else:
            # if excluding mitochondrial Pfam domains, remove them from analysis
            if args.mito_pfams == False:
                mito_pfams_file = "%s/mito_pfams.txt" % self.resource_dir
                if not os.path.isfile(mito_pfams_file):
                    sys.exit('ERROR: mitochondrial Pfams file cannot be found in the resource directory')
                with open(mito_pfams_file, "r") as rf:
                    pfam_rem = rf.read().splitlines()
                for pfam in pfam_rem:
                    try:
                        del emissions[pfam]
                    except KeyError:
                        pass
                self.excluded_string += 'm'
            
            # if excluding transposon Pfam domains, remove them from analysis
            if args.transposon_pfams == False:
                transposon_pfams_file = "%s/transposon_pfams.txt" % self.resource_dir
                if not os.path.isfile(transposon_pfams_file):
                    sys.exit('ERROR: transposon Pfams file cannot be found in the resource directory')
                with open(transposon_pfams_file, "r") as rf:
                    pfam_rem = rf.read().splitlines()
                for pfam in pfam_rem:
                    try:
                        del emissions[pfam]
                    except KeyError:
                        pass
                self.excluded_string += 't'
            
            # if excluding viral Pfam domains, remove them from analysis
            if args.viral_pfams == False:
                viral_pfams_file = "%s/viral_pfams.txt" % self.resource_dir
                if not os.path.isfile(viral_pfams_file):
                    sys.exit('ERROR: viral Pfams file cannot be found in the resource directory')
                with open(viral_pfams_file, "r") as rf:
                    pfam_rem = rf.read().splitlines()
                for pfam in pfam_rem:
                    try:
                        del emissions[pfam]
                    except KeyError:
                        pass
                self.excluded_string += 'v'
            
            # if excluding selenocysteine-containing Pfam domains, remove them from analysis
            if args.selenocysteine_pfams == False:
                seleno_pfams_file = "%s/selenocysteine_pfams.txt" % self.resource_dir
                if not os.path.isfile(seleno_pfams_file):
                    sys.exit('ERROR: selenocysteine Pfams file cannot be found in the resource directory')
                with open(seleno_pfams_file, "r") as rf:
                    pfam_rem = rf.read().splitlines()
                for pfam in pfam_rem:
                    try:
                        del emissions[pfam]
                    except KeyError:
                        pass
                self.excluded_string += 'u'
            
            # if excluding pyrrolysine-containing Pfam domains, remove them from analysis
            if args.pyrrolysine_pfams == False:
                pyrro_pfams_file = "%s/pyrrolysine_pfams.txt" % self.resource_dir
                if not os.path.isfile(pyrro_pfams_file):
                    sys.exit('ERROR: pyrrolysine Pfams file cannot be found in the resource directory')
                with open(pyrro_pfams_file, "r") as rf:
                    pfam_rem = rf.read().splitlines()
                for pfam in pfam_rem:
                    try:
                        del emissions[pfam]
                    except KeyError:
                        pass
                self.excluded_string += 'y'
        
        if len(self.excluded_string) == 0:
            self.excluded_string = 'none'
        
        # initialize some other important values
        self.n_consensus = np.zeros(64,)
        self.likelihoods = np.zeros(shape=(21, 64), dtype=np.float32)
        self.decoding_probs = np.zeros_like(self.likelihoods)
        self.npieces = None
        
        # if we're downloading a sequence and target file name is specified, use just the identifier
        if args.sequence_file == None and args.identifier != None:
            args.sequence_file = args.identifier + '.fna'
            if validate_file_path(args.sequence_file)==False:
                sys.exit('ERROR: default sequence_file is not a valid file path!')

        # check that user-provided target sequence file name is valid
        if validate_file_path(args.sequence_file)==True:
            self.genome_path = args.sequence_file
        else:
            sys.exit('ERROR: [--sequence_file] is not a valid file path!')

        # set path prefix for alignment output files and check validity 
        if args.align_output == None:
            self.align_output = self.genome_path
        elif validate_file_path(args.align_output) == True:
            self.align_output = args.align_output
        else:
            sys.exit('ERROR: [--align_output] is not a valid file path!')
        
        # set path for inference output file and check validity 
        if args.inference_output == None:
            self.inference_file = "%s.%s.%s_%s_%s_excl-%s.genetic_code.out" % (self.align_output, self.profiles, str(self.e_value_threshold), 
                                        str(self.probability_threshold), str(self.max_fraction), self.excluded_string)
        else:
            self.inference_file = args.inference_output

        if validate_file_path(self.inference_file) == False:
            sys.exit('ERROR: [--inference_output] is not a valid file path!')

        self.scratch_dir = '%s.%s.temp_files' % (self.align_output, self.profiles)
        self.alignment_output = '%s.%s.alignment_output.txt' % (self.align_output, self.profiles)
    
    def get_genome(self):
        """
        Downloads the nucleotide sequence from NCBI through genome_download function
        """
        # download file
        file_path = genome_download(self.identifier, self.download, self.resource_dir, self.genome_path)
        
        # if genome is STILL not downloaded, then quit
        if file_path == 1 or not os.path.isfile(self.genome_path):
            sys.exit('This sequence could not be downloaded from Genbank')
        
        # check that sequence file satisfies FASTA format
        if not validate_fasta(self.genome_path):
            sys.exit('ERROR: Sequence file is not in fasta format')
        
        print('Genome was downloaded from FTP url %s\nTo file %s' % (str(file_path), self.genome_path))
    
    def processing_genome(self):
        """
        Reads in input nucleotide sequence line by line and simultaneously 1) breaks nucleotide 
        sequences in to pieces <100,000 nt, 2) creates a preliminary translation and esl-sfetch index 
        and 3) writes hmmscan scripts, which are launched at the end.
        """
        
        ## Checking that genome sequence is OK
        print('Reading in sequences from %s' % self.genome_path)
        
        # if genome file does not exist
        if not os.path.isfile(self.genome_path):
            sys.exit('ERROR: Input FASTA file does not exist.')
        
        # check that sequence file satisfies FASTA format
        if not validate_fasta(self.genome_path):
            sys.exit('ERROR: Input sequence file is not in FASTA format')
                
        ## Checking that profile HMM database has been pressed
        print('Profile HMM database used is %s' % self.profiles)
        
        if not os.path.exists('%s/%s.h3m' % (self.resource_dir, self.profiles)):
            sys.exit("ERROR: profile HMM database has not been pressed with hmmpress")
        
        ## Initializing sequence pieces and preliminary translation files, temp directory
        sequence_pieces_file = '%s.sequence_pieces.fna' % self.align_output
        preliminary_translation_file = '%s.preliminary_translation.faa' % self.align_output
                
        # initialize sequence pieces file
        print('Sequences (broken into pieces <100,000 nt) will be written to %s' % sequence_pieces_file)
        try:
            with open(sequence_pieces_file, 'w') as gpf:
                pass
        except FileNotFoundError:
            sys.exit('ERROR: could not open file path %s for writing' % sequence_pieces_file)
        
        # initialize the preliminary translation file
        print('Preliminary translation will be written to %s' % preliminary_translation_file)
        try:
            with open(preliminary_translation_file, 'w') as ptf:
                pass
        except FileNotFoundError:
            sys.exit('ERROR: could not open file path %s for writing' % preliminary_translation_file)
                
        # make a temporary files directory if it does not already exist
        print('All temp files will be written to %s/' % self.scratch_dir)
        if not os.path.exists(self.scratch_dir):
            os.makedirs(self.scratch_dir)

        
        ## Initializing variables to keep track for # of hmmscan commands per script
        len_analyzed = 0
        n_hmm = 0
        shell_count = 0
        indices_to_analyze = list()
        
        ## Main business: Read input sequence file in as FASTA, and break into genome pieces to process
        with open(self.genome_path) as f:
            n_piece = 0         # keeping track of number of resulting sequences (split if too long)
            
            # read in all lines corresponding to a single input sequence
            for header, seq_lines in itertools.groupby(f, lambda l: l.startswith(">")):
                if header == True:
                    continue
                fasta_string = ''.join([l.rstrip() for l in seq_lines])
                seqs = list()
                seqs.append(fasta_string)
                
                # if input sequence is longer than 100,000 nts, split into two even pieces
                while sum(np.array([len(a) for a in seqs]) > 100000) > 0:
                    long_inds = np.where(np.array([len(a) for a in seqs]) > 100000)[0]
                    for li in long_inds:
                        firstpart, secondpart = seqs[li][:len(seqs[li])//2], seqs[li][len(seqs[li])//2:]
                        seqs[li] = firstpart
                        seqs.append(secondpart)
                
                # go through resulting sequences: write each to sequence file, preliminary translation 
                # file and chuck together for hmmscan jobs
                for seq in seqs:
                    # write to sequence file
                    with open(sequence_pieces_file, 'a') as gpf:
                        gpf.write('>piece_%i\n' % n_piece)
                        gpf.write(seq + '\n')
                    
                    ## write preliminary 6-frame translation
                    dna = [seq, seq[1:], seq[2:], 
                           reverse_complement(seq), 
                           reverse_complement(seq[:-1]), 
                           reverse_complement(seq[:-2])]
                    
                    # translate each fragment, turn in-frame stop codons into X
                    prot = [replace_stop(translate(dnas, gencode)) for dnas in dna]
                    
                    # step through each of six frames and write to file
                    for i in range(6):
                        with open(preliminary_translation_file, 'a') as ptf:
                            ptf.write('>piece_%i_%i\n'  % (n_piece, i))
                            ptf.write(prot[i] + '\n')
                        
                        ## This part groups together some number of translated sequences for hmmscan script
                        # add this index to current hmmscan script
                        indices_to_analyze.append('%i_%i' % (n_piece, i))
                        
                        # increment sequence length and n sequences of current grouping
                        seq_frame_len = len(prot[i])
                        len_analyzed += seq_frame_len
                        n_hmm += 1

                        # if length or num sequences are big enough, write an hmmscan script
                        if len_analyzed > 2000000 or n_hmm > 2000:

                            # copy template batch script and add hmmscan commands
                            shell_script = '%s/hmmscan_%s.sh' % (self.scratch_dir, shell_count)
                            seq_names_file = '%s/seq_names_%s.txt' % (self.scratch_dir, shell_count)
                            with open(seq_names_file, 'w') as snf:
                                for indices_str in indices_to_analyze:
                                    snf.write('piece_%s\n' % indices_str)
                            with open(shell_script, 'w') as batch_file:
                                batch_file.write('#!/bin/bash\n\n')
                                batch_file.write('%s/esl-sfetch -f %s %s | %s/hmmscan --nobias --textw 100000 -o %s/hmm_output_%s %s/%s -\n' % 
                                        (self.hmmer_dir, preliminary_translation_file, seq_names_file, self.hmmer_dir, self.scratch_dir, shell_count, self.resource_dir, self.profiles))
                            
                            # reset for next hmmscan script
                            indices_to_analyze = list()
                            len_analyzed = 0
                            n_hmm = 0
                            shell_count += 1
                    
                    n_piece += 1
        
        self.npieces = n_piece
        
        # create esl-sfetch SSI index for preliminary translations (hmmscan uses this)
        with open(os.devnull, "w") as f:
            prelim_index_cmd = '%s/esl-sfetch --index %s' % (self.hmmer_dir, preliminary_translation_file)
            p = Popen(prelim_index_cmd, shell=True, stdout=f, stderr=f)
            p.wait()
        
        ## checking that preliminary translation length is as expected
        # get pieces file length
        pieces_len_cmd = 'grep -v ">" %s | tr -d "\n" | wc' % (sequence_pieces_file)
        p = Popen(pieces_len_cmd, shell=True, stdout=PIPE)
        output_p = p.communicate()[0].decode()
        pieces_len = int(output_p.split()[2])   
        
        # get genome length
        transl_len_cmd = 'grep -v ">" %s | tr -d "\n" | wc' % (preliminary_translation_file)
        p = Popen(transl_len_cmd, shell=True, stdout=PIPE)
        output_t = p.communicate()[0].decode()
        transl_len = int(output_t.split()[2])     
        
        expected_len = 2*pieces_len - 4*self.npieces
        if expected_len != transl_len:
            raise TypeError('Preliminary translation is not of expected length')

        # add remaining hmmscan commands to batch script
        if len(indices_to_analyze) > 0:
            shell_script = '%s/hmmscan_%s.sh' % (self.scratch_dir, shell_count)
            seq_names_file = '%s/seq_names_%s.txt' % (self.scratch_dir, shell_count)
            with open(seq_names_file, 'w') as snf:
                for indices_str in indices_to_analyze:
                    snf.write('piece_%s\n' % indices_str)
            with open(shell_script, 'w') as batch_file:
                batch_file.write('#!/bin/bash\n\n')
                batch_file.write('%s/esl-sfetch -f %s %s | %s/hmmscan --nobias --textw 100000 -o %s/hmm_output_%s %s/%s -\n' % 
                    (self.hmmer_dir, preliminary_translation_file, seq_names_file, self.hmmer_dir, self.scratch_dir, shell_count, self.resource_dir, self.profiles))
            shell_count += 1
        
        # Run hmmscan shell scripts serially
        if self.parallelize_hmmscan == None:
            for shell_i in range(shell_count):
                print('Running hmmscan shell script %i out of %i' % (shell_i + 1, shell_count))
                shell_script = '%s/hmmscan_%i.sh' % (self.scratch_dir, shell_i)
                dum = call(["chmod", "777", shell_script])
                dum = call([shell_script])
        # Run scripts locally in parallel
        elif self.parallelize_hmmscan == 'l':
            print('Running %i hmmscan shell scripts parallelized in batches of %i jobs' % (shell_count, self.njobs))
            shell_scripts = ['%s/hmmscan_%i.sh' % (self.scratch_dir, shell_i) for shell_i in range(shell_count)]
            _exec_script_parallel(shell_scripts, self.njobs)
        # Submit jobs to SLURM array
        elif self.parallelize_hmmscan == 's':
            # remember to change the partition name in the resources/template_job_array.sh file!
            print('Submitting a SLURM job array of %i hmmscan jobs' % shell_count)
            job_array_script = '%s/hmmscan_jobarray.sh' % self.scratch_dir
            dum = call(["cp", "%s/template_jobarray.sh" % self.resource_dir, job_array_script])
            with open(job_array_script, 'a') as jf:
                jf.write('#SBATCH --array=0-%i' % shell_count)
                jf.write('\n\nchmod 777 %s/hmmscan_${SLURM_ARRAY_TASK_ID}.sh' % self.scratch_dir)
                jf.write('\n%s/hmmscan_${SLURM_ARRAY_TASK_ID}.sh' % self.scratch_dir)
            with open(job_array_script) as f:
                p = Popen(['sbatch'], stdin=f, stdout=PIPE, stderr=PIPE)
                p.wait()
        # If LSF parallelization is turned on
        #elif self.parallelize_hmmscan == 'l':
        #    ## stuff
        #    dum = 1
    
    def write_outputs(self,  gen_code_preconv):
        """
        Writes long-form output about genetic code inference to an output file and one line output to summary file
        """
        with open(self.inference_file, 'a') as of:
            # write parameters
            of.write('# Analysis arguments\n')
            of.write('alignment_prefix   %s\n' % self.align_output)
            of.write('profile_database   %s\n' % self.profiles)
            of.write('results_summary    %s\n' % self.summary_file)
            of.write('evalue_threshold   %s\n' % str(self.e_value_threshold))
            of.write('prob_threshold     %s\n' % str(self.probability_threshold))
            of.write('max_fraction       %s\n' % str(self.max_fraction))
            of.write('excluded_pfams     %s\n' % self.excluded_string)
            
            # write N consensus cols and AA for each codon
            gencode_diff = ''.join(['?' if self.gen_code[c] == '?' else 'N' if self.gen_code[c] == std_gen_code[c] else 'Y' for c in range(64)])
            of.write('#\n# Codon inferences                      Consensus columns\n')
            of.write('# codon   inference  std code  diff?    N aligned  N used\n')
            for c in range(64):
                of.write('%-10s%-11s%-10s%-9s%-11i%-10i\n' % (''.join(codons[c]), self.gen_code[c], std_gen_code[c], gencode_diff[c], \
                                                              self.original_n_cols[c], self.n_consensus[c]))
            
            # write all log decoding probabilities
            of.write('#\n# Log decoding probabilities\n# codon      ')
            of.write(''.join(['%-13s']*21) % tuple(['logP(%s)' % aa_indices[a] for a in range(21)]) + '\n')
            for c in range(64):
                of.write('%-6s    %s\n' % (''.join(codons[c]), ' '.join(['%12.4f']*21) % tuple(self.decoding_probs[:,c])))
            
            # write genetic code string
            of.write('#\n# Final genetic code inference\n%s' % self.gen_code)
        
        # WRITING TO SUMMARY FILE
        summ_line = "%s,%s,%s,%s,%s,%s,%s,%s\n" % (self.align_output, self.profiles, str(self.e_value_threshold), str(self.probability_threshold), 
            str(self.max_fraction), self.excluded_string, self.gen_code, gen_code_preconv)
        if self.summary_file:
            if not os.path.isfile(self.summary_file):
                with open(self.summary_file, 'w') as sf:
                    sf.write('prefix,profile_db,evalue_threshold,prob_threshold,max_fraction,excluded_domains,inferred_gencode,inferred_gencode_best_models\n')
            with open(self.summary_file, 'a') as sf:
                sf.write(summ_line)
    
    def process_hmmscan_results(self):
        """
        Once hmmscan jobs are complete, this function will read them in, process them 
        to figure out which profile HMM columns have aligned to what codons in the sequence, 
        and then write these results to the ialignment output file. Positions are not
        filtered at this step except for poorly aligned positions and hits with 
        evalues > threshold are excluded.
        
        For each aligned profile HMM column, write the: profile HMM name, hit e-value, 
        position within profile HMM, genome piece, frame, position, codon at that position
        """
        sequence_pieces_file = '%s.sequence_pieces.fna' % self.align_output
        if not os.path.isfile(sequence_pieces_file):
            sys.exit('ERROR: sequence_pieces file (generated by codetta_align) cannot be found. Make sure you provide the correct file name prefix (do not include file extensions)')
        
        # if processing_genome wasnt run before this, then the number of sequence pieces wasn't set
        if self.npieces == None:
            if os.path.isfile(sequence_pieces_file):
                p = Popen('grep ">" %s | wc' % sequence_pieces_file, shell=True, stdout=PIPE)
                self.npieces = int(p.communicate()[0].decode().split()[0])
        
        # create list of all anticipated hmmscan output files
        file_suffixes = ['%i_%i' % (j, i) for j in range(self.npieces) for i in range(6)]
        all_possible_hmm_outs = ['piece_%s' % suff for suff in file_suffixes]
        
        print('Reading hmmscan result files from %s/' % self.scratch_dir)
        
        if not os.path.isdir(self.scratch_dir):
            sys.exit('ERROR: scratch directory of hmmscan results (generated by codetta_align) cannot be found. Make sure you provide the correct file name prefix (do not include file extensions) and correct profile HMM database file.')
        
        # get list of all hmmscan output files that exist
        p = Popen('grep "^Query:" %s/hmm_output_* | grep -o "piece_[0-9]*_[0-9]*"' % self.scratch_dir, shell=True, stdout=PIPE)
        created_files = p.communicate()[0].decode().split()
        
        # if not all files have been created, exit
        if len(set(created_files)) != len(set(all_possible_hmm_outs)):
            sys.exit('ERROR: Not all hmmscan jobs have completed')
        
        # concantenate all hmmscan outputs into a single file
        p = Popen('find %s -type f -name "hmm_output_*"' % self.scratch_dir, shell=True, stdout=PIPE)
        n_hmm_output_files = len(p.communicate()[0].decode().split())
        concat_hmm_output_file = '%s/concatentated_hmm_output' % self.scratch_dir
        with open(concat_hmm_output_file, 'w') as concat_file:
            for ind in range(n_hmm_output_files):
                hmm_output_file = '%s/hmm_output_%i' % (self.scratch_dir, ind)
                with open(hmm_output_file, 'r') as source_file:
                    dum = concat_file.write(source_file.read()) 
        
        print('Writing alignment output file to %s' % self.alignment_output)
        
        # initialize alignment output file
        alignment_output_file = self.alignment_output + '_unsorted'
        try:
            with open(alignment_output_file, 'w') as hf:
                pass
        except FileNotFoundError:
            sys.exit('ERROR: could not open file path %s for writing' % self.alignment_output)
        
        # opening file of all hmm_outputs for reading
        concat_hmm_outputs = open(concat_hmm_output_file, 'r')
        
        with open(sequence_pieces_file, 'r') as gpf:
            for x in range(self.npieces):
                
                # read in piece DNA sequence
                seq_piece_name = gpf.readline().rstrip()
                piece = gpf.readline().rstrip()
                # check that the correct file was read in
                if int(seq_piece_name.split('piece_')[1]) != x:
                    raise TypeError('Unexpected sequence piece is being compared')
                # get all 6 frames of genome piece
                dna = [piece, piece[1:], piece[2:], 
                       reverse_complement(piece), 
                       reverse_complement(piece[:-1]), 
                       reverse_complement(piece[:-2])]
                
                lines_to_write = list()
                
                # step through each of six frames
                for i in range(0,6):
                    
                    # read in hmmscan output lines that correspond to the next output file
                    hmm_file_lines = ''
                    line = concat_hmm_outputs.readline()
                    while line[:2] != '//':
                        hmm_file_lines += line
                        line = concat_hmm_outputs.readline()
                    hmm_file_lines = hmm_file_lines.split('\n')
                    
                    # validate that not corrupted
                    piece_name = 'piece_%i_%i' % (x, i)
                    if not validate_hmm_output(hmm_file_lines, piece_name):
                        raise TypeError('Hmmscan output file %i %i does not follow expected format' % (x, i))
                    
                    conserved_regions = extract_hmmscan_output(hmm_file_lines, self.e_value_threshold)
                    
                    # now I have a list of all the conserved regions in this translated fragment and how they map to profile HMMs 
                    # step through each amino acid site, determine the codon, and accumulate the emission probabilities
                    
                    # if no profile HMM hits found, skip to next genome piece
                    if len(conserved_regions) == 0:  
                        continue
                    
                    # iterate through each profile HMM hit located in genome piece
                    for c, con in enumerate(conserved_regions):
                        
                        # extracting profile HMM name and e-value for profile HMM hit
                        data = con.split(',')
                        profile_hmm_name = data[0]
                        profile_hmm_eval = float(data[1])
                        
                        ## Extract emissions data
                        hmm_inds = [int(da) for da in data[4::2]]
                        que_inds = [int(da) for da in data[5::2]]
                        
                        # iterate through codons in profile HMM hit and add them to position-emission dictionary
                        for ind, query_index in enumerate(que_inds):
                            # extract identity of codon
                            codon = dna[i][(query_index-1)*3:(query_index-1)*3+3].upper()
                            
                            # Record: codon at that position, profile HMM name, hit e-value, 
                            # position within profile HMM, genome piece, frame, position in query seq
                            try:
                                codon_index = codon_order[codon]  # this line fails if the codon contains a non-A/T/C/G char
                            except KeyError:
                                continue
                            cod_line = '%i,%i.%i%06d,%s,%i,%i,%i,%s,%s,%i\n' % (codon_index, x, i, query_index-1, codon, x, i, 
                                query_index-1, profile_hmm_name, str(profile_hmm_eval), hmm_inds[ind]-1)
                            lines_to_write.append(cod_line)
                
                if len(lines_to_write) > 0:
                    with open(alignment_output_file, 'a') as hf:
                        for line in lines_to_write:
                            dum = hf.write(line)
        
        concat_hmm_outputs.close()
        
        # sort results file using a Unix command
        with open(self.alignment_output, 'w') as f:
            p = Popen('sort -t, -n -k1,1 -k2,2 %s' % (alignment_output_file), shell=True, stdout=f)
            p.wait()
        
        print('Cleaning up temp files.')
        
        # clean up temporary files and scratch directory
        p = Popen('rm %s' % (alignment_output_file), shell=True)
        p.wait()
        
        p = Popen('find %s -type f -delete ' % self.scratch_dir, shell=True)
        p.wait()
        
        p = Popen('rmdir %s' % self.scratch_dir, shell=True)
        p.wait()
    
    def compute_decoding_probabilities(self):
        """
        This function will step through the aligned profile HMM columns (in the alignment output 
        file) and compute for each codon the likelihood and decoding probability of each model.
        """
        print('Reading alignment output file from %s' % self.alignment_output)

        # hmmscan alignment output file
        if not os.path.isfile(self.alignment_output):
            sys.exit('ERROR: alignment output file cannot be found. Make sure you provide the correct alignment prefix (do not include file extensions) and correct profile HMM database file.')
        
        # output file
        try:
            with open(self.inference_file, 'w') as of:
                pass
        except FileNotFoundError:
            sys.exit('ERROR: could not open file path %s for writing' % self.inference_file)

        # running sum of all emissions
        totsum = np.zeros(20) - np.inf    # -inf is log(0)
        
        sum_f = open(self.alignment_output, 'r')
        info = sum_f.readline().rstrip().split(',')
        self.original_n_cols = np.zeros(64)
        # iterate through each codon
        for cod in range(64):
            # get string for this codon
            codon = ''.join(codons[cod])
            
            ## Read in all results for a given codon at a time
            codon_lines = list()         # store information of all profile HMM columns aligned to this codon
            profile_hmm_pos_counts = dict()   # count how many times each profile HMM contributes a column
            
            if len(info) > 1:
                info_codon = info[2]
            else:
                info_codon = ''
            
            while info_codon == codon:  # compare to codon on the line currently being processed
                e_value = float(info[7])
                profile_hmm_name = info[6]
                
                # filter out hits with e-values above threshold
                if e_value > self.e_value_threshold:
                    next_line = sum_f.readline()
                    if len(next_line) == 0:
                        break
                    info = next_line.rstrip().split(',')
                    info_codon = info[2]
                    continue
                
                # filter out hits that belong to excluded Pfam domains (transposon, viral, etc)
                try:
                    emiss = np.float32(emissions[profile_hmm_name])
                except KeyError:
                    next_line = sum_f.readline()
                    if len(next_line) == 0:
                        break
                    info = next_line.rstrip().split(',')
                    info_codon = info[2]
                    continue
                
                # see if previous entry in list is at the same position in the query sequence! (overlapping hits)
                try:
                    last_info = codon_lines[-1]
                    last_position = last_info[3:6]
                    curr_position = info[3:6]
                    # if so, compare e-values and only keep more significant hit
                    if curr_position == last_position:
                        # compare e-values, if new < old, then remove the old one
                        if e_value < float(last_info[7]):
                            dum = codon_lines.pop()
                            old_profile_hmm_hit = last_info[6]
                            old_hmm_pos = last_info[8]
                            old_dict_key = '%s_%s' % (old_profile_hmm_hit, old_hmm_pos)
                            profile_hmm_pos_counts[old_dict_key] = profile_hmm_pos_counts[old_dict_key] - 1
                            if profile_hmm_pos_counts[old_dict_key] == 0:  # if only observation of codon was removed
                                del profile_hmm_pos_counts[old_dict_key]
                        else:
                            next_line = sum_f.readline()
                            if len(next_line) == 0:
                                break
                            info = next_line.rstrip().split(',')
                            info_codon = info[2]
                            continue
                except IndexError:
                    pass              # this happens if list is empty, just keep going
                
                # tally how many profile HMM columns for this codon have come from this profile HMM
                hmm_position = int(info[8])
                dict_key = '%s_%s' % (profile_hmm_name, hmm_position)
                try:
                    pf_counts = profile_hmm_pos_counts[dict_key]
                except KeyError:
                    pf_counts = 0.0
                profile_hmm_pos_counts[dict_key] = pf_counts + 1
                
                # add to set of columns aligning to this codon
                codon_lines.append(info)
                
                # read in the next line...
                next_line = sum_f.readline()
                if len(next_line) == 0:
                    break
                
                # check legitimacy of this line
                if not validate_codon_line(next_line):
                    raise TypeError('alignment output file has an incorrectly formatted line')
                
                info = next_line.rstrip().split(',')
                info_codon = info[2]
            
            ## calculate how much need to reduce subsample each profile HMM consensus column due to over-representation
            # in the case that very few total consensus columns aligned such that none is < max_fraction, give
            # all profile HMM positions just one observation
            initial_profile_hmm_counts = np.array(list(profile_hmm_pos_counts.values()))
            profile_hmm_pos_list = list(profile_hmm_pos_counts.keys())
            n_profile_hmm_pos = len(profile_hmm_pos_counts)
            pruned_profile_hmm_pos = list()
            
            # first if statement-- if there are fewer consensus columns than minimum to get ANY position below max fraction if all are set to 1, make all 1
            # second if statement is important-- deals with the case if no profile HMM columns observed, then codon needs to be treated at second case
            if n_profile_hmm_pos < 1.0 / self.max_fraction and n_profile_hmm_pos > 0:
                maximum_profile_hmm_counts = dict(zip(profile_hmm_pos_list, [1]*n_profile_hmm_pos))
                pruned_profile_hmm_pos = np.where(initial_profile_hmm_counts > 1)[0]
            else:
                # find profile_hmm columns where number of observed instances exceeds the max fraction
                exceed_profile_hmms = np.where(initial_profile_hmm_counts / np.sum(initial_profile_hmm_counts) > self.max_fraction)[0]
                # iteratively "remove" observations and recalculate whether anything exceeds max fraction
                while len(exceed_profile_hmms) > 0:
                    dum = [pruned_profile_hmm_pos.append(ed) for ed in exceed_profile_hmms]
                    max_count = np.floor(self.max_fraction * np.sum(initial_profile_hmm_counts))
                    initial_profile_hmm_counts[exceed_profile_hmms] = max_count
                    exceed_profile_hmms = np.where(initial_profile_hmm_counts / np.sum(initial_profile_hmm_counts) > self.max_fraction)[0]
                
                maximum_profile_hmm_counts = dict(zip(profile_hmm_pos_list, initial_profile_hmm_counts))
            
            # shuffle list of lines from codon, so when first n lines are chosen for a pruned profile HMM position, they are random
            random.shuffle(codon_lines)
            
            # go through final list of profile_hmm columns and compute model likelihoods
            profile_hmm_counts_pruned = dict(zip(profile_hmm_pos_list, np.zeros(n_profile_hmm_pos)))
            for profile_hmm_column in codon_lines:
                profile_hmm_name = profile_hmm_column[6]
                hmm_position = int(profile_hmm_column[8])
                dict_key = '%s_%s' % (profile_hmm_name, hmm_position)
                
                # only continue analyzing profile HMM position if not seen max number of times yet
                if profile_hmm_counts_pruned[dict_key] >= maximum_profile_hmm_counts[dict_key]:
                    continue
                
                emiss = np.float32(emissions[profile_hmm_name])
                hmm_probs = emiss[hmm_position]
                
                self.likelihoods[:-1,cod] += hmm_probs
                self.n_consensus[cod] += 1
                totsum = np.logaddexp(-hmm_probs, totsum)
                profile_hmm_counts_pruned[dict_key] += 1
            
            self.original_n_cols[cod] = len(codon_lines)
            #if len(pruned_profile_hmm_pos) > 0:
            #    self.n_subsampled[cod] = len(set(pruned_profile_hmm_pos))
        
        sum_f.close()
        
        # normalize likelihoods
        for c in range(64):
            m = self.n_consensus[c]
            if m == 0:
                self.likelihoods[:,c] = np.log(1.0/21.0)
            else:
                denominator = m * totsum
                self.likelihoods[:-1,c] = np.log(1.0/21.0) - self.likelihoods[:-1,c] - denominator  # for the 20 aa models
                self.likelihoods[-1,c] = np.log(1.0/21.0) - m * np.log(sum(self.n_consensus))       # for the nonspecific model
        
        # compute decoding probabilities and determine inferred genetic code string
        post_denoms = scipy.special.logsumexp(self.likelihoods, axis=0)
        self.decoding_probs = self.likelihoods - post_denoms
        exp_posts = np.exp(self.decoding_probs)
        
        maxes = np.amax(exp_posts, axis=0)
        gen_code = np.argmax(self.decoding_probs, axis=0)
        gen_code_preconv = np.array([i for i in gen_code])
        gen_code[np.where(maxes < self.probability_threshold)[0]] = -2
        self.gen_code = ''.join([aa_indices[g] for g in gen_code])
        
        # Make a string with best amino acid for uninferred codons shown in lowercase
        gen_code_preconv[np.where(self.n_consensus == 0)[0]] = -2
        gen_code_preconv = [aa_indices[g] for g in gen_code_preconv]
        for d in np.where(gen_code == -2)[0]:
            gen_code_preconv[d]=gen_code_preconv[d].lower()
        gen_code_preconv = ''.join(gen_code_preconv)
        
        # write final genetic code to file and print to stdout
        print('Writing detailed inference output to %s\n' % self.inference_file)
        self.write_outputs(gen_code_preconv)
        print('Genetic code: %s' % self.gen_code)


def main():
    args = argument_parsing()
    
    if args.resource_directory == None:
        args.resource_directory = os.path.join(os.path.dirname(__file__), 'resources')
    args.resource_directory = os.path.normpath(args.resource_directory)
    
    if args.profiles == None:
        args.profiles = 'Pfam-A_enone.hmm'
    
    args.identifier = None
    args.download_type = None

    # initialize genetic code with command line args and download genome
    initialize_globals()
    initialize_emissions_dict(args.resource_directory, args.profiles)
    gc = GeneticCode(args)
    
    # do codetta align
    print('\nSTEP 1-- codetta_align. Aligning profile HMM database to input nucleotide sequence')
    gc.processing_genome()
    
    # do codetta summary
    print('\nSTEP 2-- codetta_summary. Collating alignments into an alignment output file')
    gc.process_hmmscan_results()

    # do codetta infer
    print('\nSTEP 3-- codetta_infer. Inferring the genetic code')
    gc.compute_decoding_probabilities()


if __name__ == "__main__":
    sys.exit(main())
