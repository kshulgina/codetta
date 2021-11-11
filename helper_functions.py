import os
import sys
from subprocess import call, Popen, PIPE
import re
import numpy as np
from ftplib import FTP
import datetime

def translate(sequence, gct):
    """
    For a DNA sequence, translate it into a sequence of one symbol per codon, following 
    a specific translation table. For example, this could be translating DNA sequence into
    an amino acid sequence following some genetic code, or it could be mapping the DNA 
    sequence to a sequence of symbols corresponding to each codon. Replaces values not in 
    the translation dict with X's
    
    Args:
        sequence (string): DNA sequence of A, T, G, C characters
        gc (dict): translation table
    Returns:
        string: new sequence of symbols, 1/3 length of DNA sequence
    
    """
    return ''.join([gct.get(sequence[3*i:3*i+3].upper(),'X') for i in range(len(sequence)//3)])

def replace_stop(sequence):
    """
    For a string, replace all '_' characters with 'X's
    
    Args:
        sequence (string)
    Returns:
        string: '_' characters replaced with 'X's
    
    """
    return sequence.replace('_', 'X')

global dna_complements
dna_complements = bytes.maketrans(b"AaTtGgCcRrYySsWwKkMmBbVvDdHhNn", b"TtAaCcGgYyRrSsWwMmKkVvBbHhDdNn")

def reverse_complement(sequence):
    """
    For a string consisting of valid IUPAC DNA characters, return the reverse complement
    """
    return sequence.translate(dna_complements)[::-1]

def validate_file_path(file_path):
    """
    Returns true if the string is a valid UNIX file path and false if not
    """
    unix_safe_name = re.sub(r'[^~/\\.\d\w-]', '_', file_path)
    if len(file_path) > 0  and unix_safe_name == file_path:
        return True
    else:
        return False

def validate_codon_line(line):
    '''
    Validates that line in hmmscan_summary file is valid
    '''
    # check line length (100 chosen somewhat arbitrarily)
    if len(line) > 100:
        return False
    
    # split into fields
    info = line.rstrip().split(',')
    
    try:
        dum = int(info[0])      # codon ind
        dum = float(info[1])    # coded position
        dum = int(info[3])
        dum = int(info[4])
        dum = int(info[5])
        dum = float(info[7])
        dum = int(info[8])
    except ValueError:
        return False
    
    if not info[2].isupper():
        return False
    
    return True

def validate_fasta(fasta_file_path):
    """
    Checks that FASTA sequence file, either downloaded genome or provided file, are in 
    correct FASTA format. Arg is file path. Returns True or False.
    """
    # get number of lines that start with >
    p = Popen('grep  "^>" %s |  wc' % (fasta_file_path), shell=True, stdout=PIPE)
    n_header = int(p.communicate()[0].decode().split()[0])
    
    # get number of lines that don't start with >
    p = Popen('grep -v "^>" %s |  wc' % (fasta_file_path), shell=True, stdout=PIPE)
    n_not_header = int(p.communicate()[0].decode().split()[0])
    
    if n_not_header < n_header:
        return False
    
    # check if any illegal characters present in non-header lines
    p = Popen('grep -v "^>" %s | grep -i [^acgturykmswbdhvn] | wc' % (fasta_file_path), shell=True, stdout=PIPE)
    n_illegal_chars = int(p.communicate()[0].decode().split()[0])
    if n_illegal_chars > 0:
        return False
    
    return True

def validate_hmm_output(hmm_output_files, hmm_file_name):
    """
    Checks that the content of an hmmscan output are not corrupted. Arg is a list of 
    line-by-line hmmscan outfile file contents. Returns True or False.
    """
    if len(hmm_output_files) < 20: # somewhat arbitrary
        return False
    # check that this is the output for the right file
    if not (hmm_file_name in hmm_output_files[7]):
        return False
    # first line always is a comment with name of program
    if hmm_output_files[0] != '# hmmscan :: search sequence(s) against a profile database':
        return False
    # last line always says // (not [ok] because I split on that)
    if hmm_output_files[-1] != '//':
        return False
    
    return True

def extract_hmmscan_output(lines, eval_threshold):
    """
    Reads a single hmmscan output and parses out relevant parts for genetic code inference.
    Filters domain hits based on E-value threshold.
    
    Args:
        lines (string): line-by-line contents of hmmscan output file
        eval_threshold (float): threshold for excluding domain hits based on e-value
    Returns:
        list of strings: a list of strings for each suitable domain found, with domain and 
                         conserved residue information encoded in a string
    """
    # parse output from hmmscan and identify which conserved domains were found
    y = 0
    conserved_regions = list()
    cr_append = conserved_regions.append  # define function locally to speed it up
    while y < len(lines):
        line = lines[y]
        domains = list()
        if line[0:2] == '>>':
            domain_name = line.split(' ')[1]
            y = y+3
            line = lines[y]
            if line[0:2] != '>>':
                while len(line)>0:
                    # only analyze domains below e-value threshold
                    splits = line.split()
                    dom_eval = splits[5]
                    if float(dom_eval) < eval_threshold:
                        domains.append(splits[0]+","+domain_name+","+dom_eval+","+splits[6]+","+splits[7]+","+splits[9]+","+splits[10])
                    y = y+1
                    line = lines[y]
                y = y+2
                line = lines[y]
            while line[0:2] != '>>':
                splits2=line.split()
                if len(splits2) > 0 and splits2[0] == '==':
                    for d, dom in enumerate(domains):
                        dom_info = dom.split(",")
                        if dom_info[0] == splits2[2]:
                            try:
                                hmm_ali = lines[y+1].split()[2]
                                que_ali = lines[y+3].split()[2]
                                post_ali = lines[y+4].split()[0]
                                y = y + 1
                                amt_add = 5
                            except IndexError:
                                try:
                                    hmm_ali = lines[y+2].split()[2]
                                    que_ali = lines[y+4].split()[2]
                                    post_ali = lines[y+5].split()[0]
                                    y = y + 2
                                    amt_add = 6
                                except IndexError:
                                    hmm_ali = lines[y+3].split()[2]
                                    que_ali = lines[y+5].split()[2]
                                    post_ali = lines[y+6].split()[0]
                                    y = y + 3
                                    amt_add = 7
                            line = lines[y]
                            while int(line.split()[3]) < int(dom_info[4]):
                                y = y + amt_add
                                line = lines[y]
                                hmm_ali = hmm_ali + line.split()[2]
                                que_ali = que_ali + lines[y+2].split()[2]
                                post_ali = post_ali + lines[y+3].split()[0]
                            # turn into numpy arrays so can mask
                            hmm_ali = np.frombuffer(hmm_ali.encode(), dtype='a1')
                            que_ali = np.frombuffer(que_ali.encode(), dtype='a1')
                            post_ali = np.frombuffer(post_ali.encode(), dtype='a1')
                            # now transform this into pairs of hmm and query index that have no gaps and pass 
                            # alignment posterior quality filter
                            mask1 = hmm_ali != b'.'
                            mask2 = que_ali != b'-'
                            mask3 = post_ali == b'*'
                            mask = mask1 & mask2 & mask3
                            
                            hmm_inds = np.ones(len(hmm_ali), dtype=int)
                            hmm_inds[~mask1] = 0
                            que_inds = np.ones(len(que_ali), dtype=int)
                            que_inds[~mask2] = 0
                            
                            hmms = np.cumsum(hmm_inds)-1+int(dom_info[3])
                            hmms = hmms[mask]
                            que = np.cumsum(que_inds)-1+int(dom_info[5])
                            que = que[mask]
                            
                            if len(hmms) == 0:
                                continue
                            region = np.empty((hmms.size + que.size,), dtype=int)
                            region[0::2] = hmms
                            region[1::2] = que
                            region = ','.join(region.astype(str))
                            cr_append("%s,%s,%s,%s,%s" % (dom_info[1], dom_info[2], dom_info[5], dom_info[6], region))
                
                y = y+1
                if y >= len(lines):
                    break
                line = lines[y]
        else:
            y = y+1
    
    return conserved_regions

def genome_download(species_id, type_download, database_dir, download_output_file):
    """
    For a given species ID, function will attempt to download a genome
    
    Args:
        species_id : species ID that will be matched to genomes in GenBank
        type_download ('a' or 'c'): specifying is ID is GenBank assembly or accession number
    Returns:
        Nothing
    """
    # database files
    genbank_file = '%s/assembly_summary_genbank.txt' % database_dir
    
    # try to download from genbank 
    genbank_res = download_genbank(species_id, genbank_file, type_download, download_output_file)
    return genbank_res

def update_genbank(genbank_file):
    '''
    Updates the assembly_summary file from GenBank, used to download assembly accessions
    '''
    # get modification time of local copy of GenBank database file (in UTC)
    UTC_OFFSET = datetime.datetime.utcnow() - datetime.datetime.now()
    if os.path.isfile(genbank_file):
        local_time = int((datetime.datetime.fromtimestamp(os.path.getmtime(genbank_file)) + 
            UTC_OFFSET).strftime("%Y%m%d%H%M"))
    else:
        local_time = 0
    
    # get modification time of database file on GenBank ftp server (in UTC)
    ftp = FTP('ftp.ncbi.nlm.nih.gov') 
    ftp.login()
    ftp.cwd('./genomes/ASSEMBLY_REPORTS')
    genbank_time = ftp.sendcmd("MDTM assembly_summary_genbank.txt")
    if genbank_time[:3] == "213":
        genbank_time = int(genbank_time[3:-2].strip())
    
    # if local version of file is younger than the GenBank version, download new version
    if local_time < genbank_time:
        print('Downloading updated version of GenBank assembly reference file')
        with open(os.devnull, "w") as f:
            dum = call('wget -O %s ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt' 
                % (genbank_file), shell=True, stdout=f, stderr=f)

def download_genbank(species_id, genbank_file, type_download, download_output_file):
    '''
    Downloads nucleotide sequence from GenBank, either assembly accession or nucleotide 
    accession. Returns either the ftp URL from which the file was downloaded OR a 1 if failed
    '''
    # deal with case where species id is a GenBank genome assembly ID
    if type_download == 'a':
        # see if newer versions of genome database file is available, if so download
        update_genbank(genbank_file)
        
        # get urls
        p = Popen("awk -F'\t' -v OFS='\t' '$1 == \"%s\" {print $20}' %s" 
            % (species_id, genbank_file), shell=True, stdout=PIPE)
        info = p.communicate()[0].decode().split('\n')[:-1]
        
        # if species id does not exist in GenBank
        if len(info) == 0:
            return 1
        
        # if multiple entries exist, just pick the first one
        u = info[0]
        prefix = u.split('/')[-1]
        download_url = '%s/%s_genomic.fna.gz' % (u, prefix)
        
        # check if sequence is being downloaded into a directory that exists
        download_dir = os.path.dirname(download_output_file)
        if download_dir != '' and not os.path.isdir(download_dir):
            sys.exit('ERROR: the path leading up to prefix has not been created')

        # download the genome! 
        genome_path = '%s.gz' % download_output_file
        with open(os.devnull, "w") as f:
            dum = call('wget -O %s %s' % (re.escape(genome_path), download_url), shell=True, stdout=f, stderr=f)
        
        # unzip
        os.system('gunzip -f %s' % (re.escape(genome_path)))
        
        return download_url
    
    # deal with case where species id is a GenBank accession number
    elif type_download == 'c':
        # check if sequence is being downloaded into a directory that exists
        download_dir = os.path.dirname(download_output_file)
        if download_dir != '' and not os.path.isdir(download_dir):
            sys.exit('ERROR: the path leading up to prefix has not been created')
        
        download_url = 'https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text' % species_id
        wget_command = "wget -q -O %s '%s'" % (download_output_file, download_url)
        with open(os.devnull, "w") as f:
            dum = call(wget_command, shell=True, stdout=f, stderr=f)
        return download_url
