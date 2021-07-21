#!~/anaconda/bin/python3
# -*- coding: utf-8 -*-
"""
Use downloaded GFF files (~/data/ncbi_genomes_bacteria/features) and FASTA
sequences (~/data/ncbi_genomes_bacteria/sequences) to extract annotated 16S
sequences and store them in a single FASTA file.


Created on Fri Jul 21 08:15:42 2017

@author: igor
"""

import gzip, os, sys
import argparse
import pandas as pd
from Bio import SeqIO
from _include import log

script_title = 'Extract 16S sequences from full genomes.'

# ass_sub = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/bacteria_assembly_summary_sub.txt'
# ass_sub = '/Users/igor/cloud/research/microbiome/genomes/data/archaea/archaea_assembly_summary_sub_2017-08-24.txt'

def main(ass_sub, fa_out_name, ass_f_dir='features/', ass_s_dir='sequences/', 
         delim=';', match_str='16S ribosomal RNA.*', len_min=500, len_max=2000,
         fa_out_exist_action='append'):
    
    """
    Load accession summary table, and extract 16S sequences into a
    single-line FASTA. Arguments:
        - ass_sub: (filtered) assembly table, loaded as a pd.dataFrame
        - fa_out_name: name of the output FASTA file
        - ass_f_dir: directory with previously downloaded feature.gz files
        - ass_s_dir: directory with previously downloaded sequence.fa.gz files
        - delim: FASTA meta-data line delimiter to store extracted seq info
        - match_str: regex to match for name column in the GFF file, when looking for 16S seq
        - len_min, len_max: keep only 16S sequences within this range (inclusive both)
        - fa_out_exist_action: allowed values are either 'append' or 'overwrite', if overwrite
            it will overwrite the fa_out_name file if it already exists. if 'append' it will
            load the FASTA file, check what's already processed and skip those. this is used
            to resume the sdcript in case of early termination / errors.
        
    """
    
    #%% first load the main assembly summary file
    ass = pd.read_csv(ass_sub, sep='\t', usecols=['# assembly_accession', 'ftp_path'],
                      dtype={'# assembly_accession': str, 'ftp_path': str})
    
    #%% load gff annotation for each entry
    # seqs = defaultdict(list)
    log('Extract 16S sequences from full genomes')
    log('* matching feature names: ' + str(match_str))
    
    #
    # If FASTA file already exists, extract all already processed entries
    # then keep them in a list so that we won't have to reopen gz files again.
    #
    
    # Keep track of already processed assembly IDs in this list
    ass_exist = set([])
    last_assid = ''
    i0 = 0
    if os.path.isfile(fa_out_name) and fa_out_exist_action == 'append':
        log('* output '+fa_out_name+' already exists: ', end='')
        if fa_out_exist_action == 'append':
            log('appending to file...', end='')
            with open(fa_out_name, 'r') as f_out:
                exist_count = 0
                for record in SeqIO.parse(f_out, 'fasta'):
                    last_assid = record.id.strip().split(delim)[0]
                    ass_exist |= set([last_assid])
                    exist_count += 1
            try:
                i0 = ass.index[ass['# assembly_accession'] == last_assid].tolist()[0]                    
            except IndexError:
                i0 = 0
            log('OK. Added % 5d sequences.' % (exist_count))
            write_mode = 'a'
        elif fa_out_exist_action == 'overwrite':
            log('overwriting!')
            write_mode = 'w'
        else:
            log('wrong exist action!')
            sys.exit('Execution stopped due to wrong fa_out_exist_action.')        
    else:
        write_mode = 'w'
    
    log('* extracting sequences with lengths: '+str(len_min)+' <= len <= '+str(len_max))
    with open(fa_out_name, write_mode) as fa_out:
        log('* writing file ' + fa_out_name + '...')
        n_wrong_len = 0
        n_missing_assids = 0 # Keep track of missing records
        n_empty_assids = 0 # Keep track of the number of empty records (either GFF.GZ or FNA.GZ)
        n_existing_assids = 0 # Number of records already in the FASTA file
        
        for i in range(i0, len(ass)): 
            
            log('\r* % 5d out of % 5d assemblies (miss: % 4d, empty % 4d, % 5d seqs diff len, exist: % 4d).' % (i, len(ass),
                    n_missing_assids, n_empty_assids, n_wrong_len, n_existing_assids), end='')
            
            ass_id = ass.iloc[i]['# assembly_accession']

            # Check if this assembly ID is already in the FASTA file
            if ass_id in ass_exist:
                n_existing_assids += 1
                continue
            
            # ass_f_path = '../data/bacteria/features/' + \
            ass_f_path = os.path.join(ass_f_dir,
             os.path.basename(ass.iloc[i]['ftp_path']) + '_feature_table.txt.gz')
            
            # ass_s_path = '../data/bacteria/sequences/' + \
            ass_s_path = os.path.join(ass_s_dir,
             os.path.basename(ass.iloc[i]['ftp_path']) + '_genomic.fna.gz')
            
            # Check if any of the files is missing
            if not os.path.isfile(ass_f_path) or not os.path.isfile(ass_s_path):
                n_missing_assids += 1
                continue
            # Check if any of the files is empty 
            if os.path.getsize(ass_f_path) == 0 or os.path.getsize(ass_s_path) == 0:
                n_empty_assids += 1
                continue
            
            
            # Now load feature, extract coordinates for '16s ribosomal RNA' annotations
            try:
                ass_f = pd.read_csv(ass_f_path, sep='\t', compression='gzip',
                                    usecols=['assembly', 'genomic_accession', 'start', 'end', 'strand', 
                                             'name', 'attributes'],
                                    dtype={'assembly': str, 'genomic_accession':str, 'start': int, 'end': int, 
                                           'strand': str, 'name': str, 'attributes': str},
                                    keep_default_na=False)
            except EOFError:
                # Error reading GFF file, skip it.
                n_missing_assids += 1
                continue
    
            # ass_f_sub = ass_f.loc[(ass_f['name'].str.match(match_str)) & \
            #                      (ass_f['attributes'] != 'partial')]
            ass_f_sub = ass_f.loc[(ass_f['name'].str.match(match_str))]
            
            # Now if we have full length 16Ss, iterate over all sequence entries
            # in FASTA sequence file to find them.            
            if not ass_f_sub.empty:
                # We found some 16S hits that are not partial (i.e. are full)
                # Now read use these coordinates to extract sequences.
                try:
                    with gzip.open(ass_s_path, 'rt') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            # Check if record.id is in ass_f_sub['genomic_accession']
                            # If yes, then extract the sequence based on start and end.
                            ass_f_sub_sub = ass_f_sub[ass_f_sub['genomic_accession'] == record.id]
                            if not ass_f_sub_sub.empty:
                                for j in range(len(ass_f_sub_sub)):
                                    start = ass_f_sub_sub.iloc[j]['start']
                                    end = ass_f_sub_sub.iloc[j]['end'] + 1
                                    length = end - start
                                    # Check length
                                    if length < len_min or length > len_max:
                                        n_wrong_len = n_wrong_len + 1
                                    else:
                                        if ass_f_sub_sub.iloc[j]['strand'] == '-':
                                            seq = str(record.seq[start:end].reverse_complement())
                                        else:
                                            seq = str(record.seq[start:end])
                                        fa_out.write('>' + ass_id + delim + record.id + delim + 'loc:' + \
                                                     str(start) + ',' + str(end-1) + delim + 'strand:' + \
                                                     ass_f_sub_sub.iloc[j]['strand'] + delim + 'length:' + \
                                                     str(length) + '\n' + \
                                                     seq + '\n')
                except:
                    # Error reading sequence data, skip it.
                    n_missing_assids += 1
                    continue
                
        # Last progress bar after we're done.        
        log('\r* % 5d out of % 5d assemblies (miss: % 4d, empty % 4d, % 5d seqs diff len, exist: % 4d).' % (i, len(ass),
            n_missing_assids, n_empty_assids, n_wrong_len, n_existing_assids))

            

def parse_input():
    parser = argparse.ArgumentParser(
        description=script_title,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-asssum', required=True, 
                        help='Input assembly summary file (either bacteria or archaea).')
    parser.add_argument('--input-features-dir', required=True, 
                        help='Folder with assembly features (typically data/features/).')
    parser.add_argument('--input-sequences-dir', required=True, 
                        help='Folder with assembly sequences (typically data/sequences/).')
    parser.add_argument('--output-fasta', required=True, 
                        help='Output FASTA with 16S sequences.')
    parser.add_argument('--action', required=False, default='append',
                        help='Whether to "append" or "overwrite" output FASTA. Used for resuming from errors.')
    args = parser.parse_args()
    return args


#%% run this is script is run directly and not imported
if __name__ == '__main__':
#    ass_f_dir = '/Users/igor/data/ncbi_genomes/archaea/features/'
#    ass_s_dir = '/Users/igor/data/ncbi_genomes/archaea/sequences/'
#    out_fa = '/Users/igor/cloud/research/microbiome/genomes/data/archaea/16s_from_genomes_2017-08-24.fasta'
#    in_file = '/Users/igor/cloud/research/microbiome/genomes/data/archaea/archaea_assembly_summary_sub_2017-08-24.txt'
#    ass_f_dir = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/features/'
#    ass_s_dir = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/sequences/'
#    out_fa = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/16s_from_genomes_2017-07-20.fasta'
#    in_file = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria/bacteria_assembly_summary_sub_2017-07-20.txt'
    # in_file = sys.argv[1]
    # ass_f_dir = sys.argv[2]
    # ass_s_dir = sys.argv[3]
    # out_fa = sys.argv[4]
    # if len(sys.argv) == 6:
    #     action = sys.argv[5]
    # else:
    #     action = 'append'
    log(script_title)
    args = parse_input()
    
    in_file = args.input_asssum
    ass_f_dir = args.input_features_dir
    ass_s_dir = args.input_sequences_dir
    out_fa = args.output_fasta
    action = args.action
    
    ass_sub = in_file
    fa_out_name = out_fa
    main(in_file, out_fa, ass_f_dir=ass_f_dir, 
         ass_s_dir=ass_s_dir, fa_out_exist_action=action)
