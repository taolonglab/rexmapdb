#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Use the tab-delimited text file (Exported from Excel) pcr_primers_table.txt
to generate FASTA file with standard DNA code for each primer pair.

We will use this, to generate a copy number of unique amplicon variants for
each of these pair, to adjust 16S abundances from our pipeline.


Created on Fri Mar 17 13:30:43 2017

@author: igor
"""

import argparse
import itertools as it
import pandas as pd
import os
from Bio.Seq import Seq

#path = os.path.expanduser('~/cloud/research/microbiome/genomes/')
#nucl_codes_file = path+'data/nucleotide_codes.txt'
#primers_table = path+'data/pcr_primers_table.txt'
#primers_out_folder = path+'data/pcr_primers'


def parse_input():
    parser = argparse.ArgumentParser(
        description='Generate FASTAs of PCR primers.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-nucl-codes', 
                        default='primers/nucleotide_codes.txt',
                        help='Input file with DNA nucleotide extended codes.')
    parser.add_argument('--input-pcr-primers-table', 
                        default='primers/pcr_primers_table.txt',
                        help='Input file with a table of PCR primers..')
    parser.add_argument('--output-pcr-primers-folder', 
                        default='data/pcr_primers',
                        help='Output folder with PCR primers files.')
    args = parser.parse_args()
    return args


def expand_extended_code (seq):
    """ Expand extended nucleic acid code according to the codes dictionary. 
    Requires codes{} dictionary loaded. """
    expanded_seq = [codes[char] for char in seq]
    seq_combs = [''.join(exp_seq) for exp_seq in it.product(*expanded_seq)]
    return(seq_combs)



if __name__ == '__main__':
    
    # Parse command line arguments
    args = parse_input()
    
    # nucl_codes_file = sys.argv[1]
    # primers_table = sys.argv[2]
    # primers_out_folder = sys.argv[3]
    nucl_codes_file = args.input_nucl_codes
    primers_table = args.input_pcr_primers_table
    primers_out_folder = args.output_pcr_primers_folder
    
    # Load nucleotide codes, parse into dictionary
    codes = {}
    with open(nucl_codes_file) as nucl_file:
        for line in nucl_file:
            if line[0:2] == 'nt': # skip header line
                continue
            else:
                key, values = line.strip().split('\t')
                values = values.split(',')
                codes[key] = values
    
    
    # Load table with PCR primer sequences and annotation
    pcr_primer_df = pd.read_csv(primers_table, sep='\t')

    # For each primer pair, we first generate FASTA file with the standard
    # nucleic acid code. 
    
    
    # Expand the primer sequences into standed code for each primer pair
    pcr_primers_expanded = {}
    for id, row in pcr_primer_df.iterrows():
        p_fwd_exp = expand_extended_code(row['Primer1_sequence_5to3'])
        p_rev_exp = [str(Seq(s).reverse_complement()) for s in expand_extended_code(row['Primer2_sequence_3to5'])]
        pair_id = row['Hypervariable_region']+'_'+row['Primer1']+'-'+row['Primer2']
        pcr_primers_expanded[pair_id] = []
        # Add all expanded forward primers
        for i, s in enumerate(p_fwd_exp):
            pcr_primers_expanded[pair_id].append(('Forward_'+row['Primer1']+'_'+str(i).zfill(2), s))
        for i, s in enumerate(p_rev_exp):
            pcr_primers_expanded[pair_id].append(('Reverse_'+row['Primer2']+'_'+str(i).zfill(2), s))
    
    # Save these to FASTA files
    if not os.path.isdir(primers_out_folder):
        os.makedirs(primers_out_folder)
    
    for pair_id, primer_list in pcr_primers_expanded.items():
        pair_filename = primers_out_folder+'/'+pair_id+'.fasta'
        with open(pair_filename, 'w') as pair_f:
            for primer_id, primer_seq in primer_list:
                pair_f.write('>'+primer_id+'\n'+primer_seq+'\n')
            