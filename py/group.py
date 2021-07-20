#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Group unique sequences in a FASTA file. Takes output from
count.py

Skip any empty sequences and filter out sequences with too many
consecutive Ns.

Created on Mon Nov 20 16:08:06 2017

@author: igor
"""
import re
import argparse
import pandas as pd
from _include import log

script_title = 'Group unique sequences.'


#%% Table processor

def table_to_dict (table_filename):
    """ Load table with the information on v4 region copy number for each
    unique strain. """
    tab = pd.read_csv(table_filename, sep='\t', usecols=[1,2])
    return tab.set_index('variant_name').to_dict()['count']

#%% Species text 
def species_from_strain (strain_name, delim='_', out_delim='_', na_list=['sp.', 'bacterium']):
    """ Extract species from strain name. Use the first
    two fields separated by delim. """
    strain_data = strain_name.split('_')
    if strain_data[1] in na_list:
        return out_delim.join(strain_data[0:3])
    else:
        return out_delim.join(strain_data[0:2])


def dict_to_fasta (d, out_fasta, out_table):
    """ Write dictionary of sequences and meta-data to FASTA.
    Enumerate FASTA files with custom meta-data containing vid
    and concatenated list of unique species. Write another tab
    delimited file with full strain names and vid. """
    with open(out_fasta, 'w') as out_fa, open(out_table, 'w') as out_tab:
        out_tab.write('variant_id\tvariant_name\tcopy_number\tseq\n')
        for i, (seq, metas) in enumerate(d.items()):
            # Write meta-data first
            meta = '>' + str(i).zfill(5) + '-'
            for cp, strains in metas.items():
                unique_species = set()
                for s in strains:
                    if len(s.split('_')) > 1:
                        unique_species.add(species_from_strain(s))
                        # Write table
                        out_tab.write('\t'.join([str(i).zfill(5), s, str(cp), seq])+'\n')
                meta = meta + str(cp) + ':(' + ','.join(sorted(unique_species)) + ');'
            # Write FASTA
            out_fa.write(meta[:-1]+'\n')
            out_fa.write(seq+'\n')


#%% Main
def main(table_filename, fasta_filename, table_output, fasta_output):
    # fasta_filename = '/Users/igor/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_337F-805R_hang22_sequences.fasta'
    # table_filename = '/Users/igor/cloud/research/microbiome/genomes/data/vregions_db/16s_from_genomes_2017-07-20_V3-V4_337F-805R_hang22_counts.txt'
    # fasta_output = '/Users/igor/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_337F-805R_hang22_unique_sequences.fasta'
    # table_output = '/Users/igor/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_337F-805R_hang22_unique_sequences_table.txt'
    #%% Load table
    counts = table_to_dict(table_filename)
    #%% FASTA processor
    def fasta_to_seq_dict (fasta_filename, maxN=1):
        """ Load fasta_filename line by line into a dictionary. The
        sequences are going to be dictionary keys and meta-name 
        dictionary values as a list. """
        ren = re.compile('.*N{'+str(maxN)+',}.*')    
        d = {}
        
        def store_entry():
            if not ren.match(seq):
                if meta in counts:
                    cp = counts[meta]
                else:
                    cp = 1
                if seq in d:
                    # Retrieve copy number
                    if cp in d[seq]:
                        d[seq][cp].append(meta)
                    else:
                        d[seq][cp] = [meta]
                else:
                    d[seq] = {cp: [meta]}            
        
        with open(fasta_filename) as fa:
            seq = ''
            for i, line in enumerate(fa):
                if line[0] == '>':
                    if i > 0 and seq != '':
                        # Store previous entry
                        store_entry()
                    meta = line.strip()[1:]
                    seq = ''
                else:
                    seq = seq + line.strip()
            # Store last entry
            if seq != '':
                store_entry()
        return d
    #%% Load FASTA
    d = fasta_to_seq_dict(fasta_filename)
    #%% Write output
    dict_to_fasta(d, fasta_output, table_output)



def parse_input():
    parser = argparse.ArgumentParser(
        description=script_title,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-table', required=True, 
                        help='Input wrefseq_table file (output from add_refseq.py).')
    parser.add_argument('--input-fasta', required=True, 
                        help='Input wrefseq FASTA file (output from add_refseq.py).')
    parser.add_argument('--output-table', required=True, 
                        help='Output unique_variants table.')
    parser.add_argument('--output-fasta', required=True,
                        help='Output unique_variants FASTA.')
    args = parser.parse_args()
    return args
    
#%% Main call
if __name__ == '__main__':
    # table_filename = sys.argv[1]
    # fasta_filename = sys.argv[2]
    # table_output = sys.argv[3]
    # fasta_output = sys.argv[4]
    log(script_title)
    args = parse_input()
    table_filename = args.input_table
    fasta_filename = args.input_fasta
    table_output = args.output_table
    fasta_output = args.output_fasta    
    main(table_filename, fasta_filename, table_output, fasta_output)
    log('Done.')
