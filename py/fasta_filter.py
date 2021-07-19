#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Filter FASTA sequences based on the sequence length. Supports multi-line
FASTA files.


Created on Wed Jun 23 15:06:44 2021

@author: igor
"""

import argparse
import re

def fasta_to_dict_of_lists (input_fa, keys=None, sep=' ', renum=True):
    """
    Parse single or multi line FASTA file into a dictionary. If keys are used,
    then meta-data line is split using 'sep', and elements are picked out
    using keys, e.g. keys=[1,2,3] picks "Meta_data_line" from
    >Meta data line something else.
    
    If renum==True, each sequence is assigned a unique meta-data key; if the
    same key is reused in the FASTA file, then suffix is added to ensure
    uniqueness.

    Parameters
    ----------
    input_fa : string
        FASTA file name.
    keys : list, optional
        List of integers to extract from the meta-data line. The default is None.
    sep : character, optional
        Separator in meta-data when keys is used. The default is ' '.
    renum: boolean
        Renumerate sequences with identical meta-data lines?

    Returns
    -------
    out : dict
        Dictionary where key is a meta-data line from the FASTA file, not
        including the prefix >. Value is a sequence corresponding to that
        meta-data line.

    """
    out = {}
    meta, seq = '', ''
    with open(input_fa) as in_fa:
        for i, line in enumerate(in_fa):
            if line[0] == '>':
                if i > 0:
                    # Store previous entry
                    if keys is None:
                        key = meta
                    else:
                        meta_l = meta.split(sep)
                        key = '_'.join([meta_l[i] for i in keys])
                    if key in out:
                        out[key].append(seq)
                    else:
                        out[key] = [seq]
                meta = line.strip()[1:]
                seq = ''
            else:
                seq = seq + line.strip()
        # Store last entry
        if keys is None:
            key = meta
        else:
            meta_l = meta.split(sep)
            key = '_'.join([meta_l[i] for i in keys])
            
        # Does key already exist in the out dict?
        if key in out:
            if renum:
                key_resub = re.sub('^.*_([0-9]+)$', '\\1', key)
                if key_resub == key:
                    # Regex sub failed. Nothing was extracted.
                    key_suffix = '_1'
                else:
                    key_suffix = '_'+str(int(key_resub) + 1)
                key = key+key_suffix
                out[key] = [seq]
            else:
                out[key].append(seq)
        else:
            out[key] = [seq]
    return out   


def fasta_dict_filter(d, min_len=0, max_len=5000):
    # Filter out sequences (values) from a dictionary d with length less than
    # min_len or more than max_len.
    d_out = {k:v for k,v in d.items() if len(v[0]) >= min_len and len(v[0]) <= max_len}
    return d_out

def dict_of_lists_to_fasta(d, fasta_file):
    """
    Iterate over dictionary d, where key is a meta-data line and
    value is a list of sequences with the same meta-data line.
    """
    with open(fasta_file, 'w') as fa_out:
        for k,v in d.items():
            fa_out.write('>'+k+'\n')
            fa_out.write(v[0]+'\n')
    

def parse_input():
    parser = argparse.ArgumentParser(
        description='FASTA Filter.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-fasta', required=True, 
                        help='Input FASTA file.')
    parser.add_argument('--output-fasta', required=True, 
                        help='Output FASTA file.')
    parser.add_argument('--min-length', required=False, default=1000,
                        help='Minimum sequence length.', type=int)
    parser.add_argument('--max-length', required=False, default=2000,
                        help='Maximum sequence length.', type=int)
    args = parser.parse_args()
    return args

if __name__ == '__main__':

    args = parse_input()

    # Load input FASTA
    fa_dict = fasta_to_dict_of_lists(args.input_fasta)
    
    # Filter
    fa_filt_dict = fasta_dict_filter(fa_dict, min_len=args.min_length,
                                     max_len=args.max_length)
    
    # Save output
    dict_of_lists_to_fasta(fa_filt_dict, args.output_fasta)