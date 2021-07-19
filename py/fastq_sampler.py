#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 09:13:02 2021

@author: igor
"""

import argparse
import sys, random

script_title = 'FASTQ sequence sampler.'

def parse_input():
    parser = argparse.ArgumentParser(
        description=script_title,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-fastq', '-i', required=True, 
                        help='Input FASTQ file.')
    parser.add_argument('--output-txt', '-o', required=False, default='[stdout]',
                        help='Output text file.')
    parser.add_argument('--nsample', '-n', required=False, default=10, 
                        help='Number of sequences to sample.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_input()
    # in_file = sys.argv[1]
    # nsamples = int(sys.argv[2])
    in_file = args.input_fastq
    out_file = args.output_txt
    nsamples = args.nsample
    seqs_all = []
    # Load all sequences
    with open(in_file, 'r') as in_f:
        for i, line in enumerate(in_f):
            if i % 4 == 1:
                seqs_all.append(line.strip())
    # Randomize
    seqs_rnd = random.sample(seqs_all, nsamples)
    # Print output
    if out_file == '[stdout]':
        for s in seqs_rnd:
            sys.stdout.write(s+'\n')
    else:
        with open(out_file, 'w') as out_f:
            for s in seqs_rnd:
                out_f.write(s+'\n')
    
        