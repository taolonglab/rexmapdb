#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: the two files that are output from count_16s_amplicon_variants.py
to which we will append the information from NCBI 16s refseq search FASTA file.


Created on Tue Feb 27 16:24:07 2018

@author: igor
"""

#%% FASTA processor
import argparse
import os, sys, re
# from subprocess import Popen, PIPE
# from io import StringIO
import pandas as pd
from count import blast_primers_vs_sequences, blast_out_vregion, get_blast_path

script_title = 'Combine hypervariable regions from full genomes and 16S sequences from NCBI RefSeq search.'


def fasta_to_dict (input_fa, keys=None, sep=' ', post_process=None,
                   post_process_seq=None, ignore_duplicates=True,
                   ignore_list=['spacer', 'unknown', 'Bacterium'
                                'Agent of', 'proteobacterium']):
    """ Load FASTA file into dictionary. If key is not None,
    then it should be a list of integer specifying which fields
    to use as a key. Fields are separated with sep. 
    If we need to post-process put function name in post_process.
    """
    def not_in_ignore_list(x):
        """ Returns number of matches of string x in the ignore list. """
        return sum(map(lambda l: len(re.findall(l, x)), ignore_list)) == 0
    seq = ''
        
    out = {}
    if post_process is not None:
        post_process_f = globals()[post_process]
    if post_process_seq is not None:
        post_process_s = globals()[post_process_seq]
        
    with open(input_fa) as in_fa:
        meta = ''
        seq = ''
        for i, line in enumerate(in_fa):
            if line[0] == '>':
                if i > 0: # Store last entry
                    # store_entry()
                    if keys is None:
                        key = meta
                    else:
                        meta_l = meta.split(sep)
                        key = '_'.join([meta_l[i] for i in keys])
                    if post_process is not None:
                        key = post_process_f(key)
                    if key in out:
                        key = key + '$'
                    if post_process_seq is not None:
                        seq = post_process_s(seq)
                    out[key] = [id, seq]
                meta = line.strip()[1:]
                try:
                    id = re.findall('^([^ ]+) .*', meta)[0]
                except IndexError:
                    id = 'NA'
                seq = ''
            else:
                seq = seq + line.strip()
        # Store final entry
        # store_entry()
        try:
            id = re.findall('^([^ ]+) .*', meta)[0]
        except IndexError:
            id = 'NA'
        if keys is None:
            key = meta
        else:
            meta_l = meta.split(sep)
            key = '_'.join([meta_l[i] for i in keys])
        if post_process is not None:
            key = post_process_f(key)
        if key in out:
            key = key + '$'
        if post_process_seq is not None:
            seq = post_process_s(seq)
        out[key] = [id, seq]
        
        if ignore_duplicates:
            out = {k:v for k,v in out.items() if k[-1] != '$'}
        if len(ignore_list) > 0:
            out = {k:v for k,v in out.items() if not_in_ignore_list(k)}
    return out

def strain_name_from_refseq_string (text):
    """ Return genus species strain name from RefSeq string of 16S sequence. """
    x = re.sub('^[\S]+ ', '', text) # This removes the accession code
    x = re.sub(r'16S (ribosomal RNA|rRNA)( \(.*\))?( gene)?( genes)?(, (partial|complete) sequence)?', '', x)
    x = re.sub(r'16S (ribosomal RNA|rRNA)( \(.*\))?( gene)?( genes)?(, (partial|complete) sequence)?', '', x)
    x = re.sub('strain ', '', re.sub('\'', '', re.sub('subsp. ', '', x)))
    x = re.sub('r.* operon, partial sequence; and', '', x)
    x = x.replace('16S ribosomal RNA rRNA', '')
    x = x.replace('[', '')
    x = x.replace(']', '')
    x = x.replace('*', '').replace('`', '')
    x = re.sub(r'ATCC[ |-]?([0-9]+)', r'ATCC \1', x)
    x = re.sub(r'isolation source.*', '', x)
    x = re.sub(r'\b(\w+)( \1\b)+', r'\1', x)
    x = re.sub(r'(nov )?(rrn )?(gene|DNA) for', '', x)
    x = re.sub(r'(complete|partial) cds', '', x)
    x = re.sub(r'culture( |-)collection.*', '', x)
    x = re.sub(r' gene$', '', x)
    x = x.replace('(', '').replace(')', '').replace(',', '').replace('strain: ', '')
    x = re.sub(r'isolate(:)?', '', x.replace('small subunit', ''))
    x = re.sub(r'(sequence|gene |partial| for )', '', x)
    x = x.replace(' and ', ' ')
    x = x.replace(' s ', ' ')
    x = x.replace(' complete ', ' ')
    x = x.replace(' forward ', ' ')
    x = x.replace(' end ', ' ')
    x = x.replace(' ribosomal ', ' ')
    x = x.replace(' RNA ', ' ')
    x = re.sub(r'[ ]{2,}', ' ', x) # Replace multiple spaces with just 1
    x = x.split(';')[0]
    return re.sub(r'\b(\w+)( \1\b)+', r'\1', x.strip()).replace(' ', '_')

def seq_to_basic_code (x):
    return re.sub('[^ACGT]', 'N', x)



def parse_input():
    parser = argparse.ArgumentParser(
        description=script_title,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-nongenome-16s-fasta', required=True, 
                        help='Input 16s_RefSeq FASTA file.')
    parser.add_argument('--input-genome-16s-counts', required=True, 
                        help='Input 16s_from_genomes counts table (output from counts.py).')
    parser.add_argument('--input-genome-16s-fasta', required=True, 
                        help='Input 16S sequences from genomes FASTA (output from counts.py).')
    parser.add_argument('--input-pcr-primers-fasta', required=True,
                        help='FASTA file with PCR primers used for full genomes in counts.py.')
    parser.add_argument('--output-table', required=True, 
                        help='Output table with combined counts.')
    parser.add_argument('--output-fasta', required=True, 
                        help='Output FASTA with combined sequences.')
    parser.add_argument('--overhang', required=False, default=0,
                        help='Extra overhang on both sides of extracted hypervariable regions. Use same value as in counts.py')
    parser.add_argument('--min-seq-len', required=False, default=200,
                        help='Minimum allowed hypervariable region sequence length. Use same value as in counts.py')
    parser.add_argument('--nthreads', required=False, default=4,
                        help='Number of parallel threads to use in the BLAST step.')
    parser.add_argument('--debug', required=False, default='False',
                        help='Debug mode (True/False).')
    args = parser.parse_args()
    return args


#%% Main function
if __name__ == '__main__':
    
    print(script_title)
    args = parse_input()
    
    ncbi_fasta = args.input_nongenome_16s_fasta
    variant_table = args.input_genome_16s_counts
    variant_fasta = args.input_genome_16s_fasta
    primer_file = args.input_pcr_primers_fasta
    out_table = args.output_table
    out_fasta = args.output_fasta
    overhang = int(args.overhang)
    seq_min_len = int(args.min_seq_len)
    nthreads = int(args.nthreads)
    debug = args.debug
    
    # ncbi_fasta = sys.argv[1]
    # variant_table = sys.argv[2]
    # variant_fasta = sys.argv[3]
    # primer_file = sys.argv[4]
    # out_table = sys.argv[5]
    # out_fasta = sys.argv[6]
    
    # # overhang = 22
    # overhang = int(sys.argv[7])
    # seq_min_len = int(sys.argv[8])
    
    blast_path = get_blast_path()
    
    print('* Load NCBI 16S search FASTA...', end='')
    ncbi16s_dict = fasta_to_dict(ncbi_fasta, post_process='strain_name_from_refseq_string', 
                                 post_process_seq='seq_to_basic_code')
    print('OK.')

    # Now load full genome table and check which strains are new
    print('* Load full genome variant table...', end='')
    vartab_df = pd.read_csv(variant_table, sep='\t')
    print('OK. (', str(len(vartab_df)), 'rows)')
    if 'strain_name' in vartab_df.columns:
        # For compatibility with older version of the script
        vartab_df.rename(columns={'strain_name':'variant_name'}, inplace=True)
    fullgen_strains = set([re.sub('_@rrn[0-9]+', '', s) for s in vartab_df['variant_name']])
    # For each strain from RefSeq search check first if the strain name exist
    # if it does, just skip it. If it does not exist, check if the sequence
    # after PCR primer trimming exist in the database.
    print('* Checking for existing strains...', end='')
    ncbi16s_new_dict = {k:v for k, v in ncbi16s_dict.items() if k not in fullgen_strains }
    # del ncbi16_dict
    ncbi_fasta_reduced = os.path.splitext(os.path.basename(ncbi_fasta))[0] + \
        '_reduced.fasta'
    with open(ncbi_fasta_reduced, 'w') as f_out:
        for meta, [id, seq] in ncbi16s_new_dict.items():
            f_out.write('>'+meta+'\n'+seq+'\n')
    print('OK.')
    
    # Now we obtained the reduced dictionary but with full length sequences. We
    # will need to do BLAST alignment to extract V3-V4 region, where it exists.
    print('* BLAST primers vs RefSeq sequences...', end='')
    blast_out = blast_primers_vs_sequences(primer_file, ncbi_fasta_reduced, blast_path, nthreads=nthreads)
    os.remove(ncbi_fasta_reduced) # Cleanup temp file
    print('OK.')

    print('* Selecting best primer hits...', end='')
    blast_out[13] = 'Forward'
    blast_out.loc[blast_out[1].str.startswith('Reverse'), 13] = 'Reverse'

    blast_out_best = blast_out.iloc[blast_out.groupby([0, 13]).apply(lambda t: t[11].idxmax())]
    blast_out_best2 = blast_out_best.groupby(0).filter(lambda t: len(t) > 1)

    blast_out_best2.columns = ['variant_name', 'primer', 'pct_sim', 'aln_len', 'mismatches',
                               'gapopen', 'seq_start', 'seq_end', 'pr_start',
                               'pr_end', 'eval', 'bitscore', 'score', 'pr_type']
    print('OK.')
    if debug == 'True':
        print(blast_out_best2.head())

    # Generate reference between strain names and NCBI RefSeq IDs
    print('* Generating strain names and RefSeq IDs...', end='')
    id_vs_strain_df = pd.DataFrame(list(ncbi16s_new_dict.keys()))
    id_vs_strain_df.columns = ['variant_name']
    id_vs_strain_df['id'] = [id for k, [id, s] in ncbi16s_new_dict.items()]
    id_vs_strain_df['seq'] = [s for k, [id, s] in ncbi16s_new_dict.items()]
    blast_out_best3 = pd.merge(blast_out_best2, id_vs_strain_df, on='variant_name')
    # if debug == 'True':
    #     print(blast_out_best3.head())
    
    # Dictionary mapping variant_name to the vregion sequence
    strain_to_vreg_dict = blast_out_best3.groupby('variant_name').apply(lambda x: blast_out_vregion(x, overhang)).to_dict()
    vartab_ncbi_df = pd.DataFrame([[id_vs_strain_df.loc[id_vs_strain_df['variant_name'] == strain, 
        'id'].iloc[0], strain+'_@rrn00', 1, vreg] for strain, vreg in strain_to_vreg_dict.items()])
    vartab_ncbi_df.columns = ['assembly_id', 'variant_name', 'count', 'sequence']
    vartab_ncbi_df_index_valid = (vartab_ncbi_df['sequence'].str.len() >= seq_min_len)
    vartab_ncbi_df = vartab_ncbi_df.loc[vartab_ncbi_df_index_valid]
    valid_strains = sum(vartab_ncbi_df_index_valid)
    total_strains = len(vartab_ncbi_df_index_valid)
    print('OK.')
    print('  Kept '+str(valid_strains)+' out of '+str(total_strains)+' total strains.')
        
    # Filter out useless strain entries ?
    print('* Filtering out unidentified strains...', end='')
    vartab_ncbi_filt_df = vartab_ncbi_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^Bacterium')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^[B|b]acteria')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^Unidentified')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('producing[^b]bacterium')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('degrading[^b]bacterium')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^[^_]+_[^_]+$')]
    # Strains with no strain designation (just genus species)
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^Endocytic')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^[A-Z]\\.')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('phytoplasma')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('mycoplasma')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^[A-Z0-9]+_')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^[^_]+-like')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^Deep-sea')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['variant_name'].str.contains('^Fe-oxidizing')]
    # Final table
    vartab_all_df = pd.concat([vartab_df, vartab_ncbi_df], ignore_index=True, sort=False)
    print('OK. (Combined: ', str(len(vartab_all_df)), 'rows)')
    if debug == 'True':
        print(vartab_all_df.head())

    # Save new table
    print('* Saving table...', end='')
    vartab_all_df.to_csv(out_table, sep='\t', index=False)
    print('OK. (', str(len(vartab_all_df)), 'rows)')
    # Save intermediate files in case we need them for debugging
    
    # Save sequences to new fasta
    print('* Saving FASTA...', end='')
    with open(out_fasta, 'w') as out_f:
        for i in range(0, len(vartab_all_df)):
            row_i = vartab_all_df.iloc[i]
            # If the sequence is too short, omit it.
            if len(row_i['sequence']) >= seq_min_len:
                out_f.write('>'+row_i['variant_name']+'\n'+row_i['sequence']+'\n')
    print('OK.')
    # For each new sequence check if its exact match to any reference, if yes
    # just add it under the same variant_id.
    print('Done.')