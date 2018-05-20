#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: the two files that are output from count_16s_amplicon_variants.py
to which we will append the information from NCBI 16s refseq search FASTA file.


Created on Tue Feb 27 16:24:07 2018

@author: igor
"""

#%% FASTA processor
import os, sys, re, pickle
# from subprocess import Popen, PIPE
# from io import StringIO
import pandas as pd
from count import blast_primers_vs_sequences, blast_out_vregion, get_blast_path

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


#%% Main function
if __name__ == '__main__':
    
    ncbi_fasta = '/Users/igor/cloud/research/microbiome/genomes/data/16s_rrna_ncbi_search_refseq_2018-03-19.fasta'
    variant_table = '/Users/igor/cloud/research/microbiome/genomes/data/vregions_db/16s_from_genomes_2018-02-23_V3-V4_341F-805R_hang22_counts.txt'
    variant_fasta = '/Users/igor/cloud/research/microbiome/genomes/data/vregions_db/V3-V4_341F-805R_hang22_sequences.fasta'
    primer_file = '/Users/igor/cloud/research/microbiome/genomes/data/pcr_primers/V3-V4_341F-805R.fasta'
    out_table = os.path.join(os.path.dirname(variant_table), 'V3-V4_341F-805R_hang22_wrefseq_table.txt')
    out_fasta = os.path.join(os.path.dirname(variant_fasta), 'V3-V4_341F-805R_hang22_wrefseq_sequences.fasta')

    ncbi_fasta = sys.argv[1]
    variant_table = sys.argv[2]
    variant_fasta = sys.argv[3]
    primer_file = sys.argv[4]
    out_table = sys.argv[5]
    out_fasta = sys.argv[6]
    
    # overhang = 22
    overhang = int(sys.argv[7])
    seq_min_len = 200
    
    #blast_path='/usr/local/ncbi/blast/bin/blastn'
    blast_path = get_blast_path()
    
    print('Add RefSeq sequences.')
    print('- load NCBI 16S search FASTA...', end='')
    ncbi16s_dict = fasta_to_dict(ncbi_fasta, post_process='strain_name_from_refseq_string', 
                                 post_process_seq='seq_to_basic_code')
    print('OK.')

    # Now load full genome table and check which strains are new
    print('- load full genome variant table...', end='')
    vartab_df = pd.read_csv(variant_table, sep='\t')
    print('OK.')
    fullgen_strains = set([re.sub('_@rrn[0-9]+', '', s) for s in vartab_df['strain_name']])

    # For each strain from RefSeq search check first if the strain name exist
    # if it does, just skip it. If it does not exist, check if the sequence
    # after PCR primer trimming exist in the database.
    ncbi16s_new_dict = {k:v for k, v in ncbi16s_dict.items() if k not in fullgen_strains }
    # del ncbi16_dict
    ncbi_fasta_reduced = os.path.splitext(os.path.basename(ncbi_fasta))[0] + \
        '_reduced.fasta'
    with open(ncbi_fasta_reduced, 'w') as f_out:
        for meta, [id, seq] in ncbi16s_new_dict.items():
            f_out.write('>'+meta+'\n'+seq+'\n')
            
    
    # Now we obtained the reduced dictionary but with full length sequences. We
    # will need to do BLAST alignment to extract V3-V4 region, where it exists.
    blast_out = blast_primers_vs_sequences(primer_file, ncbi_fasta_reduced, blast_path)
    os.remove(ncbi_fasta_reduced) # Cleanup temp file
    # Pickle these results to save progress
    with open('/Users/igor/cloud/research/microbiome/genomes/data/blast_out.pickle', 'wb') as b_out:
        pickle.dump(blast_out, b_out)

    blast_out[13] = 'Forward'
    blast_out.loc[blast_out[1].str.startswith('Reverse'), 13] = 'Reverse'

    blast_out_best = blast_out.iloc[blast_out.groupby([0, 13]).apply(lambda t: t[11].idxmax())]
    blast_out_best2 = blast_out_best.groupby(0).filter(lambda t: len(t) > 1)

    blast_out_best2.columns = ['strain_name', 'primer', 'pct_sim', 'aln_len', 'mismatches',
                               'gapopen', 'seq_start', 'seq_end', 'pr_start',
                               'pr_end', 'eval', 'bitscore', 'score', 'pr_type']

    # Generate reference between strain names and NCBI RefSeq IDs
    id_vs_strain_df = pd.DataFrame(list(ncbi16s_new_dict.keys()))
    id_vs_strain_df.columns = ['strain_name']
    id_vs_strain_df['id'] = [id for k, [id, s] in ncbi16s_new_dict.items()]
    id_vs_strain_df['seq'] = [s for k, [id, s] in ncbi16s_new_dict.items()]
    blast_out_best3 = pd.merge(blast_out_best2, id_vs_strain_df, on='strain_name')
    
    # Dictionary mapping strain name to the vregion sequence
    strain_to_vreg_dict = blast_out_best3.groupby('strain_name').apply(lambda x: blast_out_vregion(x, overhang)).to_dict()
    vartab_ncbi_df = pd.DataFrame([[id_vs_strain_df.loc[id_vs_strain_df['strain_name']==strain, 
        'id'].iloc[0], strain, 1, vreg] for strain, vreg in strain_to_vreg_dict.items()])
    vartab_ncbi_df.columns = ['assembly_id', 'strain_name', 'count', 'sequence']
    
    with open('/Users/igor/cloud/research/microbiome/genomes/data/vartab_ncbi_df.pickle', 'wb') as out:
        pickle.dump(vartab_ncbi_df, out)
    
    # Filter out useless strain entries
    vartab_ncbi_filt_df = vartab_ncbi_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^Bacterium')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^[B|b]acteria')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^Unidentified')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('producing[^b]bacterium')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('degrading[^b]bacterium')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^[^_]+_[^_]+$')]
    # Strains with no strain designation (just genus species)
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^Endocytic')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^[A-Z]\\.')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('phytoplasma')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('mycoplasma')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^[A-Z0-9]+_')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^[^_]+-like')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^Deep-sea')]
    vartab_ncbi_filt_df = vartab_ncbi_filt_df.loc[~vartab_ncbi_df['strain_name'].str.contains('^Fe-oxidizing')]
    # Final table
    vartab_all_df = pd.concat([vartab_df, vartab_ncbi_df], ignore_index=True)

    # Save new table
    vartab_all_df.to_csv(out_table, sep='\t', index=False)
    
    # Save intermediate files in case we need them for debugging
    
    # Save sequences to new fasta
    with open(out_fasta, 'w') as out_f:
        for i in range(0, len(vartab_all_df)):
            r = vartab_all_df.iloc[i]
            # If the sequence is too short, omit it.
            if len(r['sequence']) >= seq_min_len:
                out_f.write('>'+r['strain_name']+'\n'+r['sequence']+'\n')

    # For each new sequence check if its exact match to any reference, if yes
    # just add it under the same variant_id.