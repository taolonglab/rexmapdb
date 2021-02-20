#!~/anaconda/bin/python3
# -*- coding: utf-8 -*-
"""
For each strain from complete genome, count the number of identical hypervariable
regions that are experimentally amplified. 

First, load primers for each V-region, V1-V2, V3-V4, V4, V4-5. Then for each primer set
run BLAST and extract the sequence in the aligned region without a primer (make this an option).

The primer_dir is scanned for FASTA files which are loaded 1 by 1. The output
file names are automatically generated and save in the same primer_dir:
    
- V4_337F-805R_blast.txt: BLAST results (outfmt 6), no header
- V4_337F-805R_primer_miss.txt: assembly ids of 16s sequences that did not
    have good alignment with at least one of these primers.
- 16s_from_genomes_2017-07-20_V3,V4_337F-805R_counts.txt: tab-delimited table
    with columns 'assembly_id', 'strain', 'count' and 'sequence' which stores
    number of identical hypervariable regions of multiple 16S sequences for
    each strain (or assembly id)



Created on Thu Jul 27 13:27:07 2017

@author: igor
"""

import argparse
import sys, os, re, platform
import pandas as pd
from io import StringIO
from subprocess import Popen, PIPE
from threading import active_count

def get_os ():
    """ Detect operating system. """
    os_sys = platform.system()
    if os_sys == 'Linux':
        return 'linux'
    elif os_sys == 'Darwin':
        return 'macos'
    else:
        return ''
    
# Script folder sys.path[0]

def get_blast_path ():
    """ Generate an absolute BLAST path. """
    blastn = Popen('which blastn', shell=True, stdout=PIPE, stderr=PIPE)
    out, err = blastn.communicate()
    if out == b'':
        sys_os = get_os()
        return os.path.join(os.path.dirname(sys.path[0]), 'bin', 'blastn_'+sys_os)
    else:
        return out.decode().strip()


def fasta_to_df (input_fa):
    """ Load FASTA file into a DataFrame with column names 'meta' and 'seq'. """
    metas = []
    seqs = []
    meta, seq = '', ''
    with open(input_fa) as in_fa:
        for i, line in enumerate(in_fa):
            if line[0] == '>':
                if i > 0:
                    # Store previous entry
                    metas.append(meta)
                    seqs.append(seq)
                meta = line.strip()[1:]
                seq = ''
            else:
                seq = seq + line.strip()
        # Store last entry
        metas.append(meta)
        seqs.append(seq)
    return pd.DataFrame({'meta': metas, 'seq':seqs})

def fasta_to_dict_of_lists (input_fa, keys=None, sep=' '):
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
        if key in out:
            out[key].append(seq)
        else:
            out[key] = [seq]
    return out   

def assid_vs_spname_to_df (input_ass_sum):
    ass_ids = []
    ass_sps = []
    with open(input_ass_sum) as in_ass:
        for line in in_ass:
            if line[0] != '#':
                fields = line.strip().split('\t')
                # out[fields[0]] = fields[7]
                ass_ids.append(fields[0])
                ass_sps.append(fields[7])
    return pd.DataFrame({'ass_id': ass_ids, 'ass_sp': ass_sps})


#%% BLAST
    
def blast_primers_vs_sequences (primer_file, sequences_fa, blast_path='/usr/local/ncbi/blast/bin/blastn',
                                blast_ws=7, match=5, mismatch=-4, blast_go=8, blast_ge=6, nthreads=4,
                                blast_format='"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score"'):
    """ Blast primers vs blast db made out of full 16S sequences. Return the output
    as Blast outfmt 6 data frame. """
    
    # if get_os() == 'linux':
    #     try:
    #         nthreads = int(os.popen('grep -c cores /proc/cpuinfo').read())
    #     except:
    #         nthreads = 4
    # else:
    #     try:
    #         nthreads = int(active_count())
    #     except:
    #         nthreads = 4
        
    blast = Popen(' '.join([blast_path, '-subject', primer_file, '-query', sequences_fa, 
                   '-word_size', str(blast_ws), '-outfmt', blast_format, '-strand', 'both',
                   '-reward', str(match), '-penalty', str(mismatch),
                   '-gapopen', str(blast_go), '-gapextend', str(blast_ge),
                   '-num_threads', str(nthreads)
                   ]), 
                    stdout=PIPE, stderr=PIPE, shell=True)
    blast_out, blast_err = blast.communicate()
    blast_out_io = StringIO(blast_out.decode())
    blast_out_df = pd.read_csv(blast_out_io, sep='\t', header=None)
    return blast_out_df
    
def blastout_best (blast_out_df):
    """ Parse the blast output dataframe (outfmt 6) and keep only the best forward
    and best reverse primer for each sequence. """
    blast_grp = blast_out_df.groupby(0)
    blast_best_list = []
    missing_primers = []
    
    for seq_str, df in blast_grp:
        # Find the best Forward read by E-score. In case of ties, this will just select
        # the first one that comes in the list. This is fine for this purpose.
        fwd = df[1].str.startswith('Forward')
        rev = df[1].str.startswith('Reverse')
        if any(fwd) and any(rev) > 0: # Do we have any hit for each primer?      
            blast_best_list.append(df.ix[df[fwd][11].idxmax()])
            blast_best_list.append(df.ix[df[rev][11].idxmax()])
        else:
            missing_primers.append(seq_str)
    return (pd.concat(blast_best_list, ignore_index=True), missing_primers)


def blast_out_vregion (df, overhang=0):
    """ Extract the v-region from dataframe df with named columns. 
    If overhang > 0, add 'overhang' number of nts extra from each end. """
    start = df[df['pr_type']=='Forward']['seq_end'].iloc[0]
    end   = df[df['pr_type']=='Reverse']['seq_start'].iloc[0]
    seq   = df.iloc[0]['seq']
    return seq[(start+1-overhang):(end+overhang)]



#%% main call
def main(ass_fasta_a, ass_fasta_b, assembly_file_a, assembly_file_b, 
         primer_dir, fasta_out_dir, 
         blast_path='/usr/local/ncbi/blast/bin/blastn', overhang=21, 
         min_len=200, primer_file_filter='V', nthreads=4):
    
    print('Count 16S hypervariable regions')
    # Load full genome assembly 16S sequences
    
    # If we want to process only specific primer set add it here as a string
    # with exact match to the FASTA filename (output from primer_all_combinations.py)
    #  primer_file_filter = 'V3-V4_341F-805R'
    print('* Minimum sequence length:', str(min_len), 'nt')
    print('* Load files: ', end='')
    # Overhang from each V-region. Need this to make sure we align full query.
    overhang = int(overhang)
    
    primer_files = [os.path.join(primer_dir, p) for p in os.listdir(primer_dir) if 
                    os.path.splitext(p)[1] == '.fasta' and primer_file_filter in p]
    

    # Load the table between full accession id (from FASTA) and full 16S sequence.
    ass_to_seq_df = pd.concat([fasta_to_df(ass_fasta_a), fasta_to_df(ass_fasta_b)],
                               ignore_index=True)
    ass_to_seq_df['ass_id'] = [m.split(';')[0] for m in ass_to_seq_df['meta']]
    
    # Load the assembly file and process strain names
    ass_df_a = pd.read_csv(assembly_file_a, sep='\t', 
                         usecols=['# assembly_accession', 'strain_name'], 
                         dtype={'# assembly_accession': str, 'strain_name': str})

    ass_df_b = pd.read_csv(assembly_file_b, sep='\t', 
                     usecols=['# assembly_accession', 'strain_name'], 
                     dtype={'# assembly_accession': str, 'strain_name': str})

    ass_df = pd.concat([ass_df_a, ass_df_b], ignore_index=True)
    
    # Rename the assembly accession column
    ass_df = ass_df.rename(columns={'# assembly_accession': 'ass_id'})
    ass_df['strain_name'] = [s.replace(' ', '_') for s in ass_df['strain_name']]

    # After filter also filter out ass_to_seq_df by ass_id column
    ass_to_seq_df = ass_to_seq_df[ass_to_seq_df['ass_id'].isin(ass_df['ass_id'])]

    # Generate a dictionary between accession id and strain name for easy access later
    ass_to_strain_dict = ass_df[['ass_id', 'strain_name']].set_index('ass_id').to_dict()['strain_name']
    print('OK.')
    print('* BLAST path: '+blast_path)
    
    # Iterate over each primer pair combination
    for primer_file in primer_files:
        
        print('* Primer FASTA: '+primer_file)
        print('* Assembly FASTA (Archaea): '+ass_fasta_a)
        print('* Assembly FASTA (Bacteria): '+ass_fasta_b)
        print('\r* '+os.path.splitext(os.path.basename(primer_file))[0]+': blast...', end='')
        # Run BLAST with 16s sequences vs PCR primer sequences
        blast_out_a_df = blast_primers_vs_sequences(primer_file, ass_fasta_a, blast_path, nthreads=nthreads)
        blast_out_b_df = blast_primers_vs_sequences(primer_file, ass_fasta_b, blast_path, nthreads=nthreads)
        blast_out_df = pd.concat([blast_out_a_df, blast_out_b_df], ignore_index=True)
        blast_out_df[13] = 'Forward'
        blast_out_df.loc[blast_out_df[1].str.startswith('Reverse'), 13] = 'Reverse'
        
        print('OK best...', end='')
        # For each sequence (0) and primer type (12) keep only the highest hit by
        # bitscore (11)
        blast_out_best = blast_out_df.iloc[blast_out_df.groupby([0, 13]).apply(lambda t: t[11].idxmax())]
        # Now select only hits for which we have >1 (2) primer hit
        blast_out_best2 = blast_out_best.groupby(0).filter(lambda t: len(t) > 1)
        seq_miss_primer = list(set(blast_out_df[0]) - set(blast_out_best2[0]))
        
        # Extract unique strain assembly ids
        strain_miss_primer = list(set([m.group(1) for m in (re.search(r'^(GCF_[0-9]+\.[0-9]+).*', l) for l in seq_miss_primer) if m]))
        # Extract all unique strains
        # strains = list(set([m.group(1) for m in (re.search(r'^(GCF_[0-9]+\.[0-9]+).*', l) for l in blast_out_df[0]) if m]))
        
        # How many missing sequences per strain
        # ass_df.loc[ass_df['# assembly_accession'].isin(strain_miss_primer)]['organism_name'].value_counts()
        # Now generate results for each primer pair
        # --
        # First, we will need a BLAST result summary, that we will use in the next
        # function to extract actual sequences and then count them. Make a histogram
        # showing the number of mismatches (and %) for each primer to see how many
        # strains we will have difficulty amplifying with this primer set.
        blast_out_best2.columns = ['meta', 'primer', 'pct_sim', 'aln_len', 'mismatches',
                                   'gapopen', 'seq_start', 'seq_end', 'pr_start',
                                   'pr_end', 'eval', 'bitscore', 'score', 'pr_type']
        blast_out_best3 = pd.merge(blast_out_best2, ass_to_seq_df, on='meta')
        # blast_out_best3['v_seq'] = blast_out_best3.apply(
        #        lambda df: df['seq'][df['seq_start']:df['seq_end']+1], axis=1)
        print('OK v-regions...', end='')        
        
        # Extract v-regions (this takes a while) into a dictionary
        ass_to_vreg_dict = blast_out_best3.groupby('meta').apply(lambda x: blast_out_vregion(x, overhang)).to_dict()
        
        # Generate a dictionary {assid1: {seq1: count1, seq2: count2}, assid2...}
        assid_to_vreg_dict = {}
        for meta, seq in ass_to_vreg_dict.items():
            if seq == '' or len(seq) < min_len: # Empty sequence
                continue
            assid = meta.split(';')[0]
            if assid not in assid_to_vreg_dict:
                assid_to_vreg_dict[assid] = {}
            if seq not in assid_to_vreg_dict[assid]:
                assid_to_vreg_dict[assid][seq] = 1
            else:
                assid_to_vreg_dict[assid][seq] += 1
        print('OK write.', end='')

        # Now iterate over dictionary and save the counts to file
        counts_filename = os.path.join(fasta_out_dir,
             os.path.splitext(os.path.basename(ass_fasta_a.replace('archaea_', '')))[0] + '_' + \
             os.path.splitext(os.path.basename(primer_file))[0] + \
                 '_hang' + str(overhang) + '_counts.txt')
        fasta_filename = os.path.join(fasta_out_dir,
             os.path.splitext(os.path.basename(primer_file))[0].replace(',', '-') + \
                 '_hang' + str(overhang) + '_sequences.fasta')

        with open(counts_filename, 'w') as counts_file, open(fasta_filename, 'w') as fasta_file:
            counts_file.write('assembly_id\tstrain_name\tcount\tsequence\n')
            prev_strain = ''
            rrna_copy = 0
            for assid, counts_dict in assid_to_vreg_dict.items():
                for seq, count in counts_dict.items():
                    if ass_to_strain_dict[assid] == prev_strain:
                        # same strain, relabel rRNA copies
                        rrna_copy += 1
                    else:
                        rrna_copy = 0
                    strain_name = ass_to_strain_dict[assid] + '_@rrn' + str(rrna_copy).zfill(2) 
                    if seq != '':
                        counts_file.write('\t'.join([assid, strain_name, 
                                                     str(count), seq])+'\n')
                        fasta_file.write('>'+strain_name+'\n'+seq+'\n')
                    prev_strain = ass_to_strain_dict[assid]
        
        # Save BLAST results to file
        blast_out_best2.to_csv(primer_dir+'/'+os.path.splitext(os.path.basename(primer_file))[0]+'_blast.txt', 
                               sep='\t', index=False)
        
        print('.', end='')
        # Save 16S with missing 1+ primer alignment to file
        pd.DataFrame(strain_miss_primer).to_csv(
                primer_dir+'/'+os.path.splitext(os.path.basename(primer_file))[0]+'_primer_miss.txt',
                sep='\t', index=False)
        print('.OK')


def parse_input():
    parser = argparse.ArgumentParser(
        description='Count hypervariable regions per strain.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input-fasta-archaea', required=True, 
                        help='Input 16s_from_genomes_archaea FASTA file.')
    parser.add_argument('--input-fasta-bacteria', required=True, 
                        help='Input 16s_from_genomes_bacteria FASTA file.')
    parser.add_argument('--input-asssum-archaea', required=True, 
                        help='Input archaea_assembly_summary table.')
    parser.add_argument('--input-asssum-bacteria', required=True, 
                        help='Input bacteria assembly summary table.')
    parser.add_argument('--input-pcr-primers-folder', required=True,
                        help='Directory with FASTA files for PCR primers.')
    parser.add_argument('--output-dir', required=True, 
                        help='Output directory for the final FASTA and table.')
    parser.add_argument('--overhang', required=False, default=0,
                        help='Extra overhang on both sides of extracted hypervariable regions.')
    parser.add_argument('--min-seq-len', required=False, default=200,
                        help='Minimum allowed hypervariable region sequence length.')
    parser.add_argument('--hypervar-region-filter', required=False, default='V',
                        help='Text filter for running one or more specific hypervariable regions.')
    parser.add_argument('--nthreads', required=False, default=4,
                        help='Number of parallel threads to use in the BLAST step.')
    args = parser.parse_args()
    return args
    
                
#%% do this if script is ran from comamnd line
if __name__ == '__main__':

    args = parse_input()
    
    ass_fasta_a = args.input_fasta_archaea
    ass_fasta_b = args.input_fasta_bacteria
    assembly_file_a = args.input_asssum_archaea
    assembly_file_b = args.input_asssum_bacteria
    primer_dir = args.input_pcr_primers_folder
    fasta_out_dir = args.output_dir
    overhang = int(args.overhang)
    min_len = int(args.min_seq_len)
    filter = args.hypervar_region_filter
    nthreads = int(args.nthreads)
    
    # ass_fasta_a = sys.argv[1]       # Assembly FASTA for archaea
    # ass_fasta_b = sys.argv[2]       # Assembly FASTA for bacteria
    # assembly_file_a = sys.argv[3]   # Assembly summary filtered for archaea
    # assembly_file_b = sys.argv[4]   # Assembly summary filtered for bacteria
    # primer_dir = sys.argv[5]        # Folder with FASTA files for PCR primers
    # fasta_out_dir = sys.argv[6]     # Output folder for the final FASTA and table
    
    # if len(sys.argv) >= 8:          # Argument 7 is overhang
    #     overhang = sys.argv[7]
    # else:
    #     overhang = 0
    
    # if len(sys.argv) >= 9:          # Argument 8 is minimum acceptable tag length
    #     min_len = int(sys.argv[8])       # after primer alignment.
    # else:
    #     min_len = 200
    
    # if len(sys.argv) >= 10:          # Argument 9 is text filter for hypervariable region
    #     filter = sys.argv[9]
    # else:
    #     filter = 'V'
        
    
    main(ass_fasta_a, ass_fasta_b, assembly_file_a, assembly_file_b, primer_dir, 
         fasta_out_dir, blast_path = get_blast_path(), min_len=min_len,
         overhang=overhang, primer_file_filter=filter, nthreads=nthreads)
