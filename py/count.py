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
# from merge_16sntncbi_fullgenomes import strain_name_from_genome_string
import sys, os, re
import pandas as pd
from io import StringIO
# from shlex import quote
from subprocess import Popen, PIPE


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
                                blast_ws=7, match=5, mismatch=-4, blast_go=8, blast_ge=6,
                                blast_format='"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score"'):
    """ Blast primers vs blast db made out of full 16S sequences. Return the output
    as Blast outfmt 6 data frame. """
    blast = Popen(' '.join([blast_path, '-subject', primer_file, '-query', sequences_fa, 
                   '-word_size', str(blast_ws), '-outfmt', blast_format, '-strand', 'both',
                   '-reward', str(match), '-penalty', str(mismatch),
                   '-gapopen', str(blast_go), '-gapextend', str(blast_ge),
                   '-dust', 'no']), 
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


#%% strain names from assembly table

#def generate_strain_names_in_assembly_table (ass_df, verbose=False, replace_space=None):
#    """ Generate a 'strain name' column in the data frame of the assembly summary
#    table. """    
#    ass_df['infraspecific_name'] = [str(x).replace('strain=', '').replace('nan', '') \
#          for x in ass_df['infraspecific_name']]
#
#    ass_df['strain_name'] = [strain_name_from_genome_string(x) for x in ass_df['organism_name']]
#
#    # Fix bacterial full genome table annotations
#    # Check if infraspecific name is in strain_name. If not, then concat it.
#    # print('- adding missing strain names ', end='')
#    n_miss_strain = 0
#    tmp = []
#    for i in range(len(ass_df)):
#        strain_name = ass_df['strain_name'][i]
#        strain = ass_df['infraspecific_name'][i]
#        strain = re.sub('substr.*', '', strain)
#        if strain not in strain_name and \
#           strain.replace(' ', '') not in strain_name and \
#           '(' not in strain_name:
#            ass_df.loc[i, 'strain_name'] = strain_name + ' ' + strain
#            n_miss_strain = n_miss_strain + 1
#            tmp.append(i)
#        # Progress bar
#        if verbose and  (i % 1000 == 0 or i == len(ass_df)-1):
#            print('\r* adding missing strain names (% 5d out of % 5d processed)' % (i, len(ass_df)), end='')
#    if replace_space is not None:
#        ass_df['strain_name'] = [s.replace(' ', replace_space) for s in ass_df['strain_name']]
#    if verbose:
#        print(' OK.\n'+str(n_miss_strain)+' strain ids added.')
    

#%% Filter the assembly table by keeping only the best quality entries
# for each unique strain.

def filter_assembly_table (ass_df, verbose=False):
    # Add columns that describe the quality of each column we are considering
    # Run this after adding strain_name column.
    ass_df['assembly_level_int'] = 0
    ass_df.loc[ass_df['assembly_level'] == 'Complete Genome', 'assembly_level_int'] = 4
    ass_df.loc[ass_df['assembly_level'] == 'Chromosome', 'assembly_level_int'] = 3
    ass_df.loc[ass_df['assembly_level'] == 'Scaffold', 'assembly_level_int'] = 2
    ass_df.loc[ass_df['assembly_level'] == 'Contig', 'assembly_level_int'] = 1
    ass_df['refseq_category_int'] = 0
    ass_df.loc[ass_df['refseq_category'] == 'reference genome', 'refseq_category_int'] = 2
    ass_df.loc[ass_df['refseq_category'] == 'representative genome', 'refseq_category_int'] = 1
    
    # Now iterate over all unique strains   
    ass_df_grp = ass_df.groupby('strain_name')
    ass_df_filter_list = []
    # Keep track of number of every category we select when there are multiple
    # assemblies for the exact strame strain name.
    # For each category select the newest genome (by seq_rel_date).
    no_ties = 0
    no_groups = len(ass_df_grp)    
    for i, (strain_name, gr) in enumerate(ass_df_grp):        
        # Find the best ranked 
        gr_b = gr.loc[gr['assembly_level_int'] == max(gr['assembly_level_int'])]
        gr_b = gr_b.loc[gr_b['refseq_category_int'] == max(gr_b['refseq_category_int'])]
        gr_b = gr_b.loc[gr_b['seq_rel_date'] == max(gr_b['seq_rel_date'])]        
        if len(gr_b) > 1:
            # Tie, grab the first one in the list
            no_ties = no_ties + 1
            gr_b = gr_b[1:2]
        ass_df_filter_list.append(gr_b)
        # Progress bar
        if verbose and (i % 100 == 0 or i == len(ass_df)-1):
            print('\r* find unique assembly data (% 5d out of % 5d)' % (i, no_groups), end='')
    if verbose:
        print(' OK.\n')            
    return pd.concat(ass_df_filter_list, ignore_index=True)


#%% main call
def main(ass_fasta_a, ass_fasta_b, assembly_file_a, assembly_file_b, 
         primer_dir, fasta_out_dir, 
         blast_path='/usr/local/ncbi/blast/bin/blastn', overhang=21, 
         primer_file_filter='V'):
    
    print('Count 16S hypervariable regions')
    # Load full genome assembly 16S sequences
    # path = os.path.expanduser('~/cloud/research/microbiome/genomes/data/bacteria/')
    #    blast_path = ['/usr/local/ncbi/blast/bin/blastn', 
    #                  'C:/Program Files/NCBI/blast-2.7.1+/bin/'][sys.platform == 'win32']
    # ass_fasta = path+'16s_from_genomes_2017-07-20.fasta'
    
    # If we want to process only specific primer set add it here as a string
    # with exact match to the FASTA filename (output from primer_all_combinations.py)
    #  primer_file_filter = 'V3-V4_341F-805R'
    
    # For bacteria use wpartial file
    # assembly_file = path+'bacteria_assembly_summary_sub_2017-07-20.txt'
    # out_folder = os.path.expanduser('~/data/fasta_by_assid_fix')
    # ass_df = pd.read_csv(assembly_file, sep='\t')
    
    # PCR primer files
    # Just load all files from the folder into a list and iterate over it.
    # primer_dir = os.path.expanduser('~/cloud/research/microbiome/genomes/data/pcr_primers')
    # fasta_out_dir = os.path.expanduser('~/cloud/research/microbiome/genomes/data/vregions_db')    
    
    # Overhang from each V-region. Need this to make sure we align full query.
    overhang = int(overhang)
    
    primer_files = [os.path.join(primer_dir, p) for p in os.listdir(primer_dir) if 
                    os.path.splitext(p)[1] == '.fasta' and primer_file_filter in p]
    

    # Load the table between full accession id (from FASTA) and full 16S sequence.
    # ass_dict = fasta_to_dict_of_lists(ass_fasta, sep=';')
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
    # Generate strain names from organism_name and infraspecific_name entries
    # generate_strain_names_in_assembly_table(ass_df, verbose=True, replace_space='_')
    
    # Filter the table by keeping only strain names with best quality assemblies
    # ass_df = filter_assembly_table(ass_df, verbose=True)
    
    # Rename the assembly accession column
    ass_df = ass_df.rename(columns={'# assembly_accession': 'ass_id'})
    ass_df['strain_name'] = [s.replace(' ', '_') for s in ass_df['strain_name']]

    # After filter also filter out ass_to_seq_df by ass_id column
    ass_to_seq_df = ass_to_seq_df[ass_to_seq_df['ass_id'].isin(ass_df['ass_id'])]

    # Generate a dictionary between accession id and strain name for easy access later
    ass_to_strain_dict = ass_df[['ass_id', 'strain_name']].set_index('ass_id').to_dict()['strain_name']
    
    # BLAST path
    # blast_path = '/usr/local/ncbi/blast/bin/blastn'
    # blast_db = os.path.expanduser('~/cloud/research/microbiome/genomes/data/bacteria/16s_from_genomes_2017-07-20')

    # Iterate over each primer pair combination
    for primer_file in primer_files:
        
        print('\r* '+os.path.splitext(os.path.basename(primer_file))[0]+': blast...', end='')
        # Run BLAST with 16s sequences vs PCR primer sequences
        blast_out_a_df = blast_primers_vs_sequences(primer_file, ass_fasta_a, blast_path)
        blast_out_b_df = blast_primers_vs_sequences(primer_file, ass_fasta_b, blast_path)
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
            if seq == '': # Empty sequence
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
    
                
#%% do this if script is ran from comamnd line
if __name__ == '__main__':
    #%% fasta processor
    # ass_fasta_a, ass_fasta_b, assembly_file_a, assembly_file_b, primer_dir, fasta_out_dir
    ass_fasta_a = os.path.expanduser('~/cloud/research/microbiome/genomes/data/16s_from_genomes_archaea_2018-02-23.fasta')
    ass_fasta_b = os.path.expanduser('~/cloud/research/microbiome/genomes/data/16s_from_genomes_bacteria_2018-02-23.fasta')
    assembly_file_a = os.path.expanduser('~/cloud/research/microbiome/genomes/data/archaea_assembly_summary_filter_2018-02-23.txt')
    assembly_file_b = os.path.expanduser('~/cloud/research/microbiome/genomes/data/bacteria_assembly_summary_filter_2018-02-23.txt')
    primer_dir = os.path.expanduser('~/cloud/research/microbiome/genomes/data/pcr_primers')
    fasta_out_dir = os.path.expanduser('~/cloud/research/microbiome/genomes/data/vregions_db')
    #assembly_file_a = ''
    #assembly_file_b = ''
    
    ass_fasta_a = sys.argv[1]
    ass_fasta_b = sys.argv[2]
    assembly_file_a = sys.argv[3]
    assembly_file_b = sys.argv[4]
    primer_dir = sys.argv[5]
    fasta_out_dir = sys.argv[6]
    if len(sys.argv) >= 8:
        overhang = sys.argv[7]
    else:
        overhang = 0
    if len(sys.argv) >= 9:
        filter = sys.argv[8]
    else:
        filter = 'V'
    
    main(ass_fasta_a, ass_fasta_b, assembly_file_a, assembly_file_b, primer_dir, fasta_out_dir,
         overhang=overhang, primer_file_filter=filter)
