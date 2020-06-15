#!~/anaconda/bin/python3
# -*- coding: utf-8 -*-
"""
This script downloads the main metadata assembly, and from it generates
two files containing URLs of the GFF feature files ("features") and FASTA
sequence files ("sequences"). One GFF and one FASTA is downloaded for each
unique strain where better and newer assemblies are prioritized.

After this script we run download_16s_from_extracted_paths.py to actually
download the sequences and features from BASH.


Created on Thu Jul 20 12:32:39 2017

@author: igor
"""

import os, sys, datetime, re
import pandas as pd
from six.moves import urllib


def download_assembly_summary (sum_file, out_file):
    """ Download a single assembly summary. If the target folder does not
    exist, create it. """
    if not os.path.isdir(out_path):
        os.makedirs(out_path)
    urllib.request.urlretrieve(sum_file, out_file)
    
def strain_name_from_genome_string (text):
    """ Return genus species strain (if any) name from full genome species name text. 
    Sometimes, we have subspecies name exact same as species, so we remove that. """
    x = text.replace('substr. ', '').replace('str. ', '').replace('subsp. ', '').replace('\'', '')
    x = re.sub(r'\b(\w+)( \1\b)+', r'\1', x).replace('[', '').replace(']', '')
    return re.sub(r'\b(\w+)( \1\b)+', r'\1', x)

#%% main call
if __name__ == '__main__':
    
    out_path = sys.argv[1]
    # out_path = '/Users/igor/Downloads/himapdb/data'
    if len(sys.argv) == 3:
        kingdoms = sys.argv[2].strip().split(',')
    else:
        kingdoms = ['archaea', 'bacteria']
    print('Kingdoms: ', ', '.join(kingdoms))
    
    # Generate current date
    dt = datetime.datetime.now()
    date = str(dt.year).zfill(4) + '-' + str(dt.month).zfill(2) + '-' + str(dt.day).zfill(2)
    
    # Generate URLs for assembly summaries for all kingdoms and output file names
    ass_base = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/'
    good_assemblies = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
    
    #%% Download assembly meta data
    
    for k in kingdoms:
        print('Kingdom: '+k)
        
        # Download assembly summary
        print('- Downloading assembly summary...', end='')
        ass_sum = ass_base + k + '/assembly_summary.txt'
        ass_out = os.path.join(out_path, k+'_assembly_summary_'+date+'.txt')
        download_assembly_summary(ass_sum, ass_out)
        print('OK.')
        
        # Load assembly summary 
        print('- Loading assembly summary...', end='')
        ass_df = pd.read_csv(
                ass_sum, delimiter='\t', skiprows=1,
                usecols=['# assembly_accession', 'refseq_category', 'taxid',
                         'organism_name', 'infraspecific_name', 
                         'assembly_level', 'seq_rel_date', 'ftp_path'], 
                dtype={'# assembly_accession': str, 'refseq_category': str, 'taxid': str,
                       'organism_name': str, 'infraspecific_name': str,
                       'assembly_level': str, 'seq_rel_date': str, 'ftp_path': str})
        print('OK.')
                
        # For each strain select only the best assembly
        print('- Processing and fixing assembly summary table')
        ass_df['infraspecific_name'] = [re.sub('_TMP.*', '', str(x).replace('strain=', '').replace('nan', '').replace('substr. ', '')) \
                   for x in ass_df['infraspecific_name']]
        ass_df['strain_name'] = [strain_name_from_genome_string(x) for x in ass_df['organism_name']]
        ass_df['assembly_level_int'] = 0
        ass_df.loc[ass_df['assembly_level'] == 'Complete Genome', 'assembly_level_int'] = 4
        ass_df.loc[ass_df['assembly_level'] == 'Chromosome', 'assembly_level_int'] = 3
        ass_df.loc[ass_df['assembly_level'] == 'Scaffold', 'assembly_level_int'] = 2
        ass_df.loc[ass_df['assembly_level'] == 'Contig', 'assembly_level_int'] = 1
        ass_df['refseq_category_int'] = 0
        ass_df.loc[ass_df['refseq_category'] == 'reference genome', 'refseq_category_int'] = 2
        ass_df.loc[ass_df['refseq_category'] == 'representative genome', 'refseq_category_int'] = 1
        
        # Fix bacterial full genome table annotations
        # Check if infraspecific name is in strain_name. If not, then concat it.
        # print('- adding missing strain names ', end='')
        n_miss_strain = 0
        tmp = []
        for i in range(len(ass_df)):
            
            strain_name = ass_df['strain_name'][i]
            strain = ass_df['infraspecific_name'][i]
            strain = re.sub('substr.*', '', strain)
            if strain not in strain_name and \
               strain.replace(' ', '') not in strain_name and \
               '(' not in strain_name:
                ass_df.loc[i, 'strain_name'] = strain_name + ' ' + strain
                n_miss_strain = n_miss_strain + 1
                tmp.append(i)
            # Progress bar
            if i % 100 == 0 or i == len(ass_df)-1:
                print('\r- '+k+': adding missing strain names (% 5d out of % 5d processed)' % (i+1, len(ass_df)), end='')
        print(' OK.\n- '+str(n_miss_strain)+' strain ids added.')
        tmp = None
    
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
            if i % 100 == 0 or i == no_groups-1:
                print('\r- find unique assembly data (% 5d out of % 5d)' % (i+1, no_groups), end='')
    
        ass_df_filter = pd.concat(ass_df_filter_list, ignore_index=True)        
        ass_df_filter_list = None
        print('\n- filtered list from full genome list obtained')
    
        # Write filtered list to file (we will use this in other scripts)
        print('- Writing clean assembly summary table...', end='')
        ass_filter_out = os.path.join(out_path, k+'_assembly_summary_filter_'+date+'.txt')
        ass_df_filter.to_csv(ass_filter_out, sep='\t', index=False)
        print('OK.')
        
        # Save a separate list of feature and sequences for download
        fea_out = os.path.join(out_path, k+'_assembly_features_urls.txt')
        seq_out = os.path.join(out_path, k+'_assembly_sequences_urls.txt')
    
        with open(fea_out, 'w') as f_out, open(seq_out, 'w') as s_out:
            for i, row in enumerate(ass_df_filter.itertuples()):
                f_out.write(os.path.join(row[8], os.path.basename(row[8])+'_feature_table.txt.gz\n'))
                s_out.write(os.path.join(row[8], os.path.basename(row[8])+'_genomic.fna.gz\n'))
                if i % 100 == 0 or i == len(ass_df_filter)-1:
                    print('\r- Processed % 5d out of % 5d assembly records.' % (i+1, len(ass_df_filter)), end='')
        print('\n- Saved feature and sequence URLs.')
