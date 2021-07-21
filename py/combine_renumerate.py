#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 09:35:43 2021

@author: igor
"""

import re
import sys
import pandas as pd
import argparse
from _include import log




def parse_input():
  parser = argparse.ArgumentParser(
    description='Combine and renumerate multiple hypervariable regions.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--tables', required=True, 
                      help='A list of tables, separated by comma.')
  parser.add_argument('--out-table', required=True,
                      help='Output table name.')
  parser.add_argument('--out-fasta', required=True,
                      help='Output FASTA name.')
  args = parser.parse_args()
  return args



def write_df_to_fasta (df, out_fa, meta_col='meta', seq_col='sequence'):
    """ Convert data frame to a fasta file. """
    with open(out_fa, 'w') as fa:
        for i in range(df.shape[0]):
            entry_i = df.iloc[i]
            fa.write('>'+entry_i[meta_col]+'\n')
            fa.write(entry_i[seq_col]+'\n')



if __name__ == '__main__':
    log('\nCombiner and renumerator.\n')
    
    args = parse_input()
    # print(args)
    
    tables = args.tables.split(',')
    output_table = args.out_table
    output_fasta = args.out_fasta
    
    # Load tables into a list of data frames then use pd.concat
    log('  Loading INPUT tables...')
    tables_df_list = [pd.read_table(t, sep='\t') for t in tables]
    log(' OK.\n')
    log(' Concatenating...')
    tables_df = pd.concat(tables_df_list)
    log(' OK.\n')
    
    # Generate strain names
    log('  Generating strain names...')
    if 'strain_name' in tables_df.columns:
        tables_df.rename(columns={'strain_name':'variant_name'}, inplace=True)
    tables_df['strain_name'] = [re.sub('_@rrn[0-9]+$', '', x) for x in 
                                tables_df['variant_name']]
    log(' OK.\n')

    # Relabel variant names (take into account all unique sequences from all
    # hypervariable regions)
    log('  Renumerating sequences...')
    tables_df_straingroups = tables_df.groupby('strain_name')
    tables_renum_list = list()
    for i, table_group in enumerate(tables_df_straingroups):
        strain_name = table_group[0]
        table_renum = table_group[1]
        good_ids = [j for j, s in enumerate(table_renum['sequence']) 
                    if str(s) != 'nan']
        table_renum = table_renum.iloc[good_ids]
        table_renum['variant_name'] = [strain_name+'_@rrn'+str(j).zfill(2) 
                                       for j in range(len(table_renum))]
        tables_renum_list.append(table_renum)
    # tables_df['variant_name'] = 
    tables_renum_df = pd.concat(tables_renum_list)
    log(' OK.\n')

    # Write output table and FASTA for the combined sequences
    log('  Saving output... ')
    tables_renum_df.pop('strain_name')
    # tables_renum_df.rename(columns={'variant_name':'strain_name'}, inplace=True)
    write_df_to_fasta(tables_renum_df, output_fasta, meta_col='strain_name',
                      seq_col='sequence')
    tables_renum_df.to_csv(output_table, sep='\t', index=False)
    log(' OK.\n')
    # print(tables)
    # print(fastas)
