#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Count the number of unique genuses, species and strains in our
database and compare it to SILVA. This script is not used to 
generate HiMAP db, just to make statistics and compare it to
SILVA.

Created on Fri Apr 20 08:30:28 2018

@author: igor
"""

from add_refseq import fasta_to_dict
from collections import Counter
import pandas as pd
import sys, re

# ncbi_genome_file_a = '/Users/igor/cloud/research/microbiome/genomes/data/archaea_assembly_summary_filter_2018-02-23a.txt'
# ncbi_genome_file_b = '/Users/igor/cloud/research/microbiome/genomes/data/bacteria_assembly_summary_filter_2018-02-23a.txt'
# ncbi_search_file = '/Users/igor/cloud/research/microbiome/genomes/data/16s_rrna_ncbi_search_refseq_2018-03-19.fasta'

#himap_v4_file = '/Users/igor/cloud/research/microbiome/himap/inst/database/V4_515F-805R_hang22_wrefseq_table_unique_variants_R.txt'
#silva_strains_file = '/Users/igor/cloud/research/microbiome/silva/SILVA_132_SSURef_tax_silva_strains.txt'
#silva_counts_file = '/Users/igor/cloud/research/microbiome/silva/SILVA_132_SSURef_tax_silva_species_unique_counts_table.txt'
#genus_out_file = '/Users/igor/cloud/research/microbiome/silva/genuses_silva_vs_himap.txt'
#species_out_file = '/Users/igor/cloud/research/microbiome/silva/species_silva_vs_himap.txt'

# Load and process NCBI search fasta. This gives us strain names.
#ncbi_search_dict = fasta_to_dict(ncbi_search_file, 
#                                 post_process='strain_name_from_refseq_string', 
#                                 post_process_seq='seq_to_basic_code')

def main (himap_v4_file, silva_strains_file, silva_counts_file,
          genus_out_file, species_out_file):

    himap_v4_df = pd.read_csv(himap_v4_file, sep='\t')
    
    # Load and process NCBI full genome fasta.
    #ncbi_genome_dict_a = pd.read_csv(ncbi_genome_file_a, delimiter='\t')
    #ncbi_genome_dict_b = pd.read_csv(ncbi_genome_file_b, delimiter='\t')
    
    # Add all strain names to a set
    #unique_strains = set([])
    #unique_strains.update([s.replace('_', ' ') for s in ncbi_search_dict.keys()])
    #unique_strains.update(ncbi_genome_dict_a['strain_name'])
    #unique_strains.update(ncbi_genome_dict_b['strain_name'])
    
    # Now just extract unique strains
    unique_strains = [re.sub(r'_@rrn[0-9]+', '', x).replace('_', ' ') for x in himap_v4_df['strain_name']]
    uniq_strains = sorted(list(unique_strains))
    
    # Generate unique species
    species = [' '.join(s.split(' ')[0:2]) for s in uniq_strains]
    species = [s for s in species if not s.endswith('sp.')]
    uniq_species = list(set([' '.join(s.split(' ')[0:2]) for s in uniq_strains if not s.endswith('sp.')]))
    
    # 117,307 unique strains and 16,349 unique species
    uniq_genuses = list(set([s.split(' ')[0] for s in uniq_species]))
    
    # How many species and strains we have per genus
    genuses_ft_st = Counter([s.split(' ')[0] for s in uniq_strains])
    genuses_ft_sp = Counter([s.split(' ')[0] for s in uniq_species])
    
    # Load SILVA data
    silva_strains = []
    with open(silva_strains_file) as silva_f:
        for line in silva_f:
            silva_strains.append(line.strip())
    
    
    # Unique silva species
    silva_species = []
    with open(silva_counts_file) as silva_f:
        for line in silva_f:
            silva_species.append(' '.join(line.strip().split(' ')[1:3]))
            
    silva_genuses_ft_st = Counter([s.split(' ')[0] for s in silva_strains])
    silva_genuses_ft_sp = Counter([s.split(' ')[0] for s in silva_species])
    
    # Now plot silva number of species per genus vs HiMAP
    silva_genuses_uniq = set(silva_genuses_ft_st.keys())
    himap_genuses_uniq = set(genuses_ft_st.keys())
    
    common_genuses = silva_genuses_uniq.union(himap_genuses_uniq)
    
    with open(genus_out_file, 'w') as out_f:
        out_f.write('genus\tsilva_num_st\tsilva_num_sp\thimap_num_st\thimap_num_sp\n')
        for g in common_genuses:
            out_f.write('\t'.join([g, str(silva_genuses_ft_st[g]), str(silva_genuses_ft_sp[g]),
                                   str(genuses_ft_st[g]), str(genuses_ft_sp[g])])+'\n')
    
        
    silva_species_in_strains = [' '.join(s.split(' ')[0:2]) for s in silva_strains]
    silva_species_in_strains = [s for s in silva_species_in_strains if not s.endswith('sp.')]
    
    himap_species_count = Counter(species)
    silva_species_count = Counter(silva_species_in_strains)
    himap_species_uniq = set(himap_species_count.keys())
    silva_species_uniq = set(silva_species_count.keys())
    common_species = himap_species_uniq.union(silva_species_uniq)
    
        
    with open(species_out_file, 'w') as out_f:
        out_f.write('species\thimap_num_st\tsilva_num_st\n')
        for s in common_species:
            out_f.write('\t'.join([s, str(himap_species_count[s]), 
                                   str(silva_species_count[s])])+'\n')

#
    
    
if __name__ == '__main__':
    
    himap_v4_file = sys.argv[1]
    silva_strains_file = sys.argv[2]
    silva_counts_file = sys.argv[3]
    genus_out_file = sys.argv[4]
    species_out_file = sys.argv[5]
    
    main(himap_v4_file, silva_strains_file, silva_counts_file,
         genus_out_file, species_out_file)