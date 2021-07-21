#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update single RExMapDB.

Created on Sun Mar 28 17:32:48 2021

@author: igor
"""

import argparse
import os, sys, glob, re
from subprocess import Popen
from _include import log, get_python_path, get_rscript_path

script_title = 'RExMapDB update specific hypervariable region database'

def parse_input():
    parser = argparse.ArgumentParser(
        description=script_title,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser._action_groups.pop()
    parser_req = parser.add_argument_group('required arguments')
    # Required arguments
    parser_req.add_argument('-d', '--data', required=True, 
                        help='Input rexmapdb/data directory.')
    parser_req.add_argument('-r', '--hypervar-region', required=True, 
                        help='Hypervariable region database designation including primer labels (e.g. V4-V5_515F-926R).')
    # Optional arguments
    parser_opt = parser.add_argument_group('optional arguments')
    parser_opt.add_argument('--data-subfolder', action='store_true',
                            help='Create a subfolder inside --data folder to store all '+
                            'intermediate and output files.')
    parser_opt.add_argument('--primer-table', required=False, default='primers/pcr_primers_table.txt',
                        help='Full path and filename to the pcr_primers_table.txt. ')
    parser_opt.add_argument('--primer-dir', required=False, default='data/pcr_primers/',
                            help='PCR primers folder in data; output of primers_fasta.py.')
    parser_opt.add_argument('--overhang', required=False, default=22, type=int,
                        help='Hypervariable region overhang.')
    parser_opt.add_argument('--nthreads', required=False, default=10, type=int,
                        help='Number of parallel threads to run.')
    parser_opt.add_argument('--len-pct-var', required=False, default=20, type=float,
                        help='Low end variation for minimum alignment length. To obtain this'+
                        ' length, calculate reverse primer - foward primer * (1 - len_pct_var/100).'+
                        ' e.g. For V3-V4 337F/805R and 20pct default , min_aln_len = (805-337)*0.8 = 374 nt')
    # parser_opt.add_argument('--multiplexed-database', action='store_true',
    #                         help='Generate a multiplexed database using a '+
    #                         'combine_renumerate.py script. Multiple hypervariable'+
    #                         ' regions should be comma delimited with no space.')
    parser_opt.add_argument('--database', required=False, default='database/',
                            help='Path to the output database folder.')
    parser_opt.add_argument('--date', required=False, default=None,
                            help='Use this date for database name instead of the latest date '+
                            'automatically extracted from assembly summary tables.')
    parser_opt.add_argument('--show-calls', action='store_true',
                            help='Show full calls to external scripts.')
    parser_opt.add_argument('--overwrite', action='store_true',
                            help='Overwrite previously generated files.')
    parser_opt.add_argument('--variant-table-suffix', default='',
                            help='Suffix for the final text '+
                            '(Used to be: _wrefseq_table_unique_variants_R).')
    parser_opt.add_argument('--variant-blastdb-suffix', default='',
                            help='Suffix for the BLAST database files.')
    args = parser.parse_args()
    return args

def find_latest_nongenome_16s_fasta(data_dir):
    """ Find the latest data/16s_RefSeq_2018-05-20.fasta file. """
    if not os.path.exists(data_dir):
        sys.exit('\nError: Input data folder does not exist.')
    
    # Search for the latest file
    ng_16s_fa_list = glob.glob(os.path.join(data_dir, '16s_RefSeq_*.fasta'))
    if len(ng_16s_fa_list) == 0:
        return({'success': False})
    
    out = {'input-nongenome-16s-fasta': sorted(ng_16s_fa_list, reverse=True)[0],
           'success': True}
    return(out)

def find_latest_inputs(data_dir):
    """ Extract latest input files, as a dictionary output. """
    # Search through data_dir for files
    if not os.path.exists(data_dir):
        sys.exit('\nError: Input data folder does not exist.')
        
    # Search for all of the latest assembly summary filter files
    a_sum = glob.glob(os.path.join(data_dir, 'archaea_assembly_summary_filter_*.txt'))
    b_sum = glob.glob(os.path.join(data_dir, 'bacteria_assembly_summary_filter_*.txt'))
    a_16g = glob.glob(os.path.join(data_dir, '16s_from_genomes_archaea_*.fasta'))
    b_16g = glob.glob(os.path.join(data_dir, '16s_from_genomes_bacteria_*.fasta'))
    # Check these searches
    if len(a_sum) == 0 or len(b_sum) == 0 or len(a_16g) == 0 or len(b_16g) == 0:
        return({'success': False})
    
    
    # Extract dates from each of these files; they are in YYYY-MM-DD format.
    # sum_re = re.compile(r".*(archaea|bacteria)_assembly_summary_filter_([0-9]{4}-[0-9]{2}-[0-9]{2})\.txt$")
    
    out = {'input-asssum-archaea': sorted(a_sum, reverse=True)[0], 
           'input-asssum-bacteria': sorted(b_sum, reverse=True)[0],
           'input-fasta-archaea': sorted(a_16g, reverse=True)[0], 
           'input-fasta-bacteria': sorted(b_16g, reverse=True)[0],
           'success': True}
    return(out)    

def get_date_from_asssum(filename):
    sum_re = re.compile(r".*(archaea|bacteria)_assembly_summary_filter_([0-9]{4}-[0-9]{2}-[0-9]{2})\.txt$")
    sum_date = sum_re.match(filename).group(2)
    return(sum_date)

def calc_min_aln_len(hr, pct):
    """ Calculate minimum alignment length based on the hypervariable region
    string "hr". It has a format V[0-9]-V[0-9]_[0-9]+[a-zA-Z]+-[0-9]+[a-zA-Z]+"""
    hr_re = re.compile(r"^.*V[0-9][-]?[V]?[0-9]?_([0-9]+)[a-zA-Z]+-([0-9]+)[a-zA-Z]+.*$")
    hr_ml = abs(float(hr_re.match(hr).group(1)) - float(hr_re.match(hr).group(2)))*(1-pct/100)
    return(hr_ml)


if __name__ == '__main__':
    
    log('---------'+script_title+'---------')
    args = parse_input()
    data_dir = args.data
    subfolder = args.data_subfolder
    hr = args.hypervar_region
    pct = args.len_pct_var
    primer_dir = args.primer_dir
    primer_tab = args.primer_table
    overhang = args.overhang
    db_path = args.database
    nthreads = max(int(args.nthreads), 1)
    overwrite = args.overwrite
    show_calls = args.show_calls
    # multiplexed = args.multiplexed_database
    date = args.date
    
    suffix_variant_table = args.variant_table_suffix
    suffix_blastdb_files = args.variant_blastdb_suffix
    
    # Check if we have multiple valid hypervariable regions in case we are
    # running a multiplexed database genenration

    #--------------------- count.py ------------------------------------------
    # Determine inputs for count.py script
    count_inputs = find_latest_inputs(data_dir)
    if not count_inputs['success']:
        sys.exit('\nError: Some or all input files not found.')
    
    data_dir_main = data_dir
    if subfolder:
        data_dir = os.path.join(data_dir, hr)
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
            log('Using subfolder: '+hr)
            log('Created output dir: '+data_dir)        
        
    log('Output files written to: '+data_dir)
    
    latest_date_method = ''
    if date is None:
        latest_date = max([get_date_from_asssum(count_inputs['input-asssum-archaea']), 
                           get_date_from_asssum(count_inputs['input-asssum-bacteria'])])
        latest_date_method = '(auto)'
    else:
        latest_date = date
        latest_date_method = '(user)'
    
    log('Database date/version '+latest_date_method+': '+latest_date)
    # sys.stdout.write('Output files will have the HYPERVARIABLE-REGION_DATE_hangOVERHANG prefix.\n')
    log('Count')    
    log('  Input:')
    for k, v in count_inputs.items():
        if k != 'success':
            log('  --'+k+' '+v+'\n')
    
    # Check if PCR primer fasta file exists
    hr_fasta = os.path.join(primer_dir, hr+'.fasta')
    if not os.path.exists(hr_fasta):
        sys.exit('\nError: FASTA file '+hr_fasta+ ' does not exist. '+
                 'Re-run primers_fasta.py after updating pcr_primers_table.txt.')
    
    ml = int(calc_min_aln_len(hr, pct))
    
    # Generate a count system call
    python_path = get_python_path()
    if python_path is None:
        sys.exit('\nError: Python3 not found.')
        
    # Check if count.py outputs exist before running it
    # data/16s_from_genomes_2018-05-20_V3-V4_337F-805R_hang22_counts.txt
    # data/V4-V6_515F-1185mR_hang22_sequences.fasta
    count_fa_out = os.path.join(data_dir, hr+'_hang'+str(overhang)+'_sequences.fasta')      
    count_tab_out = os.path.join(
        data_dir, 
        os.path.splitext(os.path.basename(count_inputs['input-fasta-archaea'].replace('archaea_', '')))[0]+
        '_'+hr+'_hang'+str(overhang)+'_counts.txt'
    )
    log('  --input-pcr-primers-folder '+primer_dir+'\n')
    log('  --hypervar-region-filter '+hr+'\n')
    log('  Output:\n')
    log('  --output-dir '+data_dir+'\n')
    log('  Options:\n')
    log('  --min-seq-len '+str(ml)+'\n')
    log('  --overhang '+str(overhang)+'\n')
    log('  --hypervar-region-filter '+hr+'\n')
    log('  --nthreads '+str(nthreads)+'\n')
    
    if os.path.exists(count_fa_out) and os.path.exists(count_tab_out) and not overwrite:
        log('  count.py already completed.\n')
    else:
        # Run count.py if --overwrite is given or any of the files are missing
        count_path = os.path.join(os.path.dirname(sys.path[0]), 'py', 'count.py')
        count_call = ' '.join([python_path, 'py/count.py', 
                              '--input-fasta-archaea', count_inputs['input-fasta-archaea'],
                              '--input-fasta-bacteria', count_inputs['input-fasta-bacteria'],
                              '--input-asssum-archaea', count_inputs['input-asssum-archaea'],
                              '--input-asssum-bacteria', count_inputs['input-asssum-bacteria'],
                              '--input-pcr-primers-folder', primer_dir,
                              '--output-dir', data_dir,
                              '--overhang', str(overhang),
                              '--min-seq-len', str(ml),
                              '--hypervar-region-filter', hr,
                              '--nthreads', str(nthreads)])
        if show_calls:
            log('  Count call: '+count_call+'\n')
        count_cmd = Popen(count_call, shell=True)
        out, err = count_cmd.communicate()
        if out == b'':
            sys.exit('\nError: running count.py call.')
    
    #------------------------ add_refseq.py -----------------------------------
    log('Add RefSeq\n')
    ng_16s_fasta = find_latest_nongenome_16s_fasta(data_dir_main)
    if not ng_16s_fasta['success']:
        sys.exit('\nError: Non-genome 16S fasta not found.')
    log('  Input files:\n')
    log('  --input-nongenome-16s-fasta '+ng_16s_fasta['input-nongenome-16s-fasta']+'\n')
    
    addrefseq_outtab = os.path.join(data_dir, hr+'_hang'+str(overhang)+'_wrefseq_table.txt')
    addrefseq_outfa = os.path.join(data_dir, hr+'_hang'+str(overhang)+'_wrefseq_sequences.fasta')
    log('  --input-genome-16s-counts '+count_tab_out+'\n')
    log('  --input-genome-16s-fasta '+count_fa_out+'\n')
    log('  --input-pcr-primers-fasta '+hr_fasta+'\n')
    log('  Output files:\n')
    log('  --output-table '+addrefseq_outtab+'\n')
    log('  --output-fasta '+addrefseq_outfa+'\n')
    log('  Options:\n')
    log('  --min-seq-len '+str(overhang)+'\n')
    log('  --overhang '+str(ml)+'\n')    
    
    if not os.path.exists(addrefseq_outtab) or not os.path.exists(addrefseq_outfa) or overwrite:
        addrefseq_call = ' '.join([python_path, 'py/add_refseq.py', 
                              '--input-nongenome-16s-fasta', ng_16s_fasta['input-nongenome-16s-fasta'],
                              '--input-genome-16s-counts', count_tab_out,
                              '--input-genome-16s-fasta', count_fa_out,
                              '--input-pcr-primers-fasta', hr_fasta,
                              '--output-table', addrefseq_outtab,
                              '--output-fasta', addrefseq_outfa,
                              '--min-seq-len', str(ml),
                              '--overhang', str(overhang)])
        if show_calls:
            log('  Add_RefSeq call: '+addrefseq_call+'\n')
        addrefseq_cmd = Popen(addrefseq_call, shell=True)
        out, err = addrefseq_cmd.communicate()
        if out == b'':
            sys.exit('\nError: running add_refseq.py call.')
    else:
        log('  add_refseq.py already completed.\n')
    
    #---------------- Optional multiplexed database --------------------------
#     python3 /data/rexmapdb/py/combine_renumerate.py \
# 	--tables data/V1-V2_27F-338R_hang22_2021-02-13_wrefseq_table.txt,data/V3-V4_337F-805R_hang22_2021-02-13_wrefseq_table.txt \
# 	--out-table data/V1-V2_V3-V4_hang22_2021-02-13_wrefseq_table.txt \
# 	--out-fasta data/V1-V2_V3-V4_hang22_2021-02-13_wrefseq_sequences.fasta

    
    #------------------------ Group ------------------------------------------
    log('Group\n')
    grp_outtab = os.path.join(data_dir, hr+'_hang'+str(overhang)+'_wrefseq_table_unique_variants.txt')
    grp_outfa = os.path.join(data_dir, hr+'_hang'+str(overhang)+'_wrefseq_sequences_unique_variants.fasta')
    log('  Input files:\n')
    log('  --input-table '+addrefseq_outtab+'\n')
    log('  --input-fasta '+addrefseq_outfa+'\n')
    log('  Output files:\n')
    log('  --output-table '+grp_outtab+'\n')
    log('  --output-fasta '+grp_outfa+'\n')
    
    if not os.path.exists(grp_outtab) or not os.path.exists(grp_outfa) or overwrite:
        grp_call = ' '.join([python_path, 'py/group.py', 
                              '--input-table', addrefseq_outtab,
                              '--input-fasta', addrefseq_outfa,
                              '--output-table', grp_outtab,
                              '--output-fasta', grp_outfa])
        if show_calls:
            log('  Group call: '+grp_call+'\n')
        grp_cmd = Popen(grp_call, shell=True)
        out, err = grp_cmd.communicate()
        if out == b'':
            sys.exit('\nError: running group.py call.')
    else:
        log('  group.py already completed.\n')
    
    #------------------- Remove taxonomic outliers ---------------------------
    log('Remove taxonomic outliers\n')
    rscript_path = get_rscript_path()
    if rscript_path is None:
        sys.exit('\nError: Cannot find Rscript executable.')
    tax_outtab = os.path.join(db_path, hr+'_'+latest_date+'_hang'+str(overhang)+suffix_variant_table+'.txt')
    # tax_outfa = os.path.join(db_path, hr+'_'+latest_date+'_hang'+str(overhang)+suffix_blastdb_files+'.fasta')
    tax_outfa = os.path.join(data_dir, hr+'_'+latest_date+'_hang'+str(overhang)+suffix_blastdb_files+'.fasta')
    tax_outlier = os.path.join(data_dir, hr+'_excluded_outlier_strains.txt')
    log('  Input files:\n')
    log('  --input-table '+grp_outtab+'\n')
    log('  --input-fasta '+grp_outfa+'\n')
    log('  Output files:\n')
    log('  --output-table '+tax_outtab+'\n')
    log('  --output-fasta '+tax_outfa+'\n')
    log('  --output-excl '+tax_outlier+'\n')
    if not os.path.exists(db_path):
        os.makedirs(db_path)
        log('  Output database folder created: '+db_path+'\n')
    
    if not os.path.exists(tax_outtab) or not os.path.exists(tax_outfa) or overwrite:
        tax_call = ' '.join([rscript_path, '--vanilla', 'R/remove_taxonomic_outliers.R', 
                              '--input-table', grp_outtab,
                              '--input-fasta', grp_outfa,
                              '--output-table', tax_outtab,
                              '--output-fasta', tax_outfa,
                              '--output-excl', tax_outlier])
        if show_calls:
            log('  Remove taxonomic outliers call: '+tax_call+'\n')
        tax_cmd = Popen(tax_call, shell=True)
        out, err = tax_cmd.communicate()
        if out == b'':
            sys.exit('\nError: running remove_taxonomic_outliers.R call.')
    else:
        log('  remove_taxonomic_outliers.R already completed.\n')
    
    #------------------- Generate database -----------------------------------
    log('Generate BLAST database\n')
    # python3 ../rexmapdb/py/makeblastdb.py \
    # --input-fasta database/V4-V5_515F-926R_hang22_2021-02-13_wrefseq_sequences_unique_variants_R.fasta \
    # --out-db-prefix database/V4-V5_515F-926R_hang22_2021-02-13_wrefseq
    db_outpre = os.path.join(db_path, hr+'_'+latest_date+'_hang'+str(overhang))
    db_call = ' '.join([python_path, 'py/makeblastdb.py',
                        '--input-fasta', tax_outfa,
                        '--out-db-prefix', db_outpre])
    if show_calls:
        log('  Make BLAST DB call: '+db_call+'\n')
    db_cmd = Popen(db_call, shell=True)
    out, err = db_cmd.communicate()
    if out == b'':
        sys.exit('\nError: running makeblastdb.py.')
    log('  makeblastdb.py completed.\n')

    log('Done.\n')
    
    
    
