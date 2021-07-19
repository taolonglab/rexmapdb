#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update single multiplexed RExMapDB.

Created on Tue Mar 30 15:37:28 2021

@author: igor
"""

import argparse
import os, sys, glob, re
from subprocess import Popen, PIPE
from datetime import datetime

script_title = 'RExMapDB update specific MULTIPLEXED hypervariable region database'

def days_between(d1, d2):
    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)

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
    parser_req.add_argument('-r', '--hypervar-regions', required=True, 
                        help='Hypervariable regions database designations,'+
                        ' delimited by comma with no spaces, e.g.'+
                        ' V1-V2_27F-338R,V4-V5_515F-926R.')
    # Optional arguments
    parser_opt = parser.add_argument_group('optional arguments')
    parser_opt.add_argument('--data-subfolder', action='store_true',
                            help='Create a subfolder inside --data folder to store all '+
                            'intermediate and output files.')
    parser_opt.add_argument('--overhang', required=False, default=22, type=int,
                        help='Preferred hypervariable region overhang for input files '+
                        ' (in case of multiple files with the same date).')
    parser_opt.add_argument('--nthreads', required=False, default=10, type=int,
                        help='Number of parallel threads to run.')
    parser_opt.add_argument('--database', required=False, default='database/',
                            help='Path to the output database folder.')
    parser_opt.add_argument('--date', required=False, default=None,
                            help='Use this date for database name instead of the latest date '+
                            'automatically extracted from assembly summary tables.')
    parser_opt.add_argument('--show-calls', action='store_true',
                            help='Show full calls to external scripts.')
    parser_opt.add_argument('--overwrite', action='store_true',
                            help='Overwrite previously generated files.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
#     python3 /data/rexmapdb/py/combine_renumerate.py \
# 	--tables data/V1-V2_27F-338R_hang22_2021-02-13_wrefseq_table.txt,data/V3-V4_337F-805R_hang22_2021-02-13_wrefseq_table.txt \
# 	--out-table data/V1-V2_V3-V4_hang22_2021-02-13_wrefseq_table.txt \
# 	--out-fasta data/V1-V2_V3-V4_hang22_2021-02-13_wrefseq_sequences.fasta
    # Combined: V1-27F-V2-338R_
    print(script_title)
    args = parse_input()
    data_dir = args.data
    hrs = args.hypervar_regions
    subfolder = args.data_subfolder
    overhang = args.overhang
    nthreads = args.nthreads
    db_path = args.database
    date = args.date
    show_calls = args.show_calls
    overwrite = args.overwrite
    
    hrs_list = hrs.split(',')  # Split into ['V1-V2_27F-338R', 'V3-V4_337F-805R']
    nhrs = len(hrs_list)
    v_re = re.compile('^.*(V[0-9]+)-(V[0-9]+).*$')
    c_re = re.compile('^.*V[0-9]+-V[0-9]+_([0-9]+[a-zA-Z]+F)-([0-9]+[a-z][A-Z]+R).*$')
    hrs_v_list = ['-'.join([v_re.match(h).group(1),
                            c_re.match(h).group(1),
                            v_re.match(h).group(2),
                            c_re.match(h).group(2)]) for h in hrs_list]
    hrs_label = '_'.join(hrs_v_list)

    # First parse the input hypervariable regions to determine how many and
    # which regions are going to get combined.
#     python3 py/combine_renumerate.py \
# 	--tables V1-V2_27F-338R_hang22_2021-02-13_wrefseq_table.txt,V3-V4_337F-805R_hang22_2021-02-13_wrefseq_table.txt,V5-V6_805F-1185mR_hang22_2021-02-13_wrefseq_table.txt \
# 	--out-table V1-V2_V3-V4_V5-V6_hang22_2021-02-13_wrefseq_table.txt \
# 	--out-fasta V1-V2_V3-V4_V5-V6_hang22_2021-02-13_wrefseq_sequences.fasta
    
    
    # Check if all required input files exist (for each simple hypervariable
    # region).
    input_tables = []
    h_txt_date_re = re.compile(r"^.*V[0-9]-V[0-9]_[0-9]+[a-zA-Z]+-[0-9]+[a-zA-Z]+_hang[0-9]+_([0-9]{4}-[0-9]{2}-[0-9]{2})_wrefseq_table\.txt$")
    for h in hrs_list:
        h_txts = glob.glob(os.path.join(
            data_dir, h+'.*_hang'+overhang+'*_wrefseq_table.txt'))
        h_txts_dates = [h_txt_date_re.match(h_txt).group(1) for h_txt in h_txts]
        # Check these dates vs date
        if date is not None:
            if date not in h_txts_dates:
                sys.exit('\nError: user given --date '+date+' not found among files for region '+h+'.\n')
            else:
                date_index = h_txts_dates.index(date)
                h_txts_dates_add = h_txts_dates[date_index]
        else:
            max_date = max(h_txts_dates)
            date_index = h_txts_dates.index(max_date)
            h_txts_dates_add = h_txts_dates[date_index]
        input_tables.append(os.path.join(data_dir, h+'_hang'))
        
    # Start with combine_renumerate.step.
    pass
