#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 13:42:31 2018

@author: igor
"""
import sys, os
import argparse
from subprocess import Popen, PIPE
from _include import log, get_os


def get_makedb_path ():
    """ Detect path for makeblastdb binary. """
    # Get the installed version in path if exists
    makedb = Popen('which makeblastdb', shell=True, stdout=PIPE, stderr=PIPE)
    out, err = makedb.communicate()
    if out == b'':
        sys_os = get_os()
        return os.path.join(os.path.dirname(sys.path[0]), 'bin', 'makeblastdb_'+sys_os)
    else:
        return out.decode().strip()

def main(in_fa, out_prefix):
    log('Generating RExMapDB BLAST database')
    makedb_path = get_makedb_path()
    makedb = Popen([makedb_path, '-dbtype', 'nucl', '-in', in_fa,
                             '-out', out_prefix], stdout=sys.stdout, stderr=sys.stderr)
    out, err = makedb.communicate()

def parse_input():
  parser = argparse.ArgumentParser(
    description='Generate BLAST database from RExMapDB FASTA file.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--input-fasta', required=True, 
                      help='Input FASTA file (unique_variants_R.fasta).')
  parser.add_argument('--out-db-prefix', required=True,
                      help='Output prefix for BLAST database files.')
  args = parser.parse_args()
  return args

if __name__ == '__main__':
    
    args = parse_input()
    in_fa = args.input_fasta
    out_prefix = args.out_db_prefix
    # in_fa = sys.argv[1]
    # out_prefix = sys.argv[2]
    main(in_fa, out_prefix)
