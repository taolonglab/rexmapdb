#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download the 16S ribosomal RNA sequences from NCBI nucleotide
database search for bacteria and archaea.


Created on Tue Feb 27 09:15:26 2018

@author: igor
"""
import sys, subprocess, re, os
import argparse
from count import get_os

eutil_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
search_term = '16s ribosomal RNA[Title] NOT uncultured[Title] AND ' + \
              '(bacteria[Filter] OR archaea[Filter]) AND ' + \
              '(1000[SLEN] : 2000[SLEN]) AND refseq[filter]'


def detect_curl_path (debug=False):
    """ Automatically detect curl path
    usually curl_path = '/usr/bin/curl' """
    p = subprocess.Popen('which curl', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    # Check if curl is not found, then use budled binary
    out = out.decode().strip()
    if not out or debug:
        sys_os = get_os()
        out = os.path.join(
                os.path.dirname(os.path.dirname(sys.argv[0])), 
                'bin', 'curl_'+sys_os
                )
    return out


def main(fasta_out_path, curl_path=None, verbose=True, retmax=10000, debug=False,
         eutil_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/',
         search_term = '16s ribosomal RNA[Title] NOT uncultured[Title] AND ' + \
              '(bacteria[Filter] OR archaea[Filter]) AND ' + \
              '(1000[SLEN] : 2000[SLEN]) AND refseq[filter]'
):
    """ Main function that calls curl to download sequences using
    NCBI eUtils on the Web. """
    if verbose:
        print('Download 16S rRNA sequences from NCBI Web Search')
    if curl_path is None:
        if verbose:
            print('- detecting curl...', end='')
        curl_path = detect_curl_path(debug=debug)
        if verbose:
            print('OK.')
    if debug:
        print('- Debug mode: ON')
    if verbose:
        print('- using '+curl_path+' for web queries')
        print('- sending search query...', end='')
    
    # Construct full search URL query (c/p to browser in case of bugs)
    search_url = eutil_url + 'esearch.fcgi?db=nucleotide&term=' + \
        search_term.replace(' ', '+') + '&usehistory=y'

    # Retrieve first 20 results. These contain total number of results,
    # query_key and WebEnv which are the three things we will need to get
    # all the results in a for loop.
    p = subprocess.Popen([curl_path, '-g', search_url],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.decode()
    if verbose:
        print('OK.')
    
    # Extract counts, query key and web env to retrieve all results
    total = int(re.findall(r'<eSearchResult><Count>([0-9]+)</Count>', out)[0])
    query_key = re.findall(r'<QueryKey>([^<]+)</QueryKey><WebEnv>', out)[0]
    web_env = re.findall(r'<WebEnv>([^<]+)</WebEnv>', out)[0]
    if verbose:
        print('- total sequences: '+str(total))
        print('- query key: '+query_key)
        print('- WebEnv: '+web_env)
    fasta_out_string = ''
    
    # NCBI allows only up to 10k results at a time, so loop over them
    for retstart in range(0, total, retmax):
        if verbose:
            print('\r- % 6d out of % 6d sequences processed' % (retstart, total),
                  end='')
        fetch_url = eutil_url + 'efetch.fcgi?db=nucleotide&query_key=' + \
            query_key + '&WebEnv=' + web_env + '&retstart=' + str(retstart) + \
            '&retmax=' + str(retmax) + '&rettype=fasta'
        p = subprocess.Popen([curl_path, '-g', fetch_url],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        out = out.decode()
        fasta_out_string += out
    if verbose:
        print('\r- % 5d out of % 5d sequences processed' % (total, total))
    
    # Write output to FASTA
    if verbose:
        print('- writing '+fasta_out_path+' ...', end='')
    with open(fasta_out_path, 'w') as f_out:
        f_out.write(fasta_out_string)
    if verbose:
        print('OK.')

def parse_input():
  parser = argparse.ArgumentParser(
    description='Download NCBI 16S RefSeq sequences.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--output-fasta', required=True, 
                      help='Output FASTA file for 16S sequences.')
  parser.add_argument('--debug', required=False, default=False, type=bool,
                      help='Debug mode.')
  parser.add_argument('--verbose', required=False, default=True, type=bool,
                      help='Verbose print output.')
  args = parser.parse_args()
  return args


if __name__ == '__main__':
    
    # fasta_out_path = sys.argv[1]
    
    # if len(sys.argv) >= 3:          # Argument 7 is overhang
    #     debug = True
    # else:
    #     debug = False
    
    args = parse_input()
    fasta_out_path = args.output_fasta
    debug = args.debug
    verbose = args.verbose
    
    main(fasta_out_path, debug=debug, verbose=verbose)
