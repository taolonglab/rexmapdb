#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Download the 16S ribosomal RNA sequences from NCBI nucleotide
database search for bacteria and archaea.


Created on Tue Feb 27 09:15:26 2018

@author: igor
"""
import sys, subprocess, re

eutil_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
search_term = '16s ribosomal RNA[Title] NOT uncultured[Title] AND ' + \
              '(bacteria[Filter] OR archaea[Filter]) AND ' + \
              '(1000[SLEN] : 2000[SLEN]) AND refseq[filter]'


def detect_curl_path ():
    """ Automatically detect curl path
    usually curl_path = '/usr/bin/curl' """
    p = subprocess.Popen('which curl', shell=True,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out.decode().strip()


def main(fasta_out_path, curl_path=None, verbose=True, retmax=10000,
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
        curl_path = detect_curl_path()
        if verbose:
            print('OK.')
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


if __name__ == '__main__':
    
    fasta_out_path = sys.argv[1]
    main(fasta_out_path)
            
#    verbose = True
#    retmax = 10000
#    fasta_out_path = 'test.fasta'
#    if verbose:
#        print('- detecting curl...', end='')
#    curl_path = detect_curl_path()
#    if verbose:
#        print('OK.')
#        print('- using '+curl_path+' for web queries')
#        print('- sending search query...', end='')
#    search_url = eutil_url + 'esearch.fcgi?db=nucleotide&term=' + \
#        search_term.replace(' ', '+') + '&usehistory=y'
#
#    # Retrieve first 20 results. These contain total number of results,
#    # query_key and WebEnv which are the three things we will need to get
#    # all the results in a for loop.
#    p = subprocess.Popen([curl_path, '-g', search_url],
#                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#    out, err = p.communicate()
#    out = out.decode()
#    if verbose:
#        print('OK.')
#    total = int(re.findall(r'<eSearchResult><Count>([0-9]+)</Count>', out)[0])
#    query_key = re.findall(r'<QueryKey>([^<]+)</QueryKey><WebEnv>', out)[0]
#    web_env = re.findall(r'<WebEnv>([^<]+)</WebEnv>', out)[0]
#    if verbose:
#        print('- total sequences: '+str(total))
#        print('- query key: '+query_key)
#        print('- WebEnv: '+web_env)
#    fasta_out_string = ''
#
#    for retstart in range(0, total, retmax):
#        if verbose:
#            print('\r- % 6d out of % 6d sequences processed' % (retstart, total),
#                  end='')
#        fetch_url = eutil_url + 'efetch.fcgi?db=nucleotide&query_key=' + \
#            query_key + '&WebEnv=' + web_env + '&retstart=' + str(retstart) + \
#            '&retmax=' + str(retmax) + '&rettype=fasta'
#        p = subprocess.Popen([curl_path, '-g', fetch_url],
#                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#        out, err = p.communicate()
#        out = out.decode()
#        fasta_out_string += out
#    if verbose:
#        print('\r- % 5d out of % 5d sequences processed' % (total, total))
#    
#    # Write output
#    if verbose:
#        print('- writing '+fasta_out_path+' ...', end='')
#    with open(fasta_out_path, 'w') as f_out:
#        f_out.write(fasta_out_string)
#    if verbose:
#        print('OK.')