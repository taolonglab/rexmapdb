#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 20 13:42:31 2018

@author: igor
"""
import sys, os
from count import get_os
from subprocess import Popen, PIPE


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
    print('Generating HiMAP database')
    makedb_path = get_makedb_path()
    makedb = Popen([makedb_path, '-dbtype', 'nucl', '-in', in_fa,
                             '-out', out_prefix], stdout=sys.stdout, stderr=sys.stderr)
    out, err = makedb.communicate()

if __name__ == '__main__':
    in_fa = sys.argv[1]
    out_prefix = sys.argv[2]
    main(in_fa, out_prefix)