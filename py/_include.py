#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:33:28 2021

@author: igor
"""

from datetime import datetime
from sys import stdout, path
from subprocess import Popen, PIPE
import platform, os

def log(x, end=True, time_stamp=True, force_end_line=True):
    """ Log current progress. Writes x to stdout with a time_stamp if 
    set to True. Line is terminated if end_line is true. If end_line==False
    and x contains \n last character, then the \n gets removed if the
    force_end_line == True. """
    
    prefix = ''
    suffix = ''
    if time_stamp:
        date_time_obj = datetime.now()
        prefix = date_time_obj.strftime("%Y-%m-%d %H:%M:%S | ")
    if force_end_line and x[-1] == '\n':
        x = x[:-1]        
    if end:
        suffix = '\n'
    stdout.write(prefix+x+suffix)
    stdout.flush()


def get_python_path ():
    """ Generate an absolute python3 path. """
    cmd = Popen('which python3', shell=True, stdout=PIPE, stderr=PIPE)
    out, err = cmd.communicate()
    if out == b'':
        return None
    else:
        return out.decode().strip()

def get_rscript_path ():
    """ Generate an absolute Rscript path. """
    cmd = Popen('which Rscript', shell=True, stdout=PIPE, stderr=PIPE)
    out, err = cmd.communicate()
    if out == b'':
        return None
    else:
        return out.decode().strip()

def get_os ():
    """ Detect operating system. """
    os_sys = platform.system()
    if os_sys == 'Linux':
        return 'linux'
    elif os_sys == 'Darwin':
        return 'macos'
    else:
        return ''
    
# Script folder sys.path[0]

def get_blast_path ():
    """ Generate an absolute BLAST path. """
    blastn = Popen('which blastn', shell=True, stdout=PIPE, stderr=PIPE)
    out, err = blastn.communicate()
    if out == b'':
        sys_os = get_os()
        return os.path.join(os.path.dirname(path[0]), 'bin', 'blastn_'+sys_os)
    else:
        return out.decode().strip()
