#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:18:33 2019

@author: igor
"""

import requests, sys

def get_token (token_file):
    with open(token_file) as tf:
        for line in tf:
            x = line.strip()
    return(x)


def upload ():
    pass

ACCESS_TOKEN = get_token(token_file)

headers = {"Content-Type": "application/json"}
r = requests.post('https://zenodo.org/api/deposit/depositions',
                   params={'access_token': ACCESS_TOKEN}, json={},
                   headers=headers)

# Generate deposition ID
deposition_id = r.json()['id']


if __name__ == '__main__':

    # Filename with the upload token    
    token_file = sys.argv[1]
    upload()




