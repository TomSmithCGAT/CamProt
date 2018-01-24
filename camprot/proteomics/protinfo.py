'''
protinfo.py - Classes and functions for getting information on proteins
=======================================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Requirements:

Code
----

'''

import urllib3
import json
import requests, sys


def getProtName(uniprot_id):
    '''
    return the name for a uniprot_id
    '''
    ebi_api_url = 'https://www.ebi.ac.uk/proteins/api/proteins/'
    http = urllib3.PoolManager()
    r = requests.get(ebi_api_url + uniprot_id, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    text = json.loads(r.text)
    return text['protein']['recommendedName']['fullName']['value']
