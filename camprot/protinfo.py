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

def getNames(reviewed=True, taxid=10090, size=-1):

    '''
    Set size to -1 to return all names
    Set taxid to None or False to return all species
    Set reviewed to None to ignore reviewed status
    '''

    ebi_api_url = 'https://www.ebi.ac.uk/proteins/api/proteins?offset=0'
    
    ebi_api_url += '&size=%i' % size

    if not reviewed is None:
        ebi_api_url += '&reviewed=%s' % str(reviewed).lower()

    if taxid:
        ebi_api_url += '&taxid=%i' % taxid

    ebi_api_url = 'http://www.ebi.ac.uk/proteins/api/features?offset=0&size=-1&reviewed=true&taxid=10090&types=MOD_RES'

    print(ebi_api_url)

    r = requests.get(ebi_api_url, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    text = json.loads(r.text)

    #return text

    requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?&size=-1&taxid=10090"

    r = requests.get(requestURL, headers={ "Accept" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    responseBody = r.text

    return responseBody

