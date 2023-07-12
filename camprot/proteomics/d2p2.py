'''
d2p2.py - Classes and functions for extracting data from D2P2 API
=================================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python


Requirements:

Code
----

'''

import urllib3
import requests
import json
import collections

import numpy as np

from time import gmtime, strftime

d2p2_url = 'http://d2p2.pro/api/seqid'


class d2p2(object):
    '''
    holds data for d2p2 entry
    '''

    def __init__(self, d2p2_text):

        self.name = d2p2_text[0]
        annotations_dict = d2p2_text[2]
        self.disranges = annotations_dict['disorder']['disranges']
        self.consensus = annotations_dict['disorder']['consensus']
        self.structure = annotations_dict['structure']

        self.length = len(self.consensus)
        self.idrs = []

    def rebuildConsensus(self,
                         tools_whitelist=None,
                         tools_blacklist=None):
        ''' Use whitelist and blacklist to rebuild the consensus
        sequence from the disranges
        Both must be list, set or tuple'''

        if not tools_whitelist and not tools_blacklist:
            return

        accepted_tools = set([x[0] for x in self.disranges])

        if tools_blacklist:
            accepted_tools = accepted_tools.difference(set(tools_blacklist))
        if tools_whitelist:
            accepted_tools = accepted_tools.intersection(set(tools_whitelist))

        new_consensus = np.zeros(len(self.consensus), dtype=np.int8)

        for disrange in self.disranges:
            tool, start, end = disrange
            if tool in accepted_tools:
                new_consensus[int(start)-1: int(end)] += 1

        self.consensus = list(new_consensus)

    def setIDRs(self, consensus_req=3, min_block_size=20):
        ''' Identify the IDRs
        consensus_req = minimal number of agreements between tools
        min_block_size = minimal size of IDR
        Note: regions are zero-indexed, open. E.g pythonic'''

        blocks = []
        start = 0

        for position, cons in enumerate(self.consensus):
            # identify where consensus is too low
            if cons < consensus_req:
                if position - start >= min_block_size:
                    blocks.append((start, position))
                start = position + 1

        if position - start >= min_block_size:
            blocks.append((start, position + 1))

        self.idrs = blocks


def getChunks(uniprot_ids, chunk_size=250):
    '''
    return chunks of D2P2 data
    This avoids excessive 'get' requests
    Note the chunk size is limited to 250 by default to prevent excessively long urls
    which would throw an error
    '''
    http = urllib3.PoolManager()

    proteins_done = 0
    for chunk_start in range(0, len(uniprot_ids), chunk_size):

        chunk_ids = list(uniprot_ids)[chunk_start: (chunk_start + chunk_size)]

        url = '%s/["%s"]' % (d2p2_url, '","'.join(chunk_ids))
        
        #request = http.request('GET', url=url, verify=False)
        #response = json.loads(request.data.decode('utf8'))
        request = requests.get(url, verify=False)
        response = json.loads(request.text)

        yield response


def iterator(uniprot_ids, chunk_size=250):

    for chunk in getChunks(uniprot_ids, chunk_size=chunk_size):
        
        for key in chunk:

            d2p2_text = chunk[key]

            if len(d2p2_text) > 0:
                yield d2p2(chunk[key][0])
