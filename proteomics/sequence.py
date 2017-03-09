'''
sequence.py - Classes and functions for working with protein sequences
======================================================================

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


def getSequence(uniprot_id, api_first=True):
    ''' get a protein sequence from uniprot by default the api will be
    used first and the flatfile output used as a backup'''

    ebi_api_url = 'http://www.ebi.ac.uk/proteins/api/features?offset=0&accession='
    flat_url = 'http://www.uniprot.org/uniprot/'

    http = urllib3.PoolManager()

    sequence = None
    method = None

    if api_first:

        # try to extract sequence from ebi api
        response = http.request('GET', ebi_api_url + uniprot_id)
        text = json.loads(response.data.decode('utf8'))

        if isinstance(text, list) and len(text)>0:
            sequence = text[0]['sequence']
            method = "ebi_api"

    # otherwise fall back on ebi flatfile output
    if sequence is None:
        request = http.request('GET', url="%s%s.txt" % (flat_url, uniprot_id))
        uniprot_info = request.data.decode('utf8')

        if len(uniprot_info) > 0:

            output = False
            sequence = ""

            for u_line in uniprot_info.split("\n"):
                u_line = u_line.strip().split(" ")
                if u_line[0] == '//':
                    output = False
                if output:
                    sequence += "".join(u_line)
                if u_line[0] == 'SQ' and u_line[3] == 'SEQUENCE':
                    output = True

            method = "uniprot_flat"

        else:
            sequence = None

    return sequence, method


def iteratePeptides(sequence="",
                    method='trypsin',
                    missed_cleavages=2,
                    min_length=4,
                    max_length=40):
    '''
    Naive in silico digestion

    Need to add other methods
    '''

    if method=='trypsin':
        cleave_bases = ["K", "R"]
    
    peptide_ids = set((0,))
    # seq, start, end, missed, allowed_missed
    peptides = [["", 0, 0, 0, 0]]
    for missed_cleavage in range(1, missed_cleavages + 1):
        peptide_ids.add(missed_cleavage)
        peptides.append(["", 0, 0, 0, missed_cleavage])

    missed = 0
    peptide_id = missed_cleavages
    start = 0

    for position, base in enumerate(sequence):
        for peptide in peptide_ids:
            peptides[peptide][0] += base
            peptides[peptide][2] += 1

        if base in cleave_bases:

            for peptide in list(peptide_ids):
                peptides[peptide][3] += 1
                seq, start, end, missed, allowed_missed, = peptides[peptide]

                if missed > allowed_missed:
                    # remove peptide id 
                    peptide_ids.remove(peptide)

                    if len(seq) > min_length and len(seq) < max_length:
                        pep_seq, start, end, count, missed = peptides[peptide]
                        yield pep_seq, start, end, missed

            peptide_ids.add(peptide_id + 1)
            peptides.append(["", position+1, position+1, 0, 0])
            for missed_cleavage in range(1, missed_cleavages + 1):
                peptide_ids.add(peptide_id + 1 + missed_cleavage)
                peptides.append(["", position+1, position+1, 0, missed_cleavage])

            peptide_id += (missed_cleavages + 1)
