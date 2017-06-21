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


def getSequence(uniprot_id, api_first=True, return_taxid=False):
    ''' get a protein sequence from uniprot by default the api will be
    used first and the flatfile output used as a backup'''

    ebi_api_url = 'http://www.ebi.ac.uk/proteins/api/features?offset=0&accession='
    flat_url = 'http://www.uniprot.org/uniprot/'

    http = urllib3.PoolManager()

    sequence = None
    method = None
    taxid = None

    if api_first:

        # try to extract sequence from ebi api
        response = http.request('GET', ebi_api_url + uniprot_id)
        text = json.loads(response.data.decode('utf8'))

        if isinstance(text, list) and len(text)>0:
            sequence = text[0]['sequence']
            method = "ebi_api"
            return sequence, method

    # otherwise fall back on ebi flatfile output
    request = http.request('GET', url="%s%s.txt" % (flat_url, uniprot_id))
    uniprot_info = request.data.decode('utf8')

    if len(uniprot_info) > 0:

        output = False
        sequence = ""

        for u_line in uniprot_info.split("\n"):
            u_line = u_line.strip().split(" ")
            if u_line[0] == 'OX' and u_line[3].startswith('NCBI_TaxID='):
                taxid = u_line[3].replace('NCBI_TaxID=', '')[:-1]

            if u_line[0] == '//':
                output = False

            if output:
                sequence += "".join(u_line)

            if u_line[0] == 'SQ' and u_line[3] == 'SEQUENCE':
                output = True

        method = "uniprot_flat"

    if return_taxid:
        return sequence, taxid, method

    else:
        return sequence, method


def getSequences(uniprot_ids):
    ''' iterate protein sequence from uniprot. By default the api will be
    used first and the flatfile output used as a backup'''

    ebi_api_url = 'http://www.ebi.ac.uk/proteins/api/features?offset=0&accession='

    http = urllib3.PoolManager()

    sequence = None
    method = None

    returned_ids = set()

    uniprot_ids = list(uniprot_ids)

    for ix in range(0, len(uniprot_ids), 100):
        # try to extract sequence from ebi api
        response = http.request(
            'GET', ebi_api_url + ",".join(uniprot_ids[ix:ix+100]))
        text = json.loads(response.data.decode('utf8'))

        # remove proteins which cause 'is not valid' error
        if isinstance(text, dict) and text['errorMessage']: 
            remove_proteins = set()
            for error_msg in text['errorMessage']:
                if error_msg.endswith(' is not valid'):
                    remove_proteins.add(error_msg.replace(' is not valid', ''))
                    new_protein_list = [x for x in uniprot_ids[ix:ix+100]
                                        if x not in remove_proteins]
                    response = http.request(
                        'GET', ebi_api_url + ",".join(new_protein_list))
                    text = json.loads(response.data.decode('utf8'))

        for protein in text:
            if 'sequence' in protein:
                accession = protein['accession']
                taxid = protein['taxid']
                seq = protein['sequence']
                returned_ids.add(accession)
                yield(accession, taxid, seq, 'ebi_api')

    # for any remaining proteins, try the flatfile approach
    remaining_proteins = set(uniprot_ids).difference(returned_ids)
    #if len(remaining_proteins) > 0:
        #print("using flatfile approach for %s proteins" % len(remaining_proteins))
    for protein in remaining_proteins:
        seq, taxid, method = getSequence(protein, api_first=False,
                                         return_taxid=True)
        yield(protein, taxid, seq, method)


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

    if method=='gluc':
        cleave_bases = ["E", "D"]

    if method=='chromotrypsin':
        cleave_bases = ["R", "Y", "W"]
    
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


def getFreeEnergy(seq):
    ''' calculate the free energy for a linear polypeptide sequence
    going from aqeuous to organic. Uses data from White, Stephen
    (2006-06-29). "Experimentally Determined Hydrophobicity Scales" to
    determine the free energy as a protein goes from water to an
    organic compound (n-octanol)

    http://blanco.biomol.uci.edu/hydrophobicity_scales.html'''

    
    aa2freeEnergy = {"I": -1.12,
                     "L": -1.25,
                     "F": -1.71,
                     "V": -0.46,
                     "M": -0.67,
                     "P": 0.14,
                     "W": -2.09,
                     "T": 0.25,
                     "Q": 0.77,
                     "C": -0.02,
                     "Y": -0.71,
                     "A": 0.5,
                     "S": 0.46,
                     "N": 0.85,
                     "R": 1.81,
                     "G": 1.15,
                     "H": 2.33,
                     "E": 3.63,
                     "K": 2.8,
                     "D": 3.64,
                     "U": 0, # Seloncysteine was not considered in linked table
                     "X": 0}  # X = missing protein

    free_energy = 0

    for amino_acid in seq:
        free_energy += aa2freeEnergy[amino_acid]

    return free_energy


def getFractionHydrophobic(seq):
    ''' calculate the fraction of amino acids which are hydrophobic'''

    hydrophobicAA = set(("A", "V", "I", "L", "M", "F", "Y", "W"))

    count_hydrophobic = 0

    for amino_acid in seq:
        if amino_acid in hydrophobicAA:
            count_hydrophobic += 1

    return float(count_hydrophobic)/len([aa for aa in seq if aa!= "X"])

def getFractionHydrophillic(seq):
    ''' calculate the fraction of amino acids which are hydrophillic'''

    hydrophillicAA = set(("R", "N", "D", "Q", "E", "H", "K", "S", "T"))

    count_hydrophillic = 0

    for amino_acid in seq:
        if amino_acid in hydrophillicAA:
            count_hydrophillic += 1

    return float(count_hydrophillic)/len([aa for aa in seq if aa!= "X"])
