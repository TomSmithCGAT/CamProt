'''annotate_rnp - add information to the output from RNPxl
=======================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python RNP Proteomics

Purpose
-------

This script takes the xlsx output from RNPxl and annotates the table with useful information for downstream analyses.

The following columns are added:

 - master_protein(s): The master protein(s) for the peptide. See below
   for how this is derived
 - master_uniprot_id(s): The uniprot id(s) for the master protein(s)
 - protein_description(s): Description(s) for the master protein(s)
 - protein_length(s): The length(s) of the master protein(s)
 - position_in_peptide: The most-likely position of the cross-link in the peptide
 - position_in_protein(s): The most-likely position of the cross-link in the protein
 - window_13-17: 13, 15 & 17 amino acid windows around the cross-link
   position. See below for further notes on these windows
 - crap_protein: Is the protein in the cRAP database of common
   proteomics proteins, e.g keratin

If a log file is requested (--log), basic statistics are collected and
written to the log file

Fasta description format
------------------------
The source of the protein (SwissProt or TrEMBL) is derived from the
protein fasta description, with SwissProt proteins starting 'sp' and
TrEMBL 'tr'. Furthermore, the description column is derived from the
fasta description too. For this reason the fasta databases must be
correctly formatted as in the examples below. This is the standard
format for fasta files from uniprot.

format:
Three-level identifier followed by protein description:
>[sp|tr]|[Uniprot id]|[Protein name] [Description]

examples:
>sp|P27361|MK03_HUMAN Mitogen-activated protein kinase 3 OS=Homo sapiens GN=MAPK3 PE=1 SV=4
>tr|F8W1T5|F8W1T5_HUMAN GTPase RhebL1 (Fragment) OS=Homo sapiens GN=RHEBL1 PE=4 SV=1


Deriving the master proteins
----------------------------

Matching peptides to their source proteins (protein inference) is a
common task in proteomics and there are many possible
approaches. Ultimately, the aim is usually to identify the most likely
source protein since taking all possible sources makes downstream
analyses very complex. Here we use the parsimonious approach to
identify a minimal set of proteins which explains all peptides
observed. In essense, the approach is as follows:
- start with list of all peptides
- sort proteins by the number of peptides observed
- take the protein(s) with the most peptides and remove these from the peptides list
- continue through the sorted proteins, removing peptides, until the
  peptides list is empty

Additionally, we prioritise matches to SwissProt proteins over TrEMBL
proteins. SwissProt proteins have been manually curated and should not
contain any redundant proteins, truncated sequences etc. On the other
hand, the set of TrEMBL proteins will ceratainly contain proteins
which are redundant with respect to the SwissProt proteins as well as
truncated proteins. It is useful to include the TrEMBL proteins to
catch peptides which are from a protein or isoform which has not been
curated into SwissProt yet. However, where a SwissProt match is found,
any TrEMBL match can safely be ignored. Here, for all peptides with
matched to both SwissProt and TrEMBL proteins, we remove all the
TrEMBL matches.
    
In some instances, it is not possible to assign a single protein to a
peptide. In these cases, the proteins names, uniprot ids, descriptions
and lengths are ';' separated in the outfile.


13, 15 & 17 amino acid windows
------------------------------
For motif analysis, 13, 15 & 17 amino acid windows around the most
likely cross-linked protein are provided. Where the cross-link
posistion is too close to the start or end of the protein for the
window, e.g cross link is at position 6 - not possible to extract a
window from -1 -> 12, the value "protein_too_short_for_window" is
given. Where there is more than one master protein, the windows is
provided only where the amino acid sequence is the same for all master
proteins. Where the sequences diverge, the value
"different_sequences_for_the_proteins" is given. Care must be taken
with any motif analysis since the cross-link site is hard to identify
so the windows may not be well centered. Furthermore, since MS is bias
towards particular peptides, care must be taken to provide a suitable
background set of windows. For example, random windows from the fasta
could simply yield motifs which are enriched in MS analyses


Usage
-----
By default, the outfile will be created in the same directory with the
suffix annotated.xlsx. You can change the outfile name by specifying
the option --outfile

python annotate_rnp.py --infile=RNP.xlsx --fasta-db=Human_201701.fasta
--fasta-crap-db=RAP_FullIdentifiers.fasta --outfile=RNP_annotated.xlsx
--logfile=RNP_annotation.log

Command line options
--------------------

'''

import argparse
import collections 
import copy
import os
import re
import sys

import pandas as pd

####################### ------------------------- ##############################
# code below 'borrowed' from CGAT.FastaIterator (trying to limit dependencies)
class FastaRecord:
    """a :term:`fasta` record.

    Attributes
    ----------
    title: string
       the title of the sequence

    sequence: string
       the sequence

    fold : int
       the number of bases per line when writing out
    """

    def __init__(self, title, sequence, fold=False):

        self.title = title
        self.sequence = sequence
        self.fold = fold

    def __str__(self):
        ''' str method for writing out'''

        if self.fold:
            seq = [self.sequence[i:i + self.fold]
                   for i in range(0, len(self.sequence), self.fold)]
        else:
            seq = (self.sequence,)

        return ">%s\n%s" % (self.title, "\n".join(seq))


class FastaIterator:
    '''a iterator of :term:`fasta` formatted files.

    Yields
    ------
    FastaRecord

    '''

    def __init__(self, f, *args, **kwargs):
        self.iterator = iterate(f)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self.iterator)

    def next(self):
        return next(self.iterator)


def iterate(infile, comment="#", fold=False):
    '''iterate over fasta data in infile

    Lines before the first fasta record are
    ignored (starting with ``>``) as well as
    lines starting with the comment character.

    Parameters
    ----------
    infile : File
        the input file
    comment : char
        comment character
    fold : int
        the number of bases before line split when writing out

    Yields
    ------
    FastaRecord
    '''

    h = infile.readline()[:-1]

    if not h:
        raise StopIteration

    # skip everything until first fasta entry starts
    while h[0] != ">":
        h = infile.readline()[:-1]
        if not h:
            raise StopIteration
        continue

    h = h[1:]
    seq = []

    for line in infile:

        if line.startswith(comment):
            continue

        if line.startswith('>'):
            yield FastaRecord(h, ''.join(seq), fold)

            h = line[1:-1]
            seq = []
            continue

        seq.append(line[:-1])

    yield FastaRecord(h, ''.join(seq), fold)
####################### ------------------------- ##############################

def getMotifWindow(positions, proteins, length):
    ''' Extract amino acid sequences of (length) n from (proteins),
    centered at (positions)'''

    assert length % 2 == 1, 'motif window must be an odd length'
    assert len(positions) == len(proteins), "must have as many positions as proteins"

    windows = set()

    for position, protein_seq in zip(positions, proteins):
        buffer = ((length - 1) / 2)
        windows.add(protein_seq[int(position-buffer): int(position+buffer) + 1])

    if min([len(x) for x in windows]) != length:
        return "protein_too_short_for_window"
    if len(windows) == 1:
        return list(windows)[0]
    else:
        return "different_sequences_for_the_proteins"

def writeSectionHeader(logfile, section_header):
    #underliner = "".join(("-",)*len(section_header))
    section_blocker = ("======================================="
                       "=======================================")
    underliner1 = ("----------------------------------------"
                  "----------------------------------------")
    logfile.write("\n%s\n%s\n" % (section_blocker, section_header))
    logfile.write("%s\n" % underliner1)
    return section_blocker


def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        argv, usage=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--infile', dest="infile", required=True,
                        help="")

    parser.add_argument('-f', '--fasta-db', dest="fasta_db", required=True,
                        help="")

    parser.add_argument('-fc', '--fasta-crap-db', dest="fasta_crap_db",
                        required=True, help="")

    parser.add_argument('-o', '--outfile', dest="outfile", default=None,
                        help="")

    parser.add_argument('-l', '--logfile', dest="logfile", default=os.devnull,
                        help="")

    args = vars(parser.parse_args())

    if args['outfile'] is None:
        args['outfile'] = args['infile'].replace(".xlsx", "_annotated.tsv")

        
    logfile = open(args['logfile'], 'w')
    logfile.write("Logfile for annotate_rnp.py\n\n")
    section_blocker = writeSectionHeader(logfile, "Script arguments:")
    for key, value in args.items():
        logfile.write("%s: %s\n" % (key, value))
    logfile.write("%s\n\n" % section_blocker)

    # read the data into a dataframe
    df = pd.read_excel(args['infile'])

    # add some basic annotations
    df['tr_only'] = [x.count("sp|") == 0 for x in df['Proteins']]
    df['matches'] = [len(x.split(",")) for x in df['Proteins']]

    #(1) Get the mappings between peptide and proteins
    pep2pro = collections.defaultdict(lambda: collections.defaultdict(set))
    pep2allpro = collections.defaultdict(set)
    pro2pep = collections.defaultdict(set)
    top_level_proteins = set()
    initial_proteins = set()

    # (1.1) extract the initial mappings between proteins and peptides
    for row_ix, row_values in df[['Proteins', 'Peptide']].iterrows():

        proteins = row_values['Proteins'].split(",")
        peptide = row_values['Peptide']

        if peptide in pep2pro:
            assert pep2allpro[peptide] == proteins, (
                "The same peptide is observed more than once with different proteins!")
        pep2allpro[peptide] = proteins
        for protein in proteins:
            initial_proteins.add(protein)
            pro2pep[protein].add(peptide)

            if protein.split("|")[0] == "sp":
                protein_level = 1
                top_level_proteins.add(protein)
            elif protein.split("|")[0] == "tr":
                protein_level = 2
            else:
                raise ValueError("Protein does not appear to be either"
                                 "SwissProt(sp) or TrEMBL(tr)")

            pep2pro[peptide][protein_level].add(protein)

    section_blocker = writeSectionHeader(logfile, "Initial file stats")
    logfile.write("# initial peptides: %i\n" % len(pep2pro))
    logfile.write("# initial proteins: %i\n" % len(pro2pep))
    logfile.write("# initial SwissProt proteins: %i\n" % len(top_level_proteins))
    logfile.write("# initial TrEMBL proteins: %i\n" % (
        len(pro2pep)-len(top_level_proteins)))
    logfile.write("%s\n\n" % section_blocker)
            
    # (1.2) find the peptides with only TrEMBL protein matches and
    # 'upgrade' these TrEMBL proteins to being equivalent to SwissProt
    tr_only_peptides = set([x for x in pep2pro.keys() if len(pep2pro[x][1])==0])

    set_upgraded = set()
    for peptide in tr_only_peptides:
        upgraded = pep2pro[peptide][2]
        set_upgraded.update(upgraded)
        top_level_proteins.update(upgraded)
        pep2pro[peptide][2] = pep2pro[peptide][2].difference(set(upgraded))
        pep2pro[peptide][1] = pep2pro[peptide][1].union(set(upgraded))

    section_blocker = writeSectionHeader(
        logfile, "Deciding which TrEMBL proteins to retain:")
    logfile.write("# peptides with only TrEMBL matches: %i\n" % (
        len(tr_only_peptides)))
    logfile.write("# TrEMBL proteins retained as no SwissProt matches for "
                  "peptide: %i\n" % (len(set_upgraded)))
    logfile.write("%s\n\n" % section_blocker)

    # (1.3) Use a parsimonious approach to identifty the minimum number
    # of proteins required to cover all the peptides:
    # Start from the protein(s) with the most peptides and mark these as covered.
    # Continue with remaining proteins in order of peptides per protein
    # until all peptides are covered
    retained_proteins = []
    peptides = copy.deepcopy(set(pep2pro.keys()))
    peptide_counts = {}

    tmppro2pep = copy.deepcopy(pro2pep)
    new_top_level_proteins = copy.deepcopy(top_level_proteins)
    new_pep2pro = collections.defaultdict(set)

    peptide_count = max(map(len, tmppro2pep.values()))

    section_blocker = writeSectionHeader(
        logfile, ("Parsimonious method to identify minimal set of proteins"
                  " to account for all peptides"))

    while True:
        # (1.3.1) If all peptides covered or the maximum peptides per
        # protein = 0, break.
        if len(peptides) == 0 or peptide_count == 0:
            logfile.write("All peptides are now accounted for\n")
            break

        peptide_count -= 1 

        top_proteins = set()
        top_score = 0
        #(1.3.2) Find the proteins with the highest number of peptide matches
        for protein in new_top_level_proteins:
            if len(tmppro2pep[protein]) == top_score:
                top_proteins.add(protein)
            elif len(tmppro2pep[protein]) > top_score:
                top_score = len(tmppro2pep[protein])
                top_proteins = set((protein,))

        logfile.write("%i remaining protein(s) with %i peptides\n" % (
            len(top_proteins), top_score))
        # (1.3.3) Remove the top proteins and the associated peptides
        for top_protein in top_proteins:

            new_top_level_proteins.remove(top_protein)

            retained_proteins.append(top_protein)

            for peptide in pro2pep[top_protein]:
                new_pep2pro[peptide].add(top_protein)
                if peptide in peptides:
                    peptides.remove(peptide)
                for protein in pep2pro[peptide][1]:
                    if protein == top_protein:
                        continue
                    if peptide in tmppro2pep[protein]:
                        tmppro2pep[protein].remove(peptide)

    logfile.write("\n%i proteins retained\n" % len(retained_proteins))
    #logfile.write("\n".join([",".join(map(str, (x, len(tmppro2pep[x]), len(pro2pep[x]))))
    #                         for x in new_top_level_proteins]))
    logfile.write("%i SwissProt proteins retained\n" % len(
        [x for x in retained_proteins if x.split("|")[0] == 'sp']))
    logfile.write("%i TrEMBL proteins retained\n" % len(
        [x for x in retained_proteins if x.split("|")[0] == 'tr']))
    logfile.write("\nNote: If not all SwissProt proteins were retained, this means\n"
                  "these proteins only included peptides which were observed\n"
                  "in other proteins which had a greater number of peptides\n")
    logfile.write("%s\n\n" % section_blocker)

    section_blocker = writeSectionHeader(logfile, "proteins per peptide:")
    counts = collections.Counter([len(x) for x in new_pep2pro.values()])
    print(counts.most_common())
    sum_counts = sum(counts.values())
    for k, v in counts.items():
        logfile.write("%i peptide(s) (%.2f %%) have %i master protein(s)\n" % (
            v, (100 * v)/sum_counts, k))
    logfile.write("%s\n\n" % section_blocker)

    # Check all the peptides are covered
    assert set(pep2pro.keys()).difference(set(new_pep2pro.keys())) == set()


    # add the top protein and uniprot id annotations
    df['master_protein(s)'] = [";".join(new_pep2pro[protein]) for protein in df['Peptide']]
    df['master_uniprot_id(s)'] = [";".join([protein_id.split("|")[1] for protein_id in new_pep2pro[protein]])
                               for protein in df['Peptide']]

    # (1.4) Build dictionaries to map from protein id to protein
    # sequence and description using the fasta database
    crap_proteins = set()
    protein2description = {
        entry.title.split(" ")[0]: " ".join(entry.title.split(" ")[1:])
        for entry in FastaIterator(open(args['fasta_db']))}
    protein2seq = {
        entry.title.split(" ")[0]:entry.sequence
        for entry in FastaIterator(open(args['fasta_db']))}
    for entry in FastaIterator(open(args['fasta_crap_db'])):
        protein2seq[entry.title.split(" ")[0]] = entry.sequence
        protein2description[entry.title.split(" ")[0]] = entry.title.split(" ")[0]
        crap_proteins.add(entry.title.split(" ")[0])

    # (1.5) derive further annotations
    protein_lengths = []
    protein_descriptions = []
    crap_protein = []

    position_in_peptide = []
    position_in_protein = []

    motif_13 = []
    motif_15 = []
    motif_17 = []

    for ix, row in df.iterrows():
        peptide = row['Best localization(s)']
        proteins = row['master_protein(s)'].split(";")
        protein_lengths.append(";".join(map(str, [len(protein2seq[x]) for x in proteins])))
        protein_descriptions.append(";".join([protein2description[x] for x in proteins]))

        # (1.5.1) does peptide match a cRAP protein?
        crap = 0
        for protein in proteins:
            if protein in crap_proteins:
                crap = 1
                break
        crap_protein.append(crap)

        # (1.5.1) Find crosslink position in protein and extract 13,
        # 15 & 17 aa windows around the crosslink position
        if row['Best localization(s)']!='nan' and row['Best loc score']>0:

            peptide_positions = [re.search(peptide.upper(), protein2seq[x]).start() for
                                 x in proteins]
            crosslink_position = None
            for ix, aa in enumerate(peptide):
                if aa == aa.lower():
                    crosslink_position = ix
            assert crosslink_position is not None, (
                "no crosslinked position was found(!): %s" % peptide)

            position_in_peptide.append(crosslink_position + 1)

            protein_positions = [crosslink_position + x for x in peptide_positions]
            position_in_protein.append(
                ";".join(map(str, [x + 1 for x in protein_positions])))

            motif_13.append(
                getMotifWindow(protein_positions, [protein2seq[x] for x in proteins], 13))
            motif_15.append(
                getMotifWindow(protein_positions, [protein2seq[x] for x in proteins], 15))
            motif_17.append(
                getMotifWindow(protein_positions, [protein2seq[x] for x in proteins], 17))

        else:
            position_in_peptide.append("no_crosslink")
            position_in_protein.append("no_crosslink")
            motif_13.append("no_crosslink")
            motif_15.append("no_crosslink")
            motif_17.append("no_crosslink")


    df['protein_length(s)'] = protein_lengths
    df['protein_description(s)'] = protein_descriptions
    df['crap_protein'] = crap_protein

    df['position_in_peptide'] = position_in_peptide
    df['position_in_protein(s)'] = position_in_protein

    df['window_13'] = motif_13
    df['window_15'] = motif_15
    df['window_17'] = motif_17

    new_column_order = [
        "Best localization(s)",
        "RNA",
        "master_protein(s)",
        "master_uniprot_id(s)",
        'protein_description(s)',
        'protein_length(s)',
        'position_in_peptide',
        'position_in_protein(s)', 
        'window_13', 'window_15', 'window_17',
        'crap_protein',
        "Peptide",
        "Proteins"]

    new_column_order.extend([x for x in df.columns if x not in new_column_order])
    final_df = df[new_column_order]

    final_df.to_csv(args['outfile'], index=False, sep="\t")
    os.chmod(args['outfile'], 0o666)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
