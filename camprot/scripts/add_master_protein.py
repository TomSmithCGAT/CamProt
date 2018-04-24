'''annotate_rnp -  assign a unique master to each peptide across
all samples using a maximum parsimony approachadd information to the output from RNPxl
=======================================================

:Author: Tom Smith, Manasa Ramakrishna
:Release: $Id$
:Date: |today|
:Tags: Python RNP Proteomics

Purpose
-------

This script takes the xlsx output from a set of input files
(*.txt/*.xlsx) and annotates the table with unique protein information
for downstream analyses.

The following columns are added:

 - master_protein: The master protein(s) for the peptide. See below
   for how this is derived
 - master_uniprot_id: The uniprot id(s) for the master protein(s)
 - protein_description: Description(s) for the master protein(s)
 - protein_length: The length(s) of the master protein(s)
 - crap_protein: Is the protein in the cRAP database of common
   proteomics proteins, e.g keratin
 - crap_associated_protein: does the protein share peptides with a
   protein in the cRAP database of common proteomics proteins

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


Deriving master proteins
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
contain any redundant proteins, truncated sequences, miss-annotations
etc. On the other hand, the set of TrEMBL proteins will ceratainly
contain proteins which are redundant with respect to the SwissProt
proteins as well as truncated and just plain wrong(!) proteins. It is
useful to include the TrEMBL proteins to catch peptides which are from
a protein or isoform which has not been curated into SwissProt
yet. However, where a SwissProt match is found, we believe it is
preferable to ignore any TrEMBL match. Here, for all peptides with
matched to both SwissProt and TrEMBL proteins, we remove all the
TrEMBL matches.
        
In some instances, it is not possible to assign a single protein to a
peptide. In these cases, the proteins names, uniprot ids, descriptions
and lengths are ';' separated in the outfile.

In addition to the conditions above, In some cases we are looking for
master proteins that are consistent across a set of samples. This is
to ensure that for a given set of peptides, the same master protein is
assigned to all samples.


Usage
-----
By default, the outfile will be created in the same directory with the
suffix annotated.xlsx. You can change the outfile name by specifying
the option --outfile

python add_master_protein.py --infile=RNP.xlsx --fasta-db=Human_201701.fasta
--fasta-crap-db=cRAP_FullIdentifiers.fasta --outfile=master_prot_annotated.txt
--logfile=master_prot_annot.log

Command line options
--------------------
'''

import argparse
import collections 
import copy
import os
import re
import sys
import io
import gzip
import math

import pandas as pd
import numpy as np

import proteomics.fasta as fasta
from time import gmtime, strftime

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
        argv, usage=__doc__)

    optional = parser.add_argument_group('optional arguments')
    required = parser.add_argument_group('required arguments')

    required.add_argument('-i', '--infile', dest="infile",
                          required=True, nargs='+',
                          help=("Provide a single file or folder with "
                                "multiple files for processing"))

    required.add_argument('-f', '--fasta-db', dest="fasta_db",
                          required=True,
                          help=("Input a fasta file for all proteins in "
                                "the species of interest"))

    required.add_argument('-fc', '--fasta-crap-db', dest="fasta_crap_db",
                          required=True,
                          help=("Input a fasta file for all proteins that "
                                "are common contaminants in a mass-spec "
                                "experiment"))

    required.add_argument('--peptide-column', dest="pep_column",
                          required=True,
                          help=("What's the name of the column with the "
                                "peptide sequence?"))

    optional.add_argument('--matches-column', dest="matches_column",
                          default=None,
                          help=("Column with the matches already identified "
                                "for the peptide"))

    optional.add_argument('--check-crap', dest="check_crap",
                          default=False, action='store_true',
                          help=("Check each peptide against the cRAP fasta"))

    optional.add_argument('--only-swissprot', dest="strict_sw",
                          default=False, action='store_true',
                          help=("Ignore matches to non-swissprot proteins"))

    optional.add_argument('--matches-separator', dest="matches_sep",
                          default=",",
                          help=("Separator for the matches column"))

    optional.add_argument('-o', '--outfile', dest="outfile",
                          default=None,
                          help=("Enter a file name for your output"))

    optional.add_argument('-os', '--outfile-suffix', dest="outfile_suffix",
                          default=None,
                          help=("Enter a suffix to add to the output files"))

    optional.add_argument('-l', '--logfile', dest="logfile",
                          default=os.devnull,
                          help=("Enter a file name for logging program "
                                "output. Else, nothing will be printed"))

    args = vars(parser.parse_args())

    if not args['outfile'] and not args['outfile_suffix']:
        raise ValueError("must supply either --outfile or "
                         "--outfile-suffix option")

    logfile = open(args['logfile'], 'w')
    logfile.write("Logfile for annotate_rnp.py %s\n\n" % (
        strftime("%Y-%m-%d %H:%M:%S", gmtime())))

    section_blocker = writeSectionHeader(logfile, "Script arguments:")
    for key, value in args.items():
        logfile.write("%s: %s\n" % (key, value))
    logfile.write("%s\n\n" % section_blocker)

    #(1) Get the mappings between peptide and proteins

    # (1.1) Build dictionaries using the fasta database to map from:
    # 1. protein accession: protein
    # 2. protein accession: sequence
    # 3. protein accession: description e.g >sp|O43707|ACTN4_HUMAN
    # 3. protein accession: long description e.g >sp|O43707|ACTN4_HUMAN|Alpha-actinin-4

    protein2description = {}
    protein2longdescription = {}
    protein2seq = {}
    tr_proteins = set()
    sp_proteins = set()

    for fa_infile in (args['fasta_db'], args['fasta_crap_db']):
        if fa_infile.endswith(".gz"):
            fa_iterator = fasta.FastaIterator(
                io.TextIOWrapper(gzip.open(fa_infile)))
        else:
            fa_iterator = fasta.FastaIterator(open(fa_infile))

        for entry in fa_iterator:
            accession = entry.title.split(" ")[0].split("|")[1]
            protein2seq[accession] = entry.sequence.upper().replace("I", "L")
            protein2description[accession] = entry.title.split(" ")[0]
            protein2longdescription[accession] = "|".join(entry.title.split(" ")[0:2])
            if entry.title.split(" ")[0].split("|")[0] == "sp":
                sp_proteins.add(accession)
            elif entry.title.split(" ")[0].split("|")[0] == "tr":
                tr_proteins.add(accession)
            else:
                raise ValueError("Protein does not appear to be either"
                                 "SwissProt(sp) or TrEMBL(tr)")

    crap_proteins = set()
    associated_crap_proteins = set()
    if args['fasta_crap_db'].endswith(".gz"):
        fa_iterator = fasta.FastaIterator(
            io.TextIOWrapper(gzip.open(args['fasta_crap_db'])))
    else:
        fa_iterator = fasta.FastaIterator(open(args['fasta_crap_db']))
    for entry in fa_iterator:
        accession = entry.title.split(" ")[0].split("|")[1]
        crap_proteins.add(accession)

    # (1.2) Parse the infiles to obtain maps of peptides to proteins and vis versa
    pep2pro = collections.defaultdict(lambda: collections.defaultdict(set))
    pep2allpro = collections.defaultdict(set)
    pro2pep = collections.defaultdict(set)
    top_level_proteins = set()
    initial_proteins = set()

    if not args['matches_column']:
        peptides = set()

    for infile in args['infile']:

        # read the data into a dataframe
        infile_split = infile.split(".")
        if infile_split[-1] == "xlsx":
            peptide_df = pd.read_excel(infile)
        elif infile_split[-1] in ["text", "txt", "tsv"]:
            peptide_df = pd.read_table(infile, sep='\t', comment=None)
        elif infile_split[-1] == "csv":
            peptide_df = pd.read_table(infile, sep=',', comment=None)
        else:
            raise ValueError("File type must one of .xlsx, "
                             ".txt, .text, .tsv, .csv")

        # add some basic annotations
        #rnp_df['tr_only'] = [x.count("sp|") == 0 for x in rnp_df['Proteins']]
        #rnp_df['matches'] = [len(x.split(",")) for x in rnp_df['Proteins']]

        # if matches have already been made, use these
        # (1.1) extract the initial mappings between proteins and peptides
        if args['matches_column']:
            for row_ix, row_values in peptide_df[
                    [args['matches_column'], args['pep_column']]].iterrows():

                # if empty match, will be converted to NaN (type=float)
                if type(row_values[args['matches_column']]) is float:
                    # could manually search against all proteins in database?
                    continue 

                proteins = row_values[args['matches_column']].split(
                    args['matches_sep'])

                peptide = row_values[args['pep_column']]

                if args['check_crap']:
                    add_crap_proteins = []
                    for prot in protein2seq:
                        if prot in crap_proteins:
                            if peptide.replace("I", "L") in protein2seq[prot]:
                                add_crap_proteins.append(prot)
                    proteins.extend(add_crap_proteins)

                for protein in proteins:
                    if protein in crap_proteins:
                        associated_crap_proteins.update(proteins)

                    if protein not in protein2seq:
                        logfile.write(
                            "protein %s matches peptide %s but is not found "
                            "in fasta database\n" % (protein, peptide))

                # remove proteins not in fasta database
                proteins = set([prot for prot in proteins if prot in protein2seq])

                if peptide in pep2pro:
                    if not pep2allpro[peptide] == proteins:
                        current_protein_matches = ", ".join(pep2allpro[peptide])
                        new_protein_matches = ", ".join(proteins)
                        logfile.write(
                            ("The same peptide is observed more than once with "
                             "different proteins! : peptide = %(peptide)s, "
                             "matching proteins = %(current_protein_matches)s "
                             "or %(new_protein_matches)s\n" % locals()))
                        pep2allpro[peptide].update(proteins)
                else:
                    pep2allpro[peptide] = proteins

                for protein in proteins:

                    initial_proteins.add(protein)
                    pro2pep[protein].add(peptide)

                    protein_description = protein2description[protein]
                    
                    if protein in sp_proteins:
                        protein_level = 1
                        top_level_proteins.add(protein)
                    elif protein in tr_proteins:
                        protein_level = 2
                    else:
                        raise ValueError("Protein does not appear to be either"
                                         "SwissProt(sp) or TrEMBL(tr)")

                    pep2pro[peptide][protein_level].add(protein)

        else: # if we don't have a column of matches, get the set of all peptides 
            peptides.update(peptide_df[args['pep_column']])

    if not args['matches_column']: 
        # search against all proteins in the provided databases
        n = 0
        for peptide in peptides:
            n += 1
            if n % 1000 == 0:
                logfile.write("searched %i peptides against database %s\n" % (
                    n, strftime("%Y-%m-%d %H:%M:%S", gmtime())))

            proteins = [prot for prot in protein2seq if
                        peptide in protein2seq[prot]]

            for protein in proteins:

                initial_proteins.add(protein)
                pro2pep[protein].add(peptide)

                protein_description = protein2description[protein]

                if protein in sp_proteins:
                    protein_level = 1
                    top_level_proteins.add(protein)
                elif protein in tr_proteins:
                    protein_level = 2
                else:
                    raise ValueError("Protein does not appear to be either"
                                     "SwissProt(sp) or TrEMBL(tr)")


                pep2pro[peptide][protein_level].add(protein)

    section_blocker = writeSectionHeader(logfile, "Initial file(s) stats")
    logfile.write("# initial peptides: %i\n" % len(pep2pro))
    logfile.write("# initial proteins: %i\n" % len(pro2pep))
    logfile.write("# initial SwissProt proteins: %i\n" % len(top_level_proteins))
    logfile.write("# initial TrEMBL proteins: %i\n" % (
        len(pro2pep)-len(top_level_proteins)))
    logfile.write("%s\n\n" % section_blocker)

    if not args['strict_sw']:
        section_blocker = writeSectionHeader(
            logfile, "Deciding which TrEMBL proteins to retain:")

        # (1.2) find the peptides with only TrEMBL protein matches and
        # 'upgrade' these TrEMBL proteins to being equivalent to SwissProt
        # across all peptides where these TrEMBL proteins match
        tr_only_peptides = set([x for x in pep2pro.keys() if len(pep2pro[x][1]) == 0])

        logfile.write("# peptides with only TrEMBL matches: %i\n" % (
            len(tr_only_peptides)))

        set_upgraded = set()
        for peptide in tr_only_peptides:
            upgraded = pep2pro[peptide][2]
            set_upgraded.update(upgraded)
            top_level_proteins.update(upgraded)

        logfile.write("# TrEMBL proteins retained as no SwissProt matches for "
                      "peptide: %i\n" % (len(set_upgraded)))

        # 'upgrade' the selected TrEMBL proteins
        for peptide in pep2pro:
            pep2pro[peptide][2] = pep2pro[peptide][2].difference(set_upgraded)
            pep2pro[peptide][1] = pep2pro[peptide][1].union(set_upgraded)

        logfile.write("%s\n\n" % section_blocker)

    # (1.3) Use a parsimonious approach to identifty the minimum number
    # of proteins required to cover all the peptides:
    # Start from the protein(s) with the most peptides and mark these as covered.
    # Continue with remaining proteins in order of peptides per protein
    # until all peptides are covered

    section_blocker = writeSectionHeader(
        logfile, ("Parsimonious method to identify minimal set of proteins"
                  " to account for all peptides %s" % (
                      strftime("%Y-%m-%d %H:%M:%S", gmtime()))))

    retained_proteins = []
    peptides = copy.deepcopy(set(pep2pro.keys()))
    peptide_counts = {}

    tmppro2pep = copy.deepcopy(pro2pep)
    new_top_level_proteins = copy.deepcopy(top_level_proteins)
    new_pep2pro = collections.defaultdict(set)

    peptide_count = max(map(len, tmppro2pep.values()))

    while True:
        # (1.3.1) If all peptides covered or the maximum peptides per
        # protein = 0, break.
        if len(peptides) == 0 or peptide_count == 0:
            logfile.write("All peptides are now accounted for %s\n" % (
                strftime("%Y-%m-%d %H:%M:%S", gmtime())))
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

        logfile.write("%i protein(s) with %i peptides\n" % (
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

    logfile.write("%i SwissProt proteins retained\n" % len(
        set(retained_proteins).intersection(sp_proteins)))

    logfile.write("%i TrEMBL proteins retained\n" % len(
        set(retained_proteins).intersection(tr_proteins)))

    logfile.write("\nNote: If not all SwissProt proteins were retained, this means\n"
                  "these proteins only included peptides which were observed\n"
                  "in other proteins which had a greater number of peptides\n")

    logfile.write("%s\n\n" % section_blocker)

    section_blocker = writeSectionHeader(logfile, "proteins per peptide:")
    counts = collections.Counter([len(x) for x in new_pep2pro.values()])
    sum_counts = sum(counts.values())
    for k, v in counts.items():
        logfile.write("%i peptide(s) (%.2f %%) have %i master protein(s)\n" % (
            v, (100 * v)/sum_counts, k))
    logfile.write("%s\n\n" % section_blocker)

    # Check all the peptides are covered
    if not args['strict_sw']:
        assert set(pep2pro.keys()).difference(set(new_pep2pro.keys())) == set()
    else:
        missing_peptides = set(pep2pro.keys()).difference(set(new_pep2pro.keys()))
        logfile.write("%i peptide(s) (%.2f %%) have no master protein(s)\n" % (
            len(missing_peptides), (100 * len(missing_peptides))/sum_counts))


    if args['outfile']:
        outfile = open(args['outfile'], "w")

    for infile in args['infile']:

        # read the data into a dataframe
        infile_split = infile.split(".")
        if infile_split[-1] == "xlsx":
            peptide_df = pd.read_excel(infile)
        elif infile_split[-1] in ["text", "txt", "tsv"]:
            peptide_df = pd.read_table(infile, sep='\t', comment=None)
        elif infile_split[-1] == "csv":
            peptide_df = pd.read_table(infile, sep=',', comment=None)
        else:
            raise ValueError("File type must one of .xlsx, "
                             ".txt, .text, .tsv, .csv")

        # add the top protein and uniprot id annotations
        peptide_df['master_protein'] = [
            ";".join(new_pep2pro[protein]) for protein in peptide_df[args['pep_column']]]

        # (1.5) derive further annotations
        protein_lengths = []
        protein_descriptions = []
        crap_protein = []
        associated_crap_protein = []
        peptide_start = []
        peptide_end = []

        for ix, row in peptide_df.iterrows():
            proteins = row['master_protein'].split(";")
            pep_sequence = row['Sequence'].upper().replace("I", "L").replace("X", "L")

            if proteins == [""]:
                protein_lengths.append("")
                protein_descriptions.append("")
                crap_protein.append("")
                associated_crap_protein.append("")
                peptide_start.append("")
                peptide_end.append("")
            else:
                protein_lengths.append(
                    ";".join(map(str, [len(protein2seq[x]) for x in proteins])))
                protein_descriptions.append(
                    ";".join([protein2description[x] for x in proteins]))

                # (1.5.1) does peptide match a cRAP protein?
                crap = 0
                for protein in proteins:
                    if protein in crap_proteins:
                        crap = 1
                        break
                crap_protein.append(crap)

                # (1.5.2) does peptide match a protein associated with a cRAP protein?
                associated_crap = 0
                for protein in proteins:
                    if protein in associated_crap_proteins:
                        associated_crap = 1
                        break
                associated_crap_protein.append(associated_crap)

                starts = []
                ends = []
                for protein in proteins:
                    protein_sequence = protein2seq[protein]
                    all_matches = re.findall(pep_sequence, protein_sequence)
                    if len(all_matches) > 1:
                        logfile.write(
                            "peptide: %s is observed more than once in protein: %s\n" % (
                                pep_sequence, protein))
                        starts.append("NA")
                        ends.append("NA")
                    elif len(all_matches) == 0:
                        logfile.write(
                            "peptide: %s is not observed in protein: %s\n" % (
                                pep_sequence, protein))
                        starts.append("NA")
                        ends.append("NA")
                    else:
                        peptide_match = re.search(pep_sequence, protein_sequence)
                        starts.append(peptide_match.start())
                        ends.append(peptide_match.end())

                try:
                    peptide_start.append(";".join(map(str, starts)))
                    peptide_end.append(";".join(map(str, ends)))
                except:
                    print(starts)
                    print(ends)
                    raise ValueError()


        peptide_df['protein_length'] = protein_lengths
        peptide_df['protein_description'] = protein_descriptions
        peptide_df['peptide_start'] = peptide_start
        peptide_df['peptide_end'] = peptide_end
        peptide_df['crap_protein'] = crap_protein
        peptide_df['associated_crap_protein'] = associated_crap_protein
        peptide_df['unique'] = [1 if len(x.split(";"))==1 else 0
                                for x in peptide_df['master_protein']]

        if args['outfile']:
            peptide_df['filename'] = infile
            peptide_df.to_csv(outfile, index=False, sep="\t", mode="a")
            os.chmod(args['outfile'], 0o666)
        else:
            outfile = ".".join(infile.split(".")[:-1]) + args['outfile_suffix']
            peptide_df.to_csv(outfile, index=False, sep="\t")
            os.chmod(outfile, 0o666)

    if args['outfile']:
        outfile.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
