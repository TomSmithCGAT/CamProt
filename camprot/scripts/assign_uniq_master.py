'''assign_uniq_master - assign a unique master to each peptide across
all samples using a maximum parsimony approach
=====================================================================

:Author: Tom Smith, Manasa Ramakrishna
:Release: $Id$
:Date: |today|
:Tags: Python Proteomics Peptide-assignment 

Purpose
-------

This script takes the xlsx output from a set of input files
(*.txt/*.xlsx) and annotates the table with unique protein information
for downstream analyses.

The following columns are added:

 - master_protein(s): The master protein(s) for the peptide. See below
   for how this is derived
 - master_uniprot_id(s): The uniprot id(s) for the master protein(s)
 - protein_description(s): Description(s) for the master protein(s)
 - protein_length(s): The length(s) of the master protein(s)
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

python assign_uniq_master.py --infile=RNP.xlsx --fasta-db=Human_201701.fasta
--fasta-crap-db=cRAP_FullIdentifiers.fasta --outfile=master_prot_annotated.txt
--logfile=master_prot_annot.log

Command line options
--------------------

'''
#!/usr/bin/env python3

import argparse
import collections 
import copy
import errno
import glob
import os
import pandas as pd
import re
import sys
import proteomics.fasta as fasta

def writeSectionHeader(logfile, section_header):
    #underliner = "".join(("-",)*len(section_header))
    section_blocker = ("======================================="
               "=======================================")
    underliner1 = ("----------------------------------------"
               "----------------------------------------")
    logfile.write("\n%s\n%s\n" % (section_blocker, section_header))
    logfile.write("%s\n" % underliner1)
    return section_blocker

####################### ------------------------- ##############################

def main(argv=sys.argv):

    parser = argparse.ArgumentParser(argv, usage=__doc__)

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

    optional.add_argument('-o', '--outfile', dest="outfile", default=None,
                          help=("Enter a file name for your output. Else "
                                "it will be the same as your input with "
                                "the suffix 'annotated'"))

    optional.add_argument('-l', '--logfile', dest="logfile",
                          default=os.devnull,
                          help=("Enter a file name for logging program "
                                "output. Else, nothing will be printed"))

    args = vars(parser.parse_args())

        # creating an output file if no name is provided
        #if args['outfile'] is None:
        #       args['outfile'] = args.pop('infile')[1].replace(".txt", "_annotated.txt")

    # create a log file (os.devnull by default)
    logfile = open(args['logfile'], 'w')
    logfile.write("\nLogfile for assign_uniq_master.py\n")
    section_blocker = writeSectionHeader(logfile, "Script arguments:")

    # Outputs each input parameter and user defined value
    for key, value in args.items():
        logfile.write("%s: %s\n" % (key, value))
    logfile.write("%s\n\n" % section_blocker)

    #------------------------------------------------------------
    # (1.0) Build dictionaries to map from protein id to protein
    # sequence and description using the fasta database
    # protein2seq : key = protein ID, value = sequence
    # crap_proteins : set with IDs of contaminant proteins
    #------------------------------------------------------------

    crap_proteins = set()

    protein2description = {
        entry.title.split(" ")[0]: " ".join(entry.title.split(" ")[1:])
        for entry in fasta.FastaIterator(open(args['fasta_db']))}

    protein2seq = {
        entry.title.split(" ")[0]:entry.sequence
        for entry in fasta.FastaIterator(open(args['fasta_db']))}

    for entry in fasta.FastaIterator(open(args['fasta_crap_db'])):
        protein2seq[entry.title.split(" ")[0]] = entry.sequence
        protein2description[entry.title.split(" ")[0]] = entry.title.split(" ")[0]
        crap_proteins.add(entry.title.split(" ")[0])

    #------------------------------------------------
    # read each infile into a dataframe and parse rows
    #------------------------------------------------

    for f in args['infile']:
        fn = f.split(".")
        if fn[1] == "xlsx":
            #print("It is an excel file")
            rnp_df = pd.read_excel(f)
            #print(f, rnp_df.shape)
        elif fn[1] in ["text", "txt"]:
            #print("It is a txt file")
            rnp_df = pd.read_table(f,sep='\t',comment=None)
            #print(f, rnp_df.shape)

        else:
            raise ValueError("File type must one of .xlsx, .txt, .text")

        new_protein_sp_ids = set()
        new_protein_tr_ids = set()

        # Loop through each peptide and find all SwissProt and trEMBL
        # proteins that they map to

        for ix, row in rnp_df.iterrows():
            peptide = row['Sequence']

            # For each peptide, find all the protein IDs in which it is a perfect match
            l = [k for k in protein2seq.keys() if peptide in protein2seq[k]]

            # Collect all the matches that were "Swissprot" matches
            l_sp = [p for p in l if "sp|" in p]

            # Extract just the protein ids from the header line
            sp_id = [i.split("|")[1] for i in l_sp] 

            # Join the protein ids into a string a print
            join_sp_id = ";".join(sp_id)
            #print("SwissPROT:",join_sp_id)
                        
            # Collect all the lower confidence "TremBL" matches
            l_tr = [ p for p in l if "tr|" in p]
            tr_id = [i.split("|")[1] for i in l_tr]
            join_tr_id = ";".join(tr_id)
            #print("Trembl : ",join_tr_id)  
        
        
            # We dont have anything other than SwissProt or trEMBL ids so no need to look for anything else
                        
            # (1.5.1) does peptide match a cRAP protein?
            crap = 0
            for protein in sp_id:
                if protein in crap_proteins:
                    crap = 1
                    break
            
            crap_protein.update(crap
            new_protein_sp_ids.update(sp_id)
            new_protein_tr_ids.updates(tr_id)
        
            rnp_df['new_protein_sp_ids'] = new_protein_sp_ids
            rnp_df['new_protein_tr_ids'] = new_protein_tr_ids
            rnp_df['crap_protein'] = crap_protein
                
        # Now that we have matched each peptide to a protein using the
        # FASTA and cRAP fasta files We want to use the maximum
        # parsimony approach so that a set of peptides arising from
        # the same protein map back uniquely to that protein. We have
        # found that some lower quality peptides map to isoforms of
        # lower quality stored in trEMBL so we only choose those that
        # are in SwissProt database. We need to go through the
        # parsimonious approach with Swissprot entries as top
        # priorities.  If this does not return a protein, then we go
        # to trEMBL and repeat the exercise
        
        #(1) Get the mappings between peptide and proteins
        pep2pro = collections.defaultdict(lambda: collections.defaultdict(set))
        pep2allpro = collections.defaultdict(set)
        pro2pep = collections.defaultdict(set)
        initial_proteins = set()
        
        # (1.1) extract the initial mappings between proteins and peptides
        
        for row_ix, row_values in rnp_df[['new_protein_sp_ids', 'Sequence']].iterrows():
                proteins = row_values['new_protein_sp_ids'].split(";")
                peptide = row_values['Sequence']

                if peptide in pep2pro:
                        assert pep2allpro[peptide] == proteins, (
                                "The same peptide is observed more than once with different proteins!")
                
                pep2allpro[peptide] = proteins
                for protein in proteins:
                        initial_proteins.add(protein)
                        pro2pep[protein].add(peptide)
                        pep2pro[peptide][1].add(protein)
        
        #print(pep2pro.keys())
        #print(pro2pep)
        #print(initial_proteins)

        section_blocker = writeSectionHeader(logfile, "Initial file stats")
        logfile.write("# initial peptides: %i\n" % len(pep2pro))
        logfile.write("# initial proteins: %i\n" % len(pro2pep))
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
        new_top_level_proteins = copy.deepcopy(initial_proteins)
        new_pep2pro = collections.defaultdict(set)

        peptide_count = max(map(len, tmppro2pep.values()))
        
        #print(tmppro2pep)
        #print(new_top_level_proteins)
        #print(new_pep2pro)
        print("max peptide count is ", peptide_count)

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
                # Iterate through to find the protein with most number of matches. 
                # top_starts off at 0. Then add a protein, then go to the next one and see if it is higher. If yes, make this the new value. If no, continue.
                for protein in new_top_level_proteins:
                        if len(tmppro2pep[protein]) == top_score:
                                top_proteins.add(protein)
                        elif len(tmppro2pep[protein]) > top_score:
                                top_score = len(tmppro2pep[protein])
                                top_proteins = set((protein,))

                logfile.write("%i remaining protein(s) with %i peptides\n" % (
                        len(top_proteins), top_score))
'''             
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
        #                                                for x in new_top_level_proteins]))
        logfile.write("%i SwissProt proteins retained\n" % len(
                [x for x in retained_proteins if x.split("|")[0] == 'sp']))
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
        assert set(pep2pro.keys()).difference(set(new_pep2pro.keys())) == set()

        # add the top protein and uniprot id annotations
        rnp_df['master_protein(s)'] = [";".join(new_pep2pro[protein]) for protein in rnp_df['new_protein_sp_ids']]

        rnp_df.to_tsv(args['outfile'], index=False, sep="\t")
        os.chmod(args['outfile'], 0o666)

'''

if __name__ == "__main__":
        sys.exit(main(sys.argv))




































  
