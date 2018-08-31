'''get_protein_properties.py - get basic properties for proteins
================================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Proteomics

Purpose
-------

Usage
-----

Command line options
--------------------
'''

import argparse
import gzip
import io
import sys

import proteomics.fasta as fasta
import proteomics.sequence as sequence

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

    required = parser.add_argument_group('required arguments')

    required.add_argument('-i', '--infile', dest="infile",
                          required=True, default=None, 
                          help=("Fasta file with proteins"))

    required.add_argument('-o', '--outfile', dest="outfile",
                          required=True, default=None,
                          help=("File name for output"))

    required.add_argument('-l', '--logfile', dest="logfile",
                          required=True, default=None,
                          help=("File name for logfile"))


    args = vars(parser.parse_args())


    logfile = open(args['logfile'], 'w')
    logfile.write("Logfile for get_protein_properties.py %s\n\n" % (
        strftime("%Y-%m-%d %H:%M:%S", gmtime())))

    section_blocker = writeSectionHeader(logfile, "Script arguments:")
    for key, value in args.items():
        logfile.write("%s: %s\n" % (key, value))
    logfile.write("%s\n\n" % section_blocker)


    with open(args['outfile'], "w") as outf:
    
        outf.write("%s\n" % "\t".join(
            ("accession", "swissprot", "fraction_hydrophobic",
             "fraction_hydrophillic", "free_energy_change", "length")))

        if args['infile'].endswith(".gz"):
            iterator = fasta.FastaIterator(
                io.TextIOWrapper(gzip.open(args['infile'])))
        else:
            iterator = fasta.FastaIterator(open(args['infile']))
        
        for entry in iterator:
            if entry.title.startswith("sp|"):
                swissprot = True
            else:
                swissprot = False

            accession = entry.title.split("|")[1]

            fractionHydrophobic = round(
                sequence.getFractionHydrophobic(entry.sequence), 3)
            fractionHydrophillic = round(
                sequence.getFractionHydrophillic(entry.sequence), 3)

            try:
                freeEnergyChange = round(sequence.getFreeEnergy(entry.sequence), 2)
                freeEnergyChange = freeEnergyChange / len(
                    [aa for aa in entry.sequence if aa not in ["X", "U"]])
            except:
                print(entry)
                break

            outf.write(
                "%s\n" % "\t".join(map(str, (accession, str(swissprot),
                                             fractionHydrophobic, fractionHydrophillic,
                                             freeEnergyChange, str(len(entry.sequence))))))

if __name__ == "__main__":
    sys.exit(main(sys.argv))
