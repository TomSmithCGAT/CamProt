'''get_idrs - Find the IDRs for 
=======================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python Proteomics Features

Purpose
-------

This scripts uses the D2P2 database (http://d2p2.pro/) to identify the
IDRs for a given set of proteins. Currently, it requires a tax id and
it will get all IDR infor for all Swiss Prot proteins for the tax id.


Usage
-----

python get_idrs.py --tax-id=9606 --cons=4 --size=20 > idrs.tsv



To do
-----
- Add option to specify proteins
- testing!

Command line options
--------------------
'''

import collections
import argparse
import sys
import os

import datetime

from bioservices import UniProt

from proteomics import d2p2
from proteomics import protinfo


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

    required.add_argument('-t', '--tax-id', dest="tax_id",
                          required=True, type=int,
                          help=("Tax id for species"))

    optional.add_argument('-o', '--outfile', dest="outfile",
                          help=("Where to write results"))

    optional.add_argument('-cons', '--consensus', dest="consensus",
                          default=4, type=int,
                          help=("How many tools must agree"))

    optional.add_argument('--size', dest="size",
                          default=20, type=int,
                          help=("Minimum size of IDR region"))

    optional.add_argument('--whitelist', dest="whitelist",
                          nargs='+',
                          help=("Comma separated list of tools to use"))

    optional.add_argument('--blacklist', dest="blacklist",
                          nargs='+',
                          help=("Comma separated list of tools to ignore"))

    optional.add_argument('-l', '--logfile', dest="logfile",
                          default=os.devnull,
                          help=("Enter a file name for logging program "
                                "output. Else, nothing will be printed"))

    args = vars(parser.parse_args())

    if args['logfile']:
        logfile = open(args['logfile'], 'w')
    else:
        logfile = open(os.devnull,"w")

    logfile.write("Logfile for get_ids.py %s\n\n" % (
        datetime.datetime.now()))

    section_blocker = writeSectionHeader(logfile, "Script arguments:")
    for key, value in args.items():
        logfile.write("%s: %s\n" % (key, value))
    logfile.write("%s\n\n" % section_blocker)

    # 1. Get all uniprot IDs for species
    u = UniProt()
    results = u.search("organism:%s+and+reviewed:yes" % args['tax_id'], columns="id")
    uniprot_ids = [x.split()[0] for x in results.strip().split("\n")[1:]]

    # 2. Get IDRs, rebuild consensus and get IDR blocks
    section_blocker = writeSectionHeader(logfile, "Getting IDR info from D2P2...")
    idr_data_available = set()
    idr_data_unavailable = set()
    sequence_length_idr = {}

    blocks = collections.defaultdict(list)

    n = 0

    start_time = datetime.datetime.now()

    for d2p2_entry in d2p2.iterator(uniprot_ids, chunk_size=50):

        d2p2_entry.rebuildConsensus(
            tools_whitelist=args['whitelist'],
            tools_blacklist=args['blacklist'])

        d2p2_entry.setIDRs(args['consensus'], args['size'])

        blocks[d2p2_entry.name] = d2p2_entry.idrs

        idr_data_available.add(d2p2_entry.name)
        sequence_length_idr[d2p2_entry.name] = d2p2_entry.length

        n += 1
        current_time = datetime.datetime.now()
        pred_finish = start_time + ((current_time-start_time) * len(uniprot_ids)/n)
        logfile.write('proteins done: %i/%i  at %s. Pred. finish = %s\r' % (
            n, len(uniprot_ids), current_time, pred_finish))
        logfile.flush()

    logfile.write('\nfinished: %s\n' % current_time)
    logfile.flush()

    # 3. Log proteins with no IDR data
    idr_data_unavailable = set(uniprot_ids).difference(idr_data_available)

    logfile.write("Unable to detect IDR for %s / %s proteins\n" % (
        len(idr_data_unavailable), len(uniprot_ids)))
    logfile.write("No IDRs for:\n%s\n" %  "\n".join(idr_data_unavailable))

    logfile.write("%s\n\n" % section_blocker)

    # 3. Write out
    if args['outfile']:
        outf = open(args['outfile'], "w")
    else:
        outf = sys.stdout
    
    outf.write(
        "\t".join(
            ("UniprotID", "protein_length", "total_idr_length",
             "fraction_idr", "idrs")) + "\n")
    for protein in idr_data_available:
        total_idr = 0
        for block in blocks[protein]:
            if not block[1] - block[0] >= args['size']:
                raise ValueError(
                    "This block is too small, how did it get here?!: %s" % block)
            total_idr += (block[1] - block[0])
        protein_length = sequence_length_idr[protein]
        fraction_idr = total_idr/protein_length
        outf.write("\t".join(map(str, (
            protein, protein_length, total_idr,
            fraction_idr, len(blocks[protein])))) + "\n")

    #if args['outfile']:
    outf.close()

    #if args['logfile']:
    logfile.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
