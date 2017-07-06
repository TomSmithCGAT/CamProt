'''plot_tm_coverage
=======================================================

:Author: Tom Smith
:Release: $Id$
:Date: |today|
:Tags: Python Proteomics

Purpose
-------

Usage
-----

Command line options
--------------------
'''

import argparse
import collections 
import os
import re
import sys
import io
import gzip
import math

import pandas as pd
import numpy as np

import requests
import json

import proteomics.fasta as fasta
import proteomics.sequence as sequence

from time import gmtime, strftime
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def writeSectionHeader(logfile, section_header):
    #underliner = "".join(("-",)*len(section_header))
    section_blocker = ("======================================="
                       "=======================================")
    underliner1 = ("----------------------------------------"
                  "----------------------------------------")
    logfile.write("\n%s\n%s\n" % (section_blocker, section_header))
    logfile.write("%s\n" % underliner1)
    return section_blocker

def iterateProteomicsSummary(infile, sep="\t"):
    ''' iterate through a proteomics summary file and yield the data per line'''
    header = next(infile)
    for line in infile:
        try:
            uniprot_id, sequence, start, end = line.strip().split(sep)[0:4]
        except:
            print(line.strip().split(sep))
            print(line.strip().split(sep)[0:4])
            print(line)
            raise ValueError()
        yield (uniprot_id, sequence, start, end)

def getPeptidePosition(peptide_seq, protein_seq):
    ''' find the positions in the protein where the peptide matches. Returns a set of all matched positions'''
    if peptide_seq in protein_seq:
        if len(re.findall(peptide_seq, protein_seq)) > 1:
            #print("more than 1 match")
            return "more than 1 match"
        else:
            # 1 match
            span = re.search(peptide_seq, protein_seq).span()
            return set(range(span[0], span[1]))
            
    else:
        # no matches
        return None
    

def normaliseArraySize(input_array, size=10):
    ''' normalise an array to a set size, e.g [0,0,1,1] -> [0,0,0,1,1,1] or [0,1]
    Can expand or contract to an array of any length'''
    new_array = np.zeros(size)
    
    if len(input_array) == size:
        return input_array
    
    bin_size = len(input_array)/float(size)
    #print(bin_size)
    lower_edge = 0
    upper_edge = 0
    for array_ix in range(0, size):
        total_aas = (array_ix + 1) * bin_size
        upper_edge = total_aas
        
        combined_coverage = 0
        if bin_size < 1:
            if math.floor(lower_edge) == math.floor(upper_edge) or math.floor(upper_edge) == len(input_array):
                combined_coverage += (upper_edge - lower_edge) * input_array[math.floor(lower_edge)]    
            else:
                combined_coverage += (math.ceil(lower_edge) - lower_edge) * input_array[math.floor(lower_edge)]
                combined_coverage += (upper_edge - math.floor(upper_edge)) * input_array[math.floor(upper_edge)]
        else:
            combined_coverage += ((math.ceil(lower_edge) - lower_edge) *
                                  input_array[math.floor(lower_edge)])
            combined_coverage += input_array[math.ceil(lower_edge):math.floor(upper_edge)].sum()
            combined_coverage += ((upper_edge - math.floor(upper_edge)) *
                                  input_array[min(math.floor(upper_edge), len(input_array)-1)])

        new_array[array_ix] = combined_coverage / (upper_edge - lower_edge)
        lower_edge = upper_edge
        
    # the average coverage shouldn't have changed. Due to floating point arithemetic, 
    # possible there is a slight change
    if abs(new_array.mean() - input_array.mean()) > 0.001:
        print(input_array, input_array.mean())
        print(new_array, new_array.mean())
        raise ValueError()
    
    return new_array

def normalisePaddedArray(input_array, upstream_pad, downstream_pad, tm_length, size):
    ''' Normalise an array composed of three regions (upstream, middle and end), where middle = TM'''
    assert len(input_array) == sum((upstream_pad, downstream_pad, tm_length))
    
    up_array = input_array[0:upstream_pad]
    tm_array = input_array[upstream_pad:upstream_pad+tm_length]
    down_array = input_array[upstream_pad+tm_length:]
    new_array = np.zeros(3 * size)
    
    new_array[0:size] = normaliseArraySize(up_array, size=size)
    new_array[size:2*size] = normaliseArraySize(tm_array, size=size)
    new_array[2*size: 3*size] = normaliseArraySize(down_array, size=size)
    
    return(new_array)

def getTMCoverage(protein_coverage, tm_blocks, FEATURE_SIZE, sequence_dict, debug=False):

    tm_blocks_region_coverage_norm = collections.defaultdict(
        lambda: collections.defaultdict(
            lambda: collections.defaultdict(np.array)))

    for protein in tm_blocks:
        for tm_ix, tm in enumerate(tm_blocks[protein]):

            tm_start, tm_stop = tm

            tm_blocks_region_coverage_norm[protein][tm] = np.zeros(
                (len(protein_coverage[protein]), 3*FEATURE_SIZE))
            if tm_ix == 0:
                upstream_length = min(tm_start, FEATURE_SIZE)
                if len(tm_blocks[protein]) == 1:
                    downstream_length = min(len(sequence_dict[protein]) - tm_stop, FEATURE_SIZE)
                else:
                    downstream_length = min(tm_blocks[protein][tm_ix+1][0] - tm_stop, FEATURE_SIZE)
            elif tm_ix == (len(tm_blocks[protein]) - 1):
                upstream_length = min(tm_start - tm_blocks[protein][tm_ix-1][1], FEATURE_SIZE)
                downstream_length = min(len(sequence_dict[protein]) - tm_stop, FEATURE_SIZE)
            else:
                upstream_length = min(tm_start - tm_blocks[protein][tm_ix-1][1], FEATURE_SIZE)
                downstream_length = min(tm_blocks[protein][tm_ix+1][0] - tm_stop, FEATURE_SIZE)

            # print statements left in for debuggin purposes
            if upstream_length == 0:
                if tm_ix == (len(tm_blocks[protein]) - 1):
                    if debug:
                        print("skipping TM (up): ", protein, tm_blocks[protein][tm_ix],
                              len(sequence_dict[protein]))
                else:
                    if debug:
                        print("skipping TM (up): ", protein, tm_blocks[protein][tm_ix - 1],
                              tm_blocks[protein][tm_ix])
                continue

            if downstream_length == 0:            
                if tm_ix == (len(tm_blocks[protein]) - 1):
                    if debug:
                        print("skipping TM (down): ", protein, tm_blocks[protein][tm_ix],
                              len(sequence_dict[protein]))
                else:
                    if debug:
                        print("skipping TM (down): ", protein, tm_blocks[protein][tm_ix],
                              tm_blocks[protein][tm_ix+1])
                continue


            for array_ix, protein_array in enumerate(protein_coverage[protein]):
                tm_array_raw = protein_array[tm_start-upstream_length:tm_stop+downstream_length]
                #print(len(tm_array_raw))
                new_array = normalisePaddedArray(
                    tm_array_raw, upstream_pad=upstream_length, downstream_pad=downstream_length,
                    tm_length=tm_stop-tm_start, size=FEATURE_SIZE)
                tm_blocks_region_coverage_norm[protein][tm][array_ix] = new_array

    return tm_blocks_region_coverage_norm

def normalisePerStudy(study_ix, array_dict, tms, feature_size):
    tms_array = np.zeros((tms, (3 * feature_size)))
    n = 0
    for protein in array_dict:
        for tm in array_dict[protein]:
            tms_array[n] = array_dict[protein][tm][study_ix]
            n+=1
    
    return tms_array.mean(axis=0)

def makeMetaTMDF(tm_blocks_region_coverage_norm, ix2desc, tms, FEATURE_SIZE):

    bins = []
    desc = []
    ixs = []
    coverage_abs = []
    coverage_norm = []

    for ix in ix2desc:
        coverage_array = normalisePerStudy(ix, tm_blocks_region_coverage_norm, tms, FEATURE_SIZE)
        max_coverage = max(coverage_array)
        coverage_array_norm = [x/max_coverage for x in coverage_array]

        coverage_abs.extend(coverage_array)
        coverage_norm.extend(coverage_array_norm)

        bins.extend(range(0, len(coverage_array)))
        desc.extend((ix2desc[ix],)*len(coverage_array))
        ixs.extend((ix,)*len(coverage_array))

    coverage_profile_df = pd.DataFrame({"desc":desc, "ix":ixs,
                                        "coverage_abs":coverage_abs, "coverage_norm":coverage_norm,
                                        "bins":bins})
    return(coverage_profile_df)

def makeProteinCoveredDF(protein_coverage, tm_blocks, ix2descr, sequence_dict):

    rows = []

    for uniprot_id in tm_blocks:
        covered = np.zeros(len(protein_coverage[uniprot_id]))
        not_covered = np.zeros(len(protein_coverage[uniprot_id]))

        
        try:
            covered += protein_coverage[uniprot_id].sum(axis=1)
        except:
            print(uniprot_id)
            print(uniprot_id in tm_blocks)
            print(protein_coverage[uniprot_id])
            raise ValueError()

        not_covered += len(sequence_dict[uniprot_id]) - protein_coverage[uniprot_id].sum(axis=1)

        total_aas = covered + not_covered
        coverage = covered / total_aas
        #print(coverage)

        if total_aas[0]>0:
            for ix, cov in enumerate(coverage):
                rows.append([uniprot_id, ix, ix2descr[ix], 1, cov])

            #combined_coverage = []
            #for tm in sorted(tm_blocks_coverage[uniprot_id]):
            #    combined_coverage += list(tm_blocks_coverage[uniprot_id][tm].max(axis=0))
            #rows.append([uniprot_id, tms, study_ix+1, "combined", 1, np.mean(combined_coverage)])
            #print(rows)
            #raise ValueError()

    protein_coverage_df = pd.DataFrame.from_records(
        rows, columns=["uniprot_id", "ix", "desc", "length", "coverage"])

    return protein_coverage_df

def makeProteinCoveredDF(protein_coverage, tm_blocks, ix2descr):

    rows = []

    for uniprot_id in tm_blocks:
        covered = np.zeros(len(protein_coverage[uniprot_id]))
        not_covered = np.zeros(len(protein_coverage[uniprot_id]))

        
        try:
            covered += protein_coverage[uniprot_id].sum(axis=1)
        except:
            print(uniprot_id)
            print(uniprot_id in tm_blocks)
            print(protein_coverage[uniprot_id])
            raise ValueError()

        not_covered += len(sequence_dict[uniprot_id]) - protein_coverage[uniprot_id].sum(axis=1)

        total_aas = covered + not_covered
        coverage = covered / total_aas
        #print(coverage)

        if total_aas[0]>0:
            for ix, cov in enumerate(coverage):
                rows.append([uniprot_id, ix, ix2descr[ix], 1, cov])

            #combined_coverage = []
            #for tm in sorted(tm_blocks_coverage[uniprot_id]):
            #    combined_coverage += list(tm_blocks_coverage[uniprot_id][tm].max(axis=0))
            #rows.append([uniprot_id, tms, study_ix+1, "combined", 1, np.mean(combined_coverage)])
            #print(rows)
            #raise ValueError()

    protein_coverage_df = pd.DataFrame.from_records(
        rows, columns=["uniprot_id", "ix", "desc", "length", "coverage"])

    return protein_coverage_df

def makeTMCoveredDF(protein_coverage, tm_blocks, ix2descr):
    
    rows = []

    for uniprot_id in tm_blocks:
        covered = np.zeros(len(protein_coverage[uniprot_id]))
        not_covered = np.zeros(len(protein_coverage[uniprot_id]))

        tms = 0

        for tm in sorted(tm_blocks[uniprot_id]):
            try:
                tm_length = tm[1] - tm[0]
            except:
                continue

            covered += protein_coverage[uniprot_id][:,tm[0]:tm[1]].sum(axis=1)
            not_covered += tm_length - protein_coverage[uniprot_id][:,tm[0]:tm[1]].sum(axis=1)
            tms += 1

        #print(covered)
        #print(not_covered)
        #raise ValueError()
        
        total_aas = covered + not_covered
        coverage = covered / total_aas

        if total_aas[0]>0:
            for ix, cov in enumerate(coverage):
                rows.append([uniprot_id, tms, ix, ix2descr[ix], 1, cov])

            combined_coverage = []
            for tm in sorted(tm_blocks[uniprot_id]):
                combined_coverage += list(protein_coverage[uniprot_id][:,tm[0]:tm[1]].sum(axis=0))
            rows.append([uniprot_id, tms, ix+1, "combined", 1, np.mean(combined_coverage)])

    tm_coverage_df = pd.DataFrame.from_records(
        rows, columns=["uniprot_id", "tms", "ix", "desc", "length", "coverage"])

    return tm_coverage_df

def finaliseDataFrame(df, coverage_threshold=0):
    tmp_df = df.copy()
    tmp_df['coverage'] = tmp_df['coverage'] > coverage_threshold
    tmp_df.set_index("uniprot_id", inplace=True, drop=False)
    tmp_df['order_1'] = tmp_df[tmp_df['coverage'] > 0].groupby("uniprot_id")[
        'desc'].aggregate(lambda x: len(x))

    tmp_df['order_2'] = tmp_df[tmp_df['coverage'] > 0].groupby("uniprot_id")[
        'desc'].sum()

    tmp_df['order_3'] = tmp_df.groupby("uniprot_id")['coverage'].sum()
    tmp_df['order_1'] = tmp_df['order_1'].fillna(max(tmp_df['ix'])+1).astype("int")
    tmp_df['order_2'] = tmp_df['order_2'].fillna(0).astype("str")

    tmp_df = tmp_df.sort_values(
        ['order_1', 'order_2', 'order_3'], ascending=[True, False, False])
    tmp_df.reset_index(drop=True, inplace=True)

    tmp_df['cum_end'] = tmp_df.groupby("desc")['length'].apply(lambda x: np.cumsum(x))
    tmp_df['cum_start'] = tmp_df['cum_end'] - tmp_df['length']

    return tmp_df

plotCoverage = R('''
    function(coverage_profile_df, fraction_plotname, norm_plotname){
    library(ggplot2)

    # check the data structure

    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    #coverage_profile_df$desc = factor(coverage_profile_df$desc, levels=c(
    #    "krug", "ccp_qe_tl_90c", "ccp_lumos_chopin_tl", "ccp_lumos_trypsin"))

    my_theme = theme(
    aspect.ratio=1,
    text=element_text(size=20),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.border=element_blank(),
    panel.grid=element_blank(),
    legend.text=element_text(size=12))

    p = ggplot(coverage_profile_df) +
    theme_bw() + my_theme +
    scale_colour_manual(name="", values=cbPalette[1:(max(coverage_profile_df$ix) + 1)]) 

    # make relative coverage plot
    p1 = p + aes(bins, coverage_norm, colour=desc) + ylab("Normalised Coverage") +
    geom_rect(xmin=10, xmax=20, ymin=-0.1, ymax=-0.05, fill = "black", colour="black") + # add TM model
    geom_segment(x = 0, y = -0.075, xend = 30, yend = -0.075, color = "black") +  # add TM model
    annotate(geom="text", x=15, y=-0.075, label="TM", color="white") + # add TM model

    geom_segment(x = 0, y = 0, xend=30, yend = 0, color = "grey50") +  # add manual x-axis
    geom_segment(x = 0, y = 0, xend=0, yend = 1, color = "grey50") +  # add manual y-axis

    scale_y_continuous(limits=c(-0.1,1), breaks=seq(0,1,0.25)) +
    geom_line()

    max_abs_coverage <- max(coverage_profile_df$coverage_abs)

    # make absolute coverage plot
    p2 = p + aes(bins, coverage_abs, colour=desc) + ylab("Fraction Covered") +

    geom_rect(xmin=10, xmax=20, ymin=-(max_abs_coverage/20), ymax=-(max_abs_coverage/10),
              fill = "black", colour="black") +  # add TM model
    geom_segment(x = 0, y = -(max_abs_coverage/13.3), xend = 30,
                 yend = -(max_abs_coverage/13.3), color = "black") +  # add TM model
    annotate(geom="text", x=15, y=-(max_abs_coverage/13.3), label="TM", color="white") +  # add TM model

    geom_segment(x = 0, y = 0, xend=30, yend = 0, color = "grey50") + # add manual x-axis
    geom_segment(x = 0, y = 0, xend=0, yend = max_abs_coverage, color = "grey50") + # add manual x-axis

    scale_y_continuous(limits=c(-(max_abs_coverage/10),max_abs_coverage)) +
    geom_line()

    ggsave(p2, file=fraction_plotname)
    ggsave(p1, file=norm_plotname)

    }
    ''')

plotCoveredTM = R('''
    function(
        input_df,
        plotfilename,
        circular=TRUE,
        plot_only_covered=FALSE,
        colour_by_group=FALSE,
        group='desc'){

    library(ggplot2)

    tmp_df <- input_df
    tmp_df <- tmp_df[tmp_df[[group]] != "combined",]

    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    if (colour_by_group==TRUE){
        tmp_df$fill = tmp_df[[group]]
        tmp_df$fill[tmp_df$coverage == 0] <- NA
    }

    else{
        tmp_df$fill = tmp_df$coverage > 0
    }

    tmp_df$fill <- factor(tmp_df$fill)

    my_theme <- theme(
    text=element_text(size=20),
    legend.text=element_text(size=8),
    axis.title=element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank())

    my_theme_circ <- theme(axis.text=element_blank(),
    axis.ticks=element_blank())    

    plot_text <- "Samples from outside\n to inside:\n\n"

    n_studies <- max(tmp_df$ix)
    studies_step = 0.5 / (n_studies + 1)
    plot_ymax <- 1
    plot_ymin <- 0.5
    tmp_df$ymin <- plot_ymin + (tmp_df$ix  * studies_step)
    tmp_df$ymax <- plot_ymin + ((tmp_df$ix +1)  * studies_step)

    rect_colour <- NA
    if (plot_only_covered==TRUE){
        p <- ggplot(tmp_df[tmp_df$uniprot_id %in% tmp_df[tmp_df$coverage>0,"uniprot_id"],])
        if (circular==TRUE){rect_colour <- "grey93"}
    }
    else{
        p <- ggplot(tmp_df)
    } 

    p <- p  + theme_bw() + my_theme # add themes

    if (circular==TRUE) {

        p <- p + geom_rect(aes(xmin=cum_start, xmax=cum_end, ymin=ymin, ymax=ymax, fill=factor(fill)),
                           colour=rect_colour, size=0.1) +  coord_polar()  + my_theme_circ + ylim(0, 1)

    }
    else{
        p <- p + geom_rect(aes(xmin=cum_start, xmax=cum_end, ymin=ix, ymax=ix+1, fill=factor(fill)),
                           colour=rect_colour, size=0.15) #+
        scale_y_continuous(breaks=seq(0.5, max(tmp_df$ix+0.5, 0.5)))#, labels=levels(tmp_df$fill))
    }


    if (colour_by_group==TRUE){
        p <- p + scale_fill_manual(na.value="grey97",
                                   breaks=unique(tmp_df[[group]]),
                                   values=cbPalette[1:length(unique(tmp_df[[group]]))],
                                   name="")

        if (circular==TRUE){
            p <- p + theme(legend.position=c(0.5,0.5))
        }
    }
    else{
        p <- p + scale_fill_manual(values=c("grey97", "grey13")) + theme(legend.position="none")
        if (circular==TRUE){
            plot_text <- paste0(plot_text, paste(rev(unique(tmp_df[[group]])), collapse="\n"))
            p = p + annotate("text", x = 0, y = 0, label=plot_text, size=8)
        }
    }

    ggsave(p, file=plotfilename)
    }''')


def main(argv=sys.argv):

    parser = argparse.ArgumentParser(
        argv, usage=__doc__)

    optional = parser.add_argument_group('optional arguments')
    required = parser.add_argument_group('required arguments')

    required.add_argument('-i', '--infiles', dest="infiles",
                          required=True, nargs='+',
                          help=("Provide a list of files for processing"))

    required.add_argument('-o', '--outdir', dest="outdir",
                          required=True, help=("Outdir for plots"))

    optional.add_argument('-l', '--logfile', dest="logfile",
                          default=os.devnull,
                          help=("Enter a file name for logging program "
                                "output. Else, nothing will be printed"))

    optional.add_argument('--tax-id', dest="tax_id", default=None,
                          help=("Taxomony id (e.g 9606 for H.sapiens)"))

    optional.add_argument('--field-separator', dest="sep", default="\t",
                          help=("Separator for columns in infiles"))

    optional.add_argument('--coverage-threshold', dest="coverage_thresh", default=0,
                          help=("Threshold for a TM to be considered covered,"
                                " default is >0. Threshold is >= this value"))

    optional.add_argument('--subset-proteins', dest="subset_proteins", default=10000,
                          help=("Subset to the first N proteins - useful for "
                                "debugging/checking options"))

    optional.add_argument('--proteins-of-interest', dest="proteins_of_interest",
                          default=None,
                          help=("A file containing a list of proteins of interest, "
                                "e.g all inner membrane proteins etc"))

    optional.add_argument('--filename-samplename-map', dest="filename_samplename_map",
                          default=None,
                          help=("A file mapping filenames (column 1) to sample names "
                                "(column 2), with comma-seperated columns. This will "
                                "ensure samples are appropriately named in the figures"))

    # TS:
    # set constants - turn these into options later
    TM_PAD = 10 # size of 'padding' around the TM
    FEATURE_SIZE = 10 # The number of bins per feature (upstream pad, TM, downstream pad)


    args = vars(parser.parse_args())

    if not os.path.exists(args['outdir']):
        os.mkdir(args['outdir'])


    logfile = open(args['logfile'], 'w')
    logfile.write("Logfile for plot_tm_coverage.py %s\n\n" % (
        strftime("%Y-%m-%d %H:%M:%S", gmtime())))

    section_blocker = writeSectionHeader(logfile, "Script arguments:")

    for key, value in args.items():
        logfile.write("%s: %s\n" % (key, value))
    logfile.write("%s\n\n" % section_blocker)

    tm_blocks = collections.defaultdict(list)
    
    # 1. obtain the TM regions, either using the list of proteins provided or the tax id
    if args['proteins_of_interest']:

        tm_proteins = set()

        with open(args['proteins_of_interest'], "r") as inf:
            header = next(inf) # skip the header
            for line in inf:
                tm_proteins.add(line.strip())

            for ix in range(0, min(len(tm_proteins), args['subset_proteins']), 100):
                proteins = "%2C".join(list(tm_proteins)[ix:ix+100])
                tm_url = ('https://www.ebi.ac.uk/proteins/api/features?offset=0'
                          '&size=1000&accession=%s&types=TRANSMEM' % proteins)

                r = requests.get(tm_url, headers={"Accept" : "application/json"})
                text = json.loads(r.text)
                for p_ix in text:
                    for feature in p_ix['features']:
                        if feature['type'] == 'TRANSMEM':
                            # change to zero-indexed
                            tm_blocks[p_ix['accession']].append(
                                (int(feature['begin'])-1, int(feature['end'])))

        sequence_dict = {}
        for protein, species, seq, method in sequence.getSequences(tm_proteins):
            if seq is not None:
                sequence_dict[protein] = seq

        # log the number of TM proteins
        logfile.write("Found %i proteins with TM domains in Uniprot, out of "
                      "%i proteins in --proteins-of-interest infile\n" % (
                          len(tm_blocks), len(tm_proteins)))

    else:
        if not args['tax_id']:
            raise ValueError(
                "If not providing a list of proteins with --proteins-of-interest"
                " option, must provide a taxonomy id (--tax-id) for the species "
                " the data was generated from. This will be used to identify the"
                " TM-containg proteins using Uniprot")

        tm_url = ('https://www.ebi.ac.uk/proteins/api/features?offset=0&size=%s'
                  '&taxid=%s&types=TRANSMEM&reviewed=true' % (
                      args['subset_proteins'], args['tax_id']))

        r = requests.get(tm_url, headers={ "Accept" : "application/json"})
        text = json.loads(r.text)

        for p_ix in text:
            for feature in p_ix['features']:
                if feature['type'] == 'TRANSMEM':
                    tm_blocks[p_ix['accession']].append((int(feature['begin']), int(feature['end'])))

        # log the number of TM proteins
        logfile.write("Found %i proteins with TM domains in Uniprot\n" % len(tm_blocks))

        sequence_dict = {}
        for protein, species, seq, method in sequence.getSequences(tm_blocks):
            sequence_dict[protein] = seq

    # 2. Obtain protein coverage
    protein_coverage = collections.defaultdict(
        lambda: collections.defaultdict(np.array))

    for protein in sequence_dict:
        if protein not in sequence_dict:
            print(protein)
        try:
            protein_coverage[protein] = np.zeros(
                (len(args['infiles']), len(sequence_dict[protein])))
        except:
            print(protein, sequence_dict[protein])
            raise ValueError()

    ix2study = {}
    covered_proteins = set()

    if args['filename_samplename_map']:
        file2sample = {}
        with open(args['filename_samplename_map'], "r") as inf:
            for line in inf:
                line = line.strip().split(",")
                file2sample[line[0]] = line[1].replace("\\n", "\n")
    else:
        file2sample = None

    for study_ix, infile in enumerate(args['infiles']):

        # removes everything after the last "." from the filename to make the study name
        study = ".".join(os.path.basename(infile).split(".")[:-1])
        
        # if required, replace filename with a sample name
        if file2sample:
            try:
                study = file2sample[study]
            except KeyError:
                raise KeyError(
                    "The filename %s is not found in the filename-samplename-map:"
                    "file %s \n\n\n %s" % (
                        study, args['filename_samplename_map'], file2sample))

        ix2study[study_ix] = study

        for line in iterateProteomicsSummary(open(infile), args['sep']):

            logfile.write("%s\n" % ",".join(line))

            uniprot_id, peptide, start, end = line

            if uniprot_id in sequence_dict:
                match = getPeptidePosition(peptide, sequence_dict[uniprot_id])
                if match is None:
                    logfile.write("In file: %s, peptide %s is not in protein %s" % (
                        infile, peptide, uniprot_id))
                    continue

                covered_proteins.add(uniprot_id)

                if match == "more than 1 match":
                    continue

                for aa in match:
                    try:
                        protein_coverage[uniprot_id][(study_ix, aa)] = 1
                    except:
                        print(protein_coverage[uniprot_id])
                        raise ValueError()
                        
    if len(covered_proteins) == 0:
        raise ValueError("")

    tm_blocks_region_coverage_norm = getTMCoverage(
        protein_coverage, tm_blocks, FEATURE_SIZE, sequence_dict)

    tms = 0
    for protein in tm_blocks_region_coverage_norm:
        for tm in tm_blocks_region_coverage_norm[protein]:
            tms += 1

    coverage_profile_df = makeMetaTMDF(
        tm_blocks_region_coverage_norm, ix2study, tms, FEATURE_SIZE)

    plotCoverage(coverage_profile_df,
                 os.path.join(args['outdir'], "metaplot_fraction_coverage.png"),
                 os.path.join(args['outdir'], "metaplot_normalised_coverage.png"))

    tm_coverage_df = makeTMCoveredDF(protein_coverage, tm_blocks, ix2study)

    final_tm_coverage_df = finaliseDataFrame(
        tm_coverage_df, coverage_threshold=args['coverage_thresh'])
    
    covered_proteins = final_tm_coverage_df[final_tm_coverage_df['coverage']==True]

    unique_proteins = pd.DataFrame(covered_proteins.groupby("uniprot_id").aggregate(
        {"desc" : len}))

    # 2 because the dataframe inc. "combined" sample
    unique_proteins = set(unique_proteins[unique_proteins['desc']==2].index.tolist())

    covered_proteins_tally = pd.DataFrame(covered_proteins.groupby("desc").aggregate(
        {"uniprot_id" : len}))

    unique_proteins_df = covered_proteins[
        covered_proteins['uniprot_id'].isin(unique_proteins)]

    unique_proteins_tally = pd.DataFrame(unique_proteins_df.groupby("desc").aggregate(
        {"uniprot_id" : len}))

    covered_proteins_tally.to_csv(
        os.path.join(args['outdir'], "covered_proteins_count.tsv"), sep="\t")
    unique_proteins_tally.to_csv(
        os.path.join(args['outdir'], "unique_proteins_count.tsv"), sep="\t")
    covered_proteins.to_csv(
        os.path.join(args['outdir'], "covered_proteins.tsv"), sep="\t")
    unique_proteins_df.to_csv(
        os.path.join(args['outdir'], "unique_proteins.tsv"), sep="\t")

    plotCoveredTM(final_tm_coverage_df,
                  os.path.join(args['outdir'], "circ_covered.png"),
                  circular=1, plot_only_covered=1,
                  colour_by_group=1, group="desc")

    plotCoveredTM(final_tm_coverage_df,
                  os.path.join(args['outdir'], "circ_all.png"),
                  circular=1, plot_only_covered=0,
                  colour_by_group=1, group="desc")

    

if __name__ == "__main__":
    sys.exit(main(sys.argv))
