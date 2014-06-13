#!/usr/bin/env python
from __future__ import with_statement, division

__version__ = "$Revision: 0.01$"

"""Cytomod uses information on cytosine modification
locations to replace C symbols in a reference genome sequence
with the new symbols. Cytomod can incorporate data from both
single-base assays (point annotations in BED format)
and lower-resolution assays (region annotations in BED or WIG formats).
The output from this program is intended to allow for de novo
discovery of modification-sensitive motifs in known TFBSs.
"""

import warnings
import sys
import os
import glob
import datetime
import re

import numpy as np

MOD_BASES = {'5mC': 'm', '5hmC': 'h', '5fC': 'f', '5caC': 'c'}

SUPPORTED_FILE_FORMATS_REGEX = "\.(bed|wig|bed[gG]raph)$"
CHROMOSOME_EXCLUSION_REGEX = "random"
MOD_BASE_REGEX = "5.+C"

# XXX parameterize
modOrder = np.array([3, 2, 0, 1])


def v_print_timestamp(msg=""):
    """Print a timestamped message iff non-zero verbosity is specified."""
    sys.stderr.write(">> <Cytomod> %s: %s" % (
        datetime.datetime.now().isoformat(), msg + "\n")
        if args.verbose else "")

import argparse
parser = argparse.ArgumentParser()

genomeArchive = parser.add_mutually_exclusive_group(required=True)
genomeArchive.add_argument("-G", "--genomedataArchive",
                           help="The genome data archive. \
                           It must contain all needed \
                           sequence and track files. \
                           If one is not yet created, \
                           use \"-g\" and \"-t\" instead to create it.")
genomeArchive.add_argument("-d", "--archiveCompDirs", nargs=2,
                           help="Two arguments first specifying the directory containing \
                           the genome and then the directory containing all \
                           modified base tracks. The genome directory must \
                           contain (optionally gzipped) FASTA files of \
                           chromosomes and/or scaffolds. \
                           The track directory must contain \
                           (optionally gzipped) \
                           genome tracks. \
                           They must have an extension describing \
                           their format. We currently support: \
                           \".wig\", \".bed\", \
                           and \".bedGraph\". The filename of each track \
                           must specify what modified nucleobase it \
                           pertains to (i.e. \"5hmC\"). \
                           Ensure that all tracks are mapped to the same \
                           assembly and that this assembly matches the \
                           genome provided. This will \
                           create a genome data archive in an \"archive\". \
                           subdirectory of the provided track directory. \
                           Use \"-G\" instead to use an existing archive.")
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="count")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

from genomedata import Genome, load_genomedata

genomeDataArchive = ""
if (args.archiveCompDirs):
    v_print_timestamp("Creating genomedata archive.")
    genomeDataArchive = args.archiveCompDirs[1] + "/archive/"
    # Create the genome data archive
    # Load all supported track files in the tracks directory
    # Load all FASTA files in the sequences directory
    load_genomedata.load_genomedata(
        genomeDataArchive,
        tracks=[(track, args.archiveCompDirs[1] + track)
                for track in os.listdir(args.archiveCompDirs[1])
                if re.search(SUPPORTED_FILE_FORMATS_REGEX, track)],
        seqfilenames=glob.glob(args.archiveCompDirs[0] + "/*.fa*"),
        verbose=args.verbose)

else:
    v_print_timestamp("Using existing genomedata archive.")
    genomeDataArchive = args.genomedataArchive

with Genome(genomeDataArchive) as genome:
    v_print_timestamp("Genomedata archive successfully loaded.")
    warnings.simplefilter("ignore")  # Ignore supercontig warnings

    modBases = []
    for track in genome.tracknames_continuous:
        modBases.append(MOD_BASES[re.search(MOD_BASE_REGEX, track).group(0)])
    v_print_timestamp("The order of preference for base modifications is: "
                      + ','.join(modBases) + ".")

    # XXX parameterize
    start = 90000040
    end = 90000050

    for chromosome in [chromosome for chromosome in genome
                       if not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                        chromosome.name)]:
        modBasesA = np.where(np.isfinite(chromosome[start:end]), modBases, '0')
        orderedmodBasesA = modBasesA[:, modOrder]
        orderedmodBasesA = np.column_stack((
            orderedmodBasesA, list(chromosome.seq[start:end].tostring())))
        allbases = [filter(lambda x: x != '0', (bases))[0]
                    for bases in orderedmodBasesA if np.any(bases != '0')]
        print ''.join(allbases)
