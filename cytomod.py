#!/usr/bin/env python
from __future__ import with_statement, division

__version__ = "$Revision: 0.01$"

"""Cytomod uses information on modification
locations to replace the appropriate symbols in a reference genome sequence
with the new symbols. Cytomod can incorporate data from both
single-base assays (point annotations in BED format)
and lower-resolution assays (region annotations in BED or WIG formats).
The output from this program is intended to allow for de novo
discovery of modification-sensitive motifs in known TFBSs.
This program currently handles the following modifications
to the cytosine nucleobase: {5mC, 5hmC, 5fC, 5caC}.
"""

import warnings
import sys
import os
import glob
import datetime
import re

import numpy as np

DEFAULT_FASTA_FILENAME = 'modGenome.fa'

MOD_BASES = {'5mC': 'm', '5hmC': 'h', '5fC': 'f', '5caC': 'c'}

SUPPORTED_FILE_FORMATS_REGEX = '\.(bed|wig|bed[gG]raph)$'
CHROMOSOME_EXCLUSION_REGEX = 'random'
MOD_BASE_REGEX = '5.+C'
REGION_REGEX = '(chr\d+):(\d+)-(\d+)'

_MAX_REGION_AT_ONCE = 1000000

# XXX parameterize and/or automate
modOrder = np.array([3, 2, 0, 1])


def v_print_timestamp(msg="", threshold=1):
    """Print a timestamped message iff verbosity is at least threshold."""
    sys.stderr.write(">> <Cytomod> %s: %s" % (
        datetime.datetime.now().isoformat(), msg + "\n")
        if args.verbose >= threshold else "")


def ensureRegionValidity(genome, chr, start, end):
    """Ensures the validity of the given region. Dies if not valid."""
    chromosome = genome[chr]
    try:
        chromosome = genome[chr]
    except:
        sys.exit("Invalid region: invalid chromosme.")
    if (chromosome.start < 0) or (chromosome.start >= end):
        sys.exit("Invalid region: invalid start position.")
    if (end <= start) or (end > chromosome.end):
        sys.exit("Invalid region: invalid end position.")


def getModifiedGenome(genome, chr, start, end):
    """Returns the modified genome sequence, for the given genome,
    over the given input region."""
    chromosome = genome[chr]
    allbasesResult = ""
    # Only compute the modified genome in segments.
    # This prevents the creation of excessively large NumPy arrays.
    for s in range(int(chromosome.start), int(chromosome.end),
                   _MAX_REGION_AT_ONCE):
        e = s + _MAX_REGION_AT_ONCE
        v_print_timestamp("Now outputting " + chr + " for region: (" + str(s)
                          + ", " + str(e) + ")", 2)
        modBasesA = np.where(np.logical_and(np.isfinite(chromosome[s:e]),
                             chromosome[s:e] != 0), modBases, '0')
        orderedmodBasesA = modBasesA[:, modOrder]
        # XXX This re-creates the entire array.
        # Finding a more efficient way would be preferred.
        orderedmodBasesA = np.column_stack((
            orderedmodBasesA, list(chromosome.seq[s:e].tostring().upper())))
        allbases = [filter(lambda x: x != '0', (bases))[0]
                    for bases in orderedmodBasesA if np.any(bases != '0')]
        allbasesResult += ''.join(allbases)
    return allbasesResult


def generateFASTAFile(file, id, genome, chr, start, end):
    """Writes a FASTA file of the modified genome appending to the given file,
    using the given ID.
    No FASTA ID (i.e. '> ...') is written if no ID is given."""
    modGenomeFile = open(file, 'a')
    if id:
        modGenomeFile.write(">" + id + "\n")
    modGenomeFile.write(getModifiedGenome(genome, chr, start, end) + "\n")
    modGenomeFile.close()


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
parser.add_argument("-f", "--fastaFile", default=DEFAULT_FASTA_FILENAME, help="Provide a full path \
                    to a file to append the modified genome in FASTA format.")
parser.add_argument("-r", "--region", help="Only output the modified genome for the given region. \
                    The output will be printed to STDOUT unless a FASTA file \
                    name is provided via \"-f\", in which case it will be \
                    written to that file. The region must be specified in \
                    the format: chr<ID>:<start>-<end> (ex. chr1:500-510)")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

from genomedata import Genome, load_genomedata

genomeDataArchive = ""
if args.archiveCompDirs:
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
    warnings.simplefilter("ignore")  # Ignore supercontig warnings
    v_print_timestamp("Genomedata archive successfully loaded.")

    modBases = []
    for track in genome.tracknames_continuous:
        modBases.append(MOD_BASES[re.search(MOD_BASE_REGEX, track).group(0)])
    v_print_timestamp("The order of preference for base modifications is: "
                      + ','.join(modBases) + ".")

    if args.region:
        v_print_timestamp("Outputting the modified genome for: "
                          + args.region + ".")
        regionMatch = re.search(REGION_REGEX, args.region)
        if args.fastaFile != DEFAULT_FASTA_FILENAME:
            ensureRegionValidity(genome, regionMatch.group(1),
                                 int(regionMatch.group(2)),
                                 int(regionMatch.group(3)))
            generateFASTAFile(args.fastaFile, args.region, genome,
                              regionMatch.group(1), int(regionMatch.group(2)),
                              int(regionMatch.group(3)))
        else:
            print getModifiedGenome(genome, regionMatch.group(1),
                                    int(regionMatch.group(2)),
                                    int(regionMatch.group(3)))
    else:
        for chromosome in [chromosome for chromosome in genome
                           if not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                            chromosome.name)]:
            v_print_timestamp("Outputting the modified genome for: "
                              + chromosome.name)
            generateFASTAFile(args.fastaFile, chromosome.name,
                              genome, chromosome.name,
                              int(chromosome.start), int(chromosome.end))

v_print_timestamp("Program complete.")
