#!/usr/bin/env python
from __future__ import with_statement, division, print_function

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
import random

import numpy as np

MOD_BASES = {'5mC': 'm', '5hmC': 'h', '5fC': 'f', '5caC': 'c'}

AUTOSOME_ONLY_FLAG = 'u'
ALLOSOME_ONLY_FLAG = 'l'
MITOCHONDRIAL_ONLY_FLAG = 'm'
MITOCHONDRIAL_EXCLUSION_FLAG = MITOCHONDRIAL_ONLY_FLAG.upper()

SUPPORTED_FILE_FORMATS_REGEX = '\.(bed|wig|bed[gG]raph)$'
CHROMSOME_TYPE_REGEXES = {AUTOSOME_ONLY_FLAG: 'chr\d+',
                          ALLOSOME_ONLY_FLAG: 'chr[XY]',
                          MITOCHONDRIAL_ONLY_FLAG: 'chrM',
                          MITOCHONDRIAL_EXCLUSION_FLAG: 'chr(?:\d+|[XY])'}
CHROMOSOME_EXCLUSION_REGEX = '(?:random)'
MOD_BASE_REGEX = '5.+C'
REGION_REGEX = '(chr(?:\d+|[XYM]))(?::(?P<start>\d+)?-(?P<end>\d+)?)?'

_DEFAULT_FASTA_FILENAME = 'modGenome.fa'
_DEFAULT_RAN_LENGTH = 2000
_MAX_REGION_LEN = 2000000
_MAX_CONTIG_ATTEMPTS = 3

# XXX parameterize and/or automate
modOrder = np.array([3, 2, 0, 1])


def warn(*msg):
    """Emit a warning to STDERR."""
    print("Warning: ", *msg, file=sys.stderr)


def v_print_timestamp(msg="", threshold=1):
    """Print a timestamped message iff verbosity is at least threshold."""
    sys.stderr.write(">> <Cytomod> %s: %s" % (
        datetime.datetime.now().isoformat(), msg + "\n")
        if args.verbose >= threshold else "")


def _modifyChrExclusionRegex(additionalChrExclusionFlag):
    """Modify the chromosome exclusion regex accroding to the provided flags"""
    global CHROMOSOME_EXCLUSION_REGEX  # allow modification of global var
    # Modify the exclusion regex by adding the regex corresponding
    # to the flag that we wish to exclude. However, the dictionary
    # containing the regexes identify the group specified by the flag
    # (i.e. are inclusion regexes). We therefore invert the additional
    # regex via a modified anchored negative lookahead.
    CHROMOSOME_EXCLUSION_REGEX += '|(^((?!' + \
        CHROMSOME_TYPE_REGEXES[additionalChrExclusionFlag] + ').)*$)'


def _ensureRegionValidity(genome, chr, start, end):
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
    for s in range(start, end, _MAX_REGION_LEN):
        e = s + _MAX_REGION_LEN if (end - start) >= _MAX_REGION_LEN else end
        v_print_timestamp("Now outputting " + chr + " for region: (" + str(s)
                          + ", " + str(e) + ")", 2)
        modBasesA = np.where(np.logical_and(np.isfinite(chromosome[s:e]),
                             chromosome[s:e] != 0), modBases, '0')
        orderedmodBasesA = modBasesA[:, modOrder]
        # TODO This re-creates the entire array.
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
    with open(file, 'a') as modGenomeFile:
        if id:
            modGenomeFile.write(">" + id + "\n")
        modGenomeFile.write(getModifiedGenome(genome, chr, start, end) + "\n")


def selectRandomRegion(genome, length):
    """Selects a random, non-exluded, chromosome of sufficient length.
    This method attempts to ensure that the region selected is
    wholly within a supercontig."""
    selectableChromosomes = {chromosome.name: chromosome.end for
                             chromosome in genome if
                             (chromosome.end >= length and
                              not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                            chromosome.name))}
    if not selectableChromosomes:
        sys.exit(("The region length provided is too long or all "
                 "chromosomes have been excluded."))
    chr = random.choice(selectableChromosomes.keys())
    contigAttempts = 0
    while True:
        start = random.randint(genome[chr].start, genome[chr].end)
        end = start + length
        if genome[chr].supercontigs[start:end]:
            break
        elif contigAttempts >= _MAX_CONTIG_ATTEMPTS:
            warn("Attempts to procure sequence from a supercontig "
                 "were exhausted. Returning sequence that is not "
                 "wholly contained within a supercontig.")
            break
        contigAttempts += 1
    return chr, start, end


def parseRegion(genome, region):
    """Parses the provided region, ensuring its validity."""
    regionMatch = re.search(REGION_REGEX, region)
    chr = regionMatch.group(1)
    start = 0
    end = 1
    if regionMatch.group('start'):
        start = int(regionMatch.group('start'))
    else:
        start = genome[chr].start
    if regionMatch.group('end'):
        end = int(regionMatch.group('end'))
    else:
        end = genome[chr].end
    _ensureRegionValidity(genome, chr, start, end)
    return chr, start, end


def determineTrackPriority(genome):
    """Currently, an ad hoc and contrived means of determining
    which epigenetic modification has precedence. This is done by
    naively considering the resolution, and secondarily, the frequency
    of a given type of modification.
    NB: This method is not yet fully implemented."""
    # XXX cannot complete this without being able to
    # access the intervals over which the track data is defined.
    # Genomedata does not appear to support this.
    # TODO Ameliorate this.
    print(genome.num_datapoints)  # XXX
    # randomly sample 1000 bases of first chromosome to determine resolution
    for chromosome in genome:
        testRegion = random.randint(chromosome.start, chromosome.end)
        print(chromosome[testRegion:testRegion + 1000])
        break

import argparse
parser = argparse.ArgumentParser()

genomeArchive = parser.add_mutually_exclusive_group(required=True)
genomeArchive.add_argument('-G', '--genomedataArchive',
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
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action="count")
parser.add_argument('-f', '--fastaFile', nargs='?', type=str,
                    const=_DEFAULT_FASTA_FILENAME, help="Output to a file instead \
                    of STDOUT. Provide a full path to a file to append the \
                    modified genome in FASTA format. If this parameter is \
                    invoked without any arguments, a default filename will \
                    be used within the current directory.")
region = parser.add_mutually_exclusive_group()
region.add_argument('-r', '--region', help="Only output the modified genome \
                    for the given region. \
                    The region must be specified in the format: \
                    chr<ID>:<start>-<end> (ex. chr1:500-510).")
region.add_argument('-R', '--randomRegion', nargs='?',
                    const=_DEFAULT_RAN_LENGTH, type=int,
                    help="Output the modified genome for a random region. \
                    The chromsome will be randomly selected and its \
                    coordinate space will be randomly and uniformly sampled. \
                    A length for the random region can either be specified \
                    or it will otherwise be set to a reasonably \
                    small default. The length chosen may constrain the \
                    selection of a chromosome.")
parser.add_argument('-E', '--excludeChrs',
                    choices=CHROMSOME_TYPE_REGEXES.keys(),
                    help="Exclude chromosome \
                    types. '" + AUTOSOME_ONLY_FLAG + "': \
                    Use only autosomal chromosomes  (excludes chrM). \
                    '" + ALLOSOME_ONLY_FLAG + "': \
                    Use only allosomal chromosomes (excludes chrM). \
                    '" + MITOCHONDRIAL_ONLY_FLAG + "': \
                    Use only the mitochondrial chromosome. \
                    '" + MITOCHONDRIAL_EXCLUSION_FLAG + "': \
                    Exclude the mitochondrial chromosome. \
                    NB: This paprameter will be ignored if a \
                    specific genomic region is queried \
                    via '-r'.")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

if args.region and args.excludeChrs:
    warn("Exclusion regex ignored, since a specific region was specifed.")

if args.excludeChrs:
    _modifyChrExclusionRegex(args.excludeChrs)

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

    if args.region or args.randomRegion:
        if args.randomRegion:
            chr, start, end = selectRandomRegion(genome, args.randomRegion)
        else:
            chr, start, end = parseRegion(genome, args.region)
        regionStr = chr + ":" + str(start) + "-" + str(end)
        v_print_timestamp("Outputting the modified genome for: "
                          + regionStr + ".")
        if args.fastaFile:
            generateFASTAFile(args.fastaFile, regionStr, genome,
                              chr, start, end)
        else:
            print(getModifiedGenome(genome, chr, start, end))
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
