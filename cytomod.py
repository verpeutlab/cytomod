#!/usr/bin/env python
from __future__ import with_statement, division, print_function

__version__ = "$Revision: 0.02$"

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
from collections import OrderedDict
from itertools import izip

import numpy as np

# Define the "primary" modified bases and their corresponding
# one base codes, listed in their order of oxidation.
MOD_BASES = OrderedDict([('5mC', 'm'), ('5hmC', 'h'),
                        ('5fC', 'f'), ('5caC', 'c')])
# Create a dictionary mapping each "primary" modified base to
# the base it modifies.
_MODIFIES = dict.fromkeys(MOD_BASES.values(), 'C')

# All IUPAC nucleobases and their complements, plus 'X',
# which is just an additional alias for any nucleobase.
# TODO Remove redudancy by building from pairs.
COMPLEMENTS = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
               'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M',
               'W': 'W', 'S': 'S', 'B': 'V', 'D': 'H',
               'H': 'D', 'V': 'B', 'N': 'N', 'X': 'X'}
# Add all modified nucleobases.
modifiedBasesToComplements = \
    izip(MOD_BASES.values(), ''.
         join(str(i) for i in range(1, len(MOD_BASES) + 1)))
COMPLEMENTS.update(modifiedBasesToComplements)
# Add all modified nucleobase complements.
# They are ordered by the originating modification's oxidation order.
COMPLEMENTS.update([(c, m) for m, c in modifiedBasesToComplements])


# TODO use a decorator to properly return a scalar for scalar input
# and a list for list input.
def complement(bases):
    """Complements the given, potentially modified, base."""
    return [COMPLEMENTS[base] for base in bases]

# Update the dictionary mapping with every complemented modification.
_MODIFIES.update(izip(complement(_MODIFIES.keys()),
                 complement(_MODIFIES.values())))

AUTOSOME_ONLY_FLAG = 'u'
ALLOSOME_ONLY_FLAG = 'l'
MITOCHONDRIAL_ONLY_FLAG = 'm'
MITOCHONDRIAL_EXCLUSION_FLAG = MITOCHONDRIAL_ONLY_FLAG.upper()

SUPPORTED_FILE_FORMATS_REGEX = '\.(bed|wig|bed[gG]raph)$'
CHROMOSOME_TYPE_REGEXES = {AUTOSOME_ONLY_FLAG: 'chr\d+',
                           ALLOSOME_ONLY_FLAG: 'chr[XY]',
                           MITOCHONDRIAL_ONLY_FLAG: 'chrM',
                           MITOCHONDRIAL_EXCLUSION_FLAG: 'chr(?:\d+|[XY])'}
CHROMOSOME_EXCLUSION_REGEX = '(?:random)'
MOD_BASE_REGEX = '5.+C'
REGION_REGEX = '(chr(?:\d+|[XYM]))(?::(?P<start>\d+)?-(?P<end>\d+)?)?'

_MSG_PREFIX = '>> <Cytomod> '
_DEFAULT_FASTA_FILENAME = 'modGenome.fa'
_DEFAULT_BASE_PRIORITY = 'hmfc'
_DEFAULT_BASE_PRIORITY_COMMENT = """the resolution of the biological protocol
(i.e. single-base > any chemical > any DIP)"""
_DEFAULT_RAN_LENGTH = 2000
_MAX_REGION_LEN = 2000000
_MAX_CONTIG_ATTEMPTS = 3


def warn(*msg):
    """Emit a warning to STDERR."""
    print("Warning: ", *msg, file=sys.stderr)


def v_print_timestamp(msg="", threshold=1):
    """Print a timestamped message iff verbosity is at least threshold."""
    sys.stderr.write(_MSG_PREFIX + "%s: %s" % (
        datetime.datetime.now().isoformat(), msg + "\n")
        if args.verbose >= threshold else "")


def _modifychrmExclusionRegex(additionalchrmExclusionFlag):
    """Modify the chrmomosome exclusion regex,
    accroding to the provided flags."""
    global CHROMOSOME_EXCLUSION_REGEX  # allow modification of global var
    # Modify the exclusion regex by adding the regex corresponding
    # to the flag that we wish to exclude. However, the dictionary
    # containing the regexes identify the group specified by the flag
    # (i.e. are inclusion regexes). We therefore invert the additional
    # regex via a modified anchored negative lookahead.
    CHROMOSOME_EXCLUSION_REGEX += '|(^((?!' + \
        CHROMOSOME_EXCLUSION_REGEX[additionalchrmExclusionFlag] + ').)*$)'


def _ensureRegionValidity(genome, chrm, start, end):
    """Ensures the validity of the given region. Dies if not valid."""
    chrmomosome = genome[chrm]
    try:
        chrmomosome = genome[chrm]
    except:
        sys.exit("Invalid region: invalid chrmomosme.")
    if (chrmomosome.start < 0) or (chrmomosome.start >= end):
        sys.exit("Invalid region: invalid start position.")
    if (end <= start) or (end > chrmomosome.end):
        sys.exit("Invalid region: invalid end position.")


def getModifiedGenome(genome, modOrder, chrm, start, end):
    """Returns the modified genome sequence, for the given genome,
    over the given input region."""
    chrmomosome = genome[chrm]
    allbasesResult = ""
    # Only compute the modified genome in segments.
    # This prevents the creation of excessively large NumPy arrays.
    for s in range(start, end, _MAX_REGION_LEN):
        e = s + _MAX_REGION_LEN if (end - start) >= _MAX_REGION_LEN else end
        v_print_timestamp("Now outputting " + chrm + " for region: (" + str(s)
                          + ", " + str(e) + ")", 2)
        modBaseScores = chrmomosome[s:e]
        modBasesA = np.where(np.logical_and(np.isfinite(modBaseScores),
                             modBaseScores != 0), modBases, '0')
        orderedmodBasesA = modBasesA[:, modOrder]
        referenceSeq = list(chrmomosome.seq[s:e].tostring().upper())
        # Filter the bases to take the modified bases in priority order.
        x = np.transpose(np.nonzero(orderedmodBasesA != '0'))
        u, idx = np.unique(x[:, 0], return_index=True)

        def maybeGetModBase(m, r):
            """Returns the modified base corresponding to
            the given putatively modified base. The base returned
            is the input putatively modified base if the corresponding
            reference base is modifiable to the input base, or the complement
            of that base, if the complemented reference is modifiable to it,
            otherwise the reference base itself is returned."""
            if m not in _MODIFIES:
                return r
            else:
                if _MODIFIES[m] == r:
                    return m
                elif _MODIFIES[m] == complement(r)[0]:
                    return complement(m)[0]
                else:
                    return r

        # Initially the sequence is unmodified and we successively modify it.
        allModBases = np.array(referenceSeq, dtype=np.str)
        # We vectorize the function for convinience.
        # NumPy vectorized functions still execute the Python code at
        maybeGetModBase = np.vectorize(maybeGetModBase)
        # Mask the sequence, allowing only base modifications
        # that modify their 'target' base (i.e. '5fC' = 'f' only modifies 'C').
        # Return the reference base for all non-modifiable bases
        # and for unmodified bases.
        if (x.size > 0):
            np.put(allModBases, x[idx][:, 0],
                   maybeGetModBase(orderedmodBasesA[x[idx][:, 0],
                                   x[idx][:, 1]],
                                   allModBases[x[idx][:, 0]]))
        # Output the unmodified sequence at a verbosity level of at least 2,
        # if not too long, otherwise only output for a high verbosity level.
        v_print_timestamp("Corresponding unmodified reference sequence: \n"
                          + ''.join(referenceSeq), 2
                          if len(referenceSeq) < 10000 else 6)
        # Concatenate the vector together to form the (string) sequence
        allbasesResult += ''.join(allModBases)
    return allbasesResult


def generateFASTAFile(file, id, genome, modOrder, chrm, start, end):
    """Writes a FASTA file of the modified genome appending to the given file,
    using the given ID.
    No FASTA ID (i.e. '> ...') is written if no ID is given."""
    with open(file, 'a') as modGenomeFile:
        if id:
            modGenomeFile.write(">" + id + "\n")
        modGenomeFile.write(getModifiedGenome(genome, modOrder, chrm,
                            start, end) + "\n")


def selectRandomRegion(genome, length):
    """Selects a random, non-exluded, chrmomosome of sufficient length.
    This method attempts to ensure that the region selected is
    wholly within a supercontig."""
    selectablechrmomosomes = {chrmomosome.name: chrmomosome.end for
                              chrmomosome in genome if
                              (chrmomosome.end >= length and
                               not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                             chrmomosome.name))}
    if not selectablechrmomosomes:
        sys.exit(("The region length provided is too long or all "
                 "chrmomosomes have been excluded."))
    chrm = random.choice(selectablechrmomosomes.keys())
    contigAttempts = 0
    while True:
        start = random.randint(genome[chrm].start, genome[chrm].end)
        end = start + length
        if genome[chrm].supercontigs[start:end]:
            break
        elif contigAttempts >= _MAX_CONTIG_ATTEMPTS:
            warn("Attempts to procure sequence from a supercontig "
                 "were exhausted. Returning sequence that is not "
                 "wholly contained within a supercontig.")
            break
        contigAttempts += 1
    return chrm, start, end


def parseRegion(genome, region):
    """Parses the provided region, ensuring its validity."""
    regionMatch = re.search(REGION_REGEX, region)
    chrm = regionMatch.group(1)
    start = 0
    end = 1
    if regionMatch.group('start'):
        start = int(regionMatch.group('start'))
    else:
        start = genome[chrm].start
    if regionMatch.group('end'):
        end = int(regionMatch.group('end'))
    else:
        end = genome[chrm].end
    _ensureRegionValidity(genome, chrm, start, end)
    return chrm, start, end


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
    # randomly sample 1000 bases of first chrmomosome to determine resolution
    for chrmomosome in genome:
        testRegion = random.randint(chrmomosome.start, chrmomosome.end)
        print(chrmomosome[testRegion:testRegion + 1000])
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
                           chrmomosomes and/or scaffolds. \
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
                    chrm<ID>:<start>-<end> (ex. chrm1:500-510).")
region.add_argument('-R', '--randomRegion', nargs='?',
                    const=_DEFAULT_RAN_LENGTH, type=int,
                    help="Output the modified genome for a random region. \
                    The chrmomsome will be randomly selected and its \
                    coordinate space will be randomly and uniformly sampled. \
                    A length for the random region can either be specified \
                    or it will otherwise be set to a reasonably \
                    small default. The length chosen may constrain the \
                    selection of a chrmomosome.")
parser.add_argument('-E', '--excludechrms',
                    choices=CHROMOSOME_TYPE_REGEXES.keys(),
                    help="Exclude chrmomosome \
                    types. '" + AUTOSOME_ONLY_FLAG + "': \
                    Use only autosomal chrmomosomes  (excludes chrmM). \
                    '" + ALLOSOME_ONLY_FLAG + "': \
                    Use only allosomal chrmomosomes (excludes chrmM). \
                    '" + MITOCHONDRIAL_ONLY_FLAG + "': \
                    Use only the mitochondrial chrmomosome. \
                    '" + MITOCHONDRIAL_EXCLUSION_FLAG + "': \
                    Exclude the mitochondrial chrmomosome. \
                    NB: This paprameter will be ignored if a \
                    specific genomic region is queried \
                    via '-r'.")
parser.add_argument('-p', '--priority', default=_DEFAULT_BASE_PRIORITY,
                    choices=MOD_BASES.values(),
                    help="Specify the priority \
                    of modified bases. The default is:"
                    + _DEFAULT_BASE_PRIORITY + ", which is based upon"
                    + _DEFAULT_BASE_PRIORITY_COMMENT + ".")
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action="count")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

if args.region and args.excludechrms:
    warn("Exclusion regex ignored, since a specific region was specifed.")

if args.excludechrms:
    _modifychrmExclusionRegex(args.excludechrms)

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
    modOrder = [modBases.index(b) for b in list(args.priority)]
    v_print_timestamp("The order of preference for base modifications is: "
                      + ','.join(list(args.priority)) + ".")

    if args.region or args.randomRegion:
        if args.randomRegion:
            chrm, start, end = selectRandomRegion(genome, args.randomRegion)
        else:
            chrm, start, end = parseRegion(genome, args.region)
        regionStr = chrm + ":" + str(start) + "-" + str(end)
        v_print_timestamp("Outputting the modified genome for: "
                          + regionStr + ".")
        if args.fastaFile:
            generateFASTAFile(args.fastaFile, regionStr, genome, modOrder,
                              chrm, start, end)
        else:
            print(getModifiedGenome(genome, modOrder, chrm, start, end))
    else:
        for chrmomosome in [chrmomosome for chrmomosome in genome
                            if not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                             chrmomosome.name)]:
            v_print_timestamp("Outputting the modified genome for: "
                              + chrmomosome.name)
            generateFASTAFile(args.fastaFile, chrmomosome.name,
                              genome, modOrder, chrmomosome.name,
                              int(chrmomosome.start), int(chrmomosome.end))

v_print_timestamp("Program complete.")
