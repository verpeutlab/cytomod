#!/usr/bin/env python

from __future__ import with_statement, division, print_function

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

import argparse
import glob
import gzip
import math
import os
import random
import re
import warnings

from collections import OrderedDict

import numpy as np

import cytoUtils as cUtils
from cytoUtils import duplicates, indices

__version__ = cUtils.__version__

_AUTOSOME_ONLY_FLAG = 'u'
_ALLOSOME_ONLY_FLAG = 'l'
_MITOCHONDRIAL_ONLY_FLAG = 'm'
_MITOCHONDRIAL_INCLUSION_FLAG = _MITOCHONDRIAL_ONLY_FLAG.upper()

SUPPORTED_FORMATS_REGEX = '\.(bed|wig|bed[gG]raph)(.gz)?$'
CHROMOSOME_TYPE_REGEXES = {_AUTOSOME_ONLY_FLAG: 'chr\d+',
                           _ALLOSOME_ONLY_FLAG: 'chr[XY]',
                           _MITOCHONDRIAL_ONLY_FLAG: 'chrM',
                           _MITOCHONDRIAL_INCLUSION_FLAG: 'chr(?:\d+|[XYM])'}
# default chromosomal exclusions include: unmapped data, haplotypes, and chrM
CHROMOSOME_EXCLUSION_REGEX = '(?:random|hap|chrM)'
MOD_BASE_REGEX = '5(m|hm|f|ca)C'
REGION_REGEX = '(chr(?:\d+|[XYM]))(?::(?P<start>\d+)?-(?P<end>\d+)?)?'

_DEFAULT_ARCHIVE_NAME = 'archive'
_DEFAULT_FASTA_FILENAME = 'modGenome.fa'
# bases ordered by least frequent in our dataset with most ambiguous bases last
_DEFAULT_BASE_PRIORITY = 'fhmc' + ''.join(cUtils.AMBIG_MOD_BASES.keys()[::-1])
_DEFAULT_CENTRED_REGION_LENGTH = 500
_DEFAULT_BASE_PRIORITY_COMMENT = """the resolution of the biological protocol
(i.e. single-base > any chemical > any DIP)"""
_DEFAULT_RAN_LENGTH = 2000
_DEFAULT_MASK_VALUE = 0
_MAX_REGION_LEN = 2000000
_MAX_CONTIG_ATTEMPTS = 3
# tracks to be displayed densely for output UCSC browser tracks
_DENSE_TRACKS = 'ruler ensGene pubs cpgIslandExt oreganno rmsk snp128'
_MASK_TNAME = 'MASK'


def die(msg):
    cUtils.die(msg, os.path.basename(__file__))


def warn(msg):
    cUtils.warn(msg, os.path.basename(__file__))


def v_print_timestamp(verbosity, msg="", threshold=1):
    cUtils.v_print_timestamp(verbosity, msg, threshold,
                             os.path.basename(__file__))


def _modifychrmExclusionRegex(additionalchrmExclusionFlag):
    """Modify the chromosome exclusion regex,
    accroding to the provided flags."""
    global CHROMOSOME_EXCLUSION_REGEX  # allow modification of global var
    # Modify the exclusion regex by adding the regex corresponding
    # to the flag that we wish to exclude. However, the dictionary
    # containing the regexes identify the group specified by the flag
    # (i.e. are inclusion regexes). We therefore invert the additional
    # regex via a modified anchored negative lookahead.
    CHROMOSOME_EXCLUSION_REGEX += '|(^((?!' + \
        CHROMOSOME_TYPE_REGEXES[additionalchrmExclusionFlag] + ').)*$)'


def _ensureRegionValidity(genome, chrm, start, end):
    """Ensures the validity of the given region. Dies if not valid."""
    chromosome = genome[chrm]
    try:
        chromosome = genome[chrm]
    except:
        die("Invalid region: invalid chrmomosme.")
    if (chromosome.start < 0) or (chromosome.start >= end):
        die("Invalid region: invalid start position.")
    if (end <= start) or (end > chromosome.end):
        die("Invalid region: invalid end position.")


def _maybeGetAmbigMapping(base, ambigMap):
    """Maps the given base to its ambiguous base, if possible."""
    return (ambigMap.get(base) or
            (cUtils.complement(ambigMap.get(cUtils.complement(base)))
            if ambigMap.get(cUtils.complement(base)) else None) or base)


def getTrackHeader(modBase):
    """Generates and returns a valid UCSC track header,
    with an appropriate name, description, and colour
    for the the given modified nucleobase.
    The returned header is terminated by a newline.
    """
    colour = str(cUtils.getRGB256BaseCol(modBase)).translate(None, '() ')
    browserConfLines = "browser hide all\nbrowser dense " + \
        _DENSE_TRACKS + "\n"
    return browserConfLines + 'track name="Nucleobase ' + modBase + \
        '" description="Track denoting ' + cUtils.FULL_MOD_BASE_NAMES[modBase] + \
        ' covalently modified nucleobases.' + '" color=' + \
        re.sub('[() ]', '', colour + "\n")


def getModifiedGenome(genome, modOrder, chrm, start, end,
                      suppressFASTA, suppressBED, tnames, ambigMap,
                      maskRegionsFileVal, maskRegionTName,
                      maskAllUnsetRegions):
    """Returns the modified genome sequence, for the given genome,
    over the given input region.
    """
    hasModifiedBases = False
    chromosome = genome[chrm]
    allbasesResult = ""

    # Only compute the modified genome in segments.
    # This prevents the creation of excessively large NumPy arrays.
    for s in range(start, end, _MAX_REGION_LEN):
        e = s + _MAX_REGION_LEN if (end - start) >= _MAX_REGION_LEN else end

        v_print_timestamp(args.verbose, "Now outputting " + chrm +
                          " for region: (" + str(s) + ", " + str(e) + ")", 2)

        modBaseScores = chromosome[s:e]
        if _MASK_TNAME in tnames:
            maskTrack = chromosome[s:e, maskRegionTName]
            maskIndex = genome.tracknames_continuous.index(maskRegionTName)

            modBaseScores[:, maskIndex] = \
                np.where(np.logical_or(np.isnan(maskTrack),
                         maskTrack > maskRegionsFileVal),
                         np.zeros(maskTrack.shape[0]),
                         np.ones(maskTrack.shape[0]))
        else:  # may still need to remove mask track
            # TODO this could be made more efficient
            for item in genome.tracknames_continuous:
                if item.find(_MASK_TNAME) >= 0:
                    maskIndex = genome.tracknames_continuous.index(item)
                    modBaseScores = np.delete(modBaseScores, maskIndex, axis=1)

        if len(modOrder) != len(set(modOrder)):  # intersect, if duplicates
            v_print_timestamp(args.verbose, "Intersection is enabled.")

            # Convert 0s to nan s.t. below mean works correctly, as an
            # intersection. Has no other impact, so fine to leave this way.
            modBaseScores[modBaseScores == 0] = np.nan

            for base_order, idxs in indices(modOrder,
                                            duplicates(modOrder)).iteritems():
                v_print_timestamp(args.verbose, "Intersecting for base '{}'."
                                  .format([base for base in args.priority
                                           if base in modBases][base_order]),
                                  3)
                #
                # NB: mean is not really needed, could be more efficient w/0|1;
                #     use of mean is, however, more general, permitting the
                #     use of other non-binary values, although we do not
                #     currently make particular use of this.
                #
                # key is to make any NaNs remain NaN, and leave others finite
                modBaseScores[:, idxs] = np.mean(modBaseScores[:, idxs],
                                                 axis=1, keepdims=True)

        # use '0' for unmodified bases (results in the unmod. base after)
        defaultBases = '0'

        # if masking all unset regions, use mask base for those, and '0' o/w
        if maskAllUnsetRegions:
            defaultBases = np.where(np.all(np.isnan(modBaseScores), axis=1),
                                    cUtils.MASK_BASE, '0')
            defaultBases = np.tile(defaultBases, (len(modBases), 1)).T

        modBasesA = np.where(np.logical_and(np.isfinite(modBaseScores),
                             modBaseScores != 0), modBases,
                             # either the modified base or the mask base
                             # any masking applied here is only for masking
                             # bases without any data
                             defaultBases)

        modOrderBasedPermutation = np.array(modOrder).argsort()
        orderedmodBasesA = modBasesA[:, modOrderBasedPermutation]

        referenceSeq = np.array(list(chromosome.seq[s:e].
                                tostring().upper()), dtype=np.str)

        # Filter the bases to take the modified bases in priority order.
        x = np.transpose(np.nonzero(orderedmodBasesA != '0'))
        u, idx = np.unique(x[:, 0], return_index=True)

        def maybeGetModBase(m, r, ambigMap):
            """Returns the modified base corresponding to
            the given putatively modified base (m). The base returned
            is the input putatively modified base if the corresponding
            reference base (r) is modifiable to the input base, or the
            complement of that base, if the complemented reference is
            modifiable to it, otherwise the reference base (r) itself
            is returned. This function also maps modified bases to any
            applicable ambiguity codes that are provided in ambigMap.
            """
            m = ambigMap.get(m) or m  # maybe map to an ambiguous base
            if m not in cUtils.MOD_MAP:
                return ambigMap.get(r) or r
            else:
                if cUtils.MOD_MAP[m] == r:
                    return m
                elif cUtils.MOD_MAP[m] == cUtils.complement(r)[0]:
                    return cUtils.complement(m)[0]
                else:
                    return ambigMap.get(r) or r

        # Initially the sequence is unmodified and we successively modify it.
        allModBases = np.copy(referenceSeq)
        # We vectorize the function for convenience.
        # NumPy vectorized functions still execute the Python code at
        maybeGetModBase = np.vectorize(maybeGetModBase)
        # Mask the sequence, allowing only base modifications
        # that modify their 'target' base (i.e. '5fC' = 'f' only modifies 'C').
        # Return the reference base for all non-modifiable bases
        # and for unmodified bases.

        if x.size > 0:
            hasModifiedBases = True
            # Modify bases
            np.put(allModBases, x[idx][:, 0],
                   maybeGetModBase(orderedmodBasesA[x[idx][:, 0],
                                   x[idx][:, 1]],
                                   allModBases[x[idx][:, 0]], ambigMap))
            if ambigMap:  # Replace with ambiguous bases in unmodified sequence
                unmodBaseIndices = np.delete(np.indices(np.shape(allModBases)),
                                             x[idx][:, 0])
                np.put(allModBases, unmodBaseIndices, np.vectorize(lambda base:
                       _maybeGetAmbigMapping(base, ambigMap))
                       (np.delete(allModBases, x[idx][:, 0])))

            if not suppressBED:
                # Create a BED track for each modified base with track data
                # Use unique modified bases only, since we may otherwise obtain
                # redundant track lines.
                for base in set(modBases + cUtils.complement(modBases)):
                    # NB: This could be done in a more efficient manner.
                    baseModIdxs = np.flatnonzero(allModBases[x[idx][:, 0]]
                                                 == base)
                    if baseModIdxs.size > 0:
                        # Get the position of the modified bases in the
                        # sequence, adding the genome start coordinate of
                        # the sequence to operate in actual genome coordinates.
                        modBaseCoords = x[idx][:, 0][baseModIdxs] + s

                        modBaseStartEnd = np.column_stack((modBaseCoords,
                                                          modBaseCoords+1))
                        # Save the track, appending to a Gzipped BED file.
                        # TODO save as string buffer (using list joins)
                        # and gzip after (with cStringIO).
                        # That will allow actual compression to occur.
                        with gzip.open(tnames[base],
                                       'ab') as BEDTrack:
                            np.savetxt(BEDTrack, modBaseStartEnd,
                                       str(chrm) + "\t%d\t%d\t" + base)

        if not suppressFASTA:
            # Output the unmodified sequence at a verbosity level
            # of at least 2, if not too long, otherwise only output
            # for a high verbosity level.
            v_print_timestamp(args.verbose, """Corresponding unmodified
                              reference sequence: \n""" +
                              ''.join(referenceSeq), 2
                              if len(referenceSeq) < 10000 else 6)
            # Concatenate the vector together to form the (string) sequence
            allbasesResult += ''.join(allModBases)

    if (not hasModifiedBases and not suppressBED):
        warn(""""There are no modified bases within the requested
             region. Accordingly, no BED files have been output
             for this region.""")
    return allbasesResult


def generateFASTAFile(file, id, genome, modOrder, chrm, start,
                      end, suppressBED, tnames, ambigMap, maskRegionsFileVal,
                      maskRegionTName, maskAllUnsetRegions):
    """Writes an optionally Gzipped FASTA file of the modified genome
    appending to the given file, using the given ID.
    No FASTA ID (i.e. '> ...') is written if no ID is given.
    """
    # Write either a Gzipped file or not, by using the appropriate function
    with cUtils.maybe_gzip_open(file, 'ab') as modGenomeFile:
        if id:
            modGenomeFile.write(">" + id + "\n")
        modGenomeFile.write(getModifiedGenome(genome, modOrder, chrm,
                            start, end, False, suppressBED, tnames, ambigMap,
                            maskRegionsFileVal, maskRegionTName,
                            maskAllUnsetRegions) + "\n")


def selectRandomRegion(genome, length):
    """Selects a random, non-exluded, chromosome of sufficient length.
    This method attempts to ensure that the region selected is
    wholly within a supercontig.
    """
    selectablechromosomes = {chromosome.name: chromosome.end for
                             chromosome in genome if
                             (chromosome.end >= length and
                              not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                            chromosome.name))}
    if not selectablechromosomes:
        die("""The region length provided is too long or all
               chromosomes have been excluded.""")
    chrm = random.choice(selectablechromosomes.keys())
    contigAttempts = 0
    while True:
        start = random.randint(genome[chrm].start, genome[chrm].end)
        end = start + length
        if genome[chrm].supercontigs[start:end]:
            break
        elif contigAttempts >= _MAX_CONTIG_ATTEMPTS:
            warn("""Attempts to procure sequence from a supercontig
                 were exhausted. Returning sequence that is not
                 wholly contained within a supercontig.""")
            break
        contigAttempts += 1
    return chrm, start, end


def parseRegion(genome, region):
    """Parses the provided region, ensuring its validity."""
    region = re.sub('[, ]', '', region)  # remove unwanted characters
    regionMatch = re.search(REGION_REGEX, region)
    if not regionMatch:
        die("Invalid region: invalid format.")
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
    NB: This method is not yet fully implemented.
    """
    # TODO cannot complete this without being able to
    # access the intervals over which the track data is defined.
    # Genomedata does not appear to support this.
    # TODO Ameliorate this.
    # print(genome.num_datapoints)
    # randomly sample 1000 bases of first chromosome to determine resolution
    for chromosome in genome:
        testRegion = random.randint(chromosome.start, chromosome.end)
        print(chromosome[testRegion:testRegion + 1000])
        break


# TODO add a custom action to parse all directories, which ensures
#      that each directory contains a trailing slash.

parser = argparse.ArgumentParser()

genomeArchive = parser.add_mutually_exclusive_group(required=True)
genomeArchive.add_argument('-G', '--genomeDataArchiveFullname',
                           help="The genome data archive. \
                           It must contain all needed \
                           sequence and track files. \
                           If one is not yet created, \
                           use '-d' instead to create it.")
genomeArchive.add_argument("-d", "--archiveCompDirs", nargs=2,
                           metavar=('GENOME_DIR', 'TRACKS_DIR'),
                           help="Two arguments first specifying the directory \
                           containing the genome and then the directory \
                           containing all modified base tracks. The genome \
                           directory must contain (optionally Gzipped) \
                           FASTA files of chromosomes and/or scaffolds. \
                           The sequence files must end in either a \
                           \".fa\" or \".mask\" extension \
                           (with an optional \".gz\" suffix). \
                           The track directory must contain \
                           (optionally Gzipped) genome tracks. \
                           They must have an extension describing \
                           their format. We currently support: \
                           \".wig\", \".bed\", and \".bedGraph\". \
                           Provided BED files must have exactly \
                           four columns. The fourth column must be numeric. \
                           Any rows with a data value of zero will be \
                           ignored (this does not apply to masking; \
                           see '-M'), except that such positions will be \
                           considered as having evidence against that \
                           particular modification (e.g. will be considered \
                           as 'set' and will therefore not be masked by \
                           '--maskAllUnsetRegions'). \
                           The filename of each track must \
                           specify what modified nucleobase it pertains to; \
                           one of: {}. \
                           Track names can also be of ambiguity codes \
                           (e.g. \"5xC\") on the positive strand only. \
                           Such tracks directly specify ambiguous loci. \
                           If multiple tracks of the same type are \
                           provided, all such tracks will be added to the \
                           archive. The output sequence will default to the \
                           union of all of the same modification type \
                           (but see '-I'). Alternatively, the track name \
                           can contain {}, in which case masking can be used \
                           via '-M' (refer to that option for details). \
                           Instead of a track directory, a single filename \
                           that meets the aforementioned requirements may \
                           be provided if the archive is to contain only \
                           one track. \
                           Ensure that all tracks are mapped to the same \
                           assembly and that this assembly matches the \
                           genome provided. This will \
                           create a genome data archive in an \"archive\" \
                           sub-directory of the provided track directory. \
                           Use '-G' instead to use an existing \
                           archive.".format(cUtils.MOD_BASES.keys(),
                                            _MASK_TNAME))
parser.add_argument("--archiveOutDir",
                    help="Only applicable if '-d' is used. \
                    The directory in which to save the created \
                    genome data archive. If not specified, this \
                    defaults to the directory containing the tracks \
                    (i.e. the second argument provided to '-d'). \
                    This defaults to {}.".format(_DEFAULT_ARCHIVE_NAME))
parser.add_argument("--archiveOutName", default=_DEFAULT_ARCHIVE_NAME,
                    help="Only applicable if '-d' is used. \
                    The name of the archive (i.e. the name \
                    of the directory which comprises the genome data \
                    archive).")
region = parser.add_mutually_exclusive_group()
region.add_argument('-r', '--region', help="Only output the modified genome \
                    for the given region. This can either be via a file \
                    or a region specification string. In the latter case, \
                    the region must be specified in the format: \
                    chrm<ID>:<start>-<end> (ex. chr1:500-510). \
                    If a file is being provided, it can be in \
                    any BEDTools-supported file format \
                    (BED, VCF, GFF, and Gzipped versions thereof). \
                    The full path to the file should be provided \
                    (or just the file name for the current directory).")
parser.add_argument('-c', '--centeredRegion', nargs='?', type=int,
                    const=_DEFAULT_CENTRED_REGION_LENGTH,
                    help="If used in conjunction with '-r', \
                    only output the modified genome \
                    for the given base pair interval (defaults to 500 bp), \
                    centered around the given region. \
                    NB: This region does not necessarily correspond to the \
                    centre of the peak, since the region's start and end \
                    coordinates alone are used to find the centre, \
                    as opposed to any peak information \
                    (from a narrowPeaks file, for example).")
region.add_argument('-R', '--randomRegion', nargs='?',
                    const=_DEFAULT_RAN_LENGTH, type=int,
                    help="Output the modified genome for a random region. \
                    The chrmomsome will be randomly selected and its \
                    coordinate space will be randomly and uniformly sampled. \
                    A length for the random region can either be specified \
                    or it will otherwise be set to a reasonably \
                    small default. The length chosen may constrain the \
                    selection of a chromosome.")
parser.add_argument('-A', '--alterIncludedChromosomes',
                    choices=CHROMOSOME_TYPE_REGEXES.keys(),
                    help="Include or exclude chromosome \
                    types. '{}': Use only autosomal chromosomes \
                    (excludes chrM). \
                    '{}': Use only allosomal chromosomes (excludes chrM). \
                    '{}': Use only the mitochondrial chromosome. \
                    '{}': Include the mitochondrial chromosome. \
                    NB: default chromosomal exclusions include: \
                    unmapped data, haplotypes, and chrM. \
                    This parameter will be ignored if a \
                    specific genomic region is queried via '-r', \
                    but will be considered if a file of genomic \
                    regions is provided (also via '-r'). \
                    ".format(_AUTOSOME_ONLY_FLAG,
                             _ALLOSOME_ONLY_FLAG,
                             _MITOCHONDRIAL_ONLY_FLAG,
                             _MITOCHONDRIAL_INCLUSION_FLAG))
parser.add_argument('-p', '--priority', default=_DEFAULT_BASE_PRIORITY,
                    help="Specify the priority \
                    of modified bases. The default is:"
                    + _DEFAULT_BASE_PRIORITY + ", which is based upon "
                    + _DEFAULT_BASE_PRIORITY_COMMENT + ".")
BEDGeneration = parser.add_mutually_exclusive_group()
BEDGeneration.add_argument('-b', '--suppressBED', action='store_true',
                           help="Do not generate any BED tracks.")
BEDGeneration.add_argument('-B', '--onlyBED', action='store_true',
                           help="Only generate any BED tracks \
                           (i.e. do not output any sequence information). \
                           Note that generated BED files are always \
                           appended to and created in the CWD \
                           irrespective of the use of this option. \
                           This parameter is ignored if '-f' is used.")
parser.add_argument("--BEDOutDir", default='.',
                    help="Only applicable if '-b' is not used. \
                    The directory in which to save the created \
                    BED tracks. If not specified, this \
                    defaults to the current working directory.")
parser.add_argument('-f', '--fastaFile', nargs='?', type=str,
                    const=_DEFAULT_FASTA_FILENAME, help="Output to \
                    a file instead of STDOUT. Provide a full path \
                    to a file to append the modified genome in \
                    FASTA format. If this parameter is invoked \
                    without any arguments, a default filename \
                    will be used within the current directory. \
                    This will override the '-B' parameter \
                    (i.e. a FASTA file with always be produced). \
                    The output file will be Gzipped iff the \
                    path provided ends in \".gz\".")
parser.add_argument('-I', '--intersection', action='store_true',
                    help="If multiple files of the same modification \
                    type are given, take their intersection. \
                    This option is used to override the default, \
                    which is to take their union.")
ambigModUsage = \
    parser.add_argument_group(title="Ambiguous Modification",
                              description="Specify that some of the data \
                              provided for a given modified base is unable \
                              to differentiate between some number of \
                              modifications. This ensures that Cytomod \
                              outputs the correct ambiguity code such that \
                              modified genomes do not purport to convey \
                              greater information than they truly contain. \
                              The most general applicable ambiguities should \
                              be specified. Therefore, each modified \
                              nucleobase may reside in at most one specified \
                              set of ambiguities.")
ambigModUsage.add_argument('--mh', action='store_const', const='mh',
                           help="Specify that input data is not able to \
                           differentiate between 5mC and 5hmC. This would \
                           be the case if the data originated from a protocol \
                           which only included conventional bisulfite \
                           sequencing.", default='')
ambigModUsage.add_argument('--fC', action='store_const', const='fC',
                           help="Specify that input data is not able to \
                           differentiate between 5fC and C. This would be the \
                           case if the data originated from a protocol which \
                           only included oxidative bisulfite sequencing.",
                           default='')
ambigModUsage.add_argument('--fc', action='store_const', const='fc',
                           help="Specify that input data is not able to \
                           differentiate between 5fC and 5caC. This would be \
                           the case if the data originated from a protocol \
                           which only included M.SssI methylase-assisted \
                           bisulfite sequencing.",
                           default='')
ambigModUsage.add_argument('-M', '--maskRegions', nargs='?', type=float,
                           const=_DEFAULT_MASK_VALUE,
                           help="Hard mask C/G nucleobases to unknown state. \
                           Assumes that the archive contains or is \
                           being built with a trackname containing \"" +
                           _MASK_TNAME + "\". \
                           The containing loci will be interpreted as \
                           nucleobases of unknown modification state. \
                           They will be accordingly set to the appropriate \
                           (maximally) ambiguous base. This will override any \
                           other modifications at those loci. \
                           This parameter can accept an optional argument, \
                           indicating a value at and below which the locus \
                           is considered ambiguous. If not provided, \
                           this defaults to " + str(_DEFAULT_MASK_VALUE) + ". \
                           An example use case for this option would be to \
                           use a mask file containing coverage information \
                           and to mask all bases of insufficient coverage.")
ambigModUsage.add_argument('--maskAllUnsetRegions', action='store_true',
                           help="Hard mask all C/G nucleobases without \
                           any modification information to unknown state. \
                           Masked nucleobases are those lacking data, \
                           that is, bases not present in the archive \
                           nor in any files used to generate an archive. \
                           Set, but unmodified, bases can be provided, \
                           within any modified genomic interval file with \
                           a value of 0. See '-d' for further details.\
                           Unset modifiable bases will be accordingly set to \
                           the appropriate (maximally) ambiguous base. \
                           This will override any other modifications \
                           at those loci. An example use case for this option \
                           would be when using array data, for which \
                           only a subset of bases are queried.")

# TODO Implement this?
# NB: neither TRF nor dustmasker work upon modified genomes
# parser.add_argument('-M', '--hardMaskRepetitiveRegions',
#                     help="Hard mask low complexity regions. The program \
#                     first will attempt to use the Tandem Repeat Finder \
#                     (TRF) program, by Gary Benson, using a system call \
#                     to the 'trf' executable. If this is not available, \
#                     the program will attempt to use dustmasker, of the NCBI \
#                     toolkit, by Morgulis et al (executed via 'dustmasker'). \
#                     By default, a threshold of 16 is used with dustmasker, \
#                     incrasing its masking, as done by Frith et al. 2010. \
#                     An argument containing the first letter of a tool \
#                     may be passed to override this behaviour. \
#                     Additionally, passing 'B' or 'F' will use TRF \
#                     with parameters recommended by the TRF author \
#                     or from Frith et al., respectively. \
#                     Default behaviour is equivalent to providing the \
#                     'F' sub-argument.", action='store_const', const='F')
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action="count")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

if (args.region and (not os.path.isfile(args.region))
        and args.alterIncludedChromosomes):
    # '-A' is ignored iff '-r' is provided and is not a path to a file
    warn("""Chromosome inclusion or exclusion setting was ignored,
         since a single specific region was specifed.""")

if args.centeredRegion and not args.region:
    warn("""Centered region argument ignored, since no specific regions
            were specifed. Specify '-r' with '-c' to use centred regions.""")

if args.onlyBED and args.fastaFile:
    warn("""Request to only generate BED files ignored,
            since an output FASTA path was provided.""")

if not args.archiveCompDirs and (args.archiveOutDir or args.archiveOutName !=
                                 _DEFAULT_ARCHIVE_NAME):
    warn("""Neither the output directory nor the output name of a
            genome archive are applicable, because an archive
            is not being created.""")

if args.suppressBED and args.BEDOutDir:
    warn("""The directory provided for BED output has been ignored, since
            BED output has been suppressed.""")

if args.alterIncludedChromosomes:
    _modifychrmExclusionRegex(args.alterIncludedChromosomes)
# NB: Ensure to update this to include all arguments in the ambigModUsage group
# TODO It would be nice if there was an automated means of accomplishing this
# Ensure provided ambiguities are unique
if (''.join(OrderedDict.fromkeys(args.mh + args.fC + args.fc).keys()) !=
        args.mh + args.fC + args.fc):
    die("Provided ambiguity codes must be unique.")
# Create a map from ambiguity codes to the specified ambiguities
ambigMap = {}
for k in [args.fC, args.mh, args.fc]:
    ambigMap.update({b: cUtils.INVERTED_AMBIG_MOD_BASES[k]
                     for b in ''.join(k)})
v_print_timestamp(args.verbose, "Using the following ambiguity map: " +
                  str(ambigMap) + ".", 2)

from genomedata import Genome, load_genomedata

genomeDataArchiveFullname = ""
if args.archiveCompDirs:
    # Create the archive in the directory of the track files
    # unless another directory is specified
    defaultArchiveDir = os.path.dirname(args.archiveCompDirs[1]) if \
        os.path.isfile(args.archiveCompDirs[1]) else args.archiveCompDirs[1]
    genomeDataArchiveFullname = (args.archiveOutDir or defaultArchiveDir) + \
        '/' + args.archiveOutName + '/'
    v_print_timestamp(args.verbose, "Creating genomedata archive in {}.".
                      format(genomeDataArchiveFullname))
    # Create the genome data archive
    # Load all supported track files in the tracks directory
    # Load all FASTA files in the sequences directory
    # TODO this could be generalized and could be made less redundant
    FASTA_file_list = glob.glob(args.archiveCompDirs[0] + "/*.fa") or \
        glob.glob(args.archiveCompDirs[0] + "/*.fa.gz")
    if not FASTA_file_list:
        FASTA_file_list = glob.glob(args.archiveCompDirs[0] + "/*.mask") or \
            glob.glob(args.archiveCompDirs[0] + "/*.mask.gz")
        v_print_timestamp(args.verbose, """Detected that a repeat masked
                          genome is being used for archive creation.""")
    if not FASTA_file_list:
        die("Unable to locate FASTA input files for archive creation.")
    load_genomedata.load_genomedata(
        genomeDataArchiveFullname,
        # args.archiveCompDirs[1] is either a directory of tracks
        # or is the full path to a single track.
        # TODO check for a valid track path, that actually
        #      contains tracks, otherwise throw an error.
        # NB: the tracks parameter must be a *list of pair(s)*,
        #     specifying: (track name, track/directory path).
        tracks=([(os.path.basename(args.archiveCompDirs[1]),
                  args.archiveCompDirs[1])]
                if os.path.isfile(args.archiveCompDirs[1]) and
                re.search(SUPPORTED_FORMATS_REGEX, args.archiveCompDirs[1])
                else [(track, args.archiveCompDirs[1] + track)
                for track in os.listdir(args.archiveCompDirs[1])
                if re.search(SUPPORTED_FORMATS_REGEX, track)]),
        # args.archiveCompDirs[0] contains the directory with all FASTAs
        seqfilenames=FASTA_file_list, verbose=args.verbose)

else:
    v_print_timestamp(args.verbose, "Using existing genomedata archive.")
    genomeDataArchiveFullname = args.genomeDataArchiveFullname

with Genome(genomeDataArchiveFullname) as genome:
    warnings.simplefilter("ignore")  # Ignore supercontig warnings
    v_print_timestamp(args.verbose, "Genomedata archive successfully loaded.")

    maskRegionTName = ''
    modBases = []
    modOrder = []
    tnames = {}

    for track in genome.tracknames_continuous:
        if _MASK_TNAME in str(track):
            if args.maskRegions is not None:
                modBases.append(cUtils.MASK_BASE)

                tnames[_MASK_TNAME] = track
                maskRegionTName = track
            else:
                warn("""Genomedata archive contains a mask track, but
                        '-M' was not given.
                        Masking will not be performed.""")
        else:
            trackToBase = {covalent_mod_base for
                           covalent_mod, covalent_mod_base in
                           cUtils.COVALENT_MOD_BASES_TO_BASE_MAP.iteritems()
                           if covalent_mod in track}.pop()
            if trackToBase:  # add a regular or ambigous modified base track
                modBases.append(trackToBase)
            else:
                warn("Unrecognized track " + track + " has been ignored.")
    if args.maskRegions is not None and _MASK_TNAME not in tnames:
        die("""Masking of genome regions requires the generation of a
               Genomedata archive containing a mask track.""")

    # For modOrder, lowest numbers have higher priority (i.e. 0 is highest).
    for i, base in enumerate(modBases):  # get rel. ordering
        modOrder += [(args.priority.index(base) if args.intersection
                      # in else, no int., so ensure unique orders
                      #
                      # take fourth power to ensure priority always
                      # outweighs index and overcomes similar priorities
                      # to ensure unique orders
                      else i + (args.priority.index(base) + 1)**4)]

    # Get absolute (index-based) ordering from the relative ordering.
    modOrder = [sorted(modOrder).index(x) for x in modOrder]

    addtl_mask_msg = ""
    if args.maskRegions is not None:
        addtl_mask_msg = """All loci implicated by the mask will be masked
                            irrespective of any mods at those loci."""

        # Masked bases are assigned the highest priority (mask all others).
        modOrder = [order + 1 for order in modOrder]
        modOrder[modBases.index(cUtils.MASK_BASE)] = 0
    elif args.maskAllUnsetRegions:
        addtl_mask_msg = """All loci with missing data will be masked."""
        # No priority assigned, since this masking is performed separately.

    v_print_timestamp(args.verbose,
                      "Masking is enabled. {}".format(addtl_mask_msg), 2)

    v_print_timestamp(args.verbose, """The order of preference for base
                      modifications (from highest to lowest) is: """ +
                      ','.join(list(args.priority)) + ".")

    # Before computing modified bases in blocks, remove any existing BED files
    # and write the tracks' headers.
    # Also, store the tracks' names, keyed by modified base, for future use.
    if not args.suppressBED:
        trackID = os.path.splitext(os.path.basename(args.fastaFile
                                   or _DEFAULT_FASTA_FILENAME))[0]
        trackID += '-' if trackID else ''
        procModBases = []
        for base in modBases + cUtils.complement(modBases):
            trackFileName = "{}/track-{}{}.bed.gz".\
                format(args.BEDOutDir, trackID, base)
            tnames[base] = trackFileName
            # 'EAFP' way of removing any existing old tracks
            try:
                os.remove(trackFileName)
            except OSError:
                pass
            with gzip.open(trackFileName, 'ab') as BEDTrack:
                BEDTrack.write(getTrackHeader(base))

    if args.region or args.randomRegion:
        if args.randomRegion:  # Random region
            chrms, starts, ends = selectRandomRegion(genome, args.randomRegion)
            regions = np.matrix([chrms, starts, ends])
        elif os.path.isfile(args.region):  # 'BED-like' set of regions
            import pybedtools
            # NB: In pybedtools, all regions behave as if they are 0-based,
            # despite non-BED files being 1-based. Thus, we do not need to
            # alter any code in our program.
            BEDTool = pybedtools.BedTool(args.region)
            # generate regions, only including the appropriate chromosomes
            regions = np.matrix([[interval['chrom'], interval['start'],
                                interval['end']] for interval in BEDTool
                                if not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                                 interval['chrom'])])
        else:  # A single, specific, region ('genome browser-like')
            chrms, starts, ends = parseRegion(genome, args.region)
            regions = np.matrix([chrms, starts, ends])
        for i in xrange(0, len(regions.flat), 3):
            chrm, start, end = regions.flat[i], int(regions.flat[i + 1]), \
                int(regions.flat[i + 2])
            if args.centeredRegion:
                centre = int(math.floor((start + end) / 2))
                start = int(centre - args.centeredRegion / 2)
                if start < 0:
                    start = 0
                end = int(centre + args.centeredRegion / 2)
                if end > int(genome[chrm].end):
                    end = int(genome[chrm].end)
            regionStr = chrm + ":" + str(start) + "-" + str(end)
            v_print_timestamp(args.verbose, """Outputting the modified
                              genome for: """ + regionStr + ".")
            if args.fastaFile:
                generateFASTAFile(args.fastaFile, regionStr, genome, modOrder,
                                  chrm, start, end, args.suppressBED, tnames,
                                  ambigMap, args.maskRegions,
                                  maskRegionTName, args.maskAllUnsetRegions)
            else:
                print(getModifiedGenome(genome, modOrder, chrm, start, end,
                                        args.onlyBED, args.suppressBED,
                                        tnames, ambigMap,
                                        args.maskRegions, maskRegionTName,
                                        args.maskAllUnsetRegions))
    else:
        for chromosome in [chromosome for chromosome in genome
                           if not re.search(CHROMOSOME_EXCLUSION_REGEX,
                                            chromosome.name)]:
            v_print_timestamp(args.verbose, """Outputting the modified
                              genome for: """ + chromosome.name)
            generateFASTAFile(args.fastaFile or _DEFAULT_FASTA_FILENAME,
                              chromosome.name, genome, modOrder,
                              chromosome.name, int(chromosome.start),
                              int(chromosome.end), args.suppressBED, tnames,
                              ambigMap, args.maskRegions, maskRegionTName,
                              args.maskAllUnsetRegions)

v_print_timestamp(args.verbose, "Program complete.")
