#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement, division, print_function

"""Prints (to STDOUT) a minimal MEME text file, suitable for FIMO,
from an input set of sequences (one raw sequence per line)
or an input PWM, PFM, TRANSFAC matrix, or JASPAR matrix.
The background frequencies should be manually adjusted
to use case. Modified output (if applicable) is always
written to a file, while the MEME output written to STDOUT
always correpsonds to the unmodified motif."""

import os
import textwrap
import re

try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

import numpy as np
import numpy.testing as npt
from scipy.stats import itemfreq

import cUtils
__version__ = cUtils.VERSION

_DELIM = "\t"
BOTH_STRANDS = 0
STRANDS = '+ -' if BOTH_STRANDS else '+'
FROM_STR = 'custom'
# The alphabet frequencies used were (WRT Cs only):
# 2.95% 5mC, 0.055 ± 0.008% 5hmC, 0.0014 ± 0.0003% 5fC, and 0.000335% 5caC
# (Respectively: Ito et al. 2011, Booth et al. 2014, Booth et al. 2014,
# and Ito et al. 2011)
MOTIF_ALPHABET_BG_FREQUENCIES = \
    {'T': 0.292, 'A': 0.292, 'C': 0.201745991, 'G': 0.201745991,
     # (using mouse GC content of 41.6%)
     # From: NCBI Eukaryotic Genome Report File
     'm': 0.006136, '1': 0.006136, 'h': 0.0001144, '2': 0.0001144,
     'f': 0.000002912, '3': 0.000002912, 'c': 0.000000697, '4': 0.000000697}
npt.assert_allclose([sum(MOTIF_ALPHABET_BG_FREQUENCIES.itervalues())], [1])
# The MEME Suite uses ASCII ordering for custom alphabets
# This is the natural lexicographic sorting order, so no "key" is needed
MOTIF_ALPHABET = sorted(list(MOTIF_ALPHABET_BG_FREQUENCIES.keys()))
MOTIF_ALPHABET_BG_FREQUENCIES_OUTPUT = \
    ' '.join([str(k) + ' ' + str(v) for k, v
              in iter(sorted(MOTIF_ALPHABET_BG_FREQUENCIES.iteritems()))])

MEME_HEADER = """MEME version 4

ALPHABET "DNA with covalent modifications"
A "Adenine" 8510A8 ~ T "Thymine" A89610
C "Cytosine" A50026 ~ G "Guanine" 313695
m "5-Methylcytosine" D73027 ~ 1 "Guanine:5-Methylcytosine" 4575B4
h "5-Hydroxymethylcytosine" F46D43 ~ 2 "Guanine:5-Hydroxymethylcytosine" 74ADD1
f "5-Formylcytosine" FDAE61 ~ 3 "Guanine:5-Formylcytosine" ABD9E9
c "5-Carboxylcytosine" FEE090 ~ 4 "Guanine:5-Carboxylcytosine" E0F3F8
R = AG
Y = CT
K = GT
M = AC
S = CG
W = AT
B = CGT
D = GAT
H = ACT
V = ACG
N = ACGT
X = ACGT
END ALPHABET

strands: %s

Background letter frequencies (from %s):
%s

""" % (STRANDS, FROM_STR, MOTIF_ALPHABET_BG_FREQUENCIES_OUTPUT)


def die(msg):
    cUtils.die(msg, os.path.basename(__file__))


def warn(msg):
    cUtils.warn(msg, os.path.basename(__file__))

import argparse
parser = argparse.ArgumentParser()
inputFileGroupTitle = \
    parser.add_argument_group(title="Input File", description="The input set \
                              of sequences or matrix from which to create a \
                              MEME minimal text output file.")
inputFile = inputFileGroupTitle.add_mutually_exclusive_group(required=True)
inputFile.add_argument('-s', '--inSeqFile', type=str, help="File containing \
                       an input set of raw sequences.")
inputFile.add_argument('-p', '--inPWMFile', type=str, help="File containing \
                       an input whitespace-delimited PWM (i.e. frequency \
                       matrix). The file must only contain floats \
                       (but see '-A'), lexicographically ordered \
                       (i.e. A, C, G, T).")
inputFile.add_argument('-c', '--inPFMFile', type=str, help="File containing \
                       an input whitespace-delimited PFM (i.e. a count \
                       matrix; e.g. a hPDI matrix). The file must only \
                       contain floats (but see '-A'), lexicographically \
                       ordered (i.e. A, C, G, T).")
inputFile.add_argument('-t', '--inTRANSFACFile', type=str, help="File containing \
                       an input TRANSFAC count matrix. The file must exactly \
                       conform to the standard TRANSFAC matrix output format. \
                       It must neither begin nor be terminated with \"XX\" \
                       (i.e. those delimiters should be omitted).")
inputFile.add_argument('-f', '--inTRANSFACFreqFile', type=str, help="File containing \
                       an input TRANSFAC frequency matrix. The file must \
                       exactly conform to the STAMP TRANSFAC frequency matrix \
                       output format. The usual \"XX\" footer that this file \
                       is usually terminated by must be removed manually. \
                       The file should not contain any \"XX\" delimiters.")
inputFile.add_argument('-j', '--inJASPARFile', type=str, help="File containing \
                       an input JASPAR matrix. The file must exactly \
                       conform to the JASPAR matrix format.")
modBaseSpecifiers = parser.add_mutually_exclusive_group()
modBasePositions = parser.add_mutually_exclusive_group()
modBaseSpecifiers.add_argument('-M', '--baseModification', type=str,
                               help="Modify the motif to use the modified base \
                               provided in this argument at the position \
                               specified by '-P'. The resultant motif will \
                               use the given modified base with 100%% \
                               frequency at the specified positions. \
                               This will cause the program to write a file \
                               (as opposed to the usual output to STDOUT) \
                               for the given modification.")
modBasePositions.add_argument('-P', '--baseModPosition', type=int, help="The position \
                              at which to modify the motif (using the base \
                              specified by '-M'), * indexed from 1 *.")
modBaseSpecifiers.add_argument('-C', '--tryAllCModsAtPos', type=int,
                               nargs='?', const=cUtils._PARAM_A_CONST_VAL,
                               help="Modify the motif at the given position, \
                               whose consensus sequence should correspond to \
                               a cytosine at the given position. \
                               No position needs to be given if \
                               '-A' is also used. This will \
                               cause the program to write a file \
                               (as opposed to the usual output to STDOUT) \
                               for each cytosine modification and the \
                               unmodified motif will be output to STDOUT \
                               as well. The resultant motif will use the \
                               given modified base with 100%% frequency at \
                               the specified positions. Note that this is \
                               * indexed from 1 *.")
modBasePositions.add_argument('-A',
                              '--baseModificationAtAllModifiablePosFractions',
                              type=str, nargs='?',
                              const=str(cUtils._PARAM_A_CONST_VAL),
                              help="Modify the motif to use the modified base \
                              provided for all modifiable motif positions. \
                              A position is considered a modifiable one iff \
                              the nucleobase at said position contains some \
                              cytosine or guanine fraction (i.e. \
                              the resultant PWM has non-zero value for its \
                              'C' or 'G' entries at a given position). \
                              This option results in all modifiable positions \
                              being modified proportionally. That is, \
                              only the cytosine or guanine fraction \
                              will be modified to the provided base. \
                              The provided base can either be WRT \
                              cytosines on the positive strand (i.e. 'm') \
                              or a complementary modification (i.e. '1') \
                              In either case, the correct corresponding \
                              modification will be selected, depending upon \
                              whether the fraction being modified is a 'C' \
                              or a 'G' (e.g. input: '-A m'. result: \
                              'C' <-> 'm'; 'G' <-> '1'). \
                              This argument can be used concomittantly with \
                              '-C', in which case no modified nucleobase \
                              need be provided to this argument, since all \
                              possibilities will be attempted (for all \
                              modifiable positions). \
                              This assumes that the input motif did not \
                              already contain any modified nucleobases. \
                              This will cause the program to write a file \
                              (as opposed to the usual output to STDOUT) \
                              for the given modification.")
parser.add_argument('-a', '--annotated', action='store_true',
                    help="Assume that the provided matrix file contains \
                    identifiers in the first column. This option allows \
                    for them to be removed and prevents them from \
                    interfering with the processing of the input file.")
parser.add_argument('-S', '--skipMotifPortion', nargs='+', action='append',
                    help="Skip a portion of a motif. This might be used \
                    if one wishes to limit a motif to its highly \
                    informative portion. Generally, one would only \
                    remove bases from the start or end of a motif, \
                    but this option also permits removal of the middle \
                    of a motif. Specify the portion to remove via: \
                    space-delimited row indicies (from 0) \
                    or via \"start:end\". In the latter case, \
                    the given interval is interpreted as the \
                    half-open interval: [start, end).")
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action="count")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

if bool(args.baseModification) ^ bool(args.baseModPosition):
    warn("""Any base modification specification must specify both the
         particular base to be modified (via '-M') and the position
         for the modification to occur (via '-P'). Motif modification
         will not be performed. The '-P' parameter is ignored if
         '-C' is provided.""")

if (not args.tryAllCModsAtPos and
        (args.baseModificationAtAllModifiablePosFractions
         == cUtils._PARAM_A_CONST_VAL)):
    die("""You must either provide the modified base to '-A' or use
        '-C' to use all possible nucleobases.""")

if (not args.baseModificationAtAllModifiablePosFractions and
        args.tryAllCModsAtPos == cUtils._PARAM_A_CONST_VAL):
    die("""You must either provide the position to modify to '-C' or use
        '-A' to use all possible positions.""")

filename = (args.inSeqFile or args.inPWMFile or
            args.inPFMFile or args.inTRANSFACFile or
            args.inTRANSFACFreqFile or args.inJASPARFile)

if args.inSeqFile:
    # NB: min dimensionality of 1 is needed for the character view
    motifs = np.loadtxt(filename, dtype=str, ndmin=1)
    motifChars = motifs.view('S1').reshape((motifs.size, -1))
    totalNumBases = len(motifChars)

    freqMatrix = np.zeros((motifChars.shape[1],
                          len(MOTIF_ALPHABET_BG_FREQUENCIES)))
    for i in range(0, motifChars.shape[1]):
        motifCharsInts = motifChars[:, i].view(np.uint8)
        # NB: itemfreq internally uses bincount; we must map to and from ints
        f = itemfreq(motifCharsInts)
        bases = f[:, 0].view('U1')
        baseFreqs = f[:, 1]
        # Append to the letter frequency matrix
        matIn = StringIO(_DELIM.join(str(x) for x in
                                     (baseFreqs[idx][0]/totalNumBases if
                                     len(idx[0]) else 0 for idx in
                                     (np.nonzero(bases == base) for
                                      base in MOTIF_ALPHABET))))
        freqMatrix = np.vstack([freqMatrix, np.loadtxt(matIn)])
else:  # PWM or PFM
    # process the matrix for MEME format compatibility
    # If needed, skip columns or rows accordingly and/or
    # use unpack to transpose the matrix.
    ncols = 0
    with open(filename, 'rb') as inFile:
        if args.annotated:
                ncols = len(inFile.readline().split())
        elif args.inJASPARFile:
            # count number of columns and remove unwanted characters
            temp = ''
            for i, line in enumerate(inFile):
                line = re.sub('[\[\]]', '', line)
                if i == 0:
                    ncols = len(line.split())
                temp += line
            inFile = StringIO(temp)
        csvData = 0
        if args.inTRANSFACFile or args.inTRANSFACFreqFile:
            csvData = np.loadtxt(inFile, dtype=np.float,
                                 usecols=range(1, 5), skiprows=1)
        else:
            csvData = np.loadtxt(inFile,
                                 unpack=True, dtype=np.float,
                                 usecols=range(1, ncols)
                                 if args.annotated or args.inJASPARFile
                                 else None)
        if args.inPFMFile or args.inTRANSFACFile or args.inJASPARFile:
            """This function transforms a given PFM column to a PWM column.
            That is, it computes the frequency of each element of the input
            1D array (which is itself a count), and replaces the count with
            its computed frequency. The function also returns the resultant
            array, despite its having been modified in place."""
            def _computeFresFromCountSlice(countsForBase):
                sum = np.sum(countsForBase)
                for count in np.nditer(countsForBase, op_flags=['readwrite']):
                    count[...] = count / sum
                return countsForBase  # needed to use numpy.apply_along_axis
            np.apply_along_axis(_computeFresFromCountSlice, 1, csvData)

        freqMatrix = np.hstack((np.zeros((csvData.shape[0],
                               MOTIF_ALPHABET.index('A'))), csvData,
                               np.zeros((csvData.shape[0],
                                        (len(MOTIF_ALPHABET) - 1) -
                                        MOTIF_ALPHABET.index('T')))))

    if args.skipMotifPortion:
        skipRows = args.skipMotifPortion
        match = re.search('(\d*):(\d*)', str(skipRows).strip('[]'))
        if match:
            start = int(match.group(1)) if match.group(1) else 0
            end = int(match.group(2)) if match.group(2) \
                else freqMatrix.shape[0]
            if (start > end):
                die("""The start index of the motif rows to skip cannot be greater \
                than the end.""")
            if (end > freqMatrix.shape[0]):
                die("""The provided end index of the motif rows to skip exceeds \
                the motif length.""")
            if ((end - start) <= 1):
                warn("""The provided skip interval will only cause one motif row \
                to be excluded. Recall that intervals provided in colon \
                notation are half-open and that a single row to skip can \
                be specified by only providing that row number alone.""")
            skipRows = range(start, end)
        freqMatrix = np.delete(freqMatrix, skipRows, 0)

    totalNumBases = freqMatrix.shape[0]

    MEMEBody = textwrap.dedent("""MOTIF %s\nletter-probability matrix: nsites= %d\n""") \
        % (filename, totalNumBases)

    if ((args.baseModification and args.baseModPosition)
            or args.tryAllCModsAtPos or
            args.baseModificationAtAllModifiablePosFractions):
        modFreqMatrix = np.copy(freqMatrix)
        baseModPos = args.tryAllCModsAtPos or args.baseModPosition
        for b in (cUtils.MOD_BASE_NAMES.keys() if args.tryAllCModsAtPos
                  else (args.baseModification or
                  args.baseModificationAtAllModifiablePosFractions)):
            if args.baseModificationAtAllModifiablePosFractions:
                # modify cytosine fractions
                modFreqMatrix[:, MOTIF_ALPHABET.index(cUtils.getMBMaybeFromComp(b))] = \
                    modFreqMatrix[:, MOTIF_ALPHABET.index('C')]
                modFreqMatrix[:, MOTIF_ALPHABET.index('C')] = \
                    np.zeros(modFreqMatrix.shape[0])
                # modify guanine fractions
                modFreqMatrix[:, MOTIF_ALPHABET.index(cUtils.getCompMaybeFromMB(b))] = \
                    modFreqMatrix[:, MOTIF_ALPHABET.index('G')]
                modFreqMatrix[:, MOTIF_ALPHABET.index('G')] = \
                    np.zeros(modFreqMatrix.shape[0])
            else:
                # zero all entries along the frequency matrix,
                # for the given (row) position, except that corresponding
                # to the (column) index of the modified base,
                # which is set to unity.
                modFreqMatrix[(baseModPos - 1), ] = \
                    np.zeros((1, freqMatrix.shape[1]))
                modFreqMatrix[(baseModPos - 1),
                              MOTIF_ALPHABET.index(b)] = 1
            with open((os.path.splitext(filename)[0] + '-' +
                       cUtils.MOD_BASE_NAMES[b] + '.meme'), "a") as outFile:
                outFile.write(MEME_HEADER)
                outFile.write(MEMEBody)
                np.savetxt(outFile, modFreqMatrix, '%f', _DELIM)

output = StringIO()
np.savetxt(output, freqMatrix, '%f', _DELIM)
MEMEBody += output.getvalue()
output.close()
print(MEME_HEADER + MEMEBody)
