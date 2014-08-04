#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Prints (to STDOUT) a minimal MEME text file, suitable for FIMO,
from an input set of sequences (one raw sequence per line)
or an input PWM. The background frequencies should be manually
adjusted to use case.
"""

from __future__ import with_statement, division, print_function

__version__ = "$Revision: 0.04$"

import sys
import os
import textwrap
try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

import numpy as np
import numpy.testing as npt
from scipy.stats import itemfreq

BOTH_STRANDS = 0
STRANDS = '+ -' if BOTH_STRANDS else '+'
FROM_STR = 'custom'
MOD_BASE_NAMES = {'m': '5mC', 'h': '5hmC', 'f': '5fC', 'c': '5caC'}
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

_MSG_PREFIX = '>> <genMinMEME> '


def warn(msg):
    """Emit a warning to STDERR, deindenting the provided message."""
    print(_MSG_PREFIX + "Warning: " + textwrap.dedent(msg), file=sys.stderr)


import argparse
parser = argparse.ArgumentParser()
inputFile = parser.add_mutually_exclusive_group(required=True)
inputFile.add_argument('-s', '--inSeqFile', type=str, help="File containing \
                       an input set of raw sequences.")
inputFile.add_argument('-p', '--inPWMFile', type=str, help="File containing \
                       an input tab-delimited PWM. The file must only contain \
                       floats, lexicographically ordered (i.e. A, C, G, T).")
modBaseSpecifiers = parser.add_mutually_exclusive_group()
modBaseSpecifiers.add_argument('-M', '--baseModification', type=str,
                               help="Modify the motif to use the modified base \
                               provided in this argument at the position \
                               specified by '-P'. The resultant motif will \
                               use the given modified base with 100% \
                               frequency at the specified positions. \
                               This will cause the program to write a file \
                               (as opposed to the usual output to STDOUT) \
                               for the given modification.")
parser.add_argument('-P', '--baseModPosition', type=int, help="The position \
                    at which to modify the motif (using the base specified \
                    by '-M'), * indexed from 1 *.")
modBaseSpecifiers.add_argument('-C', '--tryAllCModsAtPos', type=int,
                               help="Modify the motif at the given position, \
                               whose consensus sequence should correspond to \
                               a cytosine at the given position. This will \
                               cause the program to write a file \
                               (as opposed to the usual output to STDOUT) \
                               for each cytosine modification and the \
                               unmodified motif will be output to STDOUT \
                               as well. The resultant motif will use the \
                               given modified base with 100% frequency at \
                               the specified positions. Note that this is \
                               * indexed from 1 *.")
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

filename = args.inPWMFile or args.inSeqFile

if (args.inPWMFile):
    # transpose the PWM for MEME format compatibility via unpack
    csvData = np.loadtxt(open(filename, 'rb'), delimiter="\t",
                         unpack=True, dtype=np.float)
    freqMatrix = np.hstack((np.zeros((csvData.shape[0],
                           MOTIF_ALPHABET.index('A'))), csvData,
                           np.zeros((csvData.shape[0],
                                    (len(MOTIF_ALPHABET) - 1) -
                                    MOTIF_ALPHABET.index('T')))))
    totalNumBases = csvData.shape[0]
else:
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
        matIn = StringIO("\t".join(str(x) for x in
                                   (baseFreqs[idx][0]/totalNumBases if
                                   len(idx[0]) else 0 for idx in
                                   (np.nonzero(bases == base) for
                                    base in MOTIF_ALPHABET))))
        freqMatrix = np.vstack([freqMatrix, np.loadtxt(matIn)])

MEMEBody = """MOTIF %s
letter-probability matrix: nsites= %d
""" % (filename, totalNumBases)

if ((args.baseModification and args.baseModPosition) or args.tryAllCModsAtPos):
    modFreqMatrix = np.copy(freqMatrix)
    baseModPos = args.tryAllCModsAtPos or args.baseModPosition
    for b in (MOD_BASE_NAMES.keys() if args.tryAllCModsAtPos
              else args.baseModification):
        # zero all entries along the frequency matrix, for the given (row)
        # position, except that corresponding to the (column) index of the
        # modified base, which is set to unity.
        modFreqMatrix[(baseModPos - 1), ] = \
            np.zeros((1, freqMatrix.shape[1]))
        modFreqMatrix[(baseModPos - 1),
                      MOTIF_ALPHABET.index(b)] = 1
        with open((os.path.splitext(filename)[0] + '-' + MOD_BASE_NAMES[b] +
                   '.meme'), "a") as outFile:
            outFile.write(MEME_HEADER)
            outFile.write(MEMEBody)
            np.savetxt(outFile, modFreqMatrix, '%f', "\t")

output = StringIO()
np.savetxt(output, freqMatrix, '%f', "\t")
MEMEBody += output.getvalue()
output.close()
print(MEME_HEADER + MEMEBody)
