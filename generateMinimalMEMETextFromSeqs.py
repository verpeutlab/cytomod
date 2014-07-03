#!/usr/bin/env python

"""Creates a minimal MEME text file, suitable for FIMO, from an
input set of sequences (one raw sequence per line).
Accepts a single argument: the full path to the input file.
The background frequencies should be manually adjusted to use case.
NB: This script is not intended for general use, but rather was
created prioritizing coding efficiency, for a specific use case.
As such, it is not optimized towards ease-of-use."""

from __future__ import with_statement, division, print_function

import sys

import numpy as np
import numpy.testing as npt
from scipy.stats import itemfreq

BOTH_STRANDS = 0
STRANDS = '+ -' if BOTH_STRANDS else '+'
FROM_STR = 'custom'
# The alphabet frequencies used were (WRT Cs only):
# 4% 5mC or 5hmC, 0.1% 5hmC, 0.002% 5fC, and 0.0003% 5caC
# (Song et al. 2012, Ito et al. 2011, and Ehrlich and Wang 1981)
MOTIF_ALPHABET_BG_FREQUENCIES = \
    {'T': 0.292, 'A': 0.292, 'C': 0.1879885, 'G': 0.1879885,
     # (using mouse GC content of 41.6%) <- * in use *
     # From: NCBI Eukaryotic Genome Report File
     # (assuming uniform background, excepting above) <- not in use
     # {'T': 0.25, 'A': 0.25, 'C': 0.2299885, 'G': 0.2299885,
     'm': 0.0195, '1': 0.0195, 'h': 0.0005, '2': 0.0005, 'f': 0.00001,
     '3': 0.00001, 'c': 0.0000015, '4': 0.0000015}
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

filename = sys.argv[1]

motifs = np.loadtxt(filename, dtype=str)
motifChars = motifs.view('S1').reshape((motifs.size, -1))

totalNumBases = len(motifChars)

MEMEBody = """MOTIF %s
letter-probability matrix: nsites= %d
""" % (filename, totalNumBases)

for i in range(0, motifChars.shape[1]):
    motifCharsInts = motifChars[:, i].view(np.uint8)
    # NB: itemfreq internally uses bincount, so we must map to and from ints
    f = itemfreq(motifCharsInts)
    bases = f[:, 0].view('U1')
    baseFreqs = f[:, 1]
    # Append the letter frequency matrix
    MEMEBody += "\t".join(str(x) for x in
                          (baseFreqs[idx][0]/totalNumBases if len(idx[0])
                          else 0 for idx in (np.nonzero(bases == base) for base
                                             in MOTIF_ALPHABET))) + "\n"

print(MEME_HEADER + MEMEBody)
