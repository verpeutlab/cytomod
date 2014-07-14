#!/usr/bin/env python

"""Computes the frequency of each nucleobase and outputs
a representative pie chart. Character frequency computation
code was adapted from: http://rosettacode.org/wiki/Letter_frequency#Python"""

from __future__ import print_function

__version__ = "$Revision: 0.01$"

import sys
import collections
import re
import gzip

import numpy as np

_STDIN_SPECIFIER = '-'
_DEFAULT_PLOT_NAME = 'nucleobaseFrequencies'
_EXPLODE_DISTANCE = 0.1

LINE_EXCLUSION_REGEX = '>'
CHAR_EXCLUSION_REGEX = '\n'

# Allow for substitution of some simple input to its corresponding unicode
_UNICODE_SUBS = ('---', u"\u2014"), ('--', u"\u2013")


def warn(*msg):
    """Emit a warning to STDERR."""
    print("Warning: ", *msg, file=sys.stderr)


def filecharcount(openfile, lineExclusionRegex, charExclusionRegex):
    """Counts characters characters, over the given file handle,
    returning a sorted collection. Excludes lines or characters matching
    the given respective exclusion regexes."""
    return sorted(collections.Counter
                  (c for l in openfile if
                   not re.search(lineExclusionRegex, l) for c in l if
                   not re.search(charExclusionRegex, c)).items())


def makePlot(charFreqs, plotPath, subtitle="", pie=False):
    """Creates a (bar or pie) chart from the provided, sorted,
    collection of character frequencies."""
    import matplotlib.pyplot as plt
    import prettyplotlib as ppl

    fig, ax = plt.subplots(1)
    if pie:
        explodeBases = []
        for c in zip(*charFreqs)[0]:
            explodeBases.append(_EXPLODE_DISTANCE) \
                if c in args.explodeBases else explodeBases.append(0)
        plt.pie(zip(*charFreqs)[1], explode=explodeBases,
                labels=zip(*charFreqs)[0], autopct='%1.1f%%')
    else:
        from brewer2mpl import diverging
        bmapType = diverging.Spectral
        bmapType.pop("max", None)  # remove the max pointer,
        # since if we need that we'll find that key eventually
        # Allocate a number of colours corresponding to the
        # the closest number of possible colours for the given class
        # to the number of distinct characters we have, plus one.
        # This is a naive, but sufficient, emulation of the
        # desired behaviour of obtaining the tightest
        # upper bound on the number of items in our list
        # TODO change this to read colour scheme directly
        # from the ExtreMEME mod_spec file or to use the
        # appropriate "RdYlBu" class with custom colours for A/T
        numColours = min(bmapType.keys(), key=lambda n:
                         abs(n - (len(charFreqs) + 1)))
        ppl.barh(ax, range(len(charFreqs)), np.array(zip(*charFreqs)[1]),
                 yticklabels=np.array(zip(*charFreqs)[0]),
                 color=bmapType[numColours].mpl_colors)
    # Tight plot boundaries pad for subtitle if it exists
    # Not currently used since padding for a subtitle is not functioning
    # fig.tight_layout(h_pad=100 if regionLabel else 0)
    plt.suptitle('Modified Genome Nucleobase Frenquencies',
                 bbox={'facecolor': '0.8', 'pad': 5}, fontsize=16)
    # Print the given subtitle and transform to corresponding unicode
    plt.title(reduce(lambda s, kv: s.replace(*kv), _UNICODE_SUBS, subtitle),
              fontsize=12)
    plt.savefig(plotPath + ('.png' if args.rasterize else '.pdf'))


import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file', required=True,
                    help="Input FASTA file for nucleobase frequency querying. \
                    Can either be an existing (optionally gzipped) FASTA file \
                    or '" + _STDIN_SPECIFIER + "' to read from STDIN.")
parser.add_argument('-E', '--explodeBases', default='',
                    help="Explode the given bases in the resulting pie chart \
                    output. The single character base representations to \
                    explode should be given in succession. For example, \
                    to explode all instances of 5hmC and 5fC, one should \
                    use: '-E hf'")
parser.add_argument('-R', '--rasterize', action='store_true',
                    help="Output a rasterized image, as opposed to the default, \
                    vectorized, output.")
parser.add_argument('-P', '--pie', action='store_true',
                    help="Generate a pie chart instead of the \
                    default horizontal bar chart")
parser.add_argument('-r', '--regionLabel', default='',
                    help="A label for the plotted region. Used as a subtitle.")
parser.add_argument('-o', '--outputPlotPath', default=_DEFAULT_PLOT_NAME,
                    help="A full path to where the plot should be saved. \
                    Any specified file extension will end up as a part \
                    of the filename itself and not as a genuine extension.")
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action="count")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

if args.explodeBases and not args.pie:
    warn("Explosion paramter ignored, since we are not making a pie chart.")

if args.file == _STDIN_SPECIFIER:
    f = sys.stdin
else:
    openHandler = gzip.open if args.file.endswith('.gz') else open
    f = openHandler(args.file, 'rb')
charFreqs = filecharcount(f, LINE_EXCLUSION_REGEX, CHAR_EXCLUSION_REGEX)
makePlot(charFreqs, args.outputPlotPath, args.regionLabel, args.pie)
