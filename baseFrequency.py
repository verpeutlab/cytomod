#!/usr/bin/env python

"""Computes the counts of each nucleobase and outputs
a representative horizontal bar plot chart with a logarithmic
scale, as well as another with the absolute percent frequencies
of only modified nucleobase. The character counts computation
code was adapted from: http://rosettacode.org/wiki/Letter_frequency#Python"""

from __future__ import with_statement, division, print_function

import sys
import re
import gzip
from collections import Counter

import numpy as np

import cUtils
__version__ = cUtils.VERSION

_STDIN_SPECIFIER = '-'
_DEFAULT_PLOT_NAME = 'nucleobaseFrequencies'
_PERCENTAGE_PLOT_SUFFIX = '-inclusionOnlyPercentagePlot'
_EXPLODE_DISTANCE = 0.1

LINE_EXCLUSION_REGEX = '>'
CHAR_EXCLUSION_REGEX = '\n'

# Default inclusion regex for output of modified base only plot
_DEFAULT_SELECTION_INCLUSION_REGEX = '[a-z0-9]'

# Allow for substitution of some simple input to its corresponding unicode
_UNICODE_SUBS = ('---', u"\u2014"), ('--', u"\u2013")


def warn(*msg):
    """Emit a warning to STDERR."""
    print("Warning: ", *msg, file=sys.stderr)


def orderBase(base):
    """Key function for ordering nucleobases.
    Sort by complemented modified bases, normal nucleobases, followed by
    "primary" modified nucleobases, in their natural order
    (i.e. by their order in the cytosine demethylation cascade).
    Ambiguous bases are placed after their last modified nucleobase.
    The result returned is the reverse of this ordering, since matplotlib
    barh plots in the reverse of the given key ordering."""
    # Account for ambiguious bases by mapping to the
    # last unequivocal modified base for ordering purposes
    ambigAdj = 0
    if base not in cUtils.getUnivocalModBases():
        # adjust value by one if ambiguous base, to ensure
        # that the ambigous bases come after the last
        # unambiguous modified base
        ambigAdj = 1
        base = cUtils.getLastUnivocalModBase(base)
    if base[0] in cUtils.MOD_BASES.values():
        # place modified bases in the reverse of their
        # natural order (by subtracting from # of mod bases)
        # and before everything else (so subtract large int)
        return (len(cUtils.MOD_BASES.values()) -
                cUtils.MOD_BASES.values().
                index(base[0]) - sys.maxsize - ambigAdj)
    elif cUtils.getMBMaybeFromComp(base[0]) in cUtils.MOD_BASES.values():
        # if complemented modified base, place in reverse
        # numeric order (so subtract from # of mod bases)
        return (len(cUtils.MOD_BASES.values()) -
                cUtils.MOD_BASES.values().
                index(cUtils.getMBMaybeFromComp(base[0])) - ambigAdj)
    elif cUtils.FULL_BASE_NAMES.get(base[0]):
        # if not mod base, place in ASCII order,
        # but before "primary" mod bases (so flip sign)
        return -ord(base[0])
    else:
        # if neither a canonical nucleobase (i.e. A/C/G/T),
        # nor a modified nucleobase, make sure it comes after the
        # canonical nucleobases (by subtracting the ASCII value of 'Z')
        return -ord(base[0]) - ord('Z')


def filecharcount(openfile, lineExclusionRegex, charExclusionRegex):
    """Counts characters characters, over the given file handle,
    returning a collection in natural sorted order (see orderBase()).
    Excludes lines or characters matching the given, respective,
    exclusion regexes."""
    return sorted(Counter
                  (c for l in openfile if
                   not re.search(lineExclusionRegex, l) for c in l if
                   not re.search(charExclusionRegex, c)).items(),
                  key=lambda baseAndFreq: orderBase(baseAndFreq[0]))


def makeCountPlot(charCounts, plotPath, subtitle="", pie=False, annPlot=False):
    """Creates a (bar or pie) chart from the provided, sorted,
    collection of character counts."""
    import matplotlib.pyplot as plt
    import prettyplotlib as ppl

    fig, ax = plt.subplots(1)
    # Scale in a way to allow appropriate plotting of both very large
    # and small values. This is done using a symmetrical logarithmic scale.
    # See: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.xscale
    ax.set_xscale('symlog')
    if pie:
        explodeBases = []
        for c in zip(*charCounts)[0]:
            explodeBases.append(_EXPLODE_DISTANCE) \
                if c in args.explodeBases else explodeBases.append(0)
        plt.pie(zip(*charCounts)[1], explode=explodeBases,
                labels=zip(*charCounts)[0], autopct='%1.1f%%')
    else:
        from brewer2mpl import diverging
        bmapType = diverging.Spectral
        bmapType.pop("max", None)  # Remove the max pointer,
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
                         abs(n - (len(charCounts) + 1)))
        ppl.barh(ax, range(len(charCounts)), np.array(zip(*charCounts)[1]),
                 yticklabels=np.array(zip(*charCounts)[0]),
                 color=bmapType[numColours].mpl_colors)
        plt.xlabel('Count')
        plt.ylabel('Nucleobase')
    # Tight plot boundaries pad for subtitle if it exists
    # Not currently used since padding for a subtitle is not functioning
    # fig.tight_layout(h_pad=100 if regionLabel else 0)

    plt.suptitle('Modified Genome Nucleobase Counts',
                 bbox={'facecolor': '0.8', 'pad': 5}, fontsize=16)
    # Print the given subtitle and transform to corresponding unicode
    plt.title(reduce(lambda s, kv: s.replace(*kv), _UNICODE_SUBS, subtitle),
              fontsize=12)
    plt.savefig(plotPath + ('.png' if args.rasterize else '.pdf'))


def makeSelFreqPlot(charFreqsP, plotPath, selectionInclusionRegex,
                    subtitle="", percentC=False, pie=False, annPlot=False):
    """Creates a (bar or pie) chart from the provided, sorted,
    collection of character frequencies.
    The frequencies (only if plotting a bar chart)
    are multiplied by 100 to plot them as percentages."""
    import matplotlib.pyplot as plt
    import prettyplotlib as ppl

    # TODO - refer to comment in makeCountPlot
    selectedCharFreqs = [(c, f*100) for c, f in charFreqsP
                         if re.search(selectionInclusionRegex, c)]

    fig, ax = plt.subplots(1)
    if pie:
        explodeBases = []
        for c in zip(*selectedCharFreqs)[0]:
            explodeBases.append(_EXPLODE_DISTANCE) \
                if c in args.explodeBases else explodeBases.append(0)
        plt.pie(zip(*selectedCharFreqs)[1], explode=explodeBases,
                labels=zip(*selectedCharFreqs)[0], autopct='%1.1f%%')
    else:
        import matplotlib.ticker as ticker
        from matplotlib.ticker import FormatStrFormatter
        from brewer2mpl import diverging

        bmapType = diverging.Spectral
        bmapType.pop("max", None)  # remove the max pointer
        numColours = min(bmapType.keys(), key=lambda n:
                         abs(n - (len(selectedCharFreqs) + 1)))
        ppl.barh(ax, range(len(selectedCharFreqs)),
                 np.array(zip(*selectedCharFreqs)[1]),
                 yticklabels=np.array(zip(*selectedCharFreqs)[0]),
                 color=bmapType[numColours].mpl_colors, annotate=annPlot)
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        # Add '%' symbol to each x-axis tick label
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%1.2f%%'))
        plt.xlabel('Percent of Total Cytosine Content' if percentC
                   else 'Percent of Total Genomic Content')
        plt.ylabel('Modified Nucleobase')
    plt.suptitle('Modified Nucleobase Percent Abundances',
                 bbox={'facecolor': '0.8', 'pad': 5}, fontsize=16)
    # Print the given subtitle and transform to corresponding unicode
    plt.title(reduce(lambda s, kv: s.replace(*kv), _UNICODE_SUBS, subtitle),
              fontsize=12)
    plt.savefig(plotPath +
                ('.png' if args.rasterize else '.pdf'))

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
parser.add_argument('-t', '--text', action='store_true',
                    help="Output the nucleobase frequencies (to STDOUT).")
parser.add_argument('-s', '--selectionInclusionRegex',
                    default=_DEFAULT_SELECTION_INCLUSION_REGEX,
                    help="A regular expression specifying which characters \
                    ought to be included for the absolute percentage plot. \
                    This defaults to include only modifified nucleobases \
                    (defined as any lower-case alphabetic character \
                    or any digit).")
parser.add_argument('-C', '--percentC', action='store_true',
                    help="In the frequency plot, use percent of total \
                    cytosines (guanines for complentary modified bases), \
                    instead of the default overall percentages.")
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action="count")
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

if args.explodeBases and not args.pie:
    warn("Explosion paramter ignored, since we are not making a pie chart.")

if args.file == _STDIN_SPECIFIER:
    fileH = sys.stdin
else:
    openHandler = gzip.open if args.file.endswith('.gz') else open
    fileH = openHandler(args.file, 'rb')
charCounts = filecharcount(fileH, LINE_EXCLUSION_REGEX, CHAR_EXCLUSION_REGEX)
charFreqs = 0
# TODO This could be made much more general (i.e. define a base or
# set of bases that this op occurs WRT and define a complentation
# operation for that base etc.).
if args.percentC:
    totalNumCs = sum((f for b, f in charCounts if b is 'C' or
                      cUtils.MODIFIES.get(b) == 'C'))
    totalNumGs = sum((f for b, f in charCounts if b is 'G' or
                      cUtils.MODIFIES.get(b) == 'G'))
    charFreqs = [(base, count / totalNumCs) for base, count in charCounts
                 if base is 'C' or cUtils.MODIFIES.get(base) == 'C']
    charFreqs += [(base, count / totalNumGs) for base, count in charCounts
                  if base is 'G' or cUtils.MODIFIES.get(base) == 'G']
else:
    totalNumChars = sum(zip(*charCounts)[1])
    charFreqs = [(base, count / totalNumChars) for base, count in charCounts]

if args.text:
    if args.verbose > 0:
        print("Program was invoked with:\n" + str(sys.argv[1:]) + "\n\n")
    # Print the frequencies, making sure to reverse the order for consistency
    print ('\n'.join(map(str, charFreqs[::-1])))

makeCountPlot(charCounts, args.outputPlotPath, args.regionLabel, args.pie)
makeSelFreqPlot(charFreqs, args.outputPlotPath + _PERCENTAGE_PLOT_SUFFIX,
                args.selectionInclusionRegex, args.regionLabel, args.percentC,
                args.pie)
