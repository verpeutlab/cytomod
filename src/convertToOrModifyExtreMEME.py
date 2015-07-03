#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO refactor.
# TODO make variable names consistent (i.e. x_y, not xY).

from __future__ import with_statement, division, print_function

"""Prints (to STDOUT) a minimal MEME text file, suitable for FIMO,
from an input set of sequences (one raw sequence per line)
or an input PWM, PFM, TRANSFAC matrix, or JASPAR matrix.
The background frequencies should be manually adjusted
to use case. Modified output (if applicable) is always
written to a file, while the MEME output written to STDOUT
always correpsonds to the unmodified motif.
"""

import os
import re
import textwrap

try:
    from cStringIO import StringIO
except:
    from StringIO import StringIO

import numpy as np
import numpy.testing as npt

import cUtils

__version__ = cUtils.__version__

_mESC_BG_NAME = 'mESC'
_AML_BG_NAME = 'AML'
_DEFAULT_BG = _mESC_BG_NAME
_DEFAULT_HEMIMODARG = '+-'

_CONTEXT_FREQ_THRESHOLD = 0
_DIST_FROM_CEN = 1  # distance from centre to still be considered "central"

_ALL_BASE_CONTEXTS = '*'
_DELIM = "\t"
_ARG_DELIM = ','
BOTH_STRANDS = 0
STRANDS = '+ -' if BOTH_STRANDS else '+'
FROM_STR = 'custom'
MEME_MINIMAL_BG_REGEX_G_BG = 'background'
MEME_MINIMAL_BG_REGEX = """Background\s+letter\s+frequencies\s+\(from.+\):\s+
                           (?P<{}>(?:\w\s+\d+\.\d+\s*)+)"""\
                           .format(MEME_MINIMAL_BG_REGEX_G_BG)
MEME_MINIMAL_BG_BASE_REGEX_G_BASE = 'base'
MEME_MINIMAL_BG_BASE_REGEX = """(?P<{}>\w)\s+\d+\.\d+\s*"""\
    .format(MEME_MINIMAL_BG_BASE_REGEX_G_BASE)
MEME_MIN_REGEX_G_ID = 'motif_ID'
MEME_MIN_REGEX_G_NAME = 'motif_name'
MEME_MIN_REGEX_G_ALPH_LEN = 'alph_len'
MEME_MIN_REGEX_G_WIDTH = 'width'
MEME_MIN_REGEX_G_NUM_SITES = 'num_sites'
MEME_MIN_REGEX_G_E_VALUE = 'E_value'
MEME_MIN_REGEX_G_PWM = 'PWM'
MEME_MINIMAL_MOTIF_REGEX = """MOTIF\s+(?P<{}>[a-zA-Z0-9_.]+\s+)?
                              (?P<{}>\w+\s+)  # motif ID line

                              letter-probability\s+matrix:  # start properties
                              \s+alength=\s*(?P<{}>\d+)\s+
                              w=\s*(?P<{}>\d+)\s+
                              nsites=\s*(?P<{}>\d+)\s+
                              E=\s*(?P<{}>\d+(?:\.\d+(?:e(?:\+|-)\d+)?)?)

                              (?P<{}>\s*(?:(?:[01]\.\d+)\s*)+)"""\
                               .format(MEME_MIN_REGEX_G_ID,
                                       MEME_MIN_REGEX_G_NAME,
                                       MEME_MIN_REGEX_G_ALPH_LEN,
                                       MEME_MIN_REGEX_G_WIDTH,
                                       MEME_MIN_REGEX_G_NUM_SITES,
                                       MEME_MIN_REGEX_G_E_VALUE,
                                       MEME_MIN_REGEX_G_PWM)
MEME_MINIMAL_REGEX_FLAGS = re.M | re.X

MEME_header = """MEME version 4

ALPHABET "DNA with covalent modifications"
A "Adenine" 8510A8 ~ T "Thymine" A89610
C "Cytosine" A50026 ~ G "Guanine" 313695
m "5-Methylcytosine" D73027 ~ 1 "Guanine:5-Methylcytosine" 4575B4
h "5-Hydroxymethylcytosine" F46D43 ~ 2 "Guanine:5-Hydroxymethylcytosine" 74ADD1
f "5-Formylcytosine" FDAE61 ~ 3 "Guanine:5-Formylcytosine" ABD9E9
c "5-Carboxylcytosine" FEE090 ~ 4 "Guanine:5-Carboxylcytosine" E0F3F8
z = Cmhfc
9 = G1234
y = Cfc
8 = G34
x = mh
7 = 12
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

strands: {}

Background letter frequencies (from {}):
""".format(STRANDS, FROM_STR)

INVALID_PWM_MSG = """Invalid PWM; rows do not sum to one."""

MODIFYING_INVALID_PWM_MSG = """The resultant PWM would not have rows summing to
                               to one. This is likely caused by an attempt to
                               modify an already modified motif, which cannot
                               be performed."""


def _modificationOptionType(arg):
    if (arg in cUtils.MOD_BASE_NAMES.keys() or
            arg in cUtils.complement(cUtils.MOD_BASE_NAMES.keys())):
        return str(arg)
    else:
        raise argparse.ArgumentTypeError("""The provided target nucleobase of
                                            the modification ({}) is not
                                            a recognized modified nucleobase.
                                         """.format(arg))


def _modifyPositionOptionType(arg):
    permittedCharList = ['c', 'A', 'T', 'G', 'C']
    if arg in permittedCharList:
        return str(arg)
    elif arg.isdigit():
        return int(arg)
    else:
        raise argparse.ArgumentTypeError("""The provided option must contain
                                            an integer or one of: {}.
                                         """.format(permittedCharList))


def die(msg):
    cUtils.die(msg, os.path.basename(__file__))


def warn(msg):
    if not args.noWarnings:
        cUtils.warn(msg, os.path.basename(__file__))


def checkPWMValidity(matrix, invalid_msg=INVALID_PWM_MSG):
    if not np.allclose(np.sum(matrix, axis=1), 1):
        die(invalid_msg)


def getAlteredSlice(slice_to_inc, slice_max, operation, value):
    if slice_to_inc == slice(None):
        return slice_to_inc
    else:
        slice_stop_plus_one = (operation(getattr(slice_to_inc, 'stop'), value)
                               if (operation(getattr(slice_to_inc, 'stop'),
                                   value)) <= slice_max else
                               slice_max - 1)

        if getattr(slice_to_inc, 'start') is None:
            return slice(slice_stop_plus_one)
        else:
            return slice(operation(getattr(slice_to_inc, 'start'), value),
                         slice_stop_plus_one,
                         getattr(slice_to_inc, 'step'))


def output_motif(freq_matrix, output_descriptor, motif_name,
                 motif_alphabet, numSites, EValue, motif_filename):
    if args.revcomp:
        output_descriptor += '-revcomp'
        motif_name += '-revcomp'
        freq_matrix = np.flipud(freq_matrix)  # reverse
        for base_idx, base in enumerate(motif_alphabet):  # complement
            if base_idx >= len(motif_alphabet) / 2:
                break
            complement_idx = motif_alphabet.index(cUtils.complement(base))
            freq_matrix[:, [base_idx, complement_idx]] = \
                freq_matrix[:, [complement_idx, base_idx]]

    if args.skipMotifPortions:
        output_descriptor += '-skip' + args.skipMotifPortions[0][0]
        skipRows = list()
        for skipMotifPortion in args.skipMotifPortions[0][0].split(','):
            match = re.search('(\d*):(\d*)', str(skipMotifPortion).strip('[]'))
            if match:
                start = int(match.group(1)) if match.group(1) else 0
                end = int(match.group(2)) if match.group(2) \
                    else freq_matrix.shape[0]
                if (start > freq_matrix.shape[0] - 1):
                    die("""The provided start index ({}) of the motif rows to
                    skip equals or exceeds the motif length.""".format(start))
                if (start > end):
                    die("""The start index ({}) of the motif rows to skip
                    cannot be greater than the end ({}).""".format(start, end))
                if (end > freq_matrix.shape[0]):
                    die("""The provided end index ({}) of the motif rows to
                    skip exceeds the motif length.""".format(end))
                if ((end - start) <= 1):
                    warn("""The provided skip interval ({}) will only cause one
                    motif row to be excluded. Recall that intervals provided
                    in colon notation are half-open and that a single row
                    to skip can be specified by only providing that row number
                    alone.""".format(end - start))
                skipRows.extend(range(start, end))
        freq_matrix = np.delete(freq_matrix, skipRows, 0)

    # add an extra final row to be removed later, to make computations simpler
    freq_matrix = np.vstack([freq_matrix, np.zeros(freq_matrix.shape[1])])

    totalNumBases = freq_matrix.shape[0]

    MEMEBody = textwrap.dedent("""\
               MOTIF {}\n
               letter-probability matrix: alength= {} w= {} nsites= {} E= {}\n
               """.format(motif_name, freq_matrix.shape[1],
                          totalNumBases, numSites, EValue))

    modCFracs = (args.baseModificationAtAllModifiablePosFractions or
                 args.modAllFractions or args.modCFractions)
    modGFracs = (args.baseModificationAtAllModifiablePosFractions or
                 args.modAllFractions or args.modGFractions)

    if ((args.baseModification and args.baseModPosition)
            or args.tryAllCModsAtPos or
            args.baseModificationAtAllModifiablePosFractions):
        warn("""NB: modified motifs are output to a file.
                The motif printed to STDOUT is the unmodified version.""")
        if (args.baseModPosition and isinstance(args.baseModPosition, int) and
                args.baseModPosition >= totalNumBases):
            die("""The provided modification position ({}) exceeds the
                   motif's final position ({}).
                """.format(args.baseModPosition, totalNumBases - 1))

        mod_base_context = _ALL_BASE_CONTEXTS
        baseModPos = args.tryAllCModsAtPos or args.baseModPosition
        # index is either positional, comprises all nucleobases, or contextual
        mod_base_index = slice(baseModPos - 1, baseModPos) \
            if (baseModPos and isinstance(baseModPos, int)
                and baseModPos != cUtils._PARAM_A_CONST_VAL) else slice(None)

        if args.baseModPosition and not isinstance(args.baseModPosition, int):
            if args.baseModPosition == 'c':
                mod_base_context = 'G'
                # would add one to centre since indexed from 1, but do not
                # due to the extra row accounting for this
                centre = totalNumBases//2
                central_start = centre - _DIST_FROM_CEN
                # +1 extra on stop, since slice is [start, stop)
                central_end = centre + _DIST_FROM_CEN + 1
                mod_base_index = slice(central_start if central_start >= 0
                                       else 0, central_end if central_end
                                       <= totalNumBases else totalNumBases - 1)
            else:
                mod_base_context = args.baseModPosition

        if args.baseModificationAtAllModifiablePosFractions:
            output_descriptor += '-allPos'
        elif baseModPos:
            output_descriptor += '-P' + str(baseModPos)

        for b in (cUtils.MOD_BASE_NAMES.keys() if args.tryAllCModsAtPos
                  else (args.baseModification or
                        args.baseModificationAtAllModifiablePosFractions)):
            modfreq_matrix = np.copy(freq_matrix)
            if ((b not in cUtils.COMPLEMENTS.keys()) and
                    (b not in cUtils.COMPLEMENTS.values())):
                warn("""Base {} provided in the background is not a
                        currently supported modified nucleobase. It
                        has, accordingly, been skipped.""".format(b))

            def _modifyMatrixPortion(matrix, mod_base_index,
                                     primary_base_to_mod, target_modified_base,
                                     mod_base_context, mod_fractions,
                                     hemimodifyOnly, old_matrix=None):
                if old_matrix is None:
                    old_matrix = matrix
                if not mod_fractions:
                    only_target_base_at_pos = \
                        np.zeros(matrix.shape[1])
                    only_target_base_at_pos[motif_alphabet.
                                            index(target_modified_base)] = 1
                    only_target_base_at_pos_comp = \
                        np.zeros(only_target_base_at_pos.shape)
                    only_target_base_at_pos_comp[motif_alphabet.
                                                 index(cUtils.complement
                                                       (target_modified_base))
                                                 ] = 1

                context_len = matrix[mod_base_index, motif_alphabet.
                                     index(target_modified_base)].shape
                if mod_base_context == _ALL_BASE_CONTEXTS:
                    correct_context = np.ones(context_len, dtype=bool)
                    correct_context_minus = correct_context
                else:
                    # create a mask which ensures that modifications are
                    # only permitted in their valid genomic context
                    # NB: the position of the insertion of 'False' depends upon
                    # the portion of the conjunction (i.e. [:-1] vs. [1:])
                    # such that it always inserts on the sequence's termination
                    correct_context = \
                        np.logical_and(
                            np.insert(old_matrix[mod_base_index,
                                      motif_alphabet.
                                      index(primary_base_to_mod)], -1, [0])
                            [:-1]
                            > _CONTEXT_FREQ_THRESHOLD,
                            np.insert(old_matrix[mod_base_index,
                                      motif_alphabet.
                                      index(mod_base_context)]
                                      > _CONTEXT_FREQ_THRESHOLD, -1, [0])[1:])
                    correct_context_minus = np.zeros(context_len, dtype=bool)

                    if '-' in hemimodifyOnly:
                        # shift Boolean values right, clipping them
                        correct_context_minus = \
                            np.insert(correct_context[:-1], 0, [0])
                    if '+' not in hemimodifyOnly:
                        # shift Boolean values right, clipping them
                        correct_context = \
                            np.zeros(context_len, dtype=bool)

                if mod_fractions:
                    for iter, fun in enumerate([lambda x: x,
                                               cUtils.complement]):
                        context_to_use = (correct_context if iter == 0
                                          else correct_context_minus)
                        # NB: this definition precludes further modification of
                        # modified bases, since it presumes residual
                        # frequencies are in the primary_base_to_mod column
                        matrix_cur_view = \
                            matrix[mod_base_index, motif_alphabet.
                                   index(fun(target_modified_base))].view()
                        matrix_modified_view = \
                            matrix[mod_base_index, motif_alphabet.
                                   index(fun(primary_base_to_mod))].view()

                        # transfer primary_base_to_mod frequency to the
                        # corresponding target_modified_base entry.

                        # make the modification
                        matrix[mod_base_index,
                               motif_alphabet.index(fun(target_modified_base))] = \
                            np.where(context_to_use, matrix_modified_view,
                                     matrix_cur_view)

                        # zero out other base frequencies at modification pos.
                        matrix[mod_base_index,
                               motif_alphabet.index(fun(primary_base_to_mod))] = \
                            np.where(context_to_use,
                                     (0 if baseModPos else
                                      np.zeros(matrix.shape[0])),
                                     matrix_modified_view)
                else:
                    # zero all entries along the frequency matrix,
                    # for the given (row) position, except that corresponding
                    # to the (column) index of the modified base,
                    # which is set to unity.
                    # use the context mask as a sub-array of the correct length
                    # (i.e. preceded and succeeded by 'False' to make up the
                    # length of the motif) to accomplish this.
                    if mod_base_index != slice(None):
                        def _resize_context(context):
                            return np.concatenate(
                                (np.zeros(getattr(mod_base_index,
                                          'start'), dtype=bool), context,
                                 np.zeros(getattr(mod_base_index,
                                          'stop') - context_len[0])))

                        correct_context = _resize_context(correct_context)
                        correct_context_minus = \
                            _resize_context(correct_context_minus)

                    matrix[correct_context.astype('bool'), ] = \
                        only_target_base_at_pos
                    if mod_base_context != _ALL_BASE_CONTEXTS:
                        matrix[correct_context_minus.astype('bool'), ] = \
                            only_target_base_at_pos_comp
                
                print() # XXX
                print(mod_base_context) # XXX
                print() # XXX
                print(matrix) # XXX
                print() # XXX
                
                checkPWMValidity(matrix[:-1], MODIFYING_INVALID_PWM_MSG)

            _modifyMatrixPortion(modfreq_matrix, mod_base_index, 'C', b,
                                 mod_base_context,
                                 True if modCFracs else False,
                                 args.hemimodifyOnly)
            # the modification for 'G' needs to use the unmodified matrix for
            # context computations, unless in non-fractional mode, in which
            # case it needs to use the already modified matrix to prevent
            # changing the target modification to its complement.
            # integer positional modifications only run the method once
            if not isinstance(args.baseModPosition, int):
                _modifyMatrixPortion(modfreq_matrix, mod_base_index, 'G',
                                     cUtils.complement(b),
                                     mod_base_context,
                                     True if modGFracs else False,
                                     args.hemimodifyOnly,
                                     freq_matrix if modGFracs else None)

            modfreq_matrix = modfreq_matrix[:-1]  # remove extra final row
            checkPWMValidity(modfreq_matrix)
            with open((os.path.basename(os.path.splitext(motif_filename)[0]) +
                       '-' + cUtils.MOD_BASE_NAMES[cUtils.
                                                   getMBMaybeFromComp(b)]
                       + output_descriptor +
                       ('-mCFracs' if modCFracs else '') +
                       ('-mGFracs' if modGFracs else '') +
                       '.meme'), 'a') as outFile:
                outFile.write(MEME_header)
                outFile.write(MEMEBody)
                np.savetxt(outFile, modfreq_matrix, '%f', _DELIM)

    freq_matrix = freq_matrix[:-1]  # remove extra final row
    checkPWMValidity(freq_matrix)
    output = StringIO()
    np.savetxt(output, freq_matrix, '%f', _DELIM)
    MEMEBody += output.getvalue()
    output.close()
    return MEMEBody


import argparse
parser = argparse.ArgumentParser()
inputFileGroupTitle = \
    parser.add_argument_group(title="Input File", description="Input set(s) \
                              of sequences or matrices from which to create a \
                              MEME minimal text output file. \
                              Ensure to also use arguments '-a' and '-S' \
                              as needed for the input file used.\
                              Assumes that only a single matrix \
                              is input per input argument, except if the \
                              input is itself already in minimal MEME format \
                              or if sequences are provided, \
                              in which case multiple input motifs \
                              are supported. Only a single matrix format \
                              can be provided (but a matrix format \
                              can be combined with '-s').")
inputFile = inputFileGroupTitle.add_mutually_exclusive_group(required=False)
inputFileGroupTitle.add_argument('-s', '--inSeq', type=str,
                                 help="File containing \
                                 an input set of raw sequences or a single \
                                 sequence provided as the argument.\
                                 Multiple files or sequences can be provided \
                       delimited by \"{}\".".format(_ARG_DELIM))
inputFile.add_argument('-p', '--inPWMFile', type=str, help="File containing \
                       an input whitespace-delimited PWM (i.e. frequency \
                       matrix). The file must only contain floats \
                       (but see '-A'), lexicographically ordered \
                       (i.e. A, C, G, T).")
inputFile.add_argument('-c', '--inPFMFile', type=str, help="File containing \
                       an input whitespace-delimited PFM (i.e. a count \
                       matrix; e.g. a hPDI matrix). The file must only \
                       contain floats (but see '-a'), lexicographically \
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
inputFile.add_argument('-m', '--inMEMEFile', type=str, help="File containing \
                       an input set of motifs in minimal MEME format.")
modBaseSpecifiers = parser.add_mutually_exclusive_group()
modBasePositions = parser.add_mutually_exclusive_group()
modBaseSpecifiers.add_argument('-M', '--baseModification',
                               type=_modificationOptionType,
                               choices=(cUtils.MOD_BASE_NAMES.keys() +
                                        cUtils.complement(cUtils.
                                                          MOD_BASE_NAMES.
                                                          keys())),
                               help="Modify the motif to use the modified base \
                               provided in this argument at the position \
                               specified by '-P'. The resultant motif will \
                               use the given modified base with 100%% \
                               frequency at the specified positions. \
                               This will cause the program to write a file \
                               (as opposed to the usual output to STDOUT) \
                               for the given modification.")
modBasePositions.add_argument('-P', '--baseModPosition',
                              type=_modifyPositionOptionType,
                              nargs='?', const=cUtils._PARAM_A_CONST_VAL,
                              help="The position at which to modify the motif \
                              (using the base specified by '-M'), \
                              * indexed from 1 *. \
                              Note that for explicitly specified positions, \
                              this does not check that \
                              modifications made make biological sense \
                              (i.e. a thymine could be changed into a \
                              5-methylcytosine using this option) \
                              Alternatively, 'c' can be provided \
                              to indicate that the position \
                              should be automatically determined by finding \
                              the centre CpG dinucleotide, in which case only \
                              the cytosine of the CpG will be modified. \
                              Furthermore, a single primary nuleobase \
                              (i.e. A/T/G/C), X (with complement 0), can be \
                              provided to modify only cytosines in \
                              a CpX context and guanines in a Gp0 context \
                              (unless only hemi-modification is being \
                              performed, via '-H').\
                              Non-numeric arguments (i.e. context-based) \
                              are currently experimental and may not \
                              function as intended.")
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
parser.add_argument('-S', '--skipMotifPortions', nargs='+', action='append',
                    help="Skip a portion of a motif. This might be used \
                    if one wishes to limit a motif to its highly \
                    informative portion. Generally, one would only \
                    remove bases from the start or end of a motif, \
                    but this option also permits removal of the middle \
                    of a motif. Specify the portion to remove via: \
                    space-delimited row indicies (from 0) \
                    or via \"start:end\". In the latter case, \
                    the given interval is interpreted as the \
                    half-open interval: [start, end). \
                    Intervals to skip are comma-delimited. \
                    For example, to include only positions 3 to 8, \
                    inclusive and indexed from 1, one would use: \
                    \"-S 0:2,6:\". Empty start or end intervals \
                    specify that the given interval should either start \
                    from 0 or go until the end of the motif.")
parser.add_argument('--revcomp', help="Reverse complement the input motif.",
                    action='store_true')
parser.add_argument('-N', '--notNucleobases', help="Do not use the \
                    provided (comma-delimited) list of nucleobases. \
                    These may be specified by either their abbreviation \
                    or base. However, the base specified must be a \
                    \"primary\" and not a complemented nucleobase \
                    (i.e. 'm' and not '1'). For example, \"c,5hmC\", \
                    would use neither 5-carboxylcytosine nor \
                    5-hydroxymethylcytosine.")
parser.add_argument('--modCFractions', action='store_true',
                    help="Modify fractions of cytosines instead of setting \
                    the modified base frequency to 1. Has no effect with \
                    options already specifying this behaviour (e.g. '-A').\
                    This option is silently ignored if no modifications \
                    are requested.")
parser.add_argument('--modGFractions', action='store_true',
                    help="Modify fractions of gunanines instead of setting \
                    the modified base frequency to 1. Has no effect with \
                    options already specifying this behaviour (e.g. '-A').\
                    This option is silently ignored if no modifications \
                    are requested.")
parser.add_argument('-F', '--modAllFractions', action='store_true',
                    help="Convinience option to set both \
                    '--modCFractions' and '--modGFractions'. \
                    This will modify fractions of both cytosine and guanine \
                    irrespective of their individual settings. \
                    Has no effect with \
                    options already specifying this behaviour (e.g. '-A').\
                    This option is silently ignored if no modifications \
                    are requested.")
parser.add_argument('-H', '--hemimodifyOnly', nargs='?',
                    default=_DEFAULT_HEMIMODARG,
                    const=_DEFAULT_HEMIMODARG, choices=['+-', '-+', '+', '-'],
                    help="Only perform hemi-modification. This option \
                    affects any other option which performs any \
                    modifications by only modifying a single strand \
                    (i.e. either cytosine or guanine bases, exclusively). \
                    If no argument is provided \
                    (or if either '+-' or '-+' is provided), the number of \
                    motifs output by any affected option will be doubled, \
                    with one motif being output with only cytosine bases \
                    modified and the other with only guanine bases modified. \
                    Alternatively, '+' or '-' can be provided as an \
                    argument, in which case only cytosine ('+') or guanine \
                    ('-') bases will be modified, respectively. \
                    This argument does not affect explicitly specified \
                    positional modifications (i.e. '-P' used with a \
                    specific motif position).")
parser.add_argument('-b', '--background', default=_DEFAULT_BG,
                    help="The background model to \
                    add to the output file for use by MEME. \
                    The background model must be a zero-order Markov model. \
                    Either a background model type can be specified or \
                    a file can be provided. The background model \
                    type can be one of: 'mESC' or 'AML'. \
                    It is important to note that the AML background model \
                    does not include 5caC and has its 5fC value estimated \
                    from a melanoma cell line (WM-266-4). \
                    If a file is provided, it must be a tab-delimited \
                    file, containing lines of exactly two columns, \
                    the first being the nucleobase(s) whose background \
                    frequency is being specified (e.g. 'A') \
                    and the second being its frequency (e.g. 0.292). \
                    All lines starting with '#' in the file will be ignored. \
                    If this argument is not provided, it will default to \
                    the mESC background model.")
parser.add_argument('-W', '--noWarnings', help="Disable warnings.",
                    action='store_true')
parser.add_argument('-v', '--verbose', help="increase output verbosity",
                    action='count')
parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)
args = parser.parse_args()

output_descriptor = "({})".format(args.hemimodifyOnly)

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

filename = (args.inPWMFile or args.inPFMFile or args.inTRANSFACFile or
            args.inTRANSFACFreqFile or args.inJASPARFile or args.inMEMEFile)

if not (args.inSeq or filename):
    die("""You must provide some input, by using one of the \"Input File\"
           arguments.""")

if (args.modCFractions and ('+' not in args.hemimodifyOnly)):
    warn("""Cytosine bases will not be modified, since only negative-strand
            modification is permitted in this hemi-modification mode.""")

if (args.modGFractions and ('-' not in args.hemimodifyOnly)):
    warn("""Guanine bases will not be modified, since only negative-strand
            modification is permitted in this hemi-modification mode.""")

if (isinstance(args.baseModPosition, int) and
        args.hemimodifyOnly != _DEFAULT_HEMIMODARG):
    die("""Hemi-modification cannot be specified with an explicit
           modification position (motif position {}).
        """.format(args.baseModPosition))

motifAlphBGFreqs = ()
if args.background:
    if os.path.isfile(args.background):
        with open(args.background) as bgFile:
            for bFreq in bgFile:
                (base, freq) = bFreq.split(_DELIM)
                motifAlphBGFreqs[int(base)] = freq
    else:
        if args.background == _mESC_BG_NAME:
            # The alphabet frequencies used are (WRT Cs only):
            # 2.95% 5mC, 0.055 ± 0.008% 5hmC, 0.0014 ± 0.0003% 5fC,
            # and 0.000335% 5caC
            # (Respectively: Ito et al. 2011, Booth et al. 2014,
            # Booth et al. 2014, and Ito et al. 2011)
            motifAlphBGFreqs = \
                {'T': 0.292, 'A': 0.292, 'C': 0.201745991, 'G': 0.201745991,
                 # (using mouse GC content of 41.6%)
                 # From: NCBI Eukaryotic Genome Report File
                 'm': 0.006136, '1': 0.006136, 'h': 0.0001144, '2': 0.0001144,
                 'f': 0.000002912, '3': 0.000002912,
                 'c': 0.000000697, '4': 0.000000697}
        elif args.background == _AML_BG_NAME:
            # The alphabet frequencies used were (WRT all bases):
            # 2.91 ± 0.11% 5mC and 0.039% 5hmC.
            # (Respectively: Liu et al. 2007 and Kroeze et al. 2014)
            # 5fC at 0.0021812% was estimated from the 5fC to 5hmC ratio
            # within the melanoma cell line WM-266-4,
            # which was analyzed by Liu S. et al.
            motifAlphBGFreqs = \
                {'T': 0.295, 'A': 0.295, 'C': 0.190244094, 'G': 0.190244094,
                 # (using human GC content of 41.0%)
                 # From: NCBI Eukaryotic Genome Report File
                 'm': 0.01455, '1': 0.01455, 'h': 0.000195, '2': 0.000195,
                 'f': 0.000010906, '3': 0.000010906}

# The zero-order Markov model must sum to unity to be sensical
npt.assert_allclose([sum(motifAlphBGFreqs.itervalues())], [1])

if args.notNucleobases:
    for base in args.notNucleobases.split(_ARG_DELIM):
        if base in cUtils.MOD_BASES:
            base = cUtils.MOD_BASES[base]
        elif base not in cUtils.MOD_BASE_NAMES:
            die(""""Only "primary" modified nucleobases can be
            specified for removal.""")
        motifAlphBGFreqs[base] += \
            motifAlphBGFreqs[base]
        motifAlphBGFreqs[cUtils.complement(base)] += \
            motifAlphBGFreqs[cUtils.complement(base)]
        del motifAlphBGFreqs[base]
        del motifAlphBGFreqs[cUtils.complement(base)]
        # remove excluded nucleobase from alphabet
        MEME_header = re.sub('.*' + cUtils.FULL_MOD_BASE_NAMES[base]
                             + '.*\n', '', MEME_header)
        # remove excluded nucleobase from any containing ambiguity codes
        MEME_header = re.sub('^(. = .*)(?:' + base + '|'
                             + cUtils.complement(base) + ')(.*\n)',
                             '\g<1>\g<2>', MEME_header, flags=re.MULTILINE)
    # remove any ambiguity codes that are now empty
    MEME_header = re.sub('(. = \n)', '', MEME_header, flags=re.MULTILINE)

# The MEME Suite uses ASCII ordering for custom alphabets
# This is the natural lexicographic sorting order, so no "key" is needed
motif_alphabet = sorted(list(motifAlphBGFreqs.keys()))

motifs_to_output = ''

motif_alphabet_bg_freq_output = \
    ' '.join([str(k) + ' ' + str(v) for k, v
             in iter(sorted(motifAlphBGFreqs.iteritems()))])
MEME_header += motif_alphabet_bg_freq_output + "\n\n"

if args.inSeq:
    for seqFileOrStr in args.inSeq.split(_ARG_DELIM):
        if (os.path.isfile(seqFileOrStr)):
            # NB: min dimensionality of 1 is needed for the character view
            motifs = np.loadtxt(seqFileOrStr, dtype=str, ndmin=1)
            motifChars = motifs.view('S1').reshape((motifs.size, -1))
        else:
            motifChars = np.expand_dims(np.array(list(seqFileOrStr)), axis=0)
        totalNumBases = len(motifChars)

        freq_matrix = np.zeros((motifChars.shape[1], len(motifAlphBGFreqs)))
        for i in range(0, motifChars.shape[1]):
            base = motifChars[:, i][0]
            unambig_base_s = (cUtils.IUPAC_BASES.get(base) or
                              cUtils.AMBIG_MOD_BASES.get(base) or base)
            # unambig_base_s should only contain core bases, which are those
            # that have assigned base colours
            if not all(unambig_base in cUtils.BASE_COLOURS.keys()
                       for unambig_base in unambig_base_s):
                die("""Unrecognized base "{}".""".format(base))
            bases, baseFreqs = np.unique(unambig_base_s, return_counts=True)
            # NB: /= appears to perform //= despite future import statement
            baseFreqs = baseFreqs / len(baseFreqs)
            # Append to the letter frequency matrix
            matIn = StringIO(_DELIM.join(str(x) for x in
                                         (baseFreqs[idx][0]/totalNumBases if
                                         len(idx[0]) else 0 for idx in
                                         (np.nonzero(bases == base) for
                                          base in motif_alphabet))))
            freq_matrix[i] = np.loadtxt(matIn)
        motifs_to_output += (output_motif(freq_matrix, output_descriptor,
                                          os.path.basename(seqFileOrStr),
                                          motif_alphabet, 0, '',
                                          seqFileOrStr)
                             + "\n")

if filename:  # PWM or PFM
    # process the matrix for MEME format compatibility
    # If needed, skip columns or rows accordingly and/or
    # use unpack to transpose the matrix.
    ncols = 0
    with open(filename, 'rb') as inFile:
        if args.inMEMEFile:
            input_file = inFile.read()
            bgIter = re.findall(MEME_MINIMAL_BG_BASE_REGEX,
                                re.search(MEME_MINIMAL_BG_REGEX, input_file,
                                          flags=MEME_MINIMAL_REGEX_FLAGS)
                                .group(MEME_MINIMAL_BG_REGEX_G_BG),
                                flags=MEME_MINIMAL_REGEX_FLAGS)
            motifIter = re.finditer(MEME_MINIMAL_MOTIF_REGEX, input_file,
                                    flags=MEME_MINIMAL_REGEX_FLAGS)
        if args.inMEMEFile:
            motif_alphabet = list(bgIter)

        if args.annotated:
                ncols = len(inFile.readline().split())
                inFile.seek(0)
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
                                 if (args.annotated or args.inJASPARFile)
                                 else None)
        if args.inPFMFile or args.inTRANSFACFile or args.inJASPARFile:
            """This function transforms a given PFM column to a PWM column.
            That is, it computes the frequency of each element of the input
            1D array (which is itself a count), and replaces the count with
            its computed frequency. The function also returns the resultant
            array, despite its having been modified in place."""
            def _computeFreqFromCountSlice(countsForBase):
                sum = np.sum(countsForBase)
                for count in np.nditer(countsForBase, op_flags=['readwrite']):
                    count[...] = count / sum
                return countsForBase  # needed to use numpy.apply_along_axis
            np.apply_along_axis(_computeFreqFromCountSlice, 1, csvData)

        if args.inMEMEFile:
            for match in motifIter:
                # behave as if read from CSV to minimize changes to other code
                freq_matrix = np.fromstring(match.group(MEME_MIN_REGEX_G_PWM),
                                            dtype=float, sep=' ') \
                    .reshape((int(match.group(MEME_MIN_REGEX_G_WIDTH)),
                             int(match.group(MEME_MIN_REGEX_G_ALPH_LEN))))
                motif_name = (match.group(MEME_MIN_REGEX_G_ID).strip() + ' ' +
                              match.group(MEME_MIN_REGEX_G_NAME).strip())
                numSites = int(match.group(MEME_MIN_REGEX_G_NUM_SITES))
                EValue = match.group(MEME_MIN_REGEX_G_E_VALUE)
                motifs_to_output += (output_motif(freq_matrix,
                                                  output_descriptor,
                                                  motif_name, motif_alphabet,
                                                  numSites, EValue, filename) +
                                     "\n")
        else:
            freq_matrix = np.hstack((np.zeros((csvData.shape[0],
                                    motif_alphabet.index('A'))), csvData,
                                    np.zeros((csvData.shape[0],
                                             (len(motif_alphabet) - 1) -
                                             motif_alphabet.index('T')))))
            motifs_to_output += output_motif(freq_matrix, output_descriptor,
                                             os.path.basename(filename),
                                             motif_alphabet, 0, '', filename)

print(MEME_header + motifs_to_output.strip())
