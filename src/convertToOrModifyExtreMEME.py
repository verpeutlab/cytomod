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

from collections import OrderedDict

import numpy as np
import numpy.testing as npt

import cUtils

__version__ = cUtils.__version__

_mESC_BG_NAME = 'mESC'
_AML_BG_NAME = 'AML'
_DEFAULT_BG = _mESC_BG_NAME
_DEFAULT_HEMIMODARG = '+-'

_PSEUDO_EVALUE_TO_USE = 9999999  # used when the E-value is unknown
_CONTEXT_FREQ_THRESHOLD = 0
_DIST_FROM_CEN = 1  # distance from centre to still be considered "central"
_MAX_BASE_COMP_DIFF_FOR_SIMILARITY = 0.05  # for '--onlyNonNegChange' parameter

_ALL_BASE_CONTEXTS = '*'
_DELIM = "\t"
_ARG_DELIM = ','
_LIST_OF_PRIMARY_MOD_BASES_AND_ABBREVS = \
    list(reduce(lambda abbrev, base: abbrev + base, cUtils.MOD_BASES.items()))

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
MEME_MINIMAL_REGEX_G_ALPH = 'alph'
MEME_MINIMAL_ALPH_REGEX = """(?P<{}>ALPHABET(?:.|\\n)*END\sALPHABET)"""\
    .format(MEME_MINIMAL_REGEX_G_ALPH)
MEME_MIN_REGEX_G_ID = 'motif_ID'
MEME_MIN_REGEX_G_NAME = 'motif_name'
MEME_MIN_REGEX_G_ALPH_LEN = 'alph_len'
MEME_MIN_REGEX_G_WIDTH = 'width'
MEME_MIN_REGEX_G_NUM_SITES = 'num_sites'
MEME_MIN_REGEX_G_E_VALUE = 'E_value'
MEME_MIN_REGEX_G_PWM = 'PWM'
MEME_MINIMAL_MOTIF_REGEX = """MOTIF\s+(?P<{}>[^\s]+\s+)?
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


def _omitModBaseOptionType(arg):
    bases_list = arg.split(_ARG_DELIM)
    for base in bases_list:
        if (base not in _LIST_OF_PRIMARY_MOD_BASES_AND_ABBREVS):
            raise argparse.ArgumentTypeError("""The provided primary modified
                                                nucleobase ({}) is not
                                                recognized.""".format(arg))
    return bases_list


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


def checkPWMValidity(matrix, motif_name,
                     invalid_msg=INVALID_PWM_MSG, allow_cont=False):
    """Checks the PWM validity. Currently this just checks that all rows
       sum to 1. Terminates the program if this is not the case, unless
       specified to allow_cont, in which case it returns False if invalid
       and True if valid."""
    if not np.allclose(np.sum(matrix, axis=1), 1):
        additional_msg = " Occured for motif: " + motif_name
        if allow_cont:
            warn(invalid_msg + additional_msg)
            return False
        else:
            die(invalid_msg + additional_msg)
    return True


def isMatrixSufficientlyDifferent(freq_matrix, modfreq_matrix, index_arr):
    """Outputs true iff the input matrices are considered sufficiently
       distinct. This is meant to be used to remove results that are not
       worth analyzing and will be empirically-tuned, rather than attempting
       to consider a rigorous formulation of a significant difference
       between the PWMs. It looks for changes in the composition of the
       specified bases being beyond a threshold."""
    sel_base_comp_diff = (np.linalg.norm(freq_matrix[:, index_arr], axis=0) -
                          np.linalg.norm(modfreq_matrix[:, index_arr],
                                         axis=0))
    return np.all(sel_base_comp_diff >= _MAX_BASE_COMP_DIFF_FOR_SIMILARITY)


def _createBG(backgroundString):
    """Creates a background ordered dictionary for the
       given, delimited, background input."""
    motifAlphBGFreqs = OrderedDict()
    for bFreq in backgroundString:
        if not bFreq:
            continue
        (base, freq) = bFreq.split(_DELIM)
        motifAlphBGFreqs[base] = float(freq)
    return motifAlphBGFreqs


def output_motif(freq_matrix, output_descriptor, motif_name,
                 motif_alphabet, numSites, EValue, motif_filename,
                 MEME_header):
    status = 0
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
    modFracs = args.modAllFractions

    # NB: width is (totalNumBases - 1), due to temp. additional matrix line
    MEMEBody = textwrap.dedent("""\
               MOTIF {}\n
               letter-probability matrix: alength= {} w= {} nsites= {} E= {}
               """.format(motif_name, freq_matrix.shape[1],
                          (totalNumBases - 1), numSites, EValue))

    if ((args.baseModification and args.baseModPosition)
            or args.tryAllCModsAtPos or
            args.baseModificationAtAllModifiablePos):
        warn("""NB: modified motifs are output to a file.
                The motif printed to STDOUT is the unmodified version.""")

        mod_base_context = _ALL_BASE_CONTEXTS
        baseModPos = args.tryAllCModsAtPos or args.baseModPosition

        if (baseModPos != cUtils._PARAM_A_CONST_VAL
                and isinstance(baseModPos, int)
                and baseModPos >= totalNumBases):
            die("""The provided modification position ({}) exceeds the
                   motif's final position ({}).
                """.format(args.baseModPosition, totalNumBases - 1))

        # index is either positional, comprises all nucleobases, or contextual
        mod_base_index = slice(baseModPos - 1, baseModPos) \
            if (baseModPos and isinstance(baseModPos, int)
                and baseModPos != cUtils._PARAM_A_CONST_VAL) else slice(None)

        if baseModPos and not isinstance(baseModPos, int):
            if baseModPos == 'c':
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
                mod_base_context = baseModPos

        if args.baseModificationAtAllModifiablePos:
            output_descriptor += '-allPos'
        elif baseModPos:
            output_descriptor += '-P' + str(baseModPos)
        if modFracs:
            output_descriptor += '-mFracs'

        if args.annotateMotifNames:
            motif_name += ('-' + args.annotateMotifNames +
                           '-' + output_descriptor)

        # replace the previous motif name with the modified motif name
        MEMEBody = re.sub(r"^(MOTIF ).*", r"\g<1>" + motif_name, MEMEBody)

        bases_to_iter_over = args.baseModification or \
            args.baseModificationAtAllModifiablePos
        if args.tryAllCModsAtPos:
            bases_to_iter_over = \
                [base for base in cUtils.MOD_BASE_NAMES.keys()
                 if base in motifAlphBGFreqs]
        for b in bases_to_iter_over:
            modfreq_matrix = np.copy(freq_matrix)
            if ((b not in cUtils.COMPLEMENTS.keys()) and
                    (b not in cUtils.COMPLEMENTS.values())):
                warn("""Base {} provided in the background is not a
                        currently supported modified nucleobase. It
                        has, accordingly, been skipped.""".format(b))

            def _modifyMatrixPortion(matrix, mod_base_index,
                                     primary_base_to_mod, target_modified_base,
                                     mod_base_context, mod_fractions,
                                     hemimodifyOnly, forward):
                """Modifies the provided portion of the provided matrix,
                   in the indicated fashion, under the provided genomic
                   context."""
                def _resize_context(context):
                    """Places the Boolean context vector in the correct
                       position WRT the bases that the indication of context
                       pertains to."""
                    return np.concatenate((np.zeros(getattr(mod_base_index,
                                                    'start'), dtype=bool),
                                           context,
                                           np.zeros(getattr(mod_base_index,
                                                            'stop')
                                                    - context_len[0],
                                                    dtype=bool)))

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

                # only permit the modifiable bases to be modified
                # the modifiable base is only that which was provided
                # as the primary_base_to_mod
                correct_context = \
                    matrix[mod_base_index, motif_alphabet.
                           index(primary_base_to_mod)] \
                    > _CONTEXT_FREQ_THRESHOLD

                if mod_base_context != _ALL_BASE_CONTEXTS:
                    # create a mask which ensures that modifications are
                    # only permitted in their valid genomic context
                    # NB: the position of the insertion of 'False' depends upon
                    # the portion of the conjunction (i.e. [:-1] vs. [1:])
                    # such that it always inserts on the sequence's termination
                    # and also depends upon operating on the strand
                    correct_context = \
                        np.logical_and(correct_context,
                                       np.insert(matrix[mod_base_index,
                                                        motif_alphabet.
                                                 index(mod_base_context)]
                                                 > _CONTEXT_FREQ_THRESHOLD,
                                                 (-1 if forward else 0),
                                                 [0])[slice(1, None, None)
                                                      if forward
                                                      else slice(None, -1,
                                                                 None)])

                correct_context_minus = np.zeros(context_len, dtype=bool)
                if '-' in hemimodifyOnly:
                    # shift Boolean values right, clipping them
                    correct_context_minus = \
                        np.insert(correct_context[:-1], 0, [0])
                    # also only permit modifiable bases to be mod on - strand
                    correct_context_minus_for_query = correct_context_minus
                    if correct_context_minus.size != matrix.shape[0]:
                        correct_context_minus_for_query = \
                            _resize_context(correct_context_minus)
                    bases_to_modify_minus = \
                        matrix[correct_context_minus_for_query]
                    if bases_to_modify_minus.size > 0:
                        correct_context_minus = \
                            np. \
                            logical_and(correct_context_minus,
                                        bases_to_modify_minus[0]
                                        [motif_alphabet.
                                         index(cUtils.
                                               complement(primary_base_to_mod)
                                               )] > _CONTEXT_FREQ_THRESHOLD)
                if '+' not in hemimodifyOnly:
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
                        matrix_target_view = \
                            matrix[mod_base_index, motif_alphabet.
                                   index(fun(target_modified_base))].view()
                        matrix_source_view = \
                            matrix[mod_base_index, motif_alphabet.
                                   index(fun(primary_base_to_mod))].view()

                        # transfer primary_base_to_mod frequency to the
                        # corresponding target_modified_base entry.

                        matrix_target_before_mod = np.copy(matrix_target_view)

                        # make the modification, setting target values
                        matrix[mod_base_index, motif_alphabet.
                               index(fun(target_modified_base))] = \
                            np.where(context_to_use, matrix_source_view,
                                     matrix_target_view)

                        # set other base frequencies (source values)
                        # to the previous values of the target(s)
                        matrix[mod_base_index, motif_alphabet.
                               index(fun(primary_base_to_mod))] = \
                            np.where(context_to_use, matrix_target_before_mod,
                                     matrix_source_view)
                else:
                    # zero all entries along the frequency matrix,
                    # for the given (row) position, except that corresponding
                    # to the (column) index of the modified base,
                    # which is set to unity.
                    # use the context mask as a sub-array of the correct length
                    # (i.e. preceded and succeeded by 'False' to make up the
                    # length of the motif) to accomplish this.
                    if mod_base_index != slice(None):
                        correct_context = _resize_context(correct_context)
                        correct_context_minus = \
                            _resize_context(correct_context_minus)

                    matrix[correct_context, ] = \
                        only_target_base_at_pos
                    if not isinstance(baseModPos, int):
                        matrix[correct_context_minus, ] = \
                            only_target_base_at_pos_comp

                checkPWMValidity(matrix[:-1], motif_name)

            _modifyMatrixPortion(modfreq_matrix, mod_base_index, 'C', b,
                                 mod_base_context,
                                 True if modFracs else False,
                                 args.hemimodifyOnly, True)
            # the modification for 'G' needs to use the unmodified matrix for
            # context computations, unless in non-fractional mode, in which
            # case it needs to use the already modified matrix to prevent
            # changing the target modification to its complement.
            # integer positional modifications only run the method once
            _modifyMatrixPortion(modfreq_matrix, mod_base_index, 'G',
                                 cUtils.complement(b),
                                 cUtils.complement(mod_base_context),
                                 True if modFracs else False,
                                 args.hemimodifyOnly, False)

            modfreq_matrix = modfreq_matrix[:-1]  # remove extra final row
            checkPWMValidity(modfreq_matrix, motif_name)

            if (not args.onlyNonNegChange or
                # only consider C/G frequency changes to assess the difference
                # (since any modification that wasn't previously present
                # would appear to be a substantial change)
                isMatrixSufficientlyDifferent(freq_matrix[:-1],
                                              modfreq_matrix,
                                              np.array([motif_alphabet.
                                                       index('C'),
                                                       motif_alphabet.
                                                       index('G')])
                                              )):
                with open((os.path.basename(os.path.
                                            splitext(motif_filename)[0]) +
                           '-' + cUtils.MOD_BASE_NAMES[cUtils.
                                                       getMBMaybeFromComp(b)]
                           + output_descriptor + '.meme'), 'a') as outFile:
                    outFile.write(MEME_header)
                    outFile.write(MEMEBody)
                    np.savetxt(outFile, modfreq_matrix, '%f', _DELIM)
            else:
                warn("""A modified matrix ({}) was deemed too close to the
                        original matrix and has accordingly not been output.
                     """.format(motif_name))
                status = -1

    freq_matrix = freq_matrix[:-1]  # remove extra final row
    checkPWMValidity(freq_matrix, motif_name)
    output = StringIO()
    np.savetxt(output, freq_matrix, '%f', _DELIM)
    MEMEBody += output.getvalue()
    output.close()
    return (status, MEMEBody)


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
                                 sequence provided as the argument.")
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
                              Note that this does a basic check that \
                              modifications made make some biological sense \
                              (e.g. a thymine cannot be changed into a \
                              5-methylcytosine using this option). \
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
modBaseSpecifiers.add_argument('-C', '--tryAllCModsAtPos',
                               type=_modifyPositionOptionType,
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
                              '--baseModificationAtAllModifiablePos',
                              type=str, nargs='?',
                              const=str(cUtils._PARAM_A_CONST_VAL),
                              help="Modify the motif to use the modified base \
                              provided for all modifiable motif positions. \
                              A position is considered a modifiable one iff \
                              the nucleobase at said position contains some \
                              cytosine or guanine fraction (i.e. \
                              the resultant PWM has non-zero value for its \
                              'C' or 'G' entries at a given position). \
                              This argument can be used concomitantly with \
                              '-C', in which case no modified nucleobase \
                              need be provided to this argument, since all \
                              possibilities will be attempted (for all \
                              modifiable positions). ")
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
parser.add_argument('-N', '--notNucleobases', type=_omitModBaseOptionType,
                    metavar=_LIST_OF_PRIMARY_MOD_BASES_AND_ABBREVS,
                    help="Do not use the \
                    provided (comma-delimited) list of nucleobases. \
                    These may be specified by either their abbreviation \
                    or base. However, the base specified must be a \
                    \"primary\" and not a complemented nucleobase \
                    (i.e. 'm' and not '1'). For example, \"c,5hmC\", \
                    would use neither 5-carboxylcytosine nor \
                    5-hydroxymethylcytosine.")
parser.add_argument('-F', '--modAllFractions', action='store_true',
                    help="Modify fractions of both cytosine and guanine \
                    (instead of setting the modified base frequency to 1). \
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
parser.add_argument('--annotateMotifNames', nargs='?', type=str,
                    help="Add information regarding modifications performed \
                    to motif names, optionally including the provided string.")
parser.add_argument('--ASCIICodeOrder', help="Sort output in ASCII-code \
                    order, instead of the default order (that specified and \
                    accepted by the MEME-Suite).", action='store_true')
parser.add_argument('-D', '--onlyNonNegChange', help="Only output modified \
                    matrices if they have changed non-negligibly from the \
                    original matrix. The definition of a non-negligible \
                    change may be arbitrarily modified and this option is not \
                    intended to encapsulate a rigorous comparison of the \
                    PWMs (which should be done post-hoc if desired), but \
                    merely to remove highly similar matrices.",
                    action='store_true')
parser.add_argument('--noWarnings', help="Disable warnings.",
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
        (args.baseModificationAtAllModifiablePos
         == cUtils._PARAM_A_CONST_VAL)):
    die("""You must either provide the modified base to '-A' or use
        '-C' to use all possible nucleobases.""")

if (not args.baseModificationAtAllModifiablePos and
        args.tryAllCModsAtPos == cUtils._PARAM_A_CONST_VAL):
    die("""You must either provide the position to modify to '-C' or use
        '-A' to use all possible positions.""")

filename = (args.inPWMFile or args.inPFMFile or args.inTRANSFACFile or
            args.inTRANSFACFreqFile or args.inJASPARFile or args.inMEMEFile)

if not (args.inSeq or filename):
    die("""You must provide some input, by using one of the \"Input File\"
           arguments.\nUse '-h' to view the help message.""")

if (((isinstance(args.baseModPosition, int) and
      args.baseModPosition != cUtils._PARAM_A_CONST_VAL) or
     (isinstance(args.tryAllCModsAtPos, int) and
      args.tryAllCModsAtPos != cUtils._PARAM_A_CONST_VAL))
        and args.hemimodifyOnly != _DEFAULT_HEMIMODARG):
            die("""Hemi-modification cannot be specified with an explicit
                   modification position (motif position {}).
                """.format(args.baseModPosition or args.tryAllCModsAtPos))

motifAlphBGFreqs = OrderedDict()
if args.background:
    if os.path.isfile(args.background):
        with open(args.background) as bgFile:
            motifAlphBGFreqs = _createBG(bgFile)
    else:
        if args.background == _mESC_BG_NAME:
            motifAlphBGFreqs = cUtils.MOUSE_ESC_BACKGROUND
        elif args.background == _AML_BG_NAME:
            motifAlphBGFreqs = cUtils.HUMAN_AML_BACKGROUND

motifs_to_output = ''
numSites = 0
EValue = _PSEUDO_EVALUE_TO_USE
motif_alphabet = list(motifAlphBGFreqs.keys())


if filename:  # PWM or PFM
    # process the matrix for MEME format compatibility
    # If needed, skip columns or rows accordingly and/or
    # use unpack to transpose the matrix.
    ncols = 0
    with open(filename, 'rb') as inFile:
        if args.inMEMEFile:
            input_file = inFile.read()
            MEME_parsed_alphabet = re.search(MEME_MINIMAL_ALPH_REGEX,
                                             input_file,
                                             flags=MEME_MINIMAL_REGEX_FLAGS) \
                .group(MEME_MINIMAL_REGEX_G_ALPH)
            bg_contents = re.search(MEME_MINIMAL_BG_REGEX, input_file,
                                    flags=MEME_MINIMAL_REGEX_FLAGS) \
                .group(MEME_MINIMAL_BG_REGEX_G_BG)
            bgIter = re.findall(MEME_MINIMAL_BG_BASE_REGEX, bg_contents,
                                flags=MEME_MINIMAL_REGEX_FLAGS)
            bg_contents = re.sub(r' ', r'\t', re.sub(r'(\d\d) ',
                                 r"\1\n", bg_contents)).split('\n')
            motifAlphBGFreqs = _createBG(bg_contents)
            motifIter = re.finditer(MEME_MINIMAL_MOTIF_REGEX, input_file,
                                    flags=MEME_MINIMAL_REGEX_FLAGS)
            motif_alphabet = list(bgIter)

            # replace header's alphabet with that just parsed
            MEME_header = re.sub(MEME_MINIMAL_ALPH_REGEX, MEME_parsed_alphabet,
                                 MEME_header)

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
        elif not args.inMEMEFile:
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

        if not args.inMEMEFile:
            freq_matrix = np.hstack((np.zeros((csvData.shape[0],
                                    motif_alphabet.index('A'))), csvData,
                                    np.zeros((csvData.shape[0],
                                             (len(motif_alphabet) - 1) -
                                             motif_alphabet.index('T')))))
            motif_name = os.path.basename(filename)

if args.inSeq:
    if (os.path.isfile(args.inSeq)):
        # NB: min dimensionality of 1 is needed for the character view
        motifs = np.loadtxt(args.inSeq, dtype=str, ndmin=1)
        motifChars = motifs.view('S1').reshape((motifs.size, -1))
    else:
        motifChars = np.expand_dims(np.array(list(args.inSeq)), axis=0)
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
    motif_name = os.path.basename(args.inSeq)

# The zero-order Markov model must sum to unity to be sensical
npt.assert_allclose([sum(motifAlphBGFreqs.itervalues())], [1])

matrix_indicies_to_delete = []

if args.notNucleobases:
    for base in args.notNucleobases:
        if base in cUtils.MOD_BASES:
            base = cUtils.MOD_BASES[base]
        if base not in motifAlphBGFreqs:
            warn("""A provided base to omit ({}) is not present in the
                    current input and has been skipped.""".format(base))
            continue
        motifAlphBGFreqs[cUtils.MOD_MAP[base]] += \
            motifAlphBGFreqs[base]
        motifAlphBGFreqs[cUtils.MOD_MAP[cUtils.complement(base)]] += \
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
        matrix_indicies_to_delete += [motif_alphabet.index(base),
                                      motif_alphabet.
                                      index(cUtils.complement(base))]
    # remove any ambiguity codes that are now empty
    MEME_header = re.sub('(. = \n)', '', MEME_header, flags=re.MULTILINE)

# The zero-order Markov model should still sum to unity to be sensical
npt.assert_allclose([sum(motifAlphBGFreqs.itervalues())], [1])

# used to sort the output in the correct sort order
# need to sort both the background and matrix to do so
sort_fun = ((lambda base: ord(base)) if args.ASCIICodeOrder
            else cUtils.baseSortOrder)
sorted_index = np.argsort(map(sort_fun, list(motifAlphBGFreqs.keys())))

motif_alphabet = [list(motifAlphBGFreqs.keys())[i] for i in sorted_index]

motif_alphabet_bg_freq_output = \
    ' '.join([str(k) + ' ' + str(v) for k, v
             in [motifAlphBGFreqs.items()[i] for i in sorted_index]])
MEME_header += motif_alphabet_bg_freq_output + "\n\n"


def _getMotif(freq_matrix, sorted_index, MEME_header):
    # delete necessary columns
    if matrix_indicies_to_delete:
        freq_matrix = np.delete(freq_matrix, matrix_indicies_to_delete, 1)
    # re-order the matrix by the sorted index
    identity_index = np.indices(freq_matrix.shape)
    freq_matrix = freq_matrix[identity_index[0],
                              np.tile(sorted_index,
                                      (identity_index[1].shape[0], 1))]

    return (output_motif(freq_matrix, output_descriptor, motif_name,
                         motif_alphabet, numSites, EValue,
                         filename or motif_name, MEME_header))

if args.inSeq or not args.inMEMEFile:
    motifs_to_output += _getMotif(freq_matrix, sorted_index, MEME_header)[1] \
        + "\n"
if args.inMEMEFile:
    output_header = True
    for match in motifIter:
        # behave as if read from CSV to minimize changes to other code
        meme_freq_matrix = np.fromstring(match.group(MEME_MIN_REGEX_G_PWM),
                                         dtype=float, sep=' ') \
            .reshape((int(match.group(MEME_MIN_REGEX_G_WIDTH)),
                     int(match.group(MEME_MIN_REGEX_G_ALPH_LEN))))
        motif_name = (match.group(MEME_MIN_REGEX_G_ID).strip() + ' ' +
                      match.group(MEME_MIN_REGEX_G_NAME).strip())
        numSites = int(match.group(MEME_MIN_REGEX_G_NUM_SITES))
        EValue = match.group(MEME_MIN_REGEX_G_E_VALUE)

        (status, MEME_body) = _getMotif(meme_freq_matrix, sorted_index,
                                        (MEME_header if output_header
                                         else "\n"))
        motifs_to_output += MEME_body + "\n"
        # only output a single header per file
        output_header = True if (status < 0 and output_header) else False

print(MEME_header + motifs_to_output.strip())
