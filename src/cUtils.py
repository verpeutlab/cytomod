#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Provide basic utility functions as well as definitions of
   concepts needed to work with modified nucleobases.



   Exports:


   Objects:

   AutoEnum                           - Auto-numbered Enums with docstrings.

   Constants:

   MOD_BASE_COMPLEMENT_NUM_ORDER      - Order of + strand bases, for comp. asn.
   BASE_TO_COVALENT_MODIFICATION_MAP  - Map of mod. bases to mod. names.
   NUCLEOTIDE_POSITIONS_MODIFIED      - Positions of the DNA that are mod.
   POS_STRAND_BASES_MODIFIED          - Positive strand mod. IUPAC bases.
   IUPAC_BASES                        - IUPAC DNA bases, incl. ambiguity codes
   MOD_BASES                          - Map mod. base abbreviations to base.
   MOD_BASE_NAMES                     - Inverse of MOD_BASES.
   AMBIG_MOD_BASES                    - Define mod. base ambiguity codes.
   INVERTED_AMBIG_MOD_BASES           - Inverse of AMBIG_MOD_BASES.
   MOD_MAP                            - Map each mod. base to its unmod. base.
   MODIFIABLE_BASES                   - Bases with defined modifications
   COMPLEMENTS                        - Map of all complements for all bases.
   FULL_BASE_NAMES                    - Full names of unmod. fundamental bases.
   FULL_MOD_BASE_NAMES                - Full names of modmified bases.
   MODIFIED_DINUCL_ORDER              - Order of mod. bases used in H-testing.
   MASK_BASE                          - The base used to mask other bases.
   BASE_COLOURS                       - Map of unequivocal base colours, per
                                        the MEME custom alphabet specification.
   MOUSE_ESC_BACKGROUND               - mESC background (Booth/Ito MS/MS data).
   HUMAN_AML_BACKGROUND               - Human AML bg. (Liu/Kroeze MS/MS data).

   Functions:

   Utility:
   getAlteredSlice           - Return a modified version of an existing Slice.
   makeList                  - Create list from scalar else identity.
   errorMsg                  - Emit error message.
   warn                      - Emit warning message.
   die                       - Emit error message and terminate.
   v_print_timestamp         - Print timestamped message if set
                               verbosity is sufficient.

   Modified base-interacting:
   complement                - Complement of given base(s).
   isModBase                 - Boolean indicating if base is mod.
   isUnivocal                - Boolean indicating if base is unambiguous
   isModifiable              - Boolean indicating if base can be modified
   isModifiableTo            - Boolean indicating if the given unmod. base
                               is modifiable to the given mod. base.
   getUnivocalModBases       - All unambiguous mod. bases.
   getFirstUnivocalBase      - First unambiguous mod. base
                               for given ambiguity code.
   getLastUnivocalModBase    - Last unambiguous mod. base
                               for given ambiguity code.
   getCompMaybeFromMB        - Map primary mod. base to its complement.
   getMBMaybeFromComp        - Map complement mod. base to primary mod base.
   getRGBBaseCol             - The (0-1) RGB colour of the given base.
   getRGB256BaseCol          - The (0-255) RGB colour of the given base.
   getHexBaseCol             - The hex colour of the given base.
   baseSortOrder             - Map bases to an integer rep. their sort order.
"""

from __future__ import with_statement, division, print_function

__version__ = "0.09"

import collections
import datetime
import enum
import functools
import operator
import re
import sys
import textwrap

from collections import OrderedDict
from functools import reduce
from gzip import open as gzip_open
from itertools import chain, izip
from itertools import product as CartesianProd
from os import extsep

from bidict import bidict
import chroma

_EXIT_FAILURE = 1  # exit code used upon failure

EXT_GZ = "gz"
SUFFIX_GZ = extsep + EXT_GZ

_MAX_BASE_NUM = 9
_PARAM_A_CONST_VAL = 999
# Operation used for chroma colour mixing (additive or subtractive)
_COLOUR_MIX_OP = operator.sub


class AutoEnum(enum.Enum):
    """Automatically numbers enum members starting from 1.
       Includes support for a custom docstring per member.

       Adapted from: Ethan Furman (http://stackoverflow.com/a/19330461).
    """
    __last_number__ = 0

    def __new__(cls, *args):
        """Ignores arguments (will be handled in __init__."""
        value = cls.__last_number__ + 1
        cls.__last_number__ = value
        obj = object.__new__(cls)
        obj._value_ = value
        return obj

    def __init__(self, *args):
        """Can handle 0 or 1 argument; more requires a custom __init__.

           0  = auto-number w/o docstring
           1  = auto-number w/ docstring
           2+ = needs custom __init__
        """
        if len(args) == 1 and isinstance(args[0], (str, unicode)):
            self.__doc__ = args[0]
        elif args:
            raise TypeError('{} not dealt with -- '
                            'need custom __init__'.format(args))


# Adapted from: http://stackoverflow.com/a/4029018
def _unpacked(function):
    """Decorator function to return a scalar if a scalar is input
       or a list for any non-basestring Iterable ABC input.
       This is similar to Perl's functionality WRT its context-dependency.
    """
    @functools.wraps(function)
    def _decorator(list_or_scalar, *args, **kwargs):
        result = function(list_or_scalar, *args, **kwargs)

        # NB: strings are scalars, despite having the Sequence ABC
        if (isinstance(list_or_scalar, basestring) or
                not isinstance(list_or_scalar, collections.Iterable)):
            # unpack to a scalar, if the input was a scalar
            result = result[0]

        return result

    return _decorator


def _consInvConcatBases(bases):
    """Return the swapped key, value pairs for all pairs in
       the given dictionary with lists as its values,
       concatenating the value list.
    """
    return dict((''.join(baseSubject), baseQuery)
                for baseQuery, baseSubject in bases.iteritems())


def _consAmbigModBases(mod_map, abases):
    """Return a dictionary of ambiguous modified bases."""
    # We require that the first entry in the list be a primary
    # modified base for all ambiguity codes in bases.
    if __debug__:
        for base in abases.itervalues():
            assert base[0] in mod_map, textwrap.fill(textwrap.dedent("""\
            The first value of the list for all ambiguity codes
            must be a primary modified base.
            %r is not such a base.""" % base[0]))
    # Return the dictionary for all ambiguity codes
    return {baseQuery: mod_map[baseSubject[0]]
            for baseQuery, baseSubject in abases.iteritems()}


def _consAllComps(comps, mbases, abases):
    """Return a dictionary with all modified nucleobase complements."""
    # Complements start at 1 and increment for modified bases
    modifiedBasesToComplements = \
        zip(mbases.values(),
            ''.join(str(i) for i in range(1, len(mbases) + 1)))
    # The ambiguity code complement of an ambiguous modified
    # base starts at _MAX_BASE_NUM and decrements
    modifiedBaseAmbigCodesToComplements = \
        zip(abases.keys(),
            ''.join(str(i) for i in range(_MAX_BASE_NUM,
                                          _MAX_BASE_NUM - len(abases), -1)))
    # Return all modified nucleobase and ambiguity code complements
    # They are ordered by the originating modification's oxidation order
    comps.update(modifiedBasesToComplements +
                 modifiedBaseAmbigCodesToComplements)
    return comps


def _consAmbigBaseNames(fullMbaseNames, abases):
    """Return a dictionary of all full names for ambiguity codes."""
    # Add specific names and for the remaining ambiguity codes,
    # just use their constitutive base names, concatenated with 'or'.
    tdict = dict(fullMbaseNames, **{'z': 'Possibly_modified_cytosine',
                 str(_MAX_BASE_NUM): 'Possibly_modified_guanine'})  # concat.
    tdict.update({baseQuery: '_or_'.join(filter(None, [tdict.get(ambigBase)
                                                for ambigBase in
                                                abases[baseQuery]]))
                 for baseQuery, baseSubject in abases.iteritems()
                 if baseQuery not in tdict})
    return tdict


def _consCompModBasePairs(mod_map):
    """Return zipped list of complemented modified base pairs."""
    return izip(complement(mod_map.keys()), complement(mod_map.values()))


def _consCompModBaseFullNames(fullbaseNames, fullMbaseNames, mod_map,
                              mbases, abases):
    """Return a dictionary of names for all modified base complements,
       using existing nomenclature.
    """
    for base in complement(chain(mbases.values(), abases.keys())):
        if base not in fullMbaseNames:  # add iff a name is not already present
            fullMbaseNames.update(izip(base, [fullbaseNames[mod_map[base]] +
                                              ':' +
                                              fullMbaseNames[complement(base)]]
                                       )
                                  )
    return fullMbaseNames


def _consModBasesModMap(nucPositionsMod, posStrandBasesMod,
                        baseToCovalentModMap, modBasesComplementOrder):
    """Return the "primary" modified bases and their corresponding
       one base codes, listed in their order of oxidation. Generalize for the
       base being modified and its position. Also return a dictionary mapping
       each "primary" modified base to the base it modifies.
    """
    MOD_BASES = {}
    MOD_MAP = {}
    for pos_modified in nucPositionsMod:
        for base_mod in posStrandBasesMod:
            for covalent_mod in baseToCovalentModMap.values():
                covalently_mod_base = (pos_modified + covalent_mod + base_mod)
                MOD_BASES[covalently_mod_base] = covalent_mod[:1]
                MOD_MAP.update(dict.fromkeys(MOD_BASES.values(), base_mod))
    # order the modified bases to facilitate assignment of complement numerals
    MOD_BASES = OrderedDict(sorted(MOD_BASES.items(), key=lambda t:
                                   modBasesComplementOrder.index(t[1])))
    return (MOD_BASES, MOD_MAP)


def _consCovalentMods(nucPositionsMod, posStrandBasesMod, baseToCovalentModMap,
                      modBases, ambigModBases):
    """Returns all names for covalently modified bases and ambiguously modified
       bases. This is the Cartesian product of the positions modified, the
       covalent modificationsm and the base being modified.
    """
    return [''.join(tuple) for tuple in CartesianProd(nucPositionsMod,
            [baseToCovalentModMap[mod_b] for mod_b in modBases.values()]
            + ambigModBases.keys(), posStrandBasesMod)]


def _consAllPosStrandModsToBase(covalentModToBaseMap, covalentModBases):
    """Returns a map of each covalent modification name (including ambiguous
       modifications) to its respective expanded alphabet base.
    """
    return {cov_mod_bases:
            (covalentModToBaseMap.get(cov_mod_bases[1:len(cov_mod_bases)-1]) or
             cov_mod_bases[1:len(cov_mod_bases) - 1])
            for cov_mod_bases in covalentModBases}


# defines the order for the assignment of complmentary numerals
MOD_BASE_COMPLEMENT_NUM_ORDER = ['m', 'h', 'f', 'c']

BASE_TO_COVALENT_MODIFICATION_MAP = bidict({'m': 'm', 'h': 'hm', 'f': 'f',
                                            'c': 'ca'})

NUCLEOTIDE_POSITIONS_MODIFIED = ['5']
POS_STRAND_BASES_MODIFIED = ['C']

(MOD_BASES, MOD_MAP) = _consModBasesModMap(NUCLEOTIDE_POSITIONS_MODIFIED,
                                           POS_STRAND_BASES_MODIFIED,
                                           BASE_TO_COVALENT_MODIFICATION_MAP,
                                           MOD_BASE_COMPLEMENT_NUM_ORDER)

# Define the IUPAC DNA bases, including ambiguity codes
# The list includes the additional MEME 'X' = 'N' as well as the
# four core nulceobases.
IUPAC_BASES = {'A': ['A'], 'T': ['T'], 'C': ['C'], 'G': ['G'],
               'R': ['A', 'G'], 'Y': ['C', 'T'],
               'K': ['G', 'T'], 'M': ['A', 'C'],
               'S': ['C', 'G'], 'W': ['A', 'T'],
               'B': ['C', 'G', 'T'], 'D': ['G', 'A', 'T'],
               'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
               'N': ['A', 'C', 'G', 'T'], 'X': ['A', 'C', 'G', 'T']}

# Define the modified base ambiguity codes, listed in
# their order of decreasing generality and then by the
# order of their value list in MOD_BASES.
# By convention, we place unmodified nucleobases last in
# the list of possible interpretations of an ambiguity code.
AMBIG_MOD_BASES = OrderedDict([('z', ['m', 'h', 'f',
                                      'c', 'C']),
                               ('y', ['f', 'c', 'C']),
                               ('x', ['m', 'h']),
                               ('w', ['f', 'c'])])

# all modified and ambiguously modified bases
ALL_NON_UNMOD_BASES = MOD_BASES.values() + AMBIG_MOD_BASES.keys()
COVALENT_MOD_BASES = _consCovalentMods(NUCLEOTIDE_POSITIONS_MODIFIED,
                                       POS_STRAND_BASES_MODIFIED,
                                       BASE_TO_COVALENT_MODIFICATION_MAP,
                                       MOD_BASES, AMBIG_MOD_BASES)

COVALENT_MOD_BASES_TO_BASE_MAP = \
    _consAllPosStrandModsToBase(BASE_TO_COVALENT_MODIFICATION_MAP.inv,
                                COVALENT_MOD_BASES)

# Permits ambiguity code lookup using concatenated modified bases
# NB: Assumes that values are unique (they should always be)
INVERTED_AMBIG_MOD_BASES = _consInvConcatBases(AMBIG_MOD_BASES)

MOD_MAP.update(_consAmbigModBases(MOD_MAP, AMBIG_MOD_BASES))

# All IUPAC nucleobases and their complements, plus 'X',
# which is just an additional alias for any nucleobase
# NB: Unmod. keys are expected to be those most likely to be covalently mod.,
# while each value should be the complement that is less commonly modified.
COMPLEMENTS = {'T': 'A', 'C': 'G',
               'Y': 'R', 'M': 'K',
               'W': 'W', 'S': 'S',
               'B': 'V', 'H': 'D',
               'N': 'N', 'X': 'X'}
COMPLEMENTS = bidict(_consAllComps(COMPLEMENTS, MOD_BASES,
                     AMBIG_MOD_BASES))
FULL_BASE_NAMES = {'A': 'Adenine', 'T': 'Thymine',
                   'G': 'Guanine', 'C': 'Cytosine'}

# We do not currently have any abbreviations of ambiguity codes
# We do not maintain any ordering for modified base abbreviations
MOD_BASE_NAMES = {base: abbrev for abbrev, base in dict(MOD_BASES).iteritems()}

FULL_MOD_BASE_NAMES = {'m': '5-Methylcytosine',
                       'h': '5-Hydroxymethylcytosine',
                       'f': '5-Formylcytosine',
                       'c': '5-Carboxylcytosine'}

BASE_COLOURS = {'A': '8510A8', 'T': 'A89610', 'C': 'A50026', 'G': '313695',
                'm': 'D73027', '1': '4575B4', 'h': 'F46D43', '2': '74ADD1',
                'f': 'FDAE61', '3': 'ABD9E9', 'c': 'FEE090', '4': 'E0F3F8'}
# We require that all primary bases have a defined colour
if __debug__:
    for base in FULL_BASE_NAMES.keys() + FULL_MOD_BASE_NAMES.keys():
        assert base in BASE_COLOURS, textwrap.fill(textwrap.dedent("""\
            No colour is defined for nucleobase %r.
            All primary bases must have a defined colour.
            """ % base))

FULL_MOD_BASE_NAMES = _consAmbigBaseNames(FULL_MOD_BASE_NAMES,
                                          AMBIG_MOD_BASES)

MASK_BASE = AMBIG_MOD_BASES.keys()[0]

# --------------------------------------------------------------------------
# define background based upon existing MS/MS data

# The alphabet frequencies used are (WRT Cs only):
# 2.95% 5mC, 0.055 ± 0.008% 5hmC, 0.0014 ± 0.0003% 5fC,
# and 0.000335% 5caC
# (Respectively: Ito et al. 2011, Booth et al. 2014,
# Booth et al. 2014, and Ito et al. 2011)
MOUSE_ESC_BACKGROUND = \
    OrderedDict({'T': 0.292, 'A': 0.292, 'C': 0.201745991, 'G': 0.201745991,
                 # (using mouse GC content of 41.6%)
                 # From: NCBI Eukaryotic Genome Report File
                 'm': 0.006136, '1': 0.006136, 'h': 0.0001144, '2': 0.0001144,
                 'f': 0.000002912, '3': 0.000002912,
                 'c': 0.000000697, '4': 0.000000697})
# --
# The alphabet frequencies used were (WRT all bases):
# 2.91 ± 0.11% 5mC and 0.039% 5hmC.
# (Respectively: Liu et al. 2007 and Kroeze et al. 2014)
# 5fC at 0.0021812% was estimated from the 5fC to 5hmC ratio
# within the melanoma cell line WM-266-4,
# which was analyzed by Liu S. et al.
HUMAN_AML_BACKGROUND = \
    OrderedDict({'T': 0.295, 'A': 0.295, 'C': 0.190244094, 'G': 0.190244094,
                 # (using human GC content of 41.0%)
                 # From: NCBI Eukaryotic Genome Report File
                 'm': 0.01455, '1': 0.01455, 'h': 0.000195, '2': 0.000195,
                 'f': 0.000010906, '3': 0.000010906})
# --------------------------------------------------------------------------


@_unpacked
def complement(bases):
    """Complement the given, potentially modified, base.
       Returns the input base if no complement for that base is found.
       This function returns a scalar if one is input, otherwise a list.
    """
    return [COMPLEMENTS.get(base) or COMPLEMENTS.inv.get(base) or base
            for base in bases]


# Update the dictionary mapping with every complemented modification
MOD_MAP.update(_consCompModBasePairs(MOD_MAP))

# construct a list of all modifiable bases, from unique values of MOD_MAP
MODIFIABLE_BASES = set(MOD_MAP.values())

FULL_MOD_BASE_NAMES = _consCompModBaseFullNames(FULL_BASE_NAMES,
                                                FULL_MOD_BASE_NAMES, MOD_MAP,
                                                MOD_BASES, AMBIG_MOD_BASES)


def _cons_mod_dinucl_order(mod_base_order, mod_to_unmod_map):
    """Returns a list with modified dinucleotides, in their natural ordering.
       This order is by the underlying order of the primary modified
       nucleobases, followed by their complements, ordered by:
       only positive strand hemi-modifications,
       only negative strand hemi-modifications, and complete modifications."""
    result = []
    for mod_base in mod_base_order:
        result += [mod_base + complement(mod_to_unmod_map[mod_base])]
        result += [mod_to_unmod_map[mod_base] + complement(mod_base)]
        result += [mod_base + complement(mod_base)]
    return result


MODIFIED_DINUCL_ORDER = _cons_mod_dinucl_order(MOD_BASE_COMPLEMENT_NUM_ORDER,
                                               MOD_MAP)


def isModBase(base):
    """Return a Boolean indicating if the base is modified (of any kind)."""
    return (True if base in MOD_MAP.keys() else False)


def getUnivocalModBases():
    """Return all modified bases, excluding ambiguity codes."""
    return MOD_BASE_NAMES.keys() + complement(MOD_BASE_NAMES.keys())


def isUnivocal(base):
    """Return a Boolean indicating if the base is univocal."""
    return (True if base in getUnivocalModBases() else False)


def isModifiable(unmodBase):
    """Return a Boolean indicating if the unmodified base has
       modifications defined for it.
    """
    return (True if unmodBase in MODIFIABLE_BASES else False)


def isModifiableTo(unmodBase, modBase):
    """Return a Boolean indicating if the unmodified base has
       the given modified base defined as a possible covalent
       modification.
    """
    return (True if MOD_MAP[modBase] == unmodBase else False)


def getFirstUnivocalBase(base):
    """Return the first univocal base for the given base.
    Return the first value in the definition for an ambiguous
    modified base or the base itself if it is already unambiguous.
    """
    if AMBIG_MOD_BASES.get(base):
        return AMBIG_MOD_BASES.get(base)[0]
    elif AMBIG_MOD_BASES.get(complement(base)):
        return complement(AMBIG_MOD_BASES.get(complement(base))[0])
    else:
        return base


def getLastUnivocalModBase(base):
    """Return the last univocal modified base for the given base.
    Return the last value in the definition for an ambiguous
    modified base, provided that it is a modified nucleobase,
    or the base itself if it is already unambiguous.
    """
    if AMBIG_MOD_BASES.get(base):
        for baseQuery in reversed(AMBIG_MOD_BASES.get(base)):
            if not FULL_BASE_NAMES.get(baseQuery):
                return baseQuery
    elif AMBIG_MOD_BASES.get(complement(base)):
        for baseQuery in reversed(AMBIG_MOD_BASES.get(complement(base))):
            if not FULL_BASE_NAMES.get(baseQuery):
                return complement(baseQuery)
    else:
        return base


def getCompMaybeFromMB(base):
    """Map the given modified based to the corresponding
    modified guanine nucleobase (i.e. the complemented modified base).
    Apply the identity transformation if the given base is
    already a complemented modified nucleobase.
    If an unmodified base is provided, this attempts to map it to a complement
    that is less likely to be covalently modified (e.g. C, T) or if already a
    base that is commonly modified, then to itself.
    """
    # use forward mapping
    return COMPLEMENTS.get(base) or base


def getMBMaybeFromComp(base):
    """Map the given modified based to the corresponding
    modified cytosine nucleobase (i.e. the actual modified base).
    Apply the identity transformation if the given base is
    already a modified cytosine nucleobase.
    If an unmodified base is provided, this attempts to map it to a complement
    that is more likely to be covalently modified (e.g. C, T) or if already a
    base that is commonly modified, then to itself.
    """
    # use reverse mapping (i.e. invert the bijection)
    return COMPLEMENTS.inv.get(base) or base


def _getBaseCol(base):
    """Return a chroma colour of the given modified or unmodified
    nucleobase. If the given base is ambiguous, derive the colour
    by mixing the constituent bases colours. For unmodified bases,
    only the canonical ATGCN colours are defined. Colour mixing
    is not implemented for the other IUPAC base definitions,
    i.e. ambiguous base colours can only be retrieved for modified
    nucleobases. Unrecognized bases will be coloured as N.
    """
    if base in BASE_COLOURS:
        colour = chroma.Color('#' + BASE_COLOURS[base])
    elif getMBMaybeFromComp(base) in AMBIG_MOD_BASES:
        bases = AMBIG_MOD_BASES.get(getMBMaybeFromComp(base))
        if base in COMPLEMENTS.inv:
            bases = complement(bases)
        # mix of constituent ambiguity code nucleobases
        colour = reduce(_COLOUR_MIX_OP,
                        [chroma.Color('#' + BASE_COLOURS[primary_base])
                         for primary_base in bases])
    else:
        # mix of all primary unmodified nucleobases
        colour = reduce(_COLOUR_MIX_OP,
                        [chroma.Color('#' + BASE_COLOURS[primary_base])
                         for primary_base in FULL_BASE_NAMES.keys()])
    return colour


def baseSortOrder(base):
    """Return an integer which represents the mapping of the
       nucleobase to its sort order.
       This order is specified by MEME. The order to corresponds to
       the unmodified bases, followed by primary modified bases,
       followed by their complements."""
    return (ord(base) - ord('A')*base.isupper() + ord('Z')*base.isdigit())


def getRGBBaseCol(base):
    """Return the given nucleobase's RGB (between 0 and 1) colour,
    via _getBaseCol.
    """
    return _getBaseCol(base).rgb


def getRGB256BaseCol(base):
    """Return the given nucleobase's RGB (between 0 and 255) colour,
    via _getBaseCol.
    """
    return _getBaseCol(base).rgb256


def getHexBaseCol(base):
    """Return the given nucleobase's hex (#RRGGBB) colour, via _getBaseCol.
    """
    return _getBaseCol(base).hex


def makeList(lstOrVal):
    """Return a list of a single item if the object passed is not
    already a list. This allows one to iterate over objects which
    may or may not already be lists (and therefore iterable).
    """
    return [lstOrVal] if not isinstance(lstOrVal, list) else lstOrVal


def getAlteredSlice(slice_to_alter, slice_max, operation, value):
    """Alters the provided slice by applying the provided operation, using
       the provided value.
       The returned slice will have stop value at most equal to slice_max."""
    if slice_to_alter == slice(None):
        return slice_to_alter
    else:
        slice_stop_plus_one = (operation(getattr(slice_to_alter, 'stop'),
                                         value)
                               if (operation(getattr(slice_to_alter, 'stop'),
                                   value)) <= slice_max else
                               slice_max)

        if getattr(slice_to_alter, 'start') is None:
            return slice(slice_stop_plus_one)
        else:
            return slice(operation(getattr(slice_to_alter, 'start'), value),
                         slice_stop_plus_one,
                         getattr(slice_to_alter, 'step'))


def maybe_gzip_open(filename, *args, **kwargs):
    """Open a gzipped file with the gzip open file handler and open
       a non-gzipped file with the default open file handler.
       Function from Genomedata by Michael Hoffman.
    """
    if filename.endswith(SUFFIX_GZ):
        return gzip_open(filename, *args, **kwargs)
    else:
        return open(filename, *args, **kwargs)


def errorMsg(msg, msgType, additionalPrefix):
    """Emit an error message to STDERR."""
    prefix = '>> <' + additionalPrefix + '> ' + msgType + ' '
    print(textwrap.dedent(textwrap.fill(msg, initial_indent=prefix,
          subsequent_indent=re.sub('.', ' ', prefix))), file=sys.stderr)


def warn(msg, additionalPrefix=""):
    """Emit a warning message to STDERR."""
    errorMsg(msg, 'Warning:', additionalPrefix)


def die(msg, additionalPrefix="", exit_code=_EXIT_FAILURE):
    """Emit a fatal error message to STDERR."""
    errorMsg(msg, 'Fatal:', additionalPrefix)
    exit(exit_code)


def v_print_timestamp(verbosity, msg="", threshold=1, additionalPrefix=""):
    """Print a timestamped message iff verbosity is at least threshold."""
    if verbosity >= threshold:
        errorMsg(msg, datetime.datetime.now().isoformat(), additionalPrefix)


def assert_or_die_with_msg(condition, msg=""):
    """Standard Python assertion statement, but augmented with
       invocation of the above die function if it fails."""
    try:
        assert condition
    except AssertionError as ex:
        die(msg, "assertion failed: {}".format(ex.args))
