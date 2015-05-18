"""Provide basic utility functions as well as definitions of
   concepts needed to work with modified nucleobases.



   Exports:


   Constants:

   MOD_BASES                 - Map mod. base abbreviations to base.
   MOD_BASE_NAMES            - Inverse of MOD_BASES.
   AMBIG_MOD_BASES           - Define mod. base ambiguity codes.
   INVERTED_AMBIG_MOD_BASES  - Inverse of AMBIG_MOD_BASES.
   MOD_MAP                   - Map each mod. base to its unmodified base.
   COMPLEMENTS               - Map of all complements for all bases.
   FULL_BASE_NAMES           - Full names of unmodmified fundamental bases.
   FULL_MOD_BASE_NAMES       - Full names of modmified bases.
   MASK_BASE                 - The base to be used to mask other bases.
   BASE_COLOURS              - Map of unequivocal base colours,
                               per the MEME custom alphabet specification.

   Functions:

   Utility:
   makeList                  - Create list from scalar else identity.
   errorMsg                  - Emit error message.
   warn                      - Emit warning message.
   die                       - Emit error message and terminate.
   v_print_timestamp         - Print timestamped message if set
                               verbosity is sufficient.

   Modified base-interacting:
   complement                - Return complement of given base(s).
   getUnivocalModBases       - Return all unambiguous mod. bases.
   getFirstUnivocalBase      - Return first unambiguous mod. base
                               for given ambiguity code.
   getLastUnivocalModBase    - Return last unambiguous mod. base
                               for given ambiguity code.
   getCompMaybeFromMB        - Map primary mod. base to its complement.
   getMBMaybeFromComp        - Map complement mod. base to primary mod base.
   getRGBBaseCol             - Return the (0-1) RGB colour of the given base.
   getRGB256BaseCol          - Return the (0-255) RGB colour of the given base.
"""

from __future__ import with_statement, division, print_function

__version__ = "0.09"

import datetime
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

EXT_GZ = "gz"
SUFFIX_GZ = extsep + EXT_GZ

_MAX_BASE_NUM = 9
_PARAM_A_CONST_VAL = 999
# Operation used for chroma colour mixing (additive or subtractive)
_COLOUR_MIX_OP = operator.sub


# FROM: http://stackoverflow.com/a/4029018
def _unpacked(method):
    """Decorator function to return scalar if a scalar is input
       or a list if a list input. This is similar to Perl's
       functionality WRT its context-dependency.
    """
    @functools.wraps(method)
    def _decorator(*args):
        results = method(*args)
        return results if len(results) != 1 else results[0]
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


BASE_TO_COVALENT_MODIFICATION_MAP = bidict({'m': 'm', 'h': 'hm', 'f': 'f',
                                            'c': 'ca'})

NUCLEOTIDE_POSITIONS_MODIFIED = ['5']
POS_STRAND_BASES_MODIFIED = ['C']

# Define the "primary" modified bases and their corresponding
# one base codes, listed in their order of oxidation
# Generalize for the base being modified and its position
# Also create a dictionary mapping each "primary" modified base to
# the base it modifies

# XXX NEEDS MASSIVE REFACTORING and put into functions....
# XXX RUN Flake8 and FIX ALL

# defines the order for the assignment of complmentary numerals
MOD_BASE_COMPLEMENT_NUM_ORDER = ['m', 'h', 'f', 'c']

MOD_BASES = {}
MOD_MAP = {}
for pos_modified in NUCLEOTIDE_POSITIONS_MODIFIED:
    for base_modified in POS_STRAND_BASES_MODIFIED:
        for covalent_mod in BASE_TO_COVALENT_MODIFICATION_MAP.values():
            covalently_mod_base = pos_modified + covalent_mod + base_modified
            MOD_BASES[pos_modified + covalent_mod + base_modified] = \
                covalent_mod[:1]
            MOD_MAP.update(dict.fromkeys(MOD_BASES.values(), base_modified))
# order the modified bases to facilitate assignment of complement numerals
MOD_BASES = OrderedDict(sorted(MOD_BASES.items(), key=lambda t:
                               MOD_BASE_COMPLEMENT_NUM_ORDER.index(t[1])))

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

# XXX NEEDS MASSIVE REFACTORING and put into functions....

# all modified and ambiguously modified bases
ALL_NON_UNMOD_BASES = MOD_BASES.values() + AMBIG_MOD_BASES.keys()
COVALENT_MOD_BASES = [''.join(tuple) for tuple in
                      CartesianProd(NUCLEOTIDE_POSITIONS_MODIFIED,
                                    [BASE_TO_COVALENT_MODIFICATION_MAP[mod_b]
                                     for mod_b in MOD_BASES.values()] +
                      AMBIG_MOD_BASES.keys(), POS_STRAND_BASES_MODIFIED)]
ALL_POS_STRAND_COVALENT_DNA_MODS = \
    {covalent_mod_bases:
     ((~BASE_TO_COVALENT_MODIFICATION_MAP)
      .get(covalent_mod_bases[1:len(covalent_mod_bases)-1]) or
      covalent_mod_bases[1:len(covalent_mod_bases) - 1]) for covalent_mod_bases
     in COVALENT_MOD_BASES}

# Permits ambiguity code lookup using concatenated modified bases
# NB: Assumes that values are unique (they should always be)
INVERTED_AMBIG_MOD_BASES = _consInvConcatBases(AMBIG_MOD_BASES)

MOD_MAP.update(_consAmbigModBases(MOD_MAP, AMBIG_MOD_BASES))

# All IUPAC nucleobases and their complements, plus 'X',
# which is just an additional alias for any nucleobase
COMPLEMENTS = {'A': 'T', 'G': 'C',
               'R': 'Y', 'M': 'K',
               'W': 'W', 'S': 'S',
               'B': 'V', 'D': 'H',
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


@_unpacked
def complement(bases):
    """Complement the given, potentially modified, base."""
    return [COMPLEMENTS.get(base) or (~COMPLEMENTS).get(base)
            for base in bases]


# Update the dictionary mapping with every complemented modification
MOD_MAP.update(_consCompModBasePairs(MOD_MAP))

FULL_MOD_BASE_NAMES = _consCompModBaseFullNames(FULL_BASE_NAMES,
                                                FULL_MOD_BASE_NAMES, MOD_MAP,
                                                MOD_BASES, AMBIG_MOD_BASES)


def getUnivocalModBases():
    """Return all modified bases, excluding ambiguity codes."""
    return MOD_BASE_NAMES.keys() + complement(MOD_BASE_NAMES.keys())


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


def getCompMaybeFromMB(modBase):
    """Map the given modified based to the corresponding
    modified guanine nucleobase (i.e. the complemented modified base).
    Apply the identity transformation if the given base is
    already a complemented modified nucleobase.
    """
    # use forward mapping
    return COMPLEMENTS.get(modBase) or modBase


def getMBMaybeFromComp(modBase):
    """Map the given modified based to the corresponding
    modified cytosine nucleobase (i.e. the actual modified base).
    Apply the identity transformation if the given base is
    already a modified cytosine nucleobase.
    """
    # use reverse mapping (i.e. invert the bijection)
    return (~COMPLEMENTS).get(modBase) or modBase


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
        if base in ~COMPLEMENTS:
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


def makeList(lstOrVal):
    """Return a list of a single item if the object passed is not
    already a list. This allows one to iterate over objects which
    may or may not already be lists (and therefore iterable).
    """
    return [lstOrVal] if not isinstance(lstOrVal, list) else lstOrVal


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


def die(msg, additionalPrefix=""):
    """Emit a fatal error message to STDERR."""
    errorMsg(msg, 'Fatal:', additionalPrefix)
    exit(1)


def v_print_timestamp(verbosity, msg="", threshold=1, additionalPrefix=""):
    """Print a timestamped message iff verbosity is at least threshold."""
    if verbosity >= threshold:
        errorMsg(msg, datetime.datetime.now().isoformat(), additionalPrefix)
