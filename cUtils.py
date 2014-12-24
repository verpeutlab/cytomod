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
"""

from __future__ import with_statement, division, print_function

__version__ = "0.06"


import datetime
import functools
import re
import sys
import textwrap

from collections import OrderedDict
from itertools import izip, chain

from bidict import bidict


_MAX_BASE_NUM = 9
_PARAM_A_CONST_VAL = 999


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


# Define the "primary" modified bases and their corresponding
# one base codes, listed in their order of oxidation
MOD_BASES = OrderedDict([('5mC', 'm'), ('5hmC', 'h'),
                        ('5fC', 'f'), ('5caC', 'c')])

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

# Permits ambiguity code lookup using concatenated modified bases
# NB: Assumes that values are unique (they should always be)
INVERTED_AMBIG_MOD_BASES = _consInvConcatBases(AMBIG_MOD_BASES)

# Create a dictionary mapping each "primary" modified base to
# the base it modifies
MOD_MAP = dict.fromkeys(MOD_BASES.values(), 'C')
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
FULL_MOD_BASE_NAMES = _consAmbigBaseNames(FULL_MOD_BASE_NAMES,
                                          AMBIG_MOD_BASES)


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


def makeList(lstOrVal):
    """Return a list of a single item if the object passed is not
    already a list. This allows one to iterate over objects which
    may or may not already be lists (and therefore iterable).
    """
    return [lstOrVal] if not isinstance(lstOrVal, list) else lstOrVal


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
