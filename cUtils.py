from __future__ import with_statement, division, print_function

VERSION = "$Revision: 0.06$"

import sys
import re
import datetime
import textwrap
import functools
from bidict import bidict
from collections import OrderedDict
from itertools import izip, chain


# FROM: http://stackoverflow.com/a/4029018
def unpacked(method):
    @functools.wraps(method)
    def _decorator(*args):
        results = method(*args)
        return results if len(results) != 1 else results[0]
    return _decorator

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
                               ('x', ['m', 'h'])])
# Permits ambiguity code lookup using concatenated modified bases
# NB: Assumes that values are unique (they should always be)
INVERTED_AMBIG_MOD_BASES = dict((''.join(b), a) for a, b
                                in AMBIG_MOD_BASES.iteritems())
# Create a dictionary mapping each "primary" modified base to
# the base it modifies
MODIFIES = dict.fromkeys(MOD_BASES.values(), 'C')
# We require that the first entry in the list be a primary
# modified base for all ambiguity codes in AMBIG_MOD_BASES.
for b in AMBIG_MOD_BASES.itervalues():
    assert b[0] in MODIFIES, textwrap.fill(textwrap.dedent("""\
    The first value of the list for all ambiguity codes
    must be a primary modified base.
    %r is not such a base.""" % b[0]))
# Add to the dictionary for all ambiguity codes
MODIFIES.update({a: MODIFIES[b[0]] for a, b in
                 AMBIG_MOD_BASES.iteritems()})

# All IUPAC nucleobases and their complements, plus 'X',
# which is just an additional alias for any nucleobase
COMPLEMENTS = {'A': 'T', 'G': 'C',
               'R': 'Y', 'M': 'K',
               'W': 'W', 'S': 'S',
               'B': 'V', 'D': 'H',
               'N': 'N', 'X': 'X'}
# Add all modified nucleobases
# Complements start at 1 and increment for modified bases
modifiedBasesToComplements = \
    izip(MOD_BASES.values(), ''.
         join(str(i) for i in range(1, len(MOD_BASES) + 1)))
# The ambiguity code complement of an ambiguous modified
# base starts at 9 and decrements
modifiedBaseAmbigCodesToComplements = \
    zip(AMBIG_MOD_BASES.keys(), ''.
        join(str(i) for i in range(9, 9 - len(AMBIG_MOD_BASES), -1)))
# Add all modified nucleobase and ambiguity code complements
# They are ordered by the originating modification's oxidation order
COMPLEMENTS.update(modifiedBasesToComplements)
COMPLEMENTS.update(modifiedBaseAmbigCodesToComplements)
COMPLEMENTS = bidict(COMPLEMENTS)

FULL_BASE_NAMES = {'A': 'Adenine', 'T': 'Thymine',
                   'G': 'Guanine', 'C': 'Cytosine'}
# We do not currently have any short names for ambiguity codes
MOD_BASE_NAMES = {'m': '5mC', 'h': '5hmC', 'f': '5fC', 'c': '5caC'}
FULL_MOD_BASE_NAMES = {'m': '5-Methylcytosine',
                       'h': '5-Hydroxymethylcytosine',
                       'f': '5-Formylcytosine',
                       'c': '5-Carboxylcytosine'}
# Add names for ambiguity codes:
# Add specific names
FULL_MOD_BASE_NAMES.update({'z': 'Possibly_modified_cytosine'})
FULL_MOD_BASE_NAMES.update({'9': 'Possibly_modified_guanine'})
# For the remaining ambiguity codes, just use their
# constitutive base names, concatenated with 'or'.
FULL_MOD_BASE_NAMES.update({a: '_or_'.join([FULL_MOD_BASE_NAMES.get(b)
                                            or FULL_BASE_NAMES.get(b) for
                                            b in AMBIG_MOD_BASES[a]])
                            for a, b in
                            AMBIG_MOD_BASES.iteritems() if
                            a not in FULL_MOD_BASE_NAMES})


@unpacked
def complement(bases):
    """Complements the given, potentially modified, base."""
    return [COMPLEMENTS.get(b) or (~COMPLEMENTS).get(b) for b in bases]
# Update the dictionary mapping with every complemented modification
MODIFIES.update(izip(complement(MODIFIES.keys()),
                complement(MODIFIES.values())))

# Add the names of all modified base complements, using existing nomenclature
for b in complement(chain(MOD_BASES.values(), AMBIG_MOD_BASES.keys())):
    if b not in FULL_MOD_BASE_NAMES:  # add iff a name is not already present
        FULL_MOD_BASE_NAMES.update(izip(b, [FULL_BASE_NAMES[MODIFIES[b]] +
                                   ':' + FULL_MOD_BASE_NAMES[complement(b)]]))

_PARAM_A_CONST_VAL = 999


def getUnivocalModBases():
    """Returns all modified bases, excluding ambiguity codes."""
    return MOD_BASES.values() + complement(MOD_BASES.values())


def getFirstUnivocalBase(a):
    """Returns the first univocal base for the given base.
    This is the first value in the definition for an ambiguous
    modified base or the base itself if it is already unambiguous."""
    if AMBIG_MOD_BASES.get(a):
        return AMBIG_MOD_BASES.get(a)[0]
    elif AMBIG_MOD_BASES.get(complement(a)):
        return complement(AMBIG_MOD_BASES.get(complement(a))[0])
    else:
        return a


def getLastUnivocalModBase(a):
    """Returns the last univocal modified base for the given base.
    This is the last value in the definition for an ambiguous
    modified base, provided that it is a modified nucleobase,
    or the base itself if it is already unambiguous."""
    if AMBIG_MOD_BASES.get(a):
        for b in reversed(AMBIG_MOD_BASES.get(a)):
            if not FULL_BASE_NAMES.get(b):
                return b
    elif AMBIG_MOD_BASES.get(complement(a)):
        for b in reversed(AMBIG_MOD_BASES.get(complement(a))):
            if not FULL_BASE_NAMES.get(b):
                return complement(b)
    else:
        return a


def makeList(lstOrVal):
    """Returns a list of a single item if the object passed is not
    already a list. This allows one to iterate over objects which
    may or may not already be lists (and therefore iterable)."""
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


def getCompMaybeFromMB(modBase):
    """Maps the given modified based to the corresponding
    modified guanine nucleobase (i.e. the complemented modified base).
    Applies the identity transformation if the given base is
    already a complemented modified nucleobase."""
    # use forward mapping
    return COMPLEMENTS.get(modBase) or modBase


def getMBMaybeFromComp(modBase):
    """Maps the given modified based to the corresponding
    modified cytosine nucleobase (i.e. the actual modified base).
    Applies the identity transformation if the given base is
    already a modified cytosine nucleobase."""
    # use reverse mapping (i.e. invert the bijection)
    return (~COMPLEMENTS).get(modBase) or modBase
