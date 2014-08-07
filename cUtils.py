from __future__ import with_statement, division, print_function

VERSION = "$Revision: 0.05$"

import sys
import re
import datetime
import textwrap
import functools
from bidict import bidict
from collections import OrderedDict
from itertools import izip


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
# Create a dictionary mapping each "primary" modified base to
# the base it modifies
_MODIFIES = dict.fromkeys(MOD_BASES.values(), 'C')


# All IUPAC nucleobases and their complements, plus 'X',
# which is just an additional alias for any nucleobase
COMPLEMENTS = {'A': 'T', 'G': 'C',
               'R': 'Y', 'M': 'K',
               'W': 'W', 'S': 'S',
               'B': 'V', 'D': 'H',
               'N': 'N', 'X': 'X'}
# Add all modified nucleobases
modifiedBasesToComplements = \
    izip(MOD_BASES.values(), ''.
         join(str(i) for i in range(1, len(MOD_BASES) + 1)))
# Add all modified nucleobase complements
# They are ordered by the originating modification's oxidation order
COMPLEMENTS.update(modifiedBasesToComplements)
COMPLEMENTS = bidict(COMPLEMENTS)

_FULL_BASE_NAMES = {'A': 'Adenine', 'T': 'Thymine',
                    'G': 'Guanine', 'C': 'Cytosine'}

MOD_BASE_NAMES = {'m': '5mC', 'h': '5hmC', 'f': '5fC', 'c': '5caC'}
_FULL_MOD_BASE_NAMES = {'m': '5-Methylcytosine',
                        'h': '5-Hydroxymethylcytosine',
                        'f': '5-Formylcytosine',
                        'c': '5-Carboxylcytosine'}


@unpacked
def complement(bases):
    """Complements the given, potentially modified, base."""
    return [COMPLEMENTS.get(b) or (~COMPLEMENTS).get(b) for b in bases]
# Update the dictionary mapping with every complemented modification
_MODIFIES.update(izip(complement(_MODIFIES.keys()),
                 complement(_MODIFIES.values())))

# Add the names of all modified base complements, using existing nomenclature
for b in complement(MOD_BASES.values()):
    _FULL_MOD_BASE_NAMES.update(izip(b, [_FULL_BASE_NAMES[_MODIFIES[b]] +
                                ':' + _FULL_MOD_BASE_NAMES[complement(b)]]))

_PARAM_A_CONST_VAL = 999


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
