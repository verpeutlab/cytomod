#!/usr/bin/env python

"""Generates all possible modified motifs from input
consensus unmodified motifs.
Additionally returns the unmodified input motifs.
"""

from __future__ import with_statement, division, print_function

import argparse

import cUtils

__version__ = cUtils.__version__

parser = argparse.ArgumentParser()

parser.add_argument('motifs', nargs='+',
                    help="Input unmodified motifs from which modified motifs"
                    "will be generated.")

# modifications will store a list of each provided modification
parser.add_argument('-m', '--modifications', default='mh', type=list,
                    help="The different modifications to use to modify"
                    "the unmodified motifs. By default, this is set to"
                    "5mC and 5hmC, but can be modified by providing the"
                    "desired symbols in a single string, e.g.: 'mh'."
                    "Complements will be automatically used when needed.",
                    )

parser.add_argument('-D', '--mixedDinucs', action='store_true',
                    help="Permit dinucleotides to have different"
                    "modifications, such as 'm2' or 'h1'."
                    "By default, these are not permitted within"
                    "a single modified dinucleotide.")

parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)

args = parser.parse_args()


if not all([cUtils.isUnivocal(modBase) for modBase in args.modifications]):
    cUtils.die("Only core modified base symbols can be provided to '-m'.")


def _cons_all_mods_at_pos_initial_call(motif, modsToUse, mixedDinucs):
    """Wrapper to simplify the initial call of _cons_all_mods_at_pos.
       Constructs the required zipped list of modifications from the
       provided set of primary (positive-strand) modifications.
    """

    zippedMods = zip(modsToUse, cUtils.complement(modsToUse))

    return _cons_all_mods_at_pos(motif, zippedMods, 0, None,
                                 not mixedDinucs, None)


def _cons_all_mods_at_pos(motif, modsAndCompsToUse, pos=0, _output=None,
                          singleModPerDinuc=True, prevDinucMod=None):
    """Recursively constructs all possible modifications of
       the provided motif, for the given modifications.
       The initial motif is assumed to be unmodified.
       The input modificatins (modsAndCompsToUse) must
       be pre-processed as a zipped list of the modified
       bases and their complements.
    """

    nextPos = pos + 1

    # True iff there are more motif bases to process
    hasNextPos = nextPos < len(motif)

    if _output is None:
        _output = [motif]

    if hasNextPos:
        # continue to modify the motif, without modifying this pos
        _cons_all_mods_at_pos(motif, modsAndCompsToUse,
                              nextPos, _output, singleModPerDinuc)

    for modBase, modBaseComp in modsAndCompsToUse:
        modBaseToUse = ''

        isPrimaryOfDinuc = False
        if cUtils.isModifiableTo(motif[pos], modBase):
            modBaseToUse = modBase
            isPrimaryOfDinuc = True
        elif cUtils.isModifiableTo(motif[pos], modBaseComp):
            modBaseToUse = modBaseComp

        if modBaseToUse:
            # prevent dinucleotides with different modifications, like m2
            if singleModPerDinuc and prevDinucMod and prevDinucMod != modBase:
                continue

            newModMotif = motif[:pos] + modBaseToUse + motif[nextPos:]

            _output.append(newModMotif)

            if hasNextPos:
                # continue to modify the motif, modifying this pos
                # set the previous modified base if starting a mod. dinuc.
                _cons_all_mods_at_pos(newModMotif, modsAndCompsToUse,
                                      nextPos, _output, singleModPerDinuc,
                                      modBase if isPrimaryOfDinuc else None)

    return _output


result = []

for motif in args.motifs:
    if any([cUtils.isModBase(base) for base in motif]):
        cUtils.die("All initial motifs must be unmodified.")

    result.extend(_cons_all_mods_at_pos_initial_call(motif, args.modifications,
                                                     args.mixedDinucs))

print(' '.join(result))
