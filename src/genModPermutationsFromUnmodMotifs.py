#!/usr/bin/env python

"""Generates all possible modified motifs from input
consensus unmodified motifs.
Additionally returns the unmodified input motifs.
"""

from __future__ import with_statement, division, print_function

import argparse

import cytoUtils as cUtils

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

# Negate the flag, and store its value in the listed destination,
# instead of the name of the argument. This provides the user
# with the intuitive flag, while allowing us to use its programmatically
# convenient negation.
parser.add_argument('-D', '--mixedDinucs', action='store_false',
                    dest='singleModPerDinuc',
                    help="Permit dinucleotides to have different"
                    "modifications, such as 'm2' or 'h1'."
                    "By default, these are not permitted within"
                    "a single modified dinucleotide.")

parser.add_argument('-V', '--version', action='version',
                    version="%(prog)s " + __version__)

args = parser.parse_args()

if not all([cUtils.isUnivocal(modBase) for modBase in args.modifications]):
    cUtils.die("Only core modified base symbols can be provided to '-m'.")


def _cons_all_mods_at_pos(motif, modsToUse, singleModPerDinuc):
    """Wrapper to simplify the initial call of _cons_all_mods_at_pos.
       Constructs the required zipped list of modifications from the
       provided set of primary (positive-strand) modifications.

       Refer to cUtil's docstring for the definition of primary.
    """

    zippedMods = zip(modsToUse, cUtils.complement(modsToUse))

    # the initial list of output motifs (_output) starts out with
    # the unmodified motif itself
    return _cons_all_mods_at_pos_with_zipped_mods(motif, zippedMods,
                                                  0, [motif],
                                                  singleModPerDinuc)


def _cons_all_mods_at_pos_with_zipped_mods(motif, modsAndCompsToUse,
                                           pos=0, _output=None,
                                           singleModPerDinuc=True,
                                           prevDinucMod=None):
    """Recursively constructs all possible modifications of
       the provided motif, for the given modifications.
       The initial motif is assumed to be unmodified.
       The input modifications (modsAndCompsToUse) must
       be pre-processed as a zipped list of the modified
       bases and their complements.

       This recursion proceeds as a multitree, with a top-level
       bifurcation between an unmodified versus modified current
       motif position. This bifurcation occurs for each primary
       (or equivalently complement pair of) modification(s) under
       consideration. Such bifurcations occur at each level of the
       recursion, for each modified base under consideration.

       Refer to cUtil's docstring for the definition of primary.

       _output is a list, containing the output of the recursion,
       which will eventually contain all desired combinatorial motifs.
       Its underscore prefix merely denotes this particular usage.
    """

    nextPos = pos + 1

    # True iff there are more motif bases to process
    hasNextPos = nextPos < len(motif)

    # N.B. we cannot merely return _output here, if there are no more
    #      motif bases to process. This is due to our still needing to
    #      add the modified motif with only the current position modified
    #      to _output, which is done after the below recursive step.

    if hasNextPos:
        # continue to modify the motif, without modifying this pos
        _cons_all_mods_at_pos_with_zipped_mods(motif, modsAndCompsToUse,
                                               nextPos, _output,
                                               singleModPerDinuc)

    for modBase, modBaseComp in modsAndCompsToUse:
        modBaseToUse = ''

        nextPrevDinucMod = None

        if cUtils.isModifiableTo(motif[pos], modBase):
            modBaseToUse = modBase

            # This is the primary base of a dinucleotide, so
            # prevDinucMod should be set to this modBase,
            # for the next recursive step.
            nextPrevDinucMod = modBase
        elif cUtils.isModifiableTo(motif[pos], modBaseComp):
            modBaseToUse = modBaseComp

        if modBaseToUse:
            # Below conditional prevents dinucleotides with
            # different modifications, like m2.

            # N.B. if prevDinucMod is None, the below condition should evaluate
            #      to False. Therefore, checking if prevDinucMod evaluates to
            #      True, does need to be its own clause in this conditional.
            if singleModPerDinuc and prevDinucMod and prevDinucMod != modBase:
                continue

            newModMotif = motif[:pos] + modBaseToUse + motif[nextPos:]

            _output.append(newModMotif)

            if hasNextPos:
                # continue to modify the motif, modifying this pos
                # set the previous modified base if starting a mod. dinuc.
                _cons_all_mods_at_pos_with_zipped_mods(newModMotif,
                                                       modsAndCompsToUse,
                                                       nextPos, _output,
                                                       singleModPerDinuc,
                                                       nextPrevDinucMod)

    return _output


result = []

for motif in args.motifs:
    if any([cUtils.isModBase(base) for base in motif]):
        cUtils.die("All initial motifs must be unmodified.")

    result.extend(_cons_all_mods_at_pos(motif, args.modifications,
                                        args.singleModPerDinuc))

print(' '.join(result))
