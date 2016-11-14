#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

source $(dirname $0)/commonTestUtils.sh

# NB: all test cases use non-default ASCII-code order
PROGRAM_PATH="$BASE_PROG_PATH/genModPermutationsFromUnmodMotifs.py "

echo -e "------------------------\ngenModPermutationsFromUnmodMotifs.py\n------------------------" >&2


# -------------------------------- Test 1 --------------------------------

# check that the expected number of combinations are output, when allowing all combinations

function check_num_combos {
    mods="$1"

    # exponent base is the number of mods + 1 (to account for unmod base)
    exp_base="$((${#mods} + 1))"

    CG_STRING=$(echo 'CG'{,,,,,,,,} | tr -d ' ')

    for i in {2..8}; do
        num_motifs=$($PROGRAM_PATH -D -m "$mods" ${CG_STRING:0:$i} | wc -w)

        if [[ $num_motifs -eq $(($exp_base**$i)) ]]; then
            passMsg "1${2:-$mods}-$i"
        else
            failMsgAndExit "1A-$i: incorrect number of motifs ($num_motifs)"
        fi
    done
}

# 1A) one modification (base 2)
check_num_combos 'm' 'A'

# 1B) two modifications (base 3)
check_num_combos 'mh' 'B'

# 2) manually check a simple case

result=$($PROGRAM_PATH -m 'm' 'CG' | tr ' ' '\n' | sort | tr '\n' ' ')

if [[ $result =~ 'C1 CG m1 mG' ]]; then
    passMsg "2"
else
    failMsgAndExit "2: got $result"
fi

# 3 ) verify that results of a simple test yield only unique motifs
result=$($PROGRAM_PATH -m 'mh' 'CACGTG' | tr ' ' '\n')

if diff -q <(echo "$result" | sort) <(echo "$result" | sort -u) > /dev/null  2>&1; then
    passMsg "3"
else
    failMsgAndExit "3: differences found"
fi

# 4) check that the difference between '-D' is indeed the expected omission of mixed dinucleotide modifications

if diff -q <($PROGRAM_PATH -D -m 'mh' 'CGCG' | tr ' ' '\n' | egrep -v 'm2|h1') <($PROGRAM_PATH -m 'mh' 'CGCG' | tr ' ' '\n'); then
    passMsg "4"
else
    failMsgAndExit "4: differences found"
fi

