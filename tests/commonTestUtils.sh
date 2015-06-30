#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

EXIT_SUCCESS=0
EXIT_FAILURE=64

RED_COLOUR_CODE='\e[31m'
GREEN_COLOUR_CODE='\e[32m'

BASE_PROG_PATH="../../src"
DATA_DIR="../data"

function failMsgAndExit {
    echo -e "${RED_COLOUR_CODE}FAILED test $1." >&2
    exit $EXIT_FAILURE
}


function passMsg {
    echo -e "${GREEN_COLOUR_CODE}PASSED test $1." >&2
}


function runContainsTest {
    # runTest <program path> <test ID> <program arguments> <correct result>
    # checks if the result contains <correct result>
    if [[ ! "$($1 $3)" =~ "$4" ]]; then
        failMsgAndExit "$2"
    else
        passMsg "$2"
    fi
}

work_dir=$(mktemp -d --tmpdir=.)
cd "$work_dir"


function cleanup {
    rm  -Rf "../$work_dir"
}
trap cleanup EXIT
