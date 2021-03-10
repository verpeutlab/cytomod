#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

EXIT_SUCCESS=0
EXIT_FAILURE=64

RED_COLOUR_CODE='\e[31m'
GREEN_COLOUR_CODE='\e[32m'
RESET_COLOUR_CODE='\x1b[39;49m'

BASE_PROG_PATH="../../src"
DATA_DIR="../data"


function failMsgAndExit {
    echo -e "${RED_COLOUR_CODE}FAILED test $1.${RESET_COLOUR_CODE}" >&2
    exit $EXIT_FAILURE
}


function passMsg {
    echo -e "${GREEN_COLOUR_CODE}PASSED test $1.${RESET_COLOUR_CODE}" >&2
}


function _performContainsTest {
    # _performContainsTest <test ID> <test (program) result> <correct result>
    # checks if the <test (program) result> contains <correct result>
    if [[ ! "$2" =~ "$3" ]]; then
        failMsgAndExit "$1"
    else
        passMsg "$1"
    fi
}


function runContainsTest {
    # runContainsTest <program path> <test ID> <program arguments> <correct result>
    _performContainsTest "$2" "$($1 $3 || true)" "$4"
}


function runContainsFileTest {
    # runContainsFileTest <program path> <test ID> <program arguments> <correct result> [file]
    run_errors="$($1 $3 2>&1 >/dev/null || true)"
    if [[ $(echo "$run_errors") ]]; then
        failMsgAndExit $2
    fi

    _performContainsTest "$2" "$(cat *.meme)" "$4"

    rm -f *.meme
}


work_dir=$(mktemp -d --tmpdir=.)
cd "$work_dir"


function cleanup {
    rm  -Rf "../$work_dir"
}
trap cleanup EXIT
