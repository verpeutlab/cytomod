#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

EXIT_SUCCESS=0
EXIT_FAILURE=64

RED_COLOUR_CODE='\e[31m'
GREEN_COLOUR_CODE='\e[32m'

ARCHIVE_PATH='../data/archive/'
TRACK_PREF='../data/mm9_chrY-only_'
VERBOSITY_ARG='-vvvv'

TEST_REGION_CHR='Y'
TEST_REGION_START=1858490
TEST_REGION_END=1858505
TEST_REGION_LEN=$(echo "$TEST_REGION_END - $TEST_REGION_START" | bc)
TEST_REGION="chr$TEST_REGION_CHR:$TEST_REGION_START-$TEST_REGION_END"
TEST_REGION_CORRECT_UNMASKED_RES='Tmm12TTf1AAxAxx'
TEST_REGION_CORRECT_MASKED_RES='Tmm12TTz9AAzAzz'

# "a value at and below which the locus is considered ambiguous"
DEFAULT_MASK_VALUE=$(grep -oP '_DEFAULT_MASK_VALUE = \K\d+' ../src/cytomod.py)


function failMsgAndExit {
    echo -e "${RED_COLOUR_CODE}FAILED test $1." >&2
    exit $EXIT_FAILURE
}


function passMsg {
    echo -e "${GREEN_COLOUR_CODE}PASSED test $1." >&2
}


work_dir=$(mktemp -d --tmpdir=.)
cd "$work_dir"


function cleanup {
    rm  -Rf "../$work_dir" "$ARCHIVE_PATH"
}
trap cleanup EXIT


test_to_run=${1:-0}

case $test_to_run in
0|1)
    # -------------------------------- Test 1 --------------------------------
    FASTA_file='test1A.fa'

    # 1A) check that no masking was performed
    ../../src/cytomod.py $VERBOSITY_ARG -d ../data/ ../data/ -f $FASTA_file \
        -r "$TEST_REGION"
    if [[ ! -z $(fgrep -v '>' $FASTA_file | grep '[z9]') ]]; then
        failMsgAndExit '1A'
    else
        passMsg '1A'
    fi

    # 1B) check that all expected BED files were generated
    if [[ $(find . -name "*track*" -and -name "*-*.bed*" | fgrep -f \
            ../expected_track_patterns.txt | awk 'END{print NR;}') -ne \
          $(awk 'END{print NR;}' ../expected_track_patterns.txt) ]]; then
        failMsgAndExit '1B'
    else
        passMsg '1B'
    fi
    ;&
0|2)
    # -------------------------------- Test 2 --------------------------------
    FASTA_file='test2.fa'

    ../../src/cytomod.py $VERBOSITY_ARG -d ../data/ ../data/ -f $FASTA_file -b -M

    declare -A BED_to_symbols
    BED_to_symbols=(["$TRACK_PREF"'5mC-fakeData.bedGraph.gz']='ATm1h2f3z9' \
                    ["$TRACK_PREF"'5hmC-fakeData.bedGraph.gz']='ATh2f3z9' \
                    ["$TRACK_PREF"'5fC-fakeData.bedGraph.gz']='ATf3z9' \
                    ["$TRACK_PREF"'5xC-fakeData.bedGraph.gz']='ATx7m1h2f3z9' \
                    ["$TRACK_PREF"'MASK-cov.bedGraph.gz']='ATz9')

    # check that the FASTA file generated contains
    # the expected modifications at the correct loci
    for key in "${!BED_to_symbols[@]}"; do
        additionalFilterCmd=''
        if [[ $key =~ 'MASK' ]]; then
            additionalFilterCmd=' | awk "$4 <= $DEFAULT_MASK_VALUE {print;}"'
        fi
        if [[ ! -z $(zcat "$key" $additionalFilterCmd | bedtools getfasta \
                     -fi $FASTA_file -bed stdin -fo stdout | fgrep -v '>' | \
                     grep -v "[${BED_to_symbols[$key]}]") ]]; then
            failMsgAndExit "2: $key\t${BED_to_symbols[$key]}"
    done
    passMsg '2'
    ;&
0|3)
    # -------------------------------- Test 3 --------------------------------
    cytomod_base_cmd="../../src/cytomod.py $VERBOSITY_ARG -G $ARCHIVE_PATH -b"
    test_len=47
    # test that a random region is retrieved and of the correct length
    if [[ $($cytomod_base_cmd -R $test_len | awk '{print length($0);}') \
          -ne $test_len ]]; then
        failMsgAndExit '3A'
    else
        passMsg '3A'
    fi
    
    region_test_unmasked_res=$($cytomod_base_cmd -r $TEST_REGION)
    region_test_masked_res=$($cytomod_base_cmd -r $TEST_REGION -M)
    # test that a queried region is retrieved and of the correct length
    if [[ $(echo "$region_test_unmasked_res" | awk '{print length($0);}') \
          -ne TEST_REGION_LEN ]]; then
        failMsgAndExit '3B'
    else
        passMsg '3B'
    fi
    
    # test that the region queried is actually the intended target region
    if [[ ! -f "$FASTA_file" ]]; then
        echo "Unable to find $FASTA_file. Test 3C skipped."
    else
        if [[ "$region_test_masked_res" != $(echo "$TEST_REGION" | \
                                            tr ':-' '\t' | \
                                            bedtools getfasta \
                                           -fi $FASTA_file \
                                           -bed stdin -fo stdout) ]]; then
            failMsgAndExit '3C'
        else
            passMsg '3C'
        fi
    fi
    
    # test that query region matches the ground-truth for that region
    if [[ "$region_test_unmasked_res" != \
          "$TEST_REGION_CORRECT_UNMASKED_RES" ]]; then
        failMsgAndExit '3D: unmasked'
    elif [[ "$region_test_masked_res" != \
          "$TEST_REGION_CORRECT_MASKED_RES" ]]; then
        failMsgAndExit '3D: masked'
    else
        passMsg '3D'
    fi
    ;&
esac

exit $EXIT_SUCCESS
