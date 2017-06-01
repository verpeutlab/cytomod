#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit

source $(dirname $0)/commonTestUtils.sh

PROGRAM_PATH="$BASE_PROG_PATH/cytomod.py"
ARCHIVE_PATH="$DATA_DIR/archive/"
TRACKS_BASE_PATH="$DATA_DIR/testTracks"
TRACK_PREF="$DATA_DIR/mm9_chrY-only_"
VERBOSITY_ARG='-vvvv'

TEST_REGION_CHR='Y'
TEST_REGION_START=1858490
TEST_REGION_END=1858513
TEST_4_REGION_END=1858558
TEST_REGION_LEN=$(echo "$TEST_REGION_END - $TEST_REGION_START" | bc)
TEST_REGION="chr$TEST_REGION_CHR:$TEST_REGION_START-$TEST_REGION_END"
TEST_4_REGION="chr$TEST_REGION_CHR:$TEST_REGION_START-$TEST_4_REGION_END"
TEST_REGION_CORRECT_UNMASKED_RES='Tmm12TTf1AAxAxxAGATATCT'
TEST_REGION_CORRECT_ONLY_UNSET_MASKED_RES='Tmm12TTf1AAxAxxAGATATzT'
TEST_REGION_CORRECT_MASKED_RES='Tmm12TTz9AAzAzzA9ATATzT'

# "a value at and below which the locus is considered ambiguous"
DEFAULT_MASK_VALUE=$(grep -oP '_DEFAULT_MASK_VALUE = \K\d+' \
                     $PROGRAM_PATH)
if [[ ! $DEFAULT_MASK_VALUE =~ [[:digit:]]+ ]]; then
    failMsgAndExit 'SETUP'
fi

function cleanup {
    rm  -Rf "../$work_dir" "$ARCHIVE_PATH" "$TRACKS_BASE_PATH"*
}

function maybeFilter {
    # NB: masking only occurs for entires with a value of 0
    #     whereas for all non-masked bases, a strictly positive value
    #     is needed for a modification at that position to occur.
    relation_op='>'
    if [[ $1 =~ 'MASK' ]]; then
        relation_op='<=' 
    fi
    zcat $1 | awk "\$4 $relation_op $DEFAULT_MASK_VALUE {print;}"
}


test_to_run=${1:-0}

echo -e "------------------------\nCytomod.py\n------------------------" >&2

case $test_to_run in
0|1)
    # -------------------------------- Test 1 --------------------------------
    FASTA_file='test1A.fa'
    track_out_dir="$TRACKS_BASE_PATH/"
    mkdir "$track_out_dir"
    
    # 1A) check that no masking was performed
    $PROGRAM_PATH $VERBOSITY_ARG -d "$DATA_DIR/" "$DATA_DIR/" -f "$FASTA_file" \
        -r "$TEST_REGION" --BEDOutDir "$track_out_dir"
    if [[ ! -z $(fgrep -v '>' $FASTA_file | grep '[z9]') ]]; then
        failMsgAndExit '1A'
    else
        passMsg '1A'
    fi

    # 1B) check that all expected BED files were generated
        # check that the use of a custom track directory was successful
    if [[ $(find "$track_out_dir" -name "*track*" -and -name "*-*.bed*" | fgrep -f \
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
    mkdir "$DATA_DIR/test"
    $PROGRAM_PATH $VERBOSITY_ARG -d "$DATA_DIR/" "$DATA_DIR/" --archiveOutDir "$DATA_DIR/test" --archiveOutName 'testArchive' -f "$FASTA_file" -b -M
    
    # check that the use of a custom archive directory and name was successful
    if [[ ! -d "$DATA_DIR/test/testArchive" ]]; then
        failMsgAndExit "2: Expected archive not present."
    fi
    
    declare -A BED_to_symbols
    BED_to_symbols=(["$TRACK_PREF"'5mC-fakeData.bedGraph.gz']='m1h2f3z9NAT' \
                    ["$TRACK_PREF"'5hmC-fakeData.bedGraph.gz']='h2f3z9NAT' \
                    ["$TRACK_PREF"'5fC-fakeData.bedGraph.gz']='f3z9NAT' \
                    ["$TRACK_PREF"'5xC-fakeData.bedGraph.gz']='x7m1h2f3z9NAT' \
                    ["$TRACK_PREF"'MASK-cov.bedGraph.gz']='z9NAT')

    # check that the FASTA file generated contains
    # the expected modifications at the correct loci
    
    for key in "${!BED_to_symbols[@]}"; do
        if [[ ! -z $(maybeFilter "$key" | bedtools getfasta \
                     -fi $FASTA_file -bed stdin -fo stdout | fgrep -v '>' | \
                     grep -v "[${BED_to_symbols[$key]}]") ]]; then
            failMsgAndExit "2: $key\t${BED_to_symbols[$key]}"
        fi
    done
    passMsg '2'
    ;&
0|3)
    # -------------------------------- Test 3 --------------------------------
    cytomod_base_cmd="$PROGRAM_PATH $VERBOSITY_ARG -G $ARCHIVE_PATH -b"
    test_len=47
    # test that a random region is retrieved and of the correct length
    if [[ $($cytomod_base_cmd -R $test_len | awk '{print length($0);}') \
          -ne $test_len ]]; then
        failMsgAndExit '3A'
    else
        passMsg '3A'
    fi
    
    region_test_unmasked_res=$($cytomod_base_cmd -r $TEST_REGION)

    # below three should be identical, in this case
    region_test_masked_res=$($cytomod_base_cmd -r $TEST_REGION -M)
    region_test_unset_masked_res=$($cytomod_base_cmd -r $TEST_REGION --maskAllUnsetRegions)
    region_test_masked_and_unset_masked_res=$($cytomod_base_cmd -r $TEST_REGION -M --maskAllUnsetRegions)

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
        if [[ "$region_test_masked_res" != "$(echo "$TEST_REGION" | \
                                            tr ':-' '\t' | \
                                            bedtools getfasta \
                                           -fi $FASTA_file \
                                           -bed stdin -fo stdout | \
                                           fgrep -v '>')" ]]; then
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
        failMsgAndExit '3D: masked ("-M")'
    elif [[ "$region_test_unset_masked_res" != \
          "$TEST_REGION_CORRECT_ONLY_UNSET_MASKED_RES" ]]; then
        failMsgAndExit '3D: unset masked ("--maskAllUnsetRegions")'
    elif [[ "$region_test_masked_and_unset_masked_res" != \
          "$TEST_REGION_CORRECT_MASKED_RES" ]]; then
        failMsgAndExit '3D: both masked ("-M" and "--maskAllUnsetRegions")'
    else
        passMsg '3D'
    fi
    ;&
0|4)
    # -------------------------------- Test 4 --------------------------------
    track_out_union_dir="${TRACKS_BASE_PATH}Union/"
    track_out_intersection_dir="${TRACKS_BASE_PATH}Intersection/"

    mkdir "$track_out_union_dir" "$track_out_intersection_dir"
    
    # 4)
    # create archive
    $PROGRAM_PATH $VERBOSITY_ARG -d "$DATA_DIR/" "$DATA_DIR/" -f /dev/null -b

    # query without and then with intersection
    $PROGRAM_PATH $VERBOSITY_ARG -G "$ARCHIVE_PATH" -B \
        -r "$TEST_4_REGION" --BEDOutDir "$track_out_union_dir"

    $PROGRAM_PATH $VERBOSITY_ARG -G "$ARCHIVE_PATH" -B \
        -r "$TEST_4_REGION" --BEDOutDir "$track_out_intersection_dir" -I

    # intersected version should have exactly one fewer 5xC site
    if [[ $(($(zcat "$track_out_union_dir/track-modGenome-x.bed.gz" | wc -l) - \
          $(zcat "$track_out_intersection_dir/track-modGenome-x.bed.gz" | wc -l))) != 1 ]]; then
        failMsgAndExit '4'
    else
        passMsg '4'
    fi
    ;&
esac

exit $EXIT_SUCCESS
