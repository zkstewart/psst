#!/bin/bash
# A stand-in for bbduk.sh, which records its arguments (one call per line, separated by
# '|') in $STUB_DIR/bbduk.log and copies each input read file to its output. If
# $STUB_FAIL is set, it fails after writing its outputs, as an interrupted run would.
set -e
LOG=()
for arg in "$@"; do
    LOG+=( "${arg}" )
    case ${arg} in
        in=*|in1=*) IN1=${arg#*=} ;;
        in2=*) IN2=${arg#*=} ;;
        out=*|out1=*) OUT1=${arg#*=} ;;
        out2=*) OUT2=${arg#*=} ;;
        stats=*) STATS=${arg#*=} ;;
    esac
done
( IFS="|"; echo "${LOG[*]}" ) >> "${STUB_DIR}/bbduk.log"
cp "${IN1}" "${OUT1}"
[[ -n ${IN2} ]] && cp "${IN2}" "${OUT2}"
echo "stub stats" > "${STATS}"
[[ -z ${STUB_FAIL} ]]
