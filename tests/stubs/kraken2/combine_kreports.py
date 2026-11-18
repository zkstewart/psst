#!/bin/bash
# A stand-in for combine_kreports.py, which records its arguments (separated by '|') in
# $STUB_DIR/combine_kreports.py.log and writes its -o file
set -e
( IFS="|"; echo "$*" ) >> "${STUB_DIR}/combine_kreports.py.log"
while (( $# > 0 )); do
    [[ $1 == "-o" ]] && echo "stub" > "$2"
    shift
done
