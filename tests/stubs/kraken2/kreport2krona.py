#!/bin/bash
# A stand-in for kreport2krona.py, which records its arguments (separated by '|') in
# $STUB_DIR/kreport2krona.py.log and writes its -o file
set -e
( IFS="|"; echo "$*" ) >> "${STUB_DIR}/kreport2krona.py.log"
while (( $# > 0 )); do
    [[ $1 == "-o" ]] && echo "stub" > "$2"
    shift
done
