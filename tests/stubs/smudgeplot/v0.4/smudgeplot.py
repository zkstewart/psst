#!/bin/bash
# A stand-in for smudgeplot < 0.5, run as 'smudgeplot.py'; see ../v0.5/smudgeplot
STUB_OLD=1 exec "$(dirname "$0")/../v0.5/smudgeplot" "$@"
