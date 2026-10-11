#!/usr/bin/env python3

# psst.py provides a frontend interface for submitting a variety of
# jobs to the QUT HPC.
# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import \
    validate_at, \
        validate_at_generate, \
    validate_b, \
        validate_b_run, \
            validate_b_run_mmseqs2
from modules.blast import blast_filter
from _version import __version__

def main():
    usage = """%(prog)s launches PBS jobs with flexible configuration of values.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-v", "--version",
                   action="version",
                   version="psst.py {version}".format(version=__version__))
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
    subParentParser.add_argument("-v", "--version",
                                 action="version",
                                 version="psst.py {version}".format(version=__version__))
    
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    # Annotable subparser
    atparser = subparsers.add_parser("annotable",
                                     parents=[p],
                                     add_help=False,
                                     help="Annotation table (annotable) jobs")
    atparser.set_defaults(func=apmain)
    
    subAnnotableParsers = apparser.add_subparsers(dest="annotableMode",
                                               required=True)
    
    # Annotable > generate mode
    atgenerateparser = subAnnotableParsers.add_parser("generate",
                                                      parents=[p],
                                                      add_help=False,
                                                      help="Generate an annotable from BLAST + idmapping + go.obo files")
    ### To be populated with arguments by Claude
    
    # Blast subparser
    bparser = subparsers.add_parser("blast",
                                    parents=[p],
                                    add_help=False,
                                    help="BLAST (/MMseqs2) job handling")
    bparser.set_defaults(func=bmain)
    
    subBlastParsers = bparser.add_subparsers(dest="blastMode",
                                             required=True)
    
    # Blast > run subparser
    brunparser = subBlastParsers.add_parser("run",
                                            parents=[p],
                                            add_help=False,
                                            help="Run a BLAST (/MMseqs2) search job")
    brunparser.set_defaults(func=bmain)
    
    subBlastRunParsers = brunparser.add_subparsers(dest="blastRunMode",
                                                   required=True)
    
    # Blast > run > mmseqs2 mode
    brunmparser = subBlastRunParsers.add_parser("mmseqs2",
                                                parents=[p],
                                                add_help=False,
                                                help="Run a MMseqs2 job")
    ### To be populated with arguments by Claude
    
    args = subParentParser.parse_args()
    
    # Split into mode-specific functions
    if args.mode == "annotable":
        print("## psst.py - Annotation table handling ##")
        validate_at(args)
        atmain(args)
    if args.mode == "blast":
        print("## psst.py - BLAST handling ##")
        validate_b(args)
        bmain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def atmain(args):
    # Split into sub-mode-specific functions
    if args.annotableMode == "generate":
        print("## Annotable generation ##")
        validate_at_generate(args)
        annotable_generate(args)

def bmain(args):
    # Split into sub-mode-specific functions
    if args.blastMode == "run":
        validate_b_run(args)
        if args.blastToMode == "mmseqs2":
            print("## MMseqs2 query ##")
            validate_b_run_mmseqs2(args)
            blast_run_mmseqs2(args)
    
    print("BLAST handling complete!")

if __name__ == "__main__":
    main()
