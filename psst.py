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
from modules.annotable import annotable_generate
from modules.blast import blast_run_mmseqs2
from modules.locations import MMSEQS_DIR, VARSCRIPT_DIR, REFSEQ_PLANT_FASTA, UNIREF90_FASTA, \
    IDMAPPING_FILE, OBO_FILE
from _version import __version__

def add_pbs_arguments(parser, jobName, walltime, mem, cpus):
    '''
    Adds PBS job arguments to a leaf subparser with job-specific default values.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
        jobName -- a string of the default PBS job name.
        walltime -- a string of the default walltime in 'HH:MM:SS' format.
        mem -- a string of the default memory request (e.g., '70G').
        cpus -- an int of the default number of CPUs to request.
    '''
    g = parser.add_argument_group("PBS job arguments")
    g.add_argument("--cpus", dest="cpus",
                   required=False,
                   type=int,
                   help=f"""Optionally, specify the number of CPUs to request;
                   default == {cpus}""",
                   default=cpus)
    g.add_argument("--mem", dest="mem",
                   required=False,
                   help=f"""Optionally, specify the memory to request;
                   default == '{mem}'""",
                   default=mem)
    g.add_argument("--walltime", dest="walltime",
                   required=False,
                   help=f"""Optionally, specify the walltime to request in 'H', 'H:MM',
                   or 'H:MM:SS' format (max 48 hours); default == '{walltime}'""",
                   default=walltime)
    g.add_argument("--name", dest="jobName",
                   required=False,
                   help=f"""Optionally, specify the PBS job name;
                   default == '{jobName}'""",
                   default=jobName)
    g.add_argument("--rpid", dest="rpid",
                   required=True,
                   help="Research Project ID (RPID) to attach to the job via '#PBS -P'")
    g.add_argument("--script", dest="scriptFile",
                   required=False,
                   help="""Optionally, specify where to write the job script;
                   default == '<outprefix>.<job>.sh'""",
                   default=None)
    g.add_argument("--qsub", dest="qsub",
                   required=False,
                   action="store_true",
                   help="""Optionally, submit the job script with qsub after
                   writing it""",
                   default=False)
    
    d = g.add_mutually_exclusive_group()
    d.add_argument("--afterok", dest="afterok",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more job IDs which must complete
                   successfully before this job starts""",
                   default=None)
    d.add_argument("--afternotok", dest="afternotok",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify one or more job IDs which must fail
                   before this job starts""",
                   default=None)

def add_mmseqs2_arguments(parser):
    '''
    Adds program location arguments to a leaf subparser for any job that runs
    Various_scripts/run_mmseqs2.py.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.add_argument("--mmseqs", dest="mmseqsDir",
                        required=False,
                        help=f"""Optionally, specify the directory containing the mmseqs
                        executable; default == '{MMSEQS_DIR}'""",
                        default=MMSEQS_DIR)
    parser.add_argument("--varscripts", dest="varScriptDir",
                        required=False,
                        help=f"""Optionally, specify the location of the Various_scripts
                        repository; default == '{VARSCRIPT_DIR}'""",
                        default=VARSCRIPT_DIR)

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
    atparser.set_defaults(func=atmain)
    
    subAnnotableParsers = atparser.add_subparsers(dest="annotableMode",
                                               required=True)
    
    # Annotable > generate mode
    atgenerateparser = subAnnotableParsers.add_parser("generate",
                                                      parents=[p],
                                                      add_help=False,
                                                      help="Generate an annotable from BLAST + idmapping + go.obo files")
    atgenerateparser.add_argument("-q", dest="query",
                                  required=True,
                                  help="Location of query FASTA file")
    atgenerateparser.add_argument("-o", dest="outputPrefix",
                                  required=True,
                                  help="Prefix for output files")
    atgenerateparser.add_argument("-d", dest="database",
                                  required=False,
                                  help=f"""Optionally, specify the location of the database FASTA file;
                                  default == '{REFSEQ_PLANT_FASTA}'""",
                                  default=REFSEQ_PLANT_FASTA)
    atgenerateparser.add_argument("--db", dest="databaseTag",
                                  required=False,
                                  choices=["UniRef100", "UniRef90", "UniRef50", "RefSeq", "."],
                                  help="""Optionally, specify the database type to match -d;
                                  default == 'RefSeq'""",
                                  default="RefSeq")
    atgenerateparser.add_argument("--idmapping", dest="idmappingFile",
                                  required=False,
                                  help=f"""Optionally, specify the location of the idmapping_selected.tab file;
                                  default == '{IDMAPPING_FILE}'""",
                                  default=IDMAPPING_FILE)
    atgenerateparser.add_argument("--obo", dest="oboFile",
                                  required=False,
                                  help=f"""Optionally, specify the location of the go.obo file;
                                  default == '{OBO_FILE}'""",
                                  default=OBO_FILE)
    atgenerateparser.add_argument("--numhits", dest="numHits",
                                  required=False,
                                  type=int,
                                  help="""Optionally, specify the number of hits for each sequence
                                  to report; default == 15""",
                                  default=15)
    atgenerateparser.add_argument("--iterations", dest="iterations",
                                  required=False,
                                  type=int,
                                  help="""Optionally, specify the number of MMseqs2 search iterations;
                                  default == 1""",
                                  default=1)
    atgenerateparser.add_argument("--tmp", dest="tmpDir",
                                  required=False,
                                  help="""Optionally, specify the MMseqs2 tmp directory;
                                  default == '<outprefix>_tmp'""",
                                  default=None)
    add_mmseqs2_arguments(atgenerateparser)
    add_pbs_arguments(atgenerateparser, jobName="annotable",
                      walltime="06:00:00", mem="70G", cpus=12)
    
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
    brunmparser.add_argument("-q", dest="query",
                             required=True,
                             help="Location of query FASTA file")
    brunmparser.add_argument("-o", dest="outputPrefix",
                             required=True,
                             help="Prefix for output files")
    brunmparser.add_argument("--tmp", dest="tmpDir",
                             required=True,
                             help="Location of MMseqs2 tmp directory")
    brunmparser.add_argument("-n", dest="iterations",
                             required=True,
                             type=int,
                             help="Number of MMseqs2 search iterations")
    brunmparser.add_argument("-s", dest="sensitivity",
                             required=True,
                             type=float,
                             choices=[1, 2, 3, 4, 5, 5.7, 6, 7, 7.5],
                             help="MMseqs2 search sensitivity")
    brunmparser.add_argument("-t", dest="searchType",
                             required=True,
                             choices=["auto", "blastn", "tblastx"],
                             help="MMseqs2 search type")
    brunmparser.add_argument("-d", dest="database",
                             required=False,
                             help=f"""Optionally, specify the location of the database FASTA file;
                             default == '{UNIREF90_FASTA}'""",
                             default=UNIREF90_FASTA)
    add_mmseqs2_arguments(brunmparser)
    add_pbs_arguments(brunmparser, jobName="mms2",
                      walltime="46:00:00", mem="750G", cpus=32)
    
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
    
    print("Annotable handling complete!")

def bmain(args):
    # Split into sub-mode-specific functions
    if args.blastMode == "run":
        validate_b_run(args)
        if args.blastRunMode == "mmseqs2":
            print("## MMseqs2 query ##")
            validate_b_run_mmseqs2(args)
            blast_run_mmseqs2(args)
    
    print("BLAST handling complete!")

if __name__ == "__main__":
    main()
