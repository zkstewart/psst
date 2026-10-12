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
    validate_as, \
        validate_as_hifiasm, validate_as_revcomp, validate_as_scaffolding, \
    validate_at, \
        validate_at_generate, \
    validate_b, \
        validate_b_run, \
            validate_b_run_mmseqs2, \
    validate_pbs, RPID_SHORTCUTS
from modules.annotable import annotable_generate
from modules.assembly import assembly_hifiasm, assembly_revcomp, assembly_scaffolding
from modules.blast import blast_run_mmseqs2
from modules.locations import MMSEQS_DIR, VARSCRIPT_DIR, GENSCRIPT_DIR, COMPLEASM_DIR, \
    REFSEQ_PLANT_FASTA, UNIREF90_FASTA, IDMAPPING_FILE, OBO_FILE, \
    BUSK_HAP1_FASTA, BUSK_HAP2_FASTA
from _version import __version__

# Argument conventions:
# - All default values are set in the parser; help text should reference them with
#   %(default)s so the help can never disagree with the actual default.
# - Arguments whose values are checked by validation (e.g., file or program locations)
#   should use action=UserSet. Validation only checks values the user explicitly gave,
#   so default locations on the HPC are never checked.
# - Arguments without a sensible default should default to None, which leaves an
#   '@@field@@' placeholder in the job script for the user to fill in.
# - Each leaf parser must set_defaults(func=<job function>, scriptSuffix=<str>).

class UserSet(argparse.Action):
    '''
    Stores an argument's value and records its dest in namespace.userSet, allowing
    validation to distinguish user-specified values from defaults.
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        if not hasattr(namespace, "userSet"):
            namespace.userSet = set()
        namespace.userSet.add(self.dest)

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
                   help="""Optionally, specify the number of CPUs to request;
                   default == %(default)s""",
                   default=cpus)
    g.add_argument("--mem", dest="mem",
                   required=False,
                   help="""Optionally, specify the memory to request;
                   default == '%(default)s'""",
                   default=mem)
    g.add_argument("--walltime", dest="walltime",
                   required=False,
                   help="""Optionally, specify the walltime to request in 'H', 'H:MM',
                   or 'H:MM:SS' format (max 48 hours); default == '%(default)s'""",
                   default=walltime)
    g.add_argument("--name", dest="jobName",
                   required=False,
                   help="""Optionally, specify the PBS job name;
                   default == '%(default)s'""",
                   default=jobName)
    shortcuts = ", ".join(f"{k}={v}" for k, v in RPID_SHORTCUTS.items())
    g.add_argument("--rpid", dest="rpid",
                   required=True,
                   help=f"""Research Project ID (RPID) to attach to the job via '#PBS -P';
                   accepts a full RPID or a shortcut ({shortcuts})""")
    g.add_argument("--script", dest="scriptFile",
                   required=False,
                   help="""Optionally, specify where to write the job script;
                   default == '<outprefix>.<job>.sh' if -o is given, otherwise
                   '<job>.sh' in the current directory""",
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

def add_mmseqs2_arguments(parser, database):
    '''
    Adds arguments shared by any job that runs Various_scripts/run_mmseqs2.py.

    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
        database -- a string of the default database FASTA file location.
    '''
    parser.add_argument("-q", "--query", dest="query",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the location of the query FASTA file;
                        default == left as a placeholder""",
                        default=None)
    parser.add_argument("-d", "--database", dest="database",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the location of the database FASTA file;
                        default == '%(default)s'""",
                        default=database)
    parser.add_argument("-o", "--outprefix", dest="outputPrefix",
                        required=False,
                        help="""Optionally, specify the prefix for output files;
                        default == left as a placeholder""",
                        default=None)
    parser.add_argument("-n", "--iterations", dest="iterations",
                        required=False,
                        type=int,
                        help="""Optionally, specify the number of MMseqs2 search iterations;
                        default == %(default)s""",
                        default=1)
    parser.add_argument("--tmp", dest="tmpDir",
                        required=False,
                        help="""Optionally, specify the MMseqs2 tmp directory;
                        default == '<outprefix>_tmp' if -o is given, otherwise
                        left as a placeholder""",
                        default=None)
    parser.add_argument("--mmseqs", dest="mmseqsDir",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the directory containing the mmseqs
                        executable; default == '%(default)s'""",
                        default=MMSEQS_DIR)
    add_varscripts_argument(parser)

def add_varscripts_argument(parser):
    '''
    Adds the Various_scripts location argument to a leaf subparser for any job that
    runs scripts from the Various_scripts repository.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.add_argument("--varscripts", dest="varScriptDir",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the location of the Various_scripts
                        repository; default == '%(default)s'""",
                        default=VARSCRIPT_DIR)

def add_genscripts_argument(parser):
    '''
    Adds the Genome_analysis_scripts location argument to a leaf subparser for any job
    that runs scripts from the Genome_analysis_scripts repository.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.add_argument("--genscripts", dest="genScriptDir",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the location of the Genome_analysis_scripts
                        repository; default == '%(default)s'""",
                        default=GENSCRIPT_DIR)

def add_preset_arguments(parser):
    '''
    Adds output directory and minimap2 preset arguments to a leaf subparser for
    Various_scripts/assembly pipelines that share these behavioural parameters.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.add_argument("--outdir", dest="outputDir",
                        required=False,
                        help="""Optionally, specify the output directory, relative to the
                        directory the job is submitted from; default == '%(default)s'""",
                        default="output")
    parser.add_argument("--preset", dest="preset",
                        required=False,
                        choices=["asm5", "asm10", "asm20"],
                        help="""Optionally, specify the minimap2 alignment preset;
                        default == '%(default)s'""",
                        default="asm20")

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

    subAnnotableParsers = atparser.add_subparsers(dest="annotableMode",
                                                  required=True)

    # Annotable > generate mode
    atgenerateparser = subAnnotableParsers.add_parser("generate",
                                                      parents=[p],
                                                      add_help=False,
                                                      help="Generate an annotable from BLAST + idmapping + go.obo files")
    atgenerateparser.set_defaults(func=annotable_generate, scriptSuffix="annotable")
    add_mmseqs2_arguments(atgenerateparser, database=REFSEQ_PLANT_FASTA)
    atgenerateparser.add_argument("--db", dest="databaseTag",
                                  required=False,
                                  choices=["UniRef100", "UniRef90", "UniRef50", "RefSeq", "."],
                                  help="""Optionally, specify the database type to match -d;
                                  default == '%(default)s'""",
                                  default="RefSeq")
    atgenerateparser.add_argument("--idmapping", dest="idmappingFile",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the idmapping_selected.tab file;
                                  default == '%(default)s'""",
                                  default=IDMAPPING_FILE)
    atgenerateparser.add_argument("--obo", dest="oboFile",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the go.obo file;
                                  default == '%(default)s'""",
                                  default=OBO_FILE)
    atgenerateparser.add_argument("--numhits", dest="numHits",
                                  required=False,
                                  type=int,
                                  help="""Optionally, specify the number of hits for each sequence
                                  to report; default == %(default)s""",
                                  default=15)
    add_pbs_arguments(atgenerateparser, jobName="annotable",
                      walltime="06:00:00", mem="70G", cpus=12)

    # Assembly subparser
    asparser = subparsers.add_parser("assembly",
                                     parents=[p],
                                     add_help=False,
                                     help="Genome assembly jobs")
    
    subAssemblyParsers = asparser.add_subparsers(dest="assemblyMode",
                                                 required=True)
    
    # Assembly > hifiasm mode
    ashifiasmparser = subAssemblyParsers.add_parser("hifiasm",
                                                    parents=[p],
                                                    add_help=False,
                                                    help="Assemble a genome with hifiasm (optionally with Hi-C)")
    ashifiasmparser.set_defaults(func=assembly_hifiasm, scriptSuffix="hifiasm")
    ashifiasmparser.add_argument("--hic", dest="hicFiles",
                                 required=False,
                                 nargs="*",
                                 metavar="FASTQ",
                                 help="""Optionally, run a Hi-C assembly; give the two Hi-C read
                                 files (R1 R2), or no files to leave them as placeholders;
                                 default == no Hi-C""",
                                 default=None)
    ashifiasmparser.add_argument("--ploidy", dest="ploidy",
                                 required=False,
                                 type=int,
                                 action=UserSet,
                                 help="""Optionally, specify the ploidy for a Hi-C assembly (requires
                                 --hic); default == %(default)s""",
                                 default=2)
    ashifiasmparser.add_argument("--lineage", dest="lineage",
                                 required=False,
                                 help="""Optionally, specify the compleasm lineage to search;
                                 default == '%(default)s'""",
                                 default="eudicotyledons")
    ashifiasmparser.add_argument("--compleasm", dest="compleasmDir",
                                 required=False,
                                 action=UserSet,
                                 help="""Optionally, specify the directory containing compleasm.py;
                                 default == '%(default)s'""",
                                 default=COMPLEASM_DIR)
    add_genscripts_argument(ashifiasmparser)
    add_pbs_arguments(ashifiasmparser, jobName="asm",
                      walltime="24:00:00", mem="90G", cpus=14)
    
    # Assembly > revcomp mode
    asrevcompparser = subAssemblyParsers.add_parser("revcomp",
                                                    parents=[p],
                                                    add_help=False,
                                                    help="Reverse complement haplotypes to match a reference genome")
    asrevcompparser.set_defaults(func=assembly_revcomp, scriptSuffix="revcomp")
    asrevcompparser.add_argument("-r", "--reference", dest="refGenome",
                                 required=False,
                                 action=UserSet,
                                 help="""Optionally, specify the location of the reference genome FASTA file;
                                 default == '%(default)s'""",
                                 default=BUSK_HAP1_FASTA)
    add_preset_arguments(asrevcompparser)
    add_varscripts_argument(asrevcompparser)
    add_pbs_arguments(asrevcompparser, jobName="revcomp",
                      walltime="12:00:00", mem="60G", cpus=8)
    
    # Assembly > scaffolding mode
    asscaffoldparser = subAssemblyParsers.add_parser("scaffolding",
                                                     parents=[p],
                                                     add_help=False,
                                                     help="Scaffold haplotypes against a diploid reference genome")
    asscaffoldparser.set_defaults(func=assembly_scaffolding, scriptSuffix="scaffolding")
    asscaffoldparser.add_argument("--refhap1", dest="refHap1",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the reference haplotype 1
                                  FASTA file; default == '%(default)s'""",
                                  default=BUSK_HAP1_FASTA)
    asscaffoldparser.add_argument("--refhap2", dest="refHap2",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the reference haplotype 2
                                  FASTA file; default == '%(default)s'""",
                                  default=BUSK_HAP2_FASTA)
    add_preset_arguments(asscaffoldparser)
    add_varscripts_argument(asscaffoldparser)
    add_pbs_arguments(asscaffoldparser, jobName="scaffolding",
                      walltime="12:00:00", mem="70G", cpus=8)
    
    # Blast subparser
    bparser = subparsers.add_parser("blast",
                                    parents=[p],
                                    add_help=False,
                                    help="BLAST (/MMseqs2) job handling")

    subBlastParsers = bparser.add_subparsers(dest="blastMode",
                                             required=True)

    # Blast > run subparser
    brunparser = subBlastParsers.add_parser("run",
                                            parents=[p],
                                            add_help=False,
                                            help="Run a BLAST (/MMseqs2) search job")

    subBlastRunParsers = brunparser.add_subparsers(dest="blastRunMode",
                                                   required=True)

    # Blast > run > mmseqs2 mode
    brunmparser = subBlastRunParsers.add_parser("mmseqs2",
                                                parents=[p],
                                                add_help=False,
                                                help="Run a MMseqs2 job")
    brunmparser.set_defaults(func=blast_run_mmseqs2, scriptSuffix="mms2")
    add_mmseqs2_arguments(brunmparser, database=UNIREF90_FASTA)
    brunmparser.add_argument("-s", "--sensitivity", dest="sensitivity",
                             required=False,
                             type=float,
                             choices=[1, 2, 3, 4, 5, 5.7, 6, 7, 7.5],
                             help="""Optionally, specify the MMseqs2 search sensitivity;
                             default == %(default)s""",
                             default=5.7)
    brunmparser.add_argument("-t", "--searchtype", dest="searchType",
                             required=False,
                             choices=["auto", "blastn", "tblastx"],
                             help="""Optionally, specify the MMseqs2 search type;
                             default == '%(default)s'""",
                             default="auto")
    add_pbs_arguments(brunmparser, jobName="mms2",
                      walltime="46:00:00", mem="750G", cpus=32)

    args = subParentParser.parse_args()
    if not hasattr(args, "userSet"): # only set by UserSet actions when they are used
        args.userSet = set()

    # Split into mode-specific validation functions
    if args.mode == "annotable":
        print("## psst.py - Annotation table handling ##")
        validate_at(args)
        atmain(args)
    if args.mode == "assembly":
        print("## psst.py - Genome assembly handling ##")
        validate_as(args)
        asmain(args)
    if args.mode == "blast":
        print("## psst.py - BLAST handling ##")
        validate_b(args)
        bmain(args)

    # Validate PBS arguments and write (+/- submit) the job script
    validate_pbs(args)
    args.func(args)

    # Print completion flag if we reach this point
    print("Program completed successfully!")

def atmain(args):
    # Split into sub-mode-specific validation functions
    if args.annotableMode == "generate":
        print("## Annotable generation ##")
        validate_at_generate(args)

def asmain(args):
    # Split into sub-mode-specific validation functions
    if args.assemblyMode == "hifiasm":
        print("## hifiasm genome assembly ##")
        validate_as_hifiasm(args)
    if args.assemblyMode == "revcomp":
        print("## Reverse complement haplotypes to match reference ##")
        validate_as_revcomp(args)
    if args.assemblyMode == "scaffolding":
        print("## Haplotype scaffolding ##")
        validate_as_scaffolding(args)

def bmain(args):
    # Split into sub-mode-specific validation functions
    if args.blastMode == "run":
        validate_b_run(args)
        if args.blastRunMode == "mmseqs2":
            print("## MMseqs2 query ##")
            validate_b_run_mmseqs2(args)

if __name__ == "__main__":
    main()
