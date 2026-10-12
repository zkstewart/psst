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

sys.path.append(os.path.dirname(os.path.realpath(__file__))) # realpath allows running via a symlink
from modules.validation import \
    validate_as, \
        validate_as_hifiasm, validate_as_revcomp, validate_as_scaffolding, \
    validate_at, \
        validate_at_generate, \
    validate_b, \
        validate_b_run, \
            validate_b_run_mmseqs2, \
    validate_bc, \
        validate_bc_dorado, \
    validate_bs, \
        validate_bs_download, validate_bs_prepare, validate_bs_archive, \
    validate_pbs_submission, validate_pbs_job, RPID_SHORTCUTS
from modules.annotable import annotable_generate
from modules.assembly import assembly_hifiasm, assembly_revcomp, assembly_scaffolding
from modules.basecalling import basecalling_dorado
from modules.basespace import basespace_download, basespace_prepare, basespace_archive
from modules.blast import blast_run_mmseqs2
from modules.locations import MMSEQS_DIR, VARSCRIPT_DIR, GENSCRIPT_DIR, COMPLEASM_DIR, DORADO_DIR, BASESPACE_DIR, \
    REFSEQ_PLANT_FASTA, UNIREF90_FASTA, IDMAPPING_FILE, OBO_FILE, \
    REFERENCE_HAP1_FASTA, REFERENCE_HAP2_FASTA
from _version import __version__

# Argument conventions:
# - All default values are set in the parser; help text should reference them with
#   %(default)s so the help can never disagree with the actual default.
# - Arguments whose values are checked by validation (e.g., file or program locations)
#   should use action=UserSet. Validation only checks values the user explicitly gave,
#   so default locations on the HPC are never checked.
# - Arguments without a sensible default should default to None, which leaves an
#   '@@field@@' placeholder in the job script for the user to fill in.
# - Each leaf parser must set_defaults(func=<job function>). Single-job modes use
#   add_pbs_arguments() and must also set_defaults(scriptSuffix=<str>); multi-job modes
#   use add_pbs_submission_arguments() with add_pbs_resource_arguments() per job, and
#   validate their own job resources and script locations. Modes which do not write
#   job scripts (e.g., wrappers that run a script immediately) use neither.

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

def add_pbs_resource_arguments(parser, walltime, mem, cpus, step=None, title=None):
    '''
    Adds PBS resource arguments for one job with job-specific default values.
    
    Parameters:
        parser -- an argparse.ArgumentParser (or argument group) to add arguments to.
        walltime -- a string of the default walltime in 'HH:MM:SS' format.
        mem -- a string of the default memory request (e.g., '70G').
        cpus -- an int of the default number of CPUs to request.
        step -- OPTIONAL; a string naming the job within a multi-job mode (e.g., 'p1'),
                giving arguments like '--p1-cpus' with dest 'p1Cpus'; default == None,
                giving arguments like '--cpus' with dest 'cpus'.
        title -- OPTIONAL; a string describing the job for an argument group heading;
                 if given, arguments are added to a new group within parser.
    '''
    if title != None:
        parser = parser.add_argument_group(title)
    flag = (lambda name: f"--{step}-{name}") if step != None else (lambda name: f"--{name}")
    dest = (lambda name: f"{step}{name.capitalize()}") if step != None else (lambda name: name)
    
    parser.add_argument(flag("cpus"), dest=dest("cpus"),
                        required=False,
                        type=int,
                        help="""Optionally, specify the number of CPUs to request;
                        default == %(default)s""",
                        default=cpus)
    parser.add_argument(flag("mem"), dest=dest("mem"),
                        required=False,
                        help="""Optionally, specify the memory to request;
                        default == '%(default)s'""",
                        default=mem)
    parser.add_argument(flag("walltime"), dest=dest("walltime"),
                        required=False,
                        help="""Optionally, specify the walltime to request in 'H', 'H:MM',
                        or 'H:MM:SS' format (max 48 hours); default == '%(default)s'""",
                        default=walltime)

def add_pbs_submission_arguments(parser):
    '''
    Adds PBS submission arguments common to all modes to a leaf subparser.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.set_defaults(pbsMode="multi") # add_pbs_arguments() overrides this
    g = parser.add_argument_group("PBS submission arguments")
    shortcuts = ", ".join(f"{k}={v}" for k, v in RPID_SHORTCUTS.items())
    g.add_argument("--rpid", dest="rpid",
                   required=True,
                   help=f"""Research Project ID (RPID) to attach to the job via '#PBS -P';
                   accepts a full RPID or a shortcut ({shortcuts})""")
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

def add_pbs_arguments(parser, jobName, walltime, mem, cpus):
    '''
    Adds PBS job and submission arguments to the leaf subparser of a single-job mode,
    with job-specific default values.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
        jobName -- a string of the default PBS job name.
        walltime -- a string of the default walltime in 'HH:MM:SS' format.
        mem -- a string of the default memory request (e.g., '70G').
        cpus -- an int of the default number of CPUs to request.
    '''
    g = parser.add_argument_group("PBS job arguments")
    add_pbs_resource_arguments(g, walltime, mem, cpus)
    g.add_argument("--name", dest="jobName",
                   required=False,
                   help="""Optionally, specify the PBS job name;
                   default == '%(default)s'""",
                   default=jobName)
    g.add_argument("--script", dest="scriptFile",
                   required=False,
                   help="""Optionally, specify where to write the job script;
                   default == '<outprefix>.<job>.sh' if -o is given, otherwise
                   '<job>.sh' in the current directory""",
                   default=None)
    add_pbs_submission_arguments(parser)
    parser.set_defaults(pbsMode="single")

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

    subParentParser.set_defaults(pbsMode=None) # overridden by leaf parsers with PBS arguments
    
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
                                 default=REFERENCE_HAP1_FASTA)
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
                                  default=REFERENCE_HAP1_FASTA)
    asscaffoldparser.add_argument("--refhap2", dest="refHap2",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the reference haplotype 2
                                  FASTA file; default == '%(default)s'""",
                                  default=REFERENCE_HAP2_FASTA)
    add_preset_arguments(asscaffoldparser)
    add_varscripts_argument(asscaffoldparser)
    add_pbs_arguments(asscaffoldparser, jobName="scaffolding",
                      walltime="12:00:00", mem="70G", cpus=8)
    
    # Basecalling subparser
    bcparser = subparsers.add_parser("basecalling",
                                     parents=[p],
                                     add_help=False,
                                     help="Nanopore basecalling jobs")
    
    subBasecallingParsers = bcparser.add_subparsers(dest="basecallingMode",
                                                    required=True)
    
    # Basecalling > dorado mode
    bcdoradoparser = subBasecallingParsers.add_parser("dorado",
                                                      parents=[p],
                                                      add_help=False,
                                                      help="""Basecall and correct nanopore reads with dorado
                                                      as a chain of jobs""")
    bcdoradoparser.set_defaults(func=basecalling_dorado)
    bcdoradoparser.add_argument("-p", "--prefix", dest="prefix",
                                required=False,
                                help="""Optionally, specify the prefix for output files (e.g., a sample ID);
                                default == left as a placeholder""",
                                default=None)
    bcdoradoparser.add_argument("--workdir", dest="workDir",
                                required=False,
                                help="""Optionally, specify the run directory containing the pod5
                                directory, where job scripts will be laid out; default == '%(default)s'""",
                                default=".")
    bcdoradoparser.add_argument("--model", dest="model",
                                required=False,
                                help="""Optionally, specify the dorado basecalling model;
                                default == '%(default)s'""",
                                default="sup")
    bcdoradoparser.add_argument("--minreadlen", dest="minReadLen",
                                required=False,
                                type=int,
                                help="""Optionally, specify the minimum read length to keep before
                                correction; default == %(default)s""",
                                default=5000)
    bcdoradoparser.add_argument("--concat-minimum", dest="concatMinimum",
                                required=False,
                                type=int,
                                help="""Optionally, specify the minimum length of corrected reads
                                to keep when concatenating; default == %(default)s""",
                                default=5000)
    bcdoradoparser.add_argument("--concat-multiline", dest="concatMultiline",
                                required=False,
                                type=int,
                                help="""Optionally, specify the line length of the concatenated
                                corrected reads FASTA; default == %(default)s""",
                                default=60)
    bcdoradoparser.add_argument("--gpuid", dest="gpuID",
                                required=False,
                                help="""Optionally, specify the GPU type to request for basecalling;
                                default == '%(default)s'""",
                                default="H100")
    bcdoradoparser.add_argument("--name", dest="jobName",
                                required=False,
                                help="""Optionally, specify a prefix for job names (e.g., '<name>_p1');
                                default == the value of -p, or 'dorado' if -p is not given""",
                                default=None)
    bcdoradoparser.add_argument("--dorado", dest="doradoDir",
                                required=False,
                                action=UserSet,
                                help="""Optionally, specify the dorado directory containing bin/ and lib/;
                                default == '%(default)s'""",
                                default=DORADO_DIR)
    add_varscripts_argument(bcdoradoparser)
    add_genscripts_argument(bcdoradoparser)
    add_pbs_resource_arguments(bcdoradoparser, step="p1", title="PBS resources: basecaller part 1 (GPU basecalling)",
                               walltime="32:00:00", mem="35G", cpus=1)
    add_pbs_resource_arguments(bcdoradoparser, step="p2", title="PBS resources: basecaller part 2 (BAM to FASTQ)",
                               walltime="16:00:00", mem="60G", cpus=1)
    add_pbs_resource_arguments(bcdoradoparser, step="corr", title="PBS resources: correction job submission",
                               walltime="00:10:00", mem="10G", cpus=1)
    add_pbs_resource_arguments(bcdoradoparser, step="s1", title="PBS resources: correction step 1 (CPU array)",
                               walltime="32:00:00", mem="380G", cpus=22)
    add_pbs_resource_arguments(bcdoradoparser, step="s2", title="PBS resources: correction step 2 (GPU array)",
                               walltime="24:00:00", mem="30G", cpus=4)
    add_pbs_resource_arguments(bcdoradoparser, step="concat", title="PBS resources: corrected read concatenation",
                               walltime="03:00:00", mem="30G", cpus=1)
    add_pbs_submission_arguments(bcdoradoparser)
    
    # Basespace subparser
    bsparser = subparsers.add_parser("basespace",
                                     parents=[p],
                                     add_help=False,
                                     help="Illumina BaseSpace read handling")
    
    subBasespaceParsers = bsparser.add_subparsers(dest="basespaceMode",
                                                  required=True)
    
    # Basespace > download mode
    bsdownloadparser = subBasespaceParsers.add_parser("download",
                                                      parents=[p],
                                                      add_help=False,
                                                      help="Download a BaseSpace project")
    bsdownloadparser.set_defaults(func=basespace_download, scriptSuffix="bs_download")
    bsdownloadparser.add_argument("-i", "--project", dest="projectID",
                                  required=False,
                                  help="""Optionally, specify the BaseSpace project ID to download (find
                                  this with 'bs projects list'); default == left as a placeholder""",
                                  default=None)
    bsdownloadparser.add_argument("--outdir", dest="outputDir",
                                  required=False,
                                  help="""Optionally, specify the directory to download into, relative to
                                  the directory the job is submitted from; default == '%(default)s'""",
                                  default="downloaded")
    bsdownloadparser.add_argument("--bs", dest="bsDir",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the directory containing the bs executable;
                                  default == '%(default)s'""",
                                  default=BASESPACE_DIR)
    add_pbs_arguments(bsdownloadparser, jobName="bs_download",
                      walltime="48:00:00", mem="35G", cpus=1)
    
    # Basespace > prepare mode
    bsprepareparser = subBasespaceParsers.add_parser("prepare",
                                                     parents=[p],
                                                     add_help=False,
                                                     help="""Run concatenate_trimmed_lanes.py to write a script that
                                                     merges reads split across lanes (runs immediately; not a PBS job)""")
    bsprepareparser.set_defaults(func=basespace_prepare)
    bsprepareparser.add_argument("-d", "--readsdir", dest="readsDir",
                                 required=False,
                                 help="""Optionally, specify the directory containing the downloaded reads;
                                 default == '%(default)s'""",
                                 default="reads")
    bsprepareparser.add_argument("-rs", "--suffix", dest="readsSuffix",
                                 required=False,
                                 help="""Optionally, specify the suffix which identifies read files, where the
                                 character before it is 1 or 2 for paired reads; default == '%(default)s'""",
                                 default="_001.fastq.gz")
    bsprepareparser.add_argument("-o", "--output", dest="outputFileName",
                                 required=False,
                                 help="""Optionally, specify where to write the read preparation script;
                                 its directory is created if needed; default == '%(default)s'""",
                                 default=os.path.join("prepared_reads", "run_read_prep.sh"))
    bsprepareparser.add_argument("--symlink", dest="symlink",
                                 required=False,
                                 action="store_true",
                                 help="""Optionally, symlink reads from samples sequenced on a single lane
                                 rather than copying them""",
                                 default=False)
    bsprepareparser.add_argument("--singleEnd", dest="singleEnd",
                                 required=False,
                                 action="store_true",
                                 help="""Optionally, indicate that reads are single-end rather than paired""",
                                 default=False)
    bsprepareparser.add_argument("--laneIdentifier", dest="laneIdentifier",
                                 required=False,
                                 help="""Optionally, specify the text preceding the lane number in file
                                 names; default == '%(default)s'""",
                                 default="_L00")
    add_varscripts_argument(bsprepareparser)
    
    # Basespace > archive mode
    bsarchiveparser = subBasespaceParsers.add_parser("archive",
                                                     parents=[p],
                                                     add_help=False,
                                                     help="Archive downloaded reads into a .tar file with an md5 checksum")
    bsarchiveparser.set_defaults(func=basespace_archive, scriptSuffix="archive")
    bsarchiveparser.add_argument("-p", "--prefix", dest="prefix",
                                 required=False,
                                 help="""Optionally, specify the prefix of the .tar archive file;
                                 default == left as a placeholder""",
                                 default=None)
    bsarchiveparser.add_argument("--readsdir", dest="readsDir",
                                 required=False,
                                 help="""Optionally, specify the directory to archive, relative to the
                                 directory the job is submitted from; default == '%(default)s'""",
                                 default="downloaded")
    add_pbs_arguments(bsarchiveparser, jobName="tar",
                      walltime="12:00:00", mem="30G", cpus=1)
    
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
    if args.mode == "basecalling":
        print("## psst.py - Basecalling handling ##")
        validate_bc(args)
        bcmain(args)
    if args.mode == "basespace":
        print("## psst.py - BaseSpace handling ##")
        validate_bs(args)
        bsmain(args)
    if args.mode == "blast":
        print("## psst.py - BLAST handling ##")
        validate_b(args)
        bmain(args)

    # Validate PBS arguments and write (+/- submit) the job script(s)
    if args.pbsMode != None:
        validate_pbs_submission(args)
    if args.pbsMode == "single":
        validate_pbs_job(args)
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

def bcmain(args):
    # Split into sub-mode-specific validation functions
    if args.basecallingMode == "dorado":
        print("## dorado basecalling and read correction ##")
        validate_bc_dorado(args)

def bsmain(args):
    # Split into sub-mode-specific validation functions
    if args.basespaceMode == "download":
        print("## BaseSpace project download ##")
        validate_bs_download(args)
    if args.basespaceMode == "prepare":
        print("## Prepare reads split across lanes ##")
        validate_bs_prepare(args)
    if args.basespaceMode == "archive":
        print("## Archive downloaded reads ##")
        validate_bs_archive(args)

def bmain(args):
    # Split into sub-mode-specific validation functions
    if args.blastMode == "run":
        validate_b_run(args)
        if args.blastRunMode == "mmseqs2":
            print("## MMseqs2 query ##")
            validate_b_run_mmseqs2(args)

if __name__ == "__main__":
    main()
