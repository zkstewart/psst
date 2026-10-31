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
    validate_cl, \
        validate_cl_mmseqs2, validate_cl_cdhit, \
    validate_m, \
        validate_m_wgs, validate_m_hic, \
        validate_m_index, \
            validate_m_index_bwa, validate_m_index_salmon, validate_m_index_star, \
        validate_m_rna, \
            validate_m_rna_salmon, \
    validate_pbs_submission, validate_pbs_job, RPID_SHORTCUTS
from modules.annotable import annotable_generate
from modules.assembly import assembly_hifiasm, assembly_revcomp, assembly_scaffolding
from modules.basecalling import basecalling_dorado
from modules.basespace import basespace_download, basespace_prepare, basespace_archive
from modules.blast import blast_run_mmseqs2
from modules.clustering import clustering_mmseqs2, clustering_cdhit
from modules.mapping import map_wgs, map_hic, map_index_bwa, map_index_salmon, map_index_star, map_rna_salmon
from modules.locations import MMSEQS_DIR, VARSCRIPT_DIR, GENSCRIPT_DIR, COMPLEASM_DIR, DORADO_DIR, BASESPACE_DIR, CDHIT_DIR, PAF2DOTPLOT_DIR, \
    THREEDDNA_DIR, PAIRTOOLS_CONDA_ENV, SALMON_DIR, STAR_MODULE, \
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
# - Argument dests must not reuse names used by psst itself to route arguments (e.g.,
#   'mode' and the other subparser dests, 'func', 'scriptSuffix', 'pbsMode', 'userSet').
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

def check_reserved_dests(parser):
    '''
    Raises an error if any argument's dest reuses a name psst relies on to route
    arguments (i.e., subparser dests and names set via set_defaults()), since argparse
    would otherwise silently overwrite the routing value.
    
    Parameters:
        parser -- the top-level argparse.ArgumentParser.
    '''
    reserved = {"func", "scriptSuffix", "pbsMode", "userSet"}
    def _walk(parser):
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                reserved.add(action.dest)
                for subparser in action.choices.values():
                    yield from _walk(subparser)
            else:
                yield parser.prog, action.dest
    argDests = list(_walk(parser)) # also populates reserved
    for prog, dest in argDests:
        if dest in reserved:
            raise ValueError(f"'{prog}' has an argument with dest '{dest}', which is reserved by psst")

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

def add_rpid_argument(parser):
    '''
    Adds the required RPID argument; used by add_pbs_submission_arguments() and by
    any mode which passes an RPID to a script that writes its own job scripts.
    
    Parameters:
        parser -- an argparse.ArgumentParser (or argument group) to add arguments to.
    '''
    shortcuts = ", ".join(f"{k}={v}" for k, v in RPID_SHORTCUTS.items())
    parser.add_argument("--rpid", dest="rpid",
                        required=True,
                        help=f"""Research Project ID (RPID) to attach to the job via '#PBS -P';
                        accepts a full RPID or a shortcut ({shortcuts})""")

def add_pbs_submission_arguments(parser):
    '''
    Adds PBS submission arguments common to all modes to a leaf subparser.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.set_defaults(pbsMode="multi") # add_pbs_arguments() overrides this
    g = parser.add_argument_group("PBS submission arguments")
    add_rpid_argument(g)
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
    add_mmseqs_dir_argument(parser)
    add_varscripts_argument(parser)

def add_mmseqs_dir_argument(parser):
    '''
    Adds the MMseqs2 location argument to a leaf subparser for any job that runs MMseqs2.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
    '''
    parser.add_argument("--mmseqs", dest="mmseqsDir",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the directory containing the mmseqs
                        executable; default == '%(default)s'""",
                        default=MMSEQS_DIR)

def add_clustering_arguments(parser, identity):
    '''
    Adds arguments shared by clustering jobs which run Various_scripts/run_clustering.py.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
        identity -- a float of the default identity threshold for clustering.
    '''
    parser.add_argument("-i", "--fasta", dest="fastaFile",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the location of the FASTA file to cluster;
                        default == left as a placeholder""",
                        default=None)
    parser.add_argument("-m", "--molecule", dest="molecule",
                        required=False,
                        choices=["protein", "nucleotide"],
                        help="""Optionally, specify the molecule type of the FASTA file;
                        default == left as a placeholder""",
                        default=None)
    parser.add_argument("-o", "--output", dest="outputFile",
                        required=False,
                        help="""Optionally, specify the output file name for clustering results;
                        default == left as a placeholder""",
                        default=None)
    parser.add_argument("--identity", dest="identity",
                        required=False,
                        type=float,
                        help="""Optionally, specify the identity threshold for clustering;
                        default == %(default)s""",
                        default=identity)
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

def add_bwa_arguments(parser, allowSingleEnd):
    '''
    Adds reference genome and reads arguments to a leaf subparser for "map" modes which
    run bwa mem.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
        allowSingleEnd -- a boolean indicating whether the mode accepts single-end reads.
    '''
    parser.add_argument("-r", "--reference", dest="genome",
                        required=False,
                        action=UserSet,
                        help="""Optionally, specify the location of the reference genome FASTA file,
                        which must be indexed with 'bwa index'; default == left as a placeholder""",
                        default=None)
    add_reads_arguments(parser, allowSingleEnd)

def add_reads_arguments(parser, allowSingleEnd):
    '''
    Adds reads arguments to a leaf subparser for "map" modes.
    
    Parameters:
        parser -- an argparse.ArgumentParser to add arguments to.
        allowSingleEnd -- a boolean indicating whether the mode accepts single-end reads.
    '''
    parser.add_argument("-d", "--readsdir", dest="readsDir",
                        required=False,
                        help="""Optionally, specify the directory containing the read files; if given,
                        read files are located and pairs checked; default == left as a placeholder""",
                        default=None)
    parser.add_argument("-rs", "--suffix", dest="readsSuffix",
                        required=False,
                        help="""Optionally, specify the suffix which identifies read files, where the
                        character before it is 1 or 2 for paired reads; default == '%(default)s'""",
                        default=".fq.gz")
    if allowSingleEnd:
        parser.add_argument("--singleEnd", dest="singleEnd",
                            required=False,
                            action="store_true",
                            help="""Optionally, indicate that reads are single-end rather than paired""",
                            default=False)

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
    asscaffoldparser.add_argument("--paf2dotplot", dest="paf2dotplotDir",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the directory containing paf2dotplot.r, used to
                                  make dotplots of the scaffolded haplotypes; default == '%(default)s'""",
                                  default=PAF2DOTPLOT_DIR)
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
    add_rpid_argument(bsprepareparser)
    
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

    # Clustering subparser
    clparser = subparsers.add_parser("clustering",
                                     parents=[p],
                                     add_help=False,
                                     help="Sequence clustering jobs")
    
    subClusteringParsers = clparser.add_subparsers(dest="clusteringMode",
                                                   required=True)
    
    # Clustering > mmseqs2 mode
    clmmseqsparser = subClusteringParsers.add_parser("mmseqs2",
                                                     parents=[p],
                                                     add_help=False,
                                                     help="Cluster sequences with MMseqs2")
    clmmseqsparser.set_defaults(func=clustering_mmseqs2, scriptSuffix="mms2_clust")
    add_clustering_arguments(clmmseqsparser, identity=0.85)
    clmmseqsparser.add_argument("-p", "--program", dest="program",
                                required=False,
                                choices=["mmseqs-cascade", "mmseqs-linclust"],
                                help="""Optionally, specify the MMseqs2 clustering algorithm;
                                default == '%(default)s'""",
                                default="mmseqs-cascade")
    clmmseqsparser.add_argument("--coverage", dest="coverage",
                                required=False,
                                type=float,
                                help="""Optionally, specify the coverage ratio for clustering;
                                default == %(default)s""",
                                default=0.6)
    clmmseqsparser.add_argument("--evalue", dest="evalue",
                                required=False,
                                type=float,
                                help="""Optionally, specify the E-value threshold for clustering;
                                default == %(default)s""",
                                default=1e-3)
    clmmseqsparser.add_argument("--mode", dest="clusterMode",
                                required=False,
                                choices=["set-cover", "connected-component", "greedy"],
                                help="""Optionally, specify the clustering mode;
                                default == '%(default)s'""",
                                default="connected-component")
    clmmseqsparser.add_argument("--sensitivity", dest="sensitivity",
                                required=False,
                                choices=["1", "2", "3", "4", "5", "5.7", "6", "7", "7.5"],
                                help="""Optionally, specify the sensitivity (cascade only);
                                default == %(default)s""",
                                default="4")
    clmmseqsparser.add_argument("--steps", dest="steps",
                                required=False,
                                type=int,
                                help="""Optionally, specify the number of cascaded clustering steps
                                (cascade only); default == %(default)s""",
                                default=3)
    clmmseqsparser.add_argument("--tmp", dest="tmpDir",
                                required=False,
                                help="""Optionally, specify the MMseqs2 tmp directory;
                                default == '<output without extension>_tmp' if -o is given,
                                otherwise left as a placeholder""",
                                default=None)
    add_mmseqs_dir_argument(clmmseqsparser)
    add_pbs_arguments(clmmseqsparser, jobName="mms2Clust",
                      walltime="48:00:00", mem="80G", cpus=12)
    
    # Clustering > cdhit mode
    clcdhitparser = subClusteringParsers.add_parser("cdhit",
                                                    parents=[p],
                                                    add_help=False,
                                                    help="Cluster sequences with CD-HIT")
    clcdhitparser.set_defaults(func=clustering_cdhit, scriptSuffix="cdhit_clust")
    add_clustering_arguments(clcdhitparser, identity=0.97)
    clcdhitparser.add_argument("--cdhit-mem", dest="cdhitMem",
                               required=False,
                               type=int,
                               help="""Optionally, specify the memory (in MB) for CD-HIT to use, which
                               should be less than --mem; default == %(default)s""",
                               default=25000)
    clcdhitparser.add_argument("--align", dest="align",
                               required=False,
                               choices=["local", "global"],
                               help="""Optionally, specify whether CD-HIT should use local or global
                               alignment; default == '%(default)s'""",
                               default="global")
    clcdhitparser.add_argument("--shorter-coverage", dest="shorterCov",
                               required=False,
                               type=float,
                               help="""Optionally, specify the alignment coverage ratio for the
                               shorter sequence; default == %(default)s""",
                               default=0.9)
    clcdhitparser.add_argument("--longer-coverage", dest="longerCov",
                               required=False,
                               type=float,
                               help="""Optionally, specify the alignment coverage ratio for the
                               longer sequence; default == %(default)s""",
                               default=0.0)
    clcdhitparser.add_argument("--cdhit", dest="cdhitDir",
                               required=False,
                               action=UserSet,
                               help="""Optionally, specify the directory containing the CD-HIT
                               executables; default == '%(default)s'""",
                               default=CDHIT_DIR)
    add_pbs_arguments(clcdhitparser, jobName="cdhit",
                      walltime="48:00:00", mem="30G", cpus=8)
    
    # Map subparser
    mparser = subparsers.add_parser("map",
                                    parents=[p],
                                    add_help=False,
                                    help="Read mapping jobs")
    
    subMapParsers = mparser.add_subparsers(dest="mapMode",
                                           required=True)
    
    # Map > wgs mode
    mwgsparser = subMapParsers.add_parser("wgs",
                                          parents=[p],
                                          add_help=False,
                                          help="""Map Illumina DNA reads with bwa mem, then sort and mark
                                          duplicates (an array job with one subjob per sample)""")
    mwgsparser.set_defaults(func=map_wgs, scriptSuffix="map_wgs")
    add_bwa_arguments(mwgsparser, allowSingleEnd=True)
    mwgsparser.add_argument("--outdir", dest="outputDir",
                            required=False,
                            help="""Optionally, specify the output directory, relative to the
                            directory the job is submitted from; default == '%(default)s'""",
                            default="mapped")
    mwgsparser.add_argument("--library", dest="library",
                            required=False,
                            help="""Optionally, specify the library name (LB) for read groups;
                            default == '%(default)s'""",
                            default="lib1")
    mwgsparser.add_argument("--sort-mem", dest="sortMem",
                            required=False,
                            help="""Optionally, specify the memory per samtools sort thread, which
                            multiplied by --cpus should be well below --mem; default == '%(default)s'""",
                            default="3G")
    add_pbs_arguments(mwgsparser, jobName="map_wgs",
                      walltime="24:00:00", mem="55G", cpus=8)
    
    # Map > hic mode
    mhicparser = subMapParsers.add_parser("hic",
                                          parents=[p],
                                          add_help=False,
                                          help="""Map Hi-C reads with bwa mem, then deduplicate contacts with
                                          pairtools and build a .hic file with 3D-DNA""")
    mhicparser.set_defaults(func=map_hic, scriptSuffix="map_hic")
    add_bwa_arguments(mhicparser, allowSingleEnd=False)
    mhicparser.add_argument("-p", "--prefix", dest="prefix",
                            required=False,
                            help="""Optionally, specify the prefix for output files; default == the
                            reference genome file name without its extension if -r is given,
                            otherwise left as a placeholder""",
                            default=None)
    mhicparser.add_argument("--outdir", dest="outputDir",
                            required=False,
                            help="""Optionally, specify the output directory, relative to the
                            directory the job is submitted from; default == '%(default)s'""",
                            default="hic")
    mhicparser.add_argument("--mapq", dest="mapq",
                            required=False,
                            type=int,
                            help="""Optionally, specify the minimum MAPQ of both reads for a contact
                            to be included in the .hic file; default == %(default)s""",
                            default=1)
    mhicparser.add_argument("--sort-mem", dest="sortMem",
                            required=False,
                            help="""Optionally, specify the memory for pairtools sort, which should
                            be well below --mem; default == '%(default)s'""",
                            default="20G")
    mhicparser.add_argument("--conda", dest="condaEnv",
                            required=False,
                            help="""Optionally, specify the conda environment containing pairtools;
                            default == '%(default)s'""",
                            default=PAIRTOOLS_CONDA_ENV)
    mhicparser.add_argument("--3ddna", dest="threeddnaDir",
                            required=False,
                            action=UserSet,
                            help="""Optionally, specify the 3D-DNA directory; default == '%(default)s'""",
                            default=THREEDDNA_DIR)
    add_pbs_arguments(mhicparser, jobName="map_hic",
                      walltime="48:00:00", mem="80G", cpus=16)
    
    # Map > index subparser
    mindexparser = subMapParsers.add_parser("index",
                                            parents=[p],
                                            add_help=False,
                                            help="Index a reference for read mapping")
    
    subMapIndexParsers = mindexparser.add_subparsers(dest="mapIndexMode",
                                                     required=True)
    
    # Map > index > bwa mode
    mindexbwaparser = subMapIndexParsers.add_parser("bwa",
                                                    parents=[p],
                                                    add_help=False,
                                                    help="Index a reference genome FASTA file with bwa index")
    mindexbwaparser.set_defaults(func=map_index_bwa, scriptSuffix="bwa_index")
    mindexbwaparser.add_argument("-r", "--reference", dest="genome",
                                 required=False,
                                 action=UserSet,
                                 help="""Optionally, specify the location of the reference genome FASTA file
                                 to index; index files are written alongside it; default == left as a
                                 placeholder""",
                                 default=None)
    add_pbs_arguments(mindexbwaparser, jobName="bwa_index",
                      walltime="03:00:00", mem="30G", cpus=1)
    
    # Map > index > salmon mode
    mindexsalmonparser = subMapIndexParsers.add_parser("salmon",
                                                       parents=[p],
                                                       add_help=False,
                                                       help="Index a transcripts FASTA file with salmon index")
    mindexsalmonparser.set_defaults(func=map_index_salmon, scriptSuffix="salmon_index")
    mindexsalmonparser.add_argument("-t", "--transcripts", dest="transcripts",
                                    required=False,
                                    action=UserSet,
                                    help="""Optionally, specify the location of the transcripts FASTA file
                                    to index; default == left as a placeholder""",
                                    default=None)
    mindexsalmonparser.add_argument("--genome", dest="genome",
                                    required=False,
                                    action=UserSet,
                                    help="""Optionally, specify the genome FASTA file to build a decoy-aware
                                    index, which reduces false mappings from unannotated genomic sequence
                                    (recommended); default == no decoys""",
                                    default=None)
    mindexsalmonparser.add_argument("-o", "--index", dest="indexDir",
                                    required=False,
                                    help="""Optionally, specify the index directory to create; default ==
                                    '<transcripts>_salmon_index' if -t is given, otherwise left as a
                                    placeholder""",
                                    default=None)
    mindexsalmonparser.add_argument("-k", "--kmer", dest="kmer",
                                    required=False,
                                    type=int,
                                    help="""Optionally, specify the k-mer length, which must be odd; consider
                                    lowering it for reads shorter than 75 bp; default == %(default)s""",
                                    default=31)
    mindexsalmonparser.add_argument("--ram-limit", dest="ramLimit",
                                    required=False,
                                    type=int,
                                    help="""Optionally, specify the RAM limit (in GiB) for salmon's index
                                    building, which must be less than --mem; default == 75%% of --mem""",
                                    default=None)
    mindexsalmonparser.add_argument("--salmon", dest="salmonDir",
                                    required=False,
                                    action=UserSet,
                                    help="""Optionally, specify the directory containing the salmon executable;
                                    default == '%(default)s'""",
                                    default=SALMON_DIR)
    add_pbs_arguments(mindexsalmonparser, jobName="salmon_index",
                      walltime="24:00:00", mem="65G", cpus=8)
    
    # Map > index > star mode
    mindexstarparser = subMapIndexParsers.add_parser("star",
                                                     parents=[p],
                                                     add_help=False,
                                                     help="Index a reference genome FASTA and GTF annotation with STAR")
    mindexstarparser.set_defaults(func=map_index_star, scriptSuffix="star_index")
    mindexstarparser.add_argument("-r", "--reference", dest="genome",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the (uncompressed) reference
                                  genome FASTA file to index; default == left as a placeholder""",
                                  default=None)
    mindexstarparser.add_argument("--gtf", dest="gtf",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the location of the (uncompressed) gene
                                  annotation GTF file; convert GFF3 files to GTF first; default == left
                                  as a placeholder""",
                                  default=None)
    mindexstarparser.add_argument("-o", "--index", dest="indexDir",
                                  required=False,
                                  help="""Optionally, specify the index directory to create; default ==
                                  '<reference>_STAR_index' if -r is given, otherwise left as a
                                  placeholder""",
                                  default=None)
    mindexstarparser.add_argument("--overhang", dest="overhang",
                                  required=False,
                                  type=int,
                                  help="""Optionally, specify the splice junction overhang, which should
                                  ideally be the read length minus 1; default == %(default)s""",
                                  default=100)
    mindexstarparser.add_argument("--module", dest="starModule",
                                  required=False,
                                  help="""Optionally, specify the STAR module to load;
                                  default == '%(default)s'""",
                                  default=STAR_MODULE)
    add_pbs_arguments(mindexstarparser, jobName="star_index",
                      walltime="04:00:00", mem="50G", cpus=8)
    
    # Map > rna subparser
    mrnaparser = subMapParsers.add_parser("rna",
                                          parents=[p],
                                          add_help=False,
                                          help="RNAseq read mapping and quantification")
    
    subMapRnaParsers = mrnaparser.add_subparsers(dest="mapRnaMode",
                                                 required=True)
    
    # Map > rna > salmon mode
    mrnasalmonparser = subMapRnaParsers.add_parser("salmon",
                                                   parents=[p],
                                                   add_help=False,
                                                   help="""Quantify RNAseq reads against a transcriptome with salmon
                                                   (an array job with one subjob per sample)""")
    mrnasalmonparser.set_defaults(func=map_rna_salmon, scriptSuffix="salmon_quant")
    mrnasalmonparser.add_argument("-i", "--index", dest="indexDir",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the salmon index directory (see 'map index
                                  salmon'); default == left as a placeholder""",
                                  default=None)
    add_reads_arguments(mrnasalmonparser, allowSingleEnd=True)
    mrnasalmonparser.add_argument("--outdir", dest="outputDir",
                                  required=False,
                                  help="""Optionally, specify the output directory, relative to the directory
                                  the job is submitted from, which will contain a subdirectory for each
                                  sample; default == '%(default)s'""",
                                  default="quant")
    mrnasalmonparser.add_argument("--libtype", dest="libType",
                                  required=False,
                                  help="""Optionally, specify the salmon library type (e.g., ISR for paired
                                  stranded reads), or A to detect it automatically; default == '%(default)s'""",
                                  default="A")
    mrnasalmonparser.add_argument("-g", "--genemap", dest="geneMap",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify a GTF/GFF annotation or transcript-to-gene TSV
                                  to also produce gene-level estimates (quant.genes.sf); default == none""",
                                  default=None)
    mrnasalmonparser.add_argument("--bias", dest="bias",
                                  required=False,
                                  choices=["none", "gc", "seq", "both"],
                                  help="""Optionally, specify which bias corrections to apply (fragment GC
                                  and/or sequence-specific); default == 'gc' for paired reads, otherwise
                                  'none'""",
                                  default=None)
    mrnasalmonparser.add_argument("--bootstraps", dest="bootstraps",
                                  required=False,
                                  type=int,
                                  help="""Optionally, specify the number of bootstrap replicates, which
                                  tools like sleuth or fishpond use to model quantification
                                  uncertainty; default == %(default)s""",
                                  default=0)
    mrnasalmonparser.add_argument("--salmon", dest="salmonDir",
                                  required=False,
                                  action=UserSet,
                                  help="""Optionally, specify the directory containing the salmon executable;
                                  default == '%(default)s'""",
                                  default=SALMON_DIR)
    add_pbs_arguments(mrnasalmonparser, jobName="salmon_quant",
                      walltime="04:00:00", mem="15G", cpus=5)
    
    check_reserved_dests(subParentParser)
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

    if args.mode == "clustering":
        print("## psst.py - Clustering handling ##")
        validate_cl(args)
        clmain(args)
    if args.mode == "map":
        print("## psst.py - Read mapping handling ##")
        validate_m(args)
        mmain(args)
    
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

def clmain(args):
    # Split into sub-mode-specific validation functions
    if args.clusteringMode == "mmseqs2":
        print("## MMseqs2 clustering ##")
        validate_cl_mmseqs2(args)
    if args.clusteringMode == "cdhit":
        print("## CD-HIT clustering ##")
        validate_cl_cdhit(args)

def mmain(args):
    # Split into sub-mode-specific validation functions
    if args.mapMode == "wgs":
        print("## bwa mem WGS read mapping ##")
        validate_m_wgs(args)
    if args.mapMode == "hic":
        print("## bwa mem Hi-C read mapping and .hic generation ##")
        validate_m_hic(args)
    if args.mapMode == "index":
        validate_m_index(args)
        if args.mapIndexMode == "bwa":
            print("## bwa reference indexing ##")
            validate_m_index_bwa(args)
        if args.mapIndexMode == "salmon":
            print("## salmon transcriptome indexing ##")
            validate_m_index_salmon(args)
        if args.mapIndexMode == "star":
            print("## STAR genome indexing ##")
            validate_m_index_star(args)
    if args.mapMode == "rna":
        validate_m_rna(args)
        if args.mapRnaMode == "salmon":
            print("## salmon RNAseq quantification ##")
            validate_m_rna_salmon(args)

if __name__ == "__main__":
    main()
