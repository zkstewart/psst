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

from .mapping import get_array_size
from .pbs import finalise_job, finalise_chained_jobs

def qc_library_trim(args):
    '''
    Generates (and optionally submits) a PBS array job script which trims adapters and
    low-quality ends from reads with bbduk. Each subjob trims one sample; the array size
    is a placeholder if no reads directory was given, and no array is requested if there
    is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "qc library trim" arguments.
    '''
    finalise_job(args, "qc/run_bbduk_trim.sh", {
        "condaEnv": args.condaEnv,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir,
        "adapters": args.adapters,
        "trimq": args.trimq,
        "minLength": args.minLength,
        "cpus": args.cpus,
        "xmx": args.xmx
    }, arraySize=get_array_size(args.numSamples))

def qc_library_fastqc(args):
    '''
    Generates (and optionally submits) the two PBS job scripts of "qc library fastqc":
        <name>.sh -- an array job running FastQC on each sample's reads, writing reports
                     to <outdir>/<sample>/
        <name>_tabulate.sh -- summarises every sample's reports in an Excel table with
                              annotarium's 'qc fastqc tabulate'

    With --qsub, the tabulation job is submitted to run after every FastQC subjob succeeds.

    Parameters:
        args -- an argparse.Namespace containing validated "qc library fastqc" arguments.
    '''
    reads = {
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir
    }
    finalise_chained_jobs(args, [ # step, template, template values, array size
        ("fastqc", "qc/run_fastqc.sh", {
            **reads,
            "fastqcModule": args.fastqcModule,
            "cpus": args.fastqcCpus
        }, get_array_size(args.numSamples)),
        ("tabulate", "qc/run_tabulate_fastqc.sh", {
            **reads,
            "annotariumDir": args.annotariumDir,
            "outputFile": args.outputFile
        }, None)
    ])

def qc_library_kraken2(args):
    '''
    Generates (and optionally submits) the two PBS job scripts of "qc library kraken2":
        <name>.sh -- classifies each sample's reads in turn with kraken2, writing its
                     classifications, report, (un)classified reads and a Krona chart to
                     <outdir>/<sample>/; samples run serially in one job (not an array),
                     since loading the database takes far longer than classifying a sample
        <name>_tabulate.sh -- combines every sample's report with KrakenTools'
                              combine_kreports.py, then reformats it with Various_scripts'
                              kreport_reformat.py

    With --qsub, the tabulation job is submitted to run after the kraken2 job succeeds.

    Parameters:
        args -- an argparse.Namespace containing validated "qc library kraken2" arguments.
    '''
    reads = {
        "condaEnv": args.condaEnv,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir
    }
    finalise_chained_jobs(args, [ # step, template, template values, array size
        ("classify", "qc/run_kraken2.sh", {
            **reads,
            "db": args.db,
            "cpus": args.classifyCpus
        }, None),
        ("tabulate", "qc/run_tabulate_kraken2.sh", {
            **reads,
            "varScriptDir": args.varScriptDir,
            "prefix": args.prefix
        }, None)
    ])

def qc_library_smudgeplot(args):
    '''
    Generates (and optionally submits) the two PBS job scripts of "qc library smudgeplot":
        <name>.sh -- an array job counting each sample's k-mers with FastK, then inferring
                     its ploidy with smudgeplot, writing results to <outdir>/<sample>/
        <name>_tabulate.sh -- combines every sample's most likely ploidy into one table

    With --qsub, the tabulation job is submitted to run after every smudgeplot subjob succeeds.

    Parameters:
        args -- an argparse.Namespace containing validated "qc library smudgeplot" arguments.
    '''
    reads = {
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir
    }
    finalise_chained_jobs(args, [ # step, template, template values, array size
        ("smudgeplot", "qc/run_smudgeplot.sh", {
            **reads,
            "condaEnv": args.condaEnv,
            "errorThreshold": "" if args.errorThreshold == None else args.errorThreshold,
            "genomeSize": args.genomeSize,
            "ploidy": args.ploidy,
            "keepTable": "true" if args.keepTables else "false",
            "cpus": args.smudgeplotCpus,
            "fastkMem": args.fastkMem
        }, get_array_size(args.numSamples)),
        ("tabulate", "qc/run_tabulate_smudgeplot.sh", {
            **reads,
            "outputFile": args.outputFile
        }, None)
    ])
