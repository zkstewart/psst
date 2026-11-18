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
from .pbs import PLACEHOLDER_REGEX, fill_template, format_pbs_header, write_job_script, submit_job, check_placeholders, finalise_job

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
                              Various_scripts' tabulate_fastqc.py

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
    jobs = [ # step, template, template values, extra header arguments
        ("fastqc", "qc/run_fastqc.sh", {
            **reads,
            "fastqcModule": args.fastqcModule,
            "cpus": args.fastqcCpus
        }, { "afterok": args.afterok, "afternotok": args.afternotok,
             "arraySize": get_array_size(args.numSamples) }),
        ("tabulate", "qc/run_tabulate_fastqc.sh", {
            **reads,
            "varScriptDir": args.varScriptDir,
            "outputFile": args.outputFile
        }, {})
    ]

    # Render both scripts before writing either, so placeholders block --qsub up front
    scripts = []
    missing = []
    for step, templateName, values, headerArgs in jobs:
        body, stepMissing = fill_template(templateName, values)
        jobName = args.jobName if step == "fastqc" else f"{args.jobName}_{step}"
        header = format_pbs_header(jobName,
                                   getattr(args, f"{step}Walltime"),
                                   getattr(args, f"{step}Mem"),
                                   getattr(args, f"{step}Cpus"),
                                   args.rpid, **headerArgs)
        scripts.append((args.scriptFiles[step], header, body))
        stepMissing = PLACEHOLDER_REGEX.findall(header) + stepMissing # e.g., the array size
        missing += [ m for m in stepMissing if m not in missing ]
    if args.qsub and missing != []:
        check_placeholders(args, missing) # raises an error before anything is written

    for scriptFile, header, body in scripts:
        write_job_script(scriptFile, header, body)
    check_placeholders(args, missing)

    # Submit the tabulation job to follow the FastQC job
    if args.qsub:
        fastqcID = submit_job(args.scriptFiles["fastqc"])
        submit_job(args.scriptFiles["tabulate"], afterok=fastqcID)
    else:
        print(f"Submit '{args.scriptFiles['tabulate']}' to run after the FastQC job succeeds, " +
              "i.e., with 'qsub -W depend=afterok:<FastQC job ID>'")
