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

import os

from .pbs import fill_template, format_pbs_header, write_job_script, submit_job, check_placeholders

# Layout of a dorado basecalling run directory, relative to its root
CORRECTED_DIR = "fastq_corrected"
DORADO_SCRIPTS = {
    "p1": "basecaller_part1.sh",
    "p2": "basecaller_part2.sh",
    "corr": "run_dorado_correction.sh",
    "concat": os.path.join(CORRECTED_DIR, "concat_corrected_fastq.sh")
}

def basecalling_dorado(args):
    '''
    Generates (and optionally submits) the four PBS job scripts of a dorado basecalling
    and read correction run, laid out within args.workDir as:
        basecaller_part1.sh -- GPU basecalling of pod5/ into bam_sup/
        basecaller_part2.sh -- BAM to filtered FASTQ in fastq/ and correction block count
        run_dorado_correction.sh -- submits the step 1 (CPU) and step 2 (GPU) correction
                                    jobs, then submits the concatenation job to run
                                    after step 2 succeeds
        fastq_corrected/concat_corrected_fastq.sh -- concatenates corrected reads

    With --qsub, jobs 1-3 are submitted as an afterok chain; job 4 is submitted by job 3
    since it depends on correction jobs which only exist once job 3 has run.

    Parameters:
        args -- an argparse.Namespace containing validated "basecalling dorado" arguments.
    '''
    common = {
        "doradoDir": args.doradoDir,
        "varScriptDir": args.varScriptDir,
        "prefix": args.prefix,
        "model": args.model,
        "minReadLen": args.minReadLen
    }
    jobs = [ # step, template values, extra header arguments
        ("p1", common, { "gpus": 1, "gpuID": args.gpuID,
                         "afterok": args.afterok, "afternotok": args.afternotok }),
        ("p2", common, {}),
        ("corr", {
            "prefix": args.prefix,
            "minReadLen": args.minReadLen,
            "doradoDir": args.doradoDir,
            "varScriptDir": args.varScriptDir,
            "jobPrefix": args.jobName,
            "s1Cpus": args.s1Cpus, "s1Mem": args.s1Mem, "s1Walltime": args.s1Walltime,
            "s2Cpus": args.s2Cpus, "s2Mem": args.s2Mem, "s2Walltime": args.s2Walltime
        }, {}),
        ("concat", {
            "varScriptDir": args.varScriptDir,
            "genScriptDir": args.genScriptDir,
            "prefix": args.prefix,
            "concatMinimum": args.concatMinimum,
            "concatMultiline": args.concatMultiline
        }, {})
    ]

    # Render all scripts before writing any, so placeholders block --qsub up front
    scripts = []
    missing = []
    for step, values, headerArgs in jobs:
        body, stepMissing = fill_template(os.path.basename(DORADO_SCRIPTS[step]), values)
        header = format_pbs_header(f"{args.jobName}_{step}",
                                   getattr(args, f"{step}Walltime"),
                                   getattr(args, f"{step}Mem"),
                                   getattr(args, f"{step}Cpus"),
                                   args.rpid, **headerArgs)
        scripts.append((args.scriptFiles[step], header, body))
        missing += [ m for m in stepMissing if m not in missing ]
    if args.qsub and missing != []:
        check_placeholders(args, missing) # raises an error before anything is written

    # Write scripts into the run directory layout
    os.makedirs(os.path.join(args.workDir, CORRECTED_DIR), exist_ok=True)
    for scriptFile, header, body in scripts:
        write_job_script(scriptFile, header, body)
    check_placeholders(args, missing)

    # Submit jobs 1-3 as a chain from the run directory
    if args.qsub:
        p1ID = submit_job(args.scriptFiles["p1"], cwd=args.workDir)
        p2ID = submit_job(args.scriptFiles["p2"], afterok=p1ID, cwd=args.workDir)
        submit_job(args.scriptFiles["corr"], afterok=p2ID, cwd=args.workDir)
        print(f"'{DORADO_SCRIPTS['concat']}' will be submitted by '{DORADO_SCRIPTS['corr']}' " +
              "once the correction jobs are queued")
