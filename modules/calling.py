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

import gzip, re, struct

from .mapping import get_array_size
from .pbs import fill_template, format_pbs_header, write_job_script, submit_job, check_placeholders

# Files which can index a BAM file, given as suffixes appended to its full name or
# replacing its '.bam' extension
BAM_INDEX_SUFFIXES = [".bai", ".csi"]

# bcftools mpileup -X platform profile for each "call" mode; modes without short read
# arguments (i.e., minBaseQ, hic, indelsCns) leave these settings to the profile
CALL_PROFILES = {
    "short": "",
    "nanopore": "ont-sup"
}

def get_bam_samples(bamFile):
    '''
    Reads the sample names (SM) of a BAM file's read groups from its header, which is
    at the start of the BGZF (i.e., gzip-compatible) stream.

    Parameters:
        bamFile -- a string of the location of a BAM file.
    Returns:
        samples -- a sorted list of sample names, empty if the BAM has no read groups,
                   or None if the file could not be read as a BAM.
    '''
    try:
        with gzip.open(bamFile, "rb") as fileIn:
            if fileIn.read(4) != b"BAM\x01":
                return None
            textLength = struct.unpack("<i", fileIn.read(4))[0]
            text = fileIn.read(textLength).decode(errors="replace")
    except (OSError, EOFError, struct.error):
        return None
    samples = set()
    for line in text.split("\n"):
        if line.startswith("@RG"):
            samples.update(re.findall(r"\tSM:([^\t]+)", line))
    return sorted(samples)

def call_variants(args):
    '''
    Generates (and optionally submits) the two PBS job scripts of a "call" mode:
        <name>.sh -- an array job which divides the genome into chunks of whole contigs,
                     each subjob calling variants jointly across all BAM files with
                     bcftools mpileup | call | norm
        <name>_concat.sh -- concatenates the chunks in genome order into
                            <outdir>/<prefix>.vcf.gz and indexes it

    With --qsub, the concatenation job is submitted to run after every calling subjob
    succeeds.

    Parameters:
        args -- an argparse.Namespace containing validated "call" mode arguments.
    '''
    common = {
        "outputDir": args.outputDir,
        "prefix": args.prefix,
        "chunks": args.chunks
    }
    jobs = [ # step, template, template values, extra header arguments
        ("call", "call/run_bcftools_call.sh", {
            **common,
            "profile": CALL_PROFILES[args.callMode],
            "genome": args.genome,
            "bamDir": args.bamDir,
            "bamSuffix": args.bamSuffix,
            "minMapQ": args.minMapQ,
            "minBaseQ": getattr(args, "minBaseQ", ""), # empty uses the profile's default
            "countOrphans": "true" if getattr(args, "hic", False) else "false", # Hi-C pairs aren't proper
            "maxDepth": args.maxDepth,
            "indelsCns": "true" if getattr(args, "indelsCns", False) else "false",
            "ploidy": args.ploidy,
            "cpus": args.callCpus
        }, { "afterok": args.afterok, "afternotok": args.afternotok,
             "arraySize": get_array_size(args.chunks) }),
        ("concat", "call/run_bcftools_concat.sh", {
            **common,
            "cpus": args.concatCpus
        }, {})
    ]

    # Render both scripts before writing either, so placeholders block --qsub up front
    scripts = []
    missing = []
    for step, templateName, values, headerArgs in jobs:
        body, stepMissing = fill_template(templateName, values)
        jobName = args.jobName if step == "call" else f"{args.jobName}_{step}"
        header = format_pbs_header(jobName,
                                   getattr(args, f"{step}Walltime"),
                                   getattr(args, f"{step}Mem"),
                                   getattr(args, f"{step}Cpus"),
                                   args.rpid, **headerArgs)
        scripts.append((args.scriptFiles[step], header, body))
        missing += [ m for m in stepMissing if m not in missing ]
    if args.qsub and missing != []:
        check_placeholders(args, missing) # raises an error before anything is written

    for scriptFile, header, body in scripts:
        write_job_script(scriptFile, header, body)
    check_placeholders(args, missing)

    # Submit the concatenation job to follow the calling job
    if args.qsub:
        callID = submit_job(args.scriptFiles["call"])
        submit_job(args.scriptFiles["concat"], afterok=callID)
    else:
        print(f"Submit '{args.scriptFiles['concat']}' to run after the calling job succeeds, " +
              "i.e., with 'qsub -W depend=afterok:<calling job ID>'")
