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

import os, re, shlex, subprocess

from .locations import TEMPLATE_DIR

# Template placeholders look like '@@fieldName@@'; this syntax is chosen to never clash
# with bash syntax, and an unfilled placeholder remains a valid (if meaningless) bash word
PLACEHOLDER_REGEX = re.compile(r"@@([A-Za-z_][A-Za-z0-9_]*)@@")

def q(value):
    '''
    Shell-quote a value for safe inclusion in a job script; values without
    special characters are returned unchanged so the script remains readable.

    Parameters:
        value -- any value that can be converted to a string.
    Returns:
        quoted -- a string safe to embed in a bash script.
    '''
    return shlex.quote(str(value))

def load_template(templateName):
    '''
    Parameters:
        templateName -- a string of the template file name within the templates directory.
    Returns:
        template -- a string of the template contents containing '@@field@@' placeholders.
    '''
    templateFile = os.path.join(TEMPLATE_DIR, templateName)
    if not os.path.isfile(templateFile):
        raise FileNotFoundError(f"Template file '{templateFile}' does not exist!")
    with open(templateFile, "r") as fileIn:
        return fileIn.read()

def fill_template(templateName, values, quote=True):
    '''
    Fills a template's '@@field@@' placeholders with values. Values that are None are
    intentionally left as placeholders for the user to fill in later.
    
    Every template field must be a key in values, and every key in values must be a
    template field; this catches typos between templates and the modules filling them.
    
    Parameters:
        templateName -- a string of the template file name within the templates directory.
        values -- a dict of template field names to values (or None).
        quote -- OPTIONAL; a boolean indicating whether to shell-quote values; default == True.
    Returns:
        filled -- a string of the template with values filled in.
        missing -- a list of template field names left as placeholders.
    '''
    template = load_template(templateName)
    
    # Check that the template and values agree on field names
    fields = list(dict.fromkeys(PLACEHOLDER_REGEX.findall(template))) # ordered, unique
    notGiven = [ f for f in fields if f not in values ]
    notUsed = [ k for k in values if k not in fields ]
    if notGiven != [] or notUsed != []:
        raise KeyError(f"Template '{templateName}' fields do not match the values given; " +
                       f"fields without values = {notGiven}, values without fields = {notUsed}")
    
    # Fill in the template
    def substitute(match):
        value = values[match.group(1)]
        if value == None:
            return match.group(0)
        return q(value) if quote else str(value)
    filled = PLACEHOLDER_REGEX.sub(substitute, template)
    missing = [ f for f in fields if values[f] == None ]
    return filled, missing

def format_pbs_header(jobName, walltime, mem, cpus, rpid, afterok=None, afternotok=None,
                      gpus=None, gpuID=None):
    '''
    Parameters:
        jobName -- a string of the PBS job name.
        walltime -- a string of the walltime in 'HH:MM:SS' format.
        mem -- a string of the memory request (e.g., '70G').
        cpus -- an int of the number of CPUs to request.
        rpid -- a string of the Research Project ID to attach to the job.
        afterok -- OPTIONAL; a list of job ID strings which must complete successfully
                   before this job can start.
        afternotok -- OPTIONAL; a list of job ID strings which must fail before this
                      job can start.
        gpus -- OPTIONAL; an int of the number of GPUs to request.
        gpuID -- OPTIONAL; a string of the GPU type to request (e.g., 'H100').
    Returns:
        header -- a string of the PBS directives and preamble for a job script.
    '''
    if afterok and afternotok:
        raise ValueError("afterok and afternotok dependencies are mutually exclusive")
    
    # Format optional directives
    extraDirectives = ""
    if gpus != None:
        extraDirectives += f"\n#PBS -l ngpus={gpus}"
    if gpuID != None:
        extraDirectives += f"\n#PBS -l gpu_id={gpuID}"
    if afterok:
        extraDirectives += "\n#PBS -W depend=afterok:" + ":".join(afterok)
    elif afternotok:
        extraDirectives += "\n#PBS -W depend=afternotok:" + ":".join(afternotok)
    
    header, _ = fill_template("pbs_header.sh", {
        "jobName": jobName,
        "walltime": walltime,
        "mem": mem,
        "cpus": cpus,
        "rpid": rpid,
        "extraDirectives": extraDirectives
    }, quote=False) # values are validated to be shell-safe; directives must not be quoted
    return header

def write_job_script(scriptFile, header, body):
    '''
    Writes a PBS job script consisting of the PBS header followed by the job body.
    
    Parameters:
        scriptFile -- a string indicating where to write the script.
        header -- a string of the PBS header as produced by format_pbs_header().
        body -- a string containing the job-specific script contents.
    '''
    with open(scriptFile, "w") as fileOut:
        fileOut.write(header + body)
    print(f"Job script written to '{scriptFile}'")

def submit_job(scriptFile, afterok=None, cwd=None):
    '''
    Submits a job script with qsub. The job's $PBS_O_WORKDIR is the directory qsub
    is run from.
    
    Parameters:
        scriptFile -- a string indicating the location of the job script to submit.
        afterok -- OPTIONAL; a string of a job ID which must complete successfully
                   before this job can start; overrides any dependency in the script.
        cwd -- OPTIONAL; a string of the directory to run qsub from; default == the
               current working directory.
    Returns:
        jobID -- a string of the job ID reported by qsub.
    '''
    cmd = ["qsub"]
    if afterok != None:
        cmd += ["-W", f"depend=afterok:{afterok}"]
    result = subprocess.run(cmd + [scriptFile], capture_output=True, text=True, cwd=cwd)
    if result.returncode != 0:
        raise RuntimeError(f"qsub failed with exit code {result.returncode}:\n{result.stderr}")
    
    jobID = result.stdout.strip()
    print(f"Job '{os.path.basename(scriptFile)}' submitted as '{jobID}'")
    return jobID

def check_placeholders(args, missing):
    '''
    Reports placeholders left in job scripts, raising an error if they would be
    submitted with --qsub.
    
    Parameters:
        args -- an argparse.Namespace containing args.qsub.
        missing -- a list of template field names left as placeholders.
    '''
    if missing == []:
        return
    if args.qsub:
        raise ValueError("--qsub cannot be used while job scripts have unfilled placeholders: " +
                         ", ".join(missing))
    print("Job script(s) contain placeholders to fill in before submission: " + ", ".join(missing))

def finalise_job(args, templateName, values):
    '''
    Fills a job template, writes the job script and, if requested, submits it with qsub.
    Submission is refused if any placeholders remain in the script. Used by modes which
    produce a single job script using the arguments of add_pbs_arguments() in psst.py.
    
    Parameters:
        args -- an argparse.Namespace containing validated PBS arguments and args.qsub.
        templateName -- a string of the job template file name within the templates directory.
        values -- a dict of template field names to values (or None).
    '''
    body, missing = fill_template(templateName, values)
    if args.qsub and missing != []:
        check_placeholders(args, missing) # raises an error before anything is written
    
    header = format_pbs_header(args.jobName, args.walltime, args.mem, args.cpus, args.rpid,
                               afterok=args.afterok, afternotok=args.afternotok)
    write_job_script(args.scriptFile, header, body)
    check_placeholders(args, missing)
    if args.qsub:
        submit_job(args.scriptFile)
