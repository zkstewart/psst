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

def format_pbs_header(args, afterok=None, afternotok=None):
    '''
    Parameters:
        args -- an argparse.Namespace containing validated PBS arguments
                (jobName, walltime, mem, cpus, rpid).
        afterok -- OPTIONAL; a list of job ID strings which must complete successfully
                   before this job can start.
        afternotok -- OPTIONAL; a list of job ID strings which must fail before this
                      job can start.
    Returns:
        header -- a string of the PBS directives and preamble for a job script.
    '''
    if afterok and afternotok:
        raise ValueError("afterok and afternotok dependencies are mutually exclusive")

    if afterok:
        dependLine = "\n#PBS -W depend=afterok:" + ":".join(afterok)
    elif afternotok:
        dependLine = "\n#PBS -W depend=afternotok:" + ":".join(afternotok)
    else:
        dependLine = ""

    header, _ = fill_template("pbs_header.sh", {
        "jobName": args.jobName,
        "walltime": args.walltime,
        "mem": args.mem,
        "cpus": args.cpus,
        "rpid": args.rpid,
        "dependLine": dependLine
    }, quote=False) # values are validated to be shell-safe; dependLine must not be quoted
    return header

def write_job_script(args, body):
    '''
    Writes a PBS job script consisting of the PBS header followed by the job body.

    Parameters:
        args -- an argparse.Namespace containing validated PBS arguments,
                args.afterok and args.afternotok job dependencies, and
                args.scriptFile indicating where to write the script.
        body -- a string containing the job-specific script contents.
    '''
    header = format_pbs_header(args, afterok=args.afterok, afternotok=args.afternotok)
    with open(args.scriptFile, "w") as fileOut:
        fileOut.write(header + body)
    print(f"Job script written to '{args.scriptFile}'")

def submit_job(scriptFile):
    '''
    Submits a job script with qsub from the current working directory.

    Parameters:
        scriptFile -- a string indicating the location of the job script to submit.
    Returns:
        jobID -- a string of the job ID reported by qsub.
    '''
    result = subprocess.run(["qsub", scriptFile], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"qsub failed with exit code {result.returncode}:\n{result.stderr}")

    jobID = result.stdout.strip()
    print(f"Job submitted as '{jobID}'")
    return jobID

def finalise_job(args, templateName, values):
    '''
    Fills a job template, writes the job script and, if requested, submits it with qsub.
    Submission is refused if any placeholders remain in the script.

    Parameters:
        args -- an argparse.Namespace containing validated PBS arguments and args.qsub.
        templateName -- a string of the job template file name within the templates directory.
        values -- a dict of template field names to values (or None).
    '''
    body, missing = fill_template(templateName, values)
    if args.qsub and missing != []:
        raise ValueError("--qsub cannot be used while the job script has unfilled placeholders: " +
                         ", ".join(missing))
    
    write_job_script(args, body)
    if missing != []:
        print("Job script contains placeholders to fill in before submission: " + ", ".join(missing))
    if args.qsub:
        submit_job(args.scriptFile)
