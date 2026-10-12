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

import os, re, shutil, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from locations import MMSEQS_DIR, VARSCRIPT_DIR, REFSEQ_PLANT_FASTA, UNIREF90_FASTA, \
    IDMAPPING_FILE, OBO_FILE

# Note: arguments with a default location are given default=None in the parser so that
# validation can tell whether the user specified a value. User-specified locations are
# checked to exist; default locations are substituted in without any checking.

MAX_WALLTIME_HOURS = 48
RPID_SHORTCUTS = {
    "GNGO": "XNOEYW5200",
    "SPC": "ULJQSK8720",
    "PHD": "HLCVGA3552"
}

class DirectoryNotFoundError(Exception):
    pass

class FormatError(Exception):
    pass

def parse_walltime(walltime):
    '''
    Parameters:
        walltime -- a string in format 'H', 'H:MM', or 'H:MM:SS'.
    Returns:
        hours -- an int of the number of hours.
        minutes -- an int of the number of minutes.
        seconds -- an int of the number of seconds.
    '''
    match = re.match(r"^(\d+)(?::(\d{1,2}))?(?::(\d{1,2}))?$", walltime)
    if match == None:
        raise FormatError(f"--walltime '{walltime}' is not in 'H', 'H:MM', or 'H:MM:SS' format")
    hours, minutes, seconds = [int(x) if x != None else 0 for x in match.groups()]
    if minutes >= 60 or seconds >= 60:
        raise FormatError(f"--walltime '{walltime}' has minutes or seconds >= 60")
    return hours, minutes, seconds

def validate_rpid(rpid):
    '''
    Parameters:
        rpid -- a string of a Research Project ID (6 letters followed by 4 numbers)
                or a shortcut key within RPID_SHORTCUTS.
    Returns:
        rpid -- a string of the full-length RPID.
    '''
    rpid = rpid.upper()
    if rpid in RPID_SHORTCUTS:
        return RPID_SHORTCUTS[rpid]
    if re.match(r"^[A-Z]{6}\d{4}$", rpid) == None:
        raise FormatError(f"--rpid '{rpid}' must be 6 letters followed by 4 numbers, or one of " +
                          ", ".join(RPID_SHORTCUTS))
    return rpid

def validate_input_file(fileName, flag, description, default=None):
    '''
    Parameters:
        fileName -- a string indicating the location of a file, or None if not specified.
        flag -- a string of the argument flag used to specify the file.
        description -- a string describing the file for error messages.
        default -- OPTIONAL; a string of the default file location to use, without
                   checking, if fileName is None.
    Returns:
        fileName -- the absolute path of the file, or default if not specified.
    '''
    if fileName == None:
        return default
    fileName = os.path.abspath(fileName)
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"{description} ({flag} {fileName}) does not exist!")
    return fileName

def validate_output_prefix(outputPrefix):
    '''
    Parameters:
        outputPrefix -- a string indicating the prefix for output files, or None if
                        not specified.
    Returns:
        outputPrefix -- the absolute path of the output prefix, or None if not specified.
    '''
    if outputPrefix == None:
        return None
    if os.path.basename(outputPrefix) == "":
        raise ValueError(f"Output prefix (-o {outputPrefix}) must not end with a directory separator")
    outputPrefix = os.path.abspath(outputPrefix)
    if not os.path.isdir(os.path.dirname(outputPrefix)):
        raise DirectoryNotFoundError(f"Output prefix (-o {outputPrefix}) is in a directory that does not exist!")
    return outputPrefix

def validate_program_dir(dirName, requiredFile, flag, default=None):
    '''
    Parameters:
        dirName -- a string indicating the location of a program directory, or None
                   if not specified.
        requiredFile -- a string of a file path, relative to dirName, which must exist.
        flag -- a string of the argument flag used to specify the directory.
        default -- OPTIONAL; a string of the default directory location to use, without
                   checking, if dirName is None.
    Returns:
        dirName -- the absolute path of the directory, or default if not specified.
    '''
    if dirName == None:
        return default
    dirName = os.path.abspath(dirName)
    if not os.path.isfile(os.path.join(dirName, requiredFile)):
        raise FileNotFoundError(f"{requiredFile} not found in {flag} {dirName}")
    return dirName

def validate_mmseqs2(args):
    '''
    Validation for MMseqs2 and Various_scripts locations shared by any job that
    runs Various_scripts/run_mmseqs2.py.
    '''
    args.mmseqsDir = validate_program_dir(args.mmseqsDir, "mmseqs", "--mmseqs", default=MMSEQS_DIR)

    args.varScriptDir = validate_program_dir(args.varScriptDir, "run_mmseqs2.py", "--varscripts", default=VARSCRIPT_DIR)

    if args.iterations != None and args.iterations < 1:
        raise ValueError("--iterations must be >= 1")

def validate_pbs(args, scriptSuffix):
    '''
    Validation for PBS arguments common to all job scripts. Should be called at the
    end of each job's validation function since the default script name is derived
    from args.outputPrefix.

    Parameters:
        args -- an argparse.Namespace containing PBS arguments and args.outputPrefix.
        scriptSuffix -- a string to name the job script with when deriving a default
                        script file name; either '<outputPrefix>.<scriptSuffix>.sh' or
                        '<scriptSuffix>.sh' if args.outputPrefix is None.
    '''
    # Validate resource requests
    if args.cpus < 1:
        raise ValueError("--cpus must be >= 1")

    if re.match(r"^\d+[kmgt]?b?$", args.mem, re.IGNORECASE) == None:
        raise FormatError(f"--mem '{args.mem}' is not a valid PBS memory value (e.g., '70G' or '70gb')")

    hours, minutes, seconds = parse_walltime(args.walltime)
    if (hours * 3600 + minutes * 60 + seconds) > MAX_WALLTIME_HOURS * 3600:
        raise ValueError(f"--walltime '{args.walltime}' exceeds the HPC limit of {MAX_WALLTIME_HOURS} hours")
    if (hours + minutes + seconds) == 0:
        raise ValueError("--walltime must be greater than zero")
    args.walltime = f"{hours:02d}:{minutes:02d}:{seconds:02d}"

    # Validate job identifiers
    if re.match(r"^[A-Za-z][^\s@#]*$", args.jobName) == None:
        raise FormatError(f"--name '{args.jobName}' must start with a letter and not contain whitespace, '@' or '#'")
    args.rpid = validate_rpid(args.rpid)

    # Validate job dependencies
    if args.afterok and args.afternotok:
        raise ValueError("--afterok and --afternotok are mutually exclusive")
    for flag, jobIDs in [("--afterok", args.afterok), ("--afternotok", args.afternotok)]:
        for jobID in (jobIDs or []):
            if re.match(r"^[^\s:,]+$", jobID) == None:
                raise FormatError(f"{flag} job ID '{jobID}' must not contain whitespace, ':' or ','")

    # Validate that qsub is available if submission is requested
    if args.qsub and shutil.which("qsub") == None:
        raise FileNotFoundError("--qsub was specified but qsub is not available in your PATH")
    
    # Validate script file name
    if args.scriptFile == None:
        if args.outputPrefix != None:
            args.scriptFile = f"{args.outputPrefix}.{scriptSuffix}.sh"
        else:
            args.scriptFile = f"{scriptSuffix}.sh"
    args.scriptFile = os.path.abspath(args.scriptFile)
    if os.path.exists(args.scriptFile):
        raise FileExistsError(f"Job script (--script {args.scriptFile}) already exists!")
    if not os.path.isdir(os.path.dirname(args.scriptFile)):
        raise DirectoryNotFoundError(f"Job script (--script {args.scriptFile}) is in a directory that does not exist!")

def validate_at(args):
    '''
    Validation for arguments common to all "annotable" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_at_generate(args):
    '''
    Validation for arguments used in "annotable generate" mode.
    '''
    # Validate input files
    args.query = validate_input_file(args.query, "-q", "Query FASTA file")
    args.database = validate_input_file(args.database, "-d", "Database FASTA file", default=REFSEQ_PLANT_FASTA)
    args.idmappingFile = validate_input_file(args.idmappingFile, "--idmapping", "idmapping file", default=IDMAPPING_FILE)
    args.oboFile = validate_input_file(args.oboFile, "--obo", "GO .obo file", default=OBO_FILE)

    # Validate program locations and parameters
    if args.varScriptDir != None:
        validate_program_dir(args.varScriptDir, os.path.join("annotation_table", "generate_annotation_table.py"), "--varscripts")
    validate_mmseqs2(args)

    if args.numHits < 1:
        raise ValueError("--numhits must be >= 1")

    # Validate output locations
    args.outputPrefix = validate_output_prefix(args.outputPrefix)
    if args.tmpDir == None and args.outputPrefix != None:
        args.tmpDir = f"{args.outputPrefix}_tmp"
    if args.tmpDir != None:
        args.tmpDir = os.path.abspath(args.tmpDir)

    validate_pbs(args, "annotable")

def validate_b(args):
    '''
    Validation for arguments common to all "blast" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_b_run(args):
    '''
    Validation for arguments common to all "blast run" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_b_run_mmseqs2(args):
    '''
    Validation for arguments used in "blast run mmseqs2" mode.
    '''
    # Validate input files
    args.query = validate_input_file(args.query, "-q", "Query FASTA file")
    args.database = validate_input_file(args.database, "-d", "Database FASTA file", default=UNIREF90_FASTA)

    # Validate program locations and parameters
    validate_mmseqs2(args)

    # Validate output locations
    args.outputPrefix = validate_output_prefix(args.outputPrefix)
    if args.tmpDir == None and args.outputPrefix != None:
        args.tmpDir = f"{args.outputPrefix}_tmp"
    if args.tmpDir != None:
        args.tmpDir = os.path.abspath(args.tmpDir)

    validate_pbs(args, "mms2")
