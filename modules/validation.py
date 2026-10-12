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

import os, re, shutil

from .basecalling import DORADO_SCRIPTS
from .config import CONFIG, RPID_REGEX

# Note: only values the user explicitly specified (recorded in args.userSet by the
# UserSet argparse action) are checked to exist; default locations are used unchecked.

MAX_WALLTIME_HOURS = 48
RPID_SHORTCUTS = CONFIG["rpid_shortcuts"]

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
    if RPID_REGEX.match(rpid) == None:
        raise FormatError(f"--rpid '{rpid}' must be 6 letters followed by 4 numbers, or one of " +
                          ", ".join(RPID_SHORTCUTS))
    return rpid

def validate_input_file(args, dest, flag, description):
    '''
    Checks that a user-specified file exists and makes its path absolute. Values not
    specified by the user (i.e., defaults or None placeholders) are left unchanged.
    
    Parameters:
        args -- an argparse.Namespace containing args.userSet and the dest attribute.
        dest -- a string of the argparse dest holding the file location.
        flag -- a string of the argument flag used to specify the file.
        description -- a string describing the file for error messages.
    '''
    if dest not in args.userSet:
        return
    fileName = os.path.abspath(getattr(args, dest))
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"{description} ({flag} {fileName}) does not exist!")
    setattr(args, dest, fileName)

def validate_program_dir(args, dest, requiredFile, flag):
    '''
    Checks that a user-specified program directory contains a required file and makes
    its path absolute. Values not specified by the user are left unchanged.
    
    Parameters:
        args -- an argparse.Namespace containing args.userSet and the dest attribute.
        dest -- a string of the argparse dest holding the directory location.
        requiredFile -- a string of a file path, relative to the directory, which must exist.
        flag -- a string of the argument flag used to specify the directory.
    '''
    if dest not in args.userSet:
        return
    dirName = os.path.abspath(getattr(args, dest))
    if not os.path.isfile(os.path.join(dirName, requiredFile)):
        raise FileNotFoundError(f"{requiredFile} not found in {flag} {dirName}")
    setattr(args, dest, dirName)

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

def validate_mmseqs2(args):
    '''
    Validation for arguments shared by any job that runs Various_scripts/run_mmseqs2.py
    (i.e., those added by add_mmseqs2_arguments() in psst.py).
    '''
    # Validate input files
    validate_input_file(args, "query", "-q", "Query FASTA file")
    validate_input_file(args, "database", "-d", "Database FASTA file")
    
    # Validate program locations
    validate_program_dir(args, "mmseqsDir", "mmseqs", "--mmseqs")
    validate_program_dir(args, "varScriptDir", "run_mmseqs2.py", "--varscripts")
    
    # Validate parameters
    if args.iterations < 1:
        raise ValueError("--iterations must be >= 1")
    
    # Validate output locations
    args.outputPrefix = validate_output_prefix(args.outputPrefix)
    if args.tmpDir == None and args.outputPrefix != None:
        args.tmpDir = f"{args.outputPrefix}_tmp"
    if args.tmpDir != None:
        args.tmpDir = os.path.abspath(args.tmpDir)

def validate_pbs_resources(cpus, mem, walltime, flagPrefix=""):
    '''
    Validates the resource requests for one PBS job.
    
    Parameters:
        cpus -- an int of the number of CPUs to request.
        mem -- a string of the memory to request (e.g., '70G').
        walltime -- a string of the walltime in 'H', 'H:MM', or 'H:MM:SS' format.
        flagPrefix -- OPTIONAL; a string prefixing the argument flags for error messages
                      (e.g., 'p1-' for '--p1-cpus'); default == ''.
    Returns:
        walltime -- a string of the walltime normalised to 'HH:MM:SS' format.
    '''
    if cpus < 1:
        raise ValueError(f"--{flagPrefix}cpus must be >= 1")
    
    if re.match(r"^\d+[kmgt]?b?$", mem, re.IGNORECASE) == None:
        raise FormatError(f"--{flagPrefix}mem '{mem}' is not a valid PBS memory value (e.g., '70G' or '70gb')")
    
    hours, minutes, seconds = parse_walltime(walltime)
    if (hours * 3600 + minutes * 60 + seconds) > MAX_WALLTIME_HOURS * 3600:
        raise ValueError(f"--{flagPrefix}walltime '{walltime}' exceeds the HPC limit of {MAX_WALLTIME_HOURS} hours")
    if (hours + minutes + seconds) == 0:
        raise ValueError(f"--{flagPrefix}walltime must be greater than zero")
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

def validate_job_name(jobName, flag="--name"):
    '''
    Parameters:
        jobName -- a string of a PBS job name.
        flag -- OPTIONAL; a string of the argument flag for error messages.
    '''
    if re.match(r"^[A-Za-z0-9_+][^\s@#]*$", jobName) == None:
        raise FormatError(f"{flag} '{jobName}' must start with a letter, number, '_' or '+' and " +
                          "not contain whitespace, '@' or '#'")

def validate_pbs_submission(args):
    '''
    Validation for PBS submission arguments common to all modes (i.e., those added
    by add_pbs_submission_arguments() in psst.py). Called once by psst.py after all
    mode-specific validation.
    '''
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

def validate_new_script(scriptFile, flag="--script"):
    '''
    Parameters:
        scriptFile -- a string indicating where a job script will be written.
        flag -- OPTIONAL; a string of the argument flag for error messages.
    Returns:
        scriptFile -- the absolute path of the script file.
    '''
    scriptFile = os.path.abspath(scriptFile)
    if os.path.exists(scriptFile):
        raise FileExistsError(f"Job script ({flag} {scriptFile}) already exists!")
    return scriptFile

def validate_pbs_job(args):
    '''
    Validation for the PBS arguments of modes which write a single job script (i.e.,
    those added by add_pbs_arguments() in psst.py). Called once by psst.py after all
    mode-specific validation, since the default script name is derived from the
    validated args.outputPrefix.
    
    Parameters:
        args -- an argparse.Namespace containing PBS arguments, args.outputPrefix, and
                args.scriptSuffix (set by each leaf parser) to name the script with as
                either '<outputPrefix>.<scriptSuffix>.sh' or '<scriptSuffix>.sh' if
                args.outputPrefix is None or is not an argument of the mode.
    '''
    args.walltime = validate_pbs_resources(args.cpus, args.mem, args.walltime)
    validate_job_name(args.jobName)
    
    # Validate script file name
    if args.scriptFile == None:
        if getattr(args, "outputPrefix", None) != None: # not all modes have an output prefix
            args.scriptFile = f"{args.outputPrefix}.{args.scriptSuffix}.sh"
        else:
            args.scriptFile = f"{args.scriptSuffix}.sh"
    args.scriptFile = validate_new_script(args.scriptFile)
    if not os.path.isdir(os.path.dirname(args.scriptFile)):
        raise DirectoryNotFoundError(f"Job script (--script {args.scriptFile}) is in a directory that does not exist!")

def validate_as(args):
    '''
    Validation for arguments common to all "assembly" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_as_hifiasm(args):
    '''
    Validation for arguments used in "assembly hifiasm" mode.
    '''
    # Validate program locations
    validate_program_dir(args, "genScriptDir", "genome_stats.py", "--genscripts")
    validate_program_dir(args, "compleasmDir", "compleasm.py", "--compleasm")
    
    # Validate Hi-C arguments
    if args.hicFiles == None:
        if "ploidy" in args.userSet:
            raise ValueError("--ploidy is only used for Hi-C assembly; specify --hic to use it")
    else:
        if len(args.hicFiles) not in [0, 2]:
            raise ValueError(f"--hic must be given either 0 or 2 files, not {len(args.hicFiles)}")
        for i, hicFile in enumerate(args.hicFiles):
            args.hicFiles[i] = os.path.abspath(hicFile)
            if not os.path.isfile(args.hicFiles[i]):
                raise FileNotFoundError(f"Hi-C file (--hic {args.hicFiles[i]}) does not exist!")
    
    # Validate parameters
    if args.ploidy < 1:
        raise ValueError("--ploidy must be >= 1")
    if args.lineage == "":
        raise ValueError("--lineage must be given a non-empty value")

def validate_as_revcomp(args):
    '''
    Validation for arguments used in "assembly revcomp" mode.
    '''
    validate_program_dir(args, "varScriptDir", os.path.join("assembly", "revcomp_to_match_reference.py"), "--varscripts")
    validate_input_file(args, "refGenome", "-r", "Reference genome FASTA file")
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")

def validate_as_scaffolding(args):
    '''
    Validation for arguments used in "assembly scaffolding" mode.
    '''
    validate_program_dir(args, "varScriptDir", os.path.join("assembly", "haplotype_scaffolding_pipeline.py"), "--varscripts")
    validate_input_file(args, "refHap1", "--refhap1", "Reference haplotype 1 FASTA file")
    validate_input_file(args, "refHap2", "--refhap2", "Reference haplotype 2 FASTA file")
    validate_program_dir(args, "paf2dotplotDir", "paf2dotplot.r", "--paf2dotplot")
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")

def validate_at(args):
    '''
    Validation for arguments common to all "annotable" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_at_generate(args):
    '''
    Validation for arguments used in "annotable generate" mode.
    '''
    validate_mmseqs2(args)
    
    # Validate annotation files
    validate_input_file(args, "idmappingFile", "--idmapping", "idmapping file")
    validate_input_file(args, "oboFile", "--obo", "GO .obo file")
    
    # Validate program locations and parameters
    validate_program_dir(args, "varScriptDir", os.path.join("annotation_table", "generate_annotation_table.py"), "--varscripts")
    if args.numHits < 1:
        raise ValueError("--numhits must be >= 1")

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
    validate_mmseqs2(args)

def validate_bc(args):
    '''
    Validation for arguments common to all "basecalling" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_bc_dorado(args):
    '''
    Validation for arguments used in "basecalling dorado" mode.
    '''
    # Validate the run directory
    args.workDir = os.path.abspath(args.workDir)
    if not os.path.isdir(args.workDir):
        raise DirectoryNotFoundError(f"Run directory (--workdir {args.workDir}) does not exist!")
    pod5Dir = os.path.join(args.workDir, "pod5")
    if not os.path.isdir(pod5Dir):
        if args.qsub:
            raise DirectoryNotFoundError(f"pod5 directory '{pod5Dir}' does not exist; cannot --qsub basecalling")
        print(f"WARNING: pod5 directory '{pod5Dir}' does not exist yet; it must exist before submission")
    
    # Validate program locations
    validate_program_dir(args, "doradoDir", os.path.join("bin", "dorado"), "--dorado")
    for requiredFile in [os.path.join("basecalling", "dorado_correct_pipeline.py"),
                         "fasta_handling_master_code.py", "fasta_concat.py"]:
        validate_program_dir(args, "varScriptDir", requiredFile, "--varscripts")
    validate_program_dir(args, "genScriptDir", "genome_stats.py", "--genscripts")
    
    # Validate behavioural parameters
    if args.prefix != None and re.match(r"^[^\s/]+$", args.prefix) == None:
        raise FormatError(f"-p '{args.prefix}' must not contain whitespace or '/'")
    if re.match(r"^\S+$", args.model) == None:
        raise FormatError(f"--model '{args.model}' must be non-empty and not contain whitespace")
    if re.match(r"^\w+$", args.gpuID) == None:
        raise FormatError(f"--gpuid '{args.gpuID}' must only contain letters, numbers or '_'")
    if args.minReadLen < 1:
        raise ValueError("--minreadlen must be >= 1")
    if args.concatMinimum < 0:
        raise ValueError("--concat-minimum must be >= 0")
    if args.concatMultiline < 0:
        raise ValueError("--concat-multiline must be >= 0")
    
    # Validate job names, which are derived from -p if --name is not given
    if args.jobName == None:
        args.jobName = args.prefix if args.prefix != None else "dorado"
    validate_job_name(args.jobName)
    
    # Validate resources of each job
    for step in ["p1", "p2", "corr", "s1", "s2", "concat"]:
        walltime = validate_pbs_resources(getattr(args, f"{step}Cpus"), getattr(args, f"{step}Mem"),
                                          getattr(args, f"{step}Walltime"), flagPrefix=f"{step}-")
        setattr(args, f"{step}Walltime", walltime)
    for step in ["s1", "s2"]: # dorado_correct_pipeline.py only accepts memory like '30G'
        if re.match(r"^\d+G$", getattr(args, f"{step}Mem")) == None:
            raise FormatError(f"--{step}-mem '{getattr(args, f'{step}Mem')}' must be in the format '<number>G' (e.g., '30G')")
    
    # Validate job script locations
    args.scriptFiles = {}
    for step, scriptName in DORADO_SCRIPTS.items():
        args.scriptFiles[step] = validate_new_script(os.path.join(args.workDir, scriptName), "--workdir")

def validate_bs(args):
    '''
    Validation for arguments common to all "basespace" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_bs_download(args):
    '''
    Validation for arguments used in "basespace download" mode.
    '''
    validate_program_dir(args, "bsDir", "bs", "--bs")
    if args.projectID != None and re.match(r"^\d+$", args.projectID) == None:
        raise FormatError(f"-i '{args.projectID}' must be a numeric BaseSpace project ID")
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")

def validate_bs_prepare(args):
    '''
    Validation for arguments used in "basespace prepare" mode. Since this mode runs a
    script immediately rather than writing a job script, default locations are checked
    as well as user-specified ones.
    '''
    # Validate program location
    args.varScriptDir = os.path.abspath(args.varScriptDir)
    requiredFile = os.path.join(args.varScriptDir, "Fasta_related", "concatenate_trimmed_lanes.py")
    if not os.path.isfile(requiredFile):
        raise FileNotFoundError(f"Fasta_related/concatenate_trimmed_lanes.py not found in --varscripts {args.varScriptDir}")
    
    # Validate input and output locations
    args.readsDir = os.path.abspath(args.readsDir)
    if not os.path.isdir(args.readsDir):
        raise DirectoryNotFoundError(f"Reads directory (-d {args.readsDir}) does not exist!")
    args.outputFileName = validate_new_script(args.outputFileName, "-o")
    
    # Validate parameters
    args.rpid = validate_rpid(args.rpid)
    if args.readsSuffix == "":
        raise ValueError("-rs must be given a non-empty value")
    if args.laneIdentifier == "":
        raise ValueError("--laneIdentifier must be given a non-empty value")

def validate_bs_archive(args):
    '''
    Validation for arguments used in "basespace archive" mode.
    '''
    if args.prefix != None and re.match(r"^\S+$", args.prefix) == None:
        raise FormatError(f"-p '{args.prefix}' must not contain whitespace")
    if args.readsDir == "":
        raise ValueError("--readsdir must be given a non-empty value")

def validate_ratio(value, flag):
    '''
    Parameters:
        value -- a float which should be in the range 0.0 -> 1.0.
        flag -- a string of the argument flag for error messages.
    '''
    if not 0 <= value <= 1:
        raise ValueError(f"{flag} must be a ratio in the range 0.0 -> 1.0")

def validate_cl(args):
    '''
    Validation for arguments common to all "clustering" mode commands (i.e., those
    added by add_clustering_arguments() in psst.py).
    '''
    validate_input_file(args, "fastaFile", "-i", "Input FASTA file")
    validate_program_dir(args, "varScriptDir", "run_clustering.py", "--varscripts")
    validate_ratio(args.identity, "--identity")
    
    # Validate output file, which run_clustering.py will not overwrite
    if args.outputFile != None:
        args.outputFile = os.path.abspath(args.outputFile)
        if os.path.exists(args.outputFile):
            raise FileExistsError(f"Output file (-o {args.outputFile}) already exists!")
        if not os.path.isdir(os.path.dirname(args.outputFile)):
            raise DirectoryNotFoundError(f"Output file (-o {args.outputFile}) is in a directory that does not exist!")

def validate_cl_mmseqs2(args):
    '''
    Validation for arguments used in "clustering mmseqs2" mode.
    '''
    validate_program_dir(args, "mmseqsDir", "mmseqs", "--mmseqs")
    validate_ratio(args.coverage, "--coverage")
    if args.evalue < 0:
        raise ValueError("--evalue must be >= 0")
    if args.steps < 1:
        raise ValueError("--steps must be >= 1")
    
    # Validate tmp directory, derived from the output file name if not given
    if args.tmpDir == None and args.outputFile != None:
        args.tmpDir = os.path.splitext(args.outputFile)[0] + "_tmp"
    if args.tmpDir != None:
        args.tmpDir = os.path.abspath(args.tmpDir)

def validate_cl_cdhit(args):
    '''
    Validation for arguments used in "clustering cdhit" mode.
    '''
    validate_program_dir(args, "cdhitDir", "cd-hit", "--cdhit")
    if args.cdhitMem < 1000:
        raise ValueError("--cdhit-mem must be >= 1000 (MB)")
    validate_ratio(args.shorterCov, "--shorter-coverage")
    validate_ratio(args.longerCov, "--longer-coverage")
    
    # Validate the .clstr output file, which run_clustering.py will not overwrite
    if args.outputFile != None and os.path.exists(args.outputFile + ".clstr"):
        raise FileExistsError(f"Output file ({args.outputFile}.clstr) already exists!")
