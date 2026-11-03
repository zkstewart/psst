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
from .calling import BAM_INDEX_SUFFIXES, get_bam_samples
from .mapping import BWA_INDEX_SUFFIXES
from .config import CONFIG, RPID_REGEX
from .reads import get_read_files

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

def parse_mem_gb(mem):
    '''
    Parameters:
        mem -- a string of a validated PBS memory value (e.g., '70G' or '70gb').
    Returns:
        memGB -- a float of the memory in gigabytes (GiB).
    '''
    match = re.match(r"^(\d+)([kmgt]?)b?$", mem, re.IGNORECASE)
    return int(match.group(1)) / {"": 1024**3, "k": 1024**2, "m": 1024, "g": 1, "t": 1/1024}[match.group(2).lower()]

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

def validate_h(args):
    '''
    Validation for arguments common to all "hic" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_h_contacts(args):
    '''
    Validation for arguments used in "hic contacts" mode.
    '''
    samples = validate_bwa_inputs(args, False)
    if samples != None:
        print(f"Found {len(samples)} read pair(s) to map as one Hi-C library")
    validate_program_dir(args, "threeddnaDir", os.path.join("visualize", "run-assembly-visualizer.sh"), "--3ddna")
    
    # Validate the output prefix, derived from the reference genome if not given
    if args.prefix == None and args.genome != None:
        args.prefix = re.sub(r"\.(fasta|fa|fna|fas)$", "", os.path.basename(args.genome), flags=re.IGNORECASE)
    if args.prefix != None and re.match(r"^[^\s/]+$", args.prefix) == None:
        raise FormatError(f"-p '{args.prefix}' must not contain whitespace or '/'")
    
    # Validate behavioural parameters
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")
    if args.mapq < 0:
        raise ValueError("--mapq must be >= 0")
    if re.match(r"^\S+$", args.condaEnv) == None:
        raise FormatError(f"--conda '{args.condaEnv}' must be non-empty and not contain whitespace")
    validate_sort_mem(args.sortMem)
    
    # 3D-DNA's juicebox_tools.sh requests a fixed 48 GB Java heap
    args.walltime = validate_pbs_resources(args.cpus, args.mem, args.walltime) # validates --mem early
    if parse_mem_gb(args.mem) < 56:
        print(f"WARNING: --mem {args.mem} is likely too little for 3D-DNA, which runs Java with a 48 GB heap")

def validate_m(args):
    '''
    Validation for arguments common to all "map" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_reads(args, singleEnd):
    '''
    Validation for the reads arguments shared by "map" modes (i.e., those added by
    add_reads_arguments() in psst.py). If a reads directory was given, its read files
    are located and pairs checked.
    
    Parameters:
        args -- an argparse.Namespace containing reads arguments.
        singleEnd -- a boolean indicating whether reads are single-end rather than paired.
    Returns:
        samples -- a list of samples as returned by get_read_files(), or None if no
                   reads directory was given.
    '''
    if args.readsSuffix == "":
        raise ValueError("-rs must be given a non-empty value")
    if args.readsDir == None:
        return None
    args.readsDir = os.path.abspath(args.readsDir)
    if not os.path.isdir(args.readsDir):
        raise DirectoryNotFoundError(f"Reads directory (-d {args.readsDir}) does not exist!")
    return get_read_files(args.readsDir, args.readsSuffix, singleEnd)

def validate_bwa_inputs(args, singleEnd):
    '''
    Validation for the reference genome and reads arguments shared by modes which
    run bwa mem (i.e., those added by add_bwa_arguments() in psst.py).
    
    Parameters:
        args -- an argparse.Namespace containing bwa input arguments.
        singleEnd -- a boolean indicating whether reads are single-end rather than paired.
    Returns:
        samples -- a list of samples as returned by validate_reads().
    '''
    # Validate the reference genome and its bwa index, which may be built by a job
    # this one depends on
    validate_input_file(args, "genome", "-r", "Reference genome FASTA file")
    if "genome" in args.userSet:
        missing = [ s for s in BWA_INDEX_SUFFIXES if not os.path.isfile(args.genome + s) ]
        if missing != [] and args.afterok:
            print(f"WARNING: reference genome (-r {args.genome}) is not indexed yet; " +
                  "the --afterok job(s) must index it")
        elif missing != []:
            raise FileNotFoundError(f"Reference genome (-r {args.genome}) is not indexed; " +
                                    f"run 'bwa index' on it first ({args.genome}{missing[0]} not found)")
    return validate_reads(args, singleEnd)

def validate_sort_mem(sortMem):
    '''
    Parameters:
        sortMem -- a string of the memory to give a sort command (e.g., '2G').
    '''
    if re.match(r"^\d+[KMG]$", sortMem) == None:
        raise FormatError(f"--sort-mem '{sortMem}' must be in the format '<number>K/M/G' (e.g., '2G')")

def validate_m_short(args):
    '''
    Validation for arguments used in "map short" mode. If a reads directory was given,
    its samples are counted to determine the array job size.
    '''
    if args.hic and args.singleEnd:
        raise ValueError("--hic and --singleEnd are mutually exclusive, since Hi-C reads are paired")
    samples = validate_bwa_inputs(args, args.singleEnd)
    args.numSamples = None
    if samples != None:
        args.numSamples = len(samples)
        print(f"Found {args.numSamples} {'single-end' if args.singleEnd else 'paired'} sample(s) to map")
    
    # Validate behavioural parameters
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")
    if re.match(r"^[^\s@]+$", args.library) == None:
        raise FormatError(f"--library '{args.library}' must be non-empty and not contain whitespace or '@'")
    validate_sort_mem(args.sortMem)

def validate_m_index(args):
    '''
    Validation for arguments common to all "map index" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_m_index_bwa(args):
    '''
    Validation for arguments used in "map index bwa" mode.
    '''
    validate_input_file(args, "genome", "-r", "Reference genome FASTA file")
    if "genome" in args.userSet:
        existing = [ args.genome + s for s in BWA_INDEX_SUFFIXES if os.path.exists(args.genome + s) ]
        if existing != []:
            raise FileExistsError(f"Reference genome (-r {args.genome}) is already indexed; " +
                                  f"remove the existing index files (e.g., {existing[0]}) to rebuild it")

def validate_m_index_salmon(args):
    '''
    Validation for arguments used in "map index salmon" mode.
    '''
    validate_program_dir(args, "salmonDir", "salmon", "--salmon")
    validate_input_file(args, "transcripts", "-t", "Transcripts FASTA file")
    validate_input_file(args, "genome", "--genome", "Decoy genome FASTA file")
    
    # Validate the index directory, derived from the transcripts file if not given
    if args.indexDir == None and args.transcripts != None:
        args.indexDir = f"{args.transcripts}_salmon_index"
    if args.indexDir != None:
        args.indexDir = os.path.abspath(args.indexDir)
        if os.path.exists(args.indexDir):
            raise FileExistsError(f"Salmon index (-o {args.indexDir}) already exists!")
        if not os.path.isdir(os.path.dirname(args.indexDir)):
            raise DirectoryNotFoundError(f"Salmon index (-o {args.indexDir}) is in a directory that does not exist!")
    
    # Validate parameters; salmon requires an odd k-mer length no greater than 31
    if args.kmer < 1 or args.kmer > 31 or args.kmer % 2 == 0:
        raise ValueError("--kmer must be an odd number from 1 to 31")
    
    # Validate the RAM limit, which salmon otherwise sizes from the whole node's memory
    args.walltime = validate_pbs_resources(args.cpus, args.mem, args.walltime) # validates --mem early
    memGB = parse_mem_gb(args.mem)
    if args.ramLimit == None:
        args.ramLimit = max(1, int(memGB * 0.75))
    if args.ramLimit < 1:
        raise ValueError("--ram-limit must be >= 1")
    if args.ramLimit >= memGB:
        raise ValueError(f"--ram-limit {args.ramLimit} must be less than --mem {args.mem}")

def validate_m_index_star(args):
    '''
    Validation for arguments used in "map index star" mode.
    '''
    # Validate input files, which STAR cannot read if compressed
    validate_input_file(args, "genome", "-r", "Reference genome FASTA file")
    validate_input_file(args, "gtf", "--gtf", "Gene annotation GTF file")
    for flag, fileName in [("-r", args.genome), ("--gtf", args.gtf)]:
        if fileName != None and fileName.endswith(".gz"):
            raise ValueError(f"{flag} {fileName} must be uncompressed for STAR")
    if args.gtf != None and re.search(r"\.gff3?$", args.gtf, re.IGNORECASE) != None:
        raise ValueError(f"--gtf {args.gtf} appears to be a GFF3; convert it to GTF first " +
                         "(e.g., agat_convert_sp_gff2gtf.pl --gff input.gff3 -o output.gtf)")
    
    # Validate the index directory, derived from the genome file if not given
    if args.indexDir == None and args.genome != None:
        args.indexDir = f"{args.genome}_STAR_index"
    if args.indexDir != None:
        args.indexDir = os.path.abspath(args.indexDir)
        if os.path.exists(args.indexDir):
            raise FileExistsError(f"STAR index (-o {args.indexDir}) already exists!")
        if not os.path.isdir(os.path.dirname(args.indexDir)):
            raise DirectoryNotFoundError(f"STAR index (-o {args.indexDir}) is in a directory that does not exist!")
    
    # Validate parameters
    if args.overhang < 1:
        raise ValueError("--overhang must be >= 1")
    if re.match(r"^\S+$", args.starModule) == None:
        raise FormatError(f"--module '{args.starModule}' must be non-empty and not contain whitespace")
    
    # Set STAR's RAM limit from --mem, leaving headroom for the rest of the process
    args.walltime = validate_pbs_resources(args.cpus, args.mem, args.walltime) # validates --mem early
    args.ramLimit = int(parse_mem_gb(args.mem) * 0.9 * 1024**3)

def validate_m_rna(args):
    '''
    Validation for arguments common to all "map rna" mode commands.
    '''
    pass # no specific validation needed for this mode

def validate_m_rna_salmon(args):
    '''
    Validation for arguments used in "map rna salmon" mode. If a reads directory was
    given, its samples are counted to determine the array job size.
    '''
    validate_program_dir(args, "salmonDir", "salmon", "--salmon")
    validate_input_file(args, "geneMap", "-g", "Gene map file")
    
    # Validate the salmon index, which may be built by a job this one depends on
    if "indexDir" in args.userSet:
        args.indexDir = os.path.abspath(args.indexDir)
        if not os.path.isfile(os.path.join(args.indexDir, "info.json")) and args.afterok:
            print(f"WARNING: salmon index (-i {args.indexDir}) does not exist yet; " +
                  "the --afterok job(s) must build it")
        elif not os.path.isfile(os.path.join(args.indexDir, "info.json")):
            raise FileNotFoundError(f"Salmon index (-i {args.indexDir}) does not exist or is incomplete!")
    
    # Validate reads and count samples
    samples = validate_reads(args, args.singleEnd)
    args.numSamples = None
    if samples != None:
        args.numSamples = len(samples)
        print(f"Found {args.numSamples} {'single-end' if args.singleEnd else 'paired'} sample(s) to quantify")
    
    # Validate the library type; paired types start with I, O or M (e.g., ISR), and
    # single-end types do not (e.g., SR)
    if args.libType != "A":
        if re.match(r"^[IOM]?(U|S[FR])$", args.libType) == None:
            raise FormatError(f"--libtype '{args.libType}' is not a valid salmon library type (e.g., A, IU, ISR, U, SR)")
        if args.singleEnd == (args.libType[0] in "IOM"):
            raise ValueError(f"--libtype '{args.libType}' is a {'paired' if args.libType[0] in 'IOM' else 'single-end'} " +
                             f"library type, but reads are {'single-end' if args.singleEnd else 'paired'}")
    
    # Validate bias correction, defaulting to GC bias correction for paired reads
    if args.bias == None:
        args.bias = "none" if args.singleEnd else "gc"
    if args.singleEnd and args.bias in ["gc", "both"]:
        raise ValueError("--bias gc requires paired reads, since it models fragment GC content")
    
    # Validate behavioural parameters
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")
    if args.bootstraps < 0:
        raise ValueError("--bootstraps must be >= 0")

def validate_m_rna_star(args):
    '''
    Validation for arguments used in "map rna star" mode. If a reads directory was
    given, its samples are counted to determine the array job size.
    '''
    if re.match(r"^\S+$", args.starModule) == None:
        raise FormatError(f"--module '{args.starModule}' must be non-empty and not contain whitespace")
    
    # Validate the STAR index, which may be built by a job this one depends on
    if "indexDir" in args.userSet:
        args.indexDir = os.path.abspath(args.indexDir)
        if not os.path.isfile(os.path.join(args.indexDir, "SAindex")) and args.afterok:
            print(f"WARNING: STAR index (-i {args.indexDir}) does not exist yet; " +
                  "the --afterok job(s) must build it")
        elif not os.path.isfile(os.path.join(args.indexDir, "SAindex")):
            raise FileNotFoundError(f"STAR index (-i {args.indexDir}) does not exist or is incomplete!")
    
    # Validate reads and count samples
    samples = validate_reads(args, args.singleEnd)
    args.numSamples = None
    if samples != None:
        args.numSamples = len(samples)
        print(f"Found {args.numSamples} {'single-end' if args.singleEnd else 'paired'} sample(s) to map")
    
    # Validate behavioural parameters
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")

def validate_m_nanopore(args):
    '''
    Validation for arguments used in "map nanopore" mode. If a reads directory was
    given, its samples are counted to determine the array job size.
    '''
    validate_program_dir(args, "minimap2Dir", "minimap2", "--minimap2")
    validate_program_dir(args, "chopperDir", "chopper", "--chopper")
    validate_input_file(args, "genome", "-r", "Reference genome FASTA file")
    
    # Validate reads and count samples
    samples = validate_reads(args, True)
    args.numSamples = None
    if samples != None:
        args.numSamples = len(samples)
        print(f"Found {args.numSamples} sample(s) to map")
    
    # Validate read filtering parameters
    for flag, value in [("--minlength", args.minLength), ("--maxlength", args.maxLength), ("--minqual", args.minQual)]:
        if value < 0:
            raise ValueError(f"{flag} must be >= 0")
    if args.maxLength > 0 and args.minLength > args.maxLength:
        raise ValueError("--minlength must not be greater than --maxlength")
    
    # Validate behavioural parameters
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")
    if re.match(r"^[^\s@]+$", args.library) == None:
        raise FormatError(f"--library '{args.library}' must be non-empty and not contain whitespace or '@'")
    validate_sort_mem(args.sortMem)

def validate_ratio(value, flag):
    '''
    Parameters:
        value -- a float which should be in the range 0.0 -> 1.0.
        flag -- a string of the argument flag for error messages.
    '''
    if not 0 <= value <= 1:
        raise ValueError(f"{flag} must be a ratio in the range 0.0 -> 1.0")

def validate_bams(args):
    '''
    Validation for the BAM directory of "call" modes. If given, its BAM files must exist
    and be indexed, unless the --afterok job(s) will produce them; samples shared by
    multiple BAM files are reported, since bcftools treats them as one sample.

    Parameters:
        args -- an argparse.Namespace containing args.bamDir and args.bamSuffix.
    '''
    if args.bamSuffix == "":
        raise ValueError("-bs must be given a non-empty value")
    if args.bamDir == None:
        return
    args.bamDir = os.path.abspath(args.bamDir)
    bamFiles = []
    if os.path.isdir(args.bamDir):
        bamFiles = sorted([ os.path.join(args.bamDir, f) for f in os.listdir(args.bamDir) if f.endswith(args.bamSuffix) ])
    if bamFiles == []:
        if args.afterok:
            print(f"WARNING: no BAM files ending with '{args.bamSuffix}' are in -b {args.bamDir} yet; " +
                  "the --afterok job(s) must create them")
            return
        raise FileNotFoundError(f"No BAM files ending with '{args.bamSuffix}' were found in -b {args.bamDir}")

    # Check that BAMs are indexed, which may be done by a job this one depends on
    unindexed = [ f for f in bamFiles
                  if not any(os.path.isfile(f + s) or os.path.isfile(re.sub(r"\.bam$", s, f))
                             for s in BAM_INDEX_SUFFIXES) ]
    if unindexed != [] and args.afterok:
        print(f"WARNING: {len(unindexed)} BAM file(s) are not indexed yet; the --afterok job(s) must index them")
    elif unindexed != []:
        raise FileNotFoundError(f"{len(unindexed)} BAM file(s) are not indexed; run 'samtools index' on " +
                                f"them first (e.g., {unindexed[0]})")

    # Report samples, warning about those which would be merged
    fileSamples = { f: get_bam_samples(f) for f in bamFiles }
    sampleFiles = {}
    for bamFile, samples in fileSamples.items():
        for sample in (samples or []):
            sampleFiles.setdefault(sample, []).append(os.path.basename(bamFile))
    print(f"Found {len(bamFiles)} BAM file(s) with {len(sampleFiles)} sample name(s) in their read groups")
    for sample, files in sampleFiles.items():
        if len(files) > 1:
            print(f"WARNING: sample '{sample}' is named in {len(files)} BAM files ({', '.join(files)}), " +
                  "which bcftools will call as one sample")
    noRG = [ os.path.basename(f) for f, samples in fileSamples.items() if samples == [] ]
    if noRG != []:
        print(f"WARNING: {len(noRG)} BAM file(s) have no read group samples, so bcftools will name " +
              f"samples by file name (e.g., {noRG[0]})")

def validate_c(args):
    '''
    Validation for arguments common to all "call" mode commands (i.e., those added by
    add_call_arguments() in psst.py).
    '''
    # Validate the reference genome, whose .fai index determines the contigs to call
    validate_input_file(args, "genome", "-r", "Reference genome FASTA file")
    numContigs = None
    if "genome" in args.userSet:
        if not os.path.isfile(args.genome + ".fai"):
            raise FileNotFoundError(f"Reference genome (-r {args.genome}) is not indexed; " +
                                    f"run 'samtools faidx' on it first ({args.genome}.fai not found)")
        with open(args.genome + ".fai", "r") as fileIn:
            numContigs = sum(1 for line in fileIn if line.strip() != "")
    validate_bams(args)

    # Validate the number of chunks, which can't usefully exceed the number of contigs
    if args.chunks < 1:
        raise ValueError("--chunks must be >= 1")
    if numContigs != None and args.chunks > numContigs:
        print(f"Reducing --chunks from {args.chunks} to the genome's {numContigs} contig(s)")
        args.chunks = numContigs

    # Validate the output prefix, derived from the reference genome if not given
    if args.prefix == None and args.genome != None:
        args.prefix = re.sub(r"\.(fasta|fa|fna|fas)$", "", os.path.basename(args.genome), flags=re.IGNORECASE)
    if args.prefix != None and re.match(r"^[^\s/]+$", args.prefix) == None:
        raise FormatError(f"-p '{args.prefix}' must not contain whitespace or '/'")
    if args.outputDir == "":
        raise ValueError("--outdir must be given a non-empty value")
    if args.prefix != None:
        outputFile = os.path.join(args.outputDir, f"{args.prefix}.vcf.gz")
        if os.path.exists(outputFile):
            raise FileExistsError(f"Output file '{os.path.abspath(outputFile)}' already exists!")

    # Validate behavioural parameters
    if args.minMapQ < 0:
        raise ValueError("-q must be >= 0")
    if getattr(args, "minBaseQ", 0) < 0: # not all modes have this argument
        raise ValueError("-Q must be >= 0")
    if args.maxDepth < 1:
        raise ValueError("--max-depth must be >= 1")

    # Validate job names and resources
    validate_job_name(args.jobName)
    for step in ["call", "concat"]:
        walltime = validate_pbs_resources(getattr(args, f"{step}Cpus"), getattr(args, f"{step}Mem"),
                                          getattr(args, f"{step}Walltime"), flagPrefix=f"{step}-")
        setattr(args, f"{step}Walltime", walltime)

    # Validate job script locations
    args.scriptFiles = {
        "call": validate_new_script(f"{args.jobName}.sh", "--name"),
        "concat": validate_new_script(f"{args.jobName}_concat.sh", "--name")
    }

def validate_c_short(args):
    '''
    Validation for arguments used in "call short" mode.
    '''
    pass # all arguments are validated by validate_c()

def validate_c_nanopore(args):
    '''
    Validation for arguments used in "call nanopore" mode.
    '''
    pass # all arguments are validated by validate_c()

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
