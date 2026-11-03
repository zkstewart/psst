# Shared fixtures and helpers for psst's tests. Tests marked 'integration' run generated
# job scripts with real tools on simulated data; they are skipped by default (see
# pytest.ini) and run with 'pytest -m integration'.

import os, re, shutil, subprocess, sys

import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
PSST = os.path.join(os.path.dirname(TESTS_DIR), "psst.py")
STUBS_DIR = os.path.join(TESTS_DIR, "stubs")
RPID = "ABCDEF1234"

def run_psst(args, cwd, env=None, check=True):
    '''
    Runs psst.py with arguments from within a directory.

    Parameters:
        args -- a list of argument strings (e.g., ["call", "short"]); --rpid is added.
        cwd -- a string of the directory to run psst from.
        env -- OPTIONAL; a dict of environment variables; default == os.environ.
        check -- OPTIONAL; a boolean indicating whether to fail if psst fails.
    Returns:
        result -- a subprocess.CompletedProcess with text stdout and stderr.
    '''
    result = subprocess.run([sys.executable, PSST, *args, "--rpid", RPID], cwd=cwd, env=env,
                            capture_output=True, text=True)
    if check and result.returncode != 0:
        raise AssertionError(f"psst {' '.join(args)} failed:\n{result.stdout}\n{result.stderr}")
    return result

def array_size(scriptFile):
    '''
    Returns:
        size -- an int of the job's '#PBS -J 1-N' array size, or 1 if it is not an array.
    '''
    with open(scriptFile) as fileIn:
        match = re.search(r"^#PBS -J 1-(\d+)$", fileIn.read(), re.MULTILINE)
    return int(match.group(1)) if match else 1

def run_job(scriptFile, workDir, index=None, env=None, check=True):
    '''
    Runs a job script as PBS would, from workDir, for one array index.

    Returns:
        result -- a subprocess.CompletedProcess with text stdout and stderr combined.
    '''
    jobEnv = dict(env if env != None else os.environ)
    jobEnv["PBS_O_WORKDIR"] = workDir # the job header runs 'cd $PBS_O_WORKDIR'
    if index != None:
        jobEnv["PBS_ARRAY_INDEX"] = str(index)
    else:
        jobEnv.pop("PBS_ARRAY_INDEX", None)
    result = subprocess.run(["bash", scriptFile], cwd=workDir, env=jobEnv,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if check and result.returncode != 0:
        raise AssertionError(f"{scriptFile} (index {index}) failed:\n{result.stdout}")
    return result

def run_array_job(scriptFile, workDir, env=None):
    '''
    Runs every subjob of a job script, as sized by its '#PBS -J' directive.
    '''
    size = array_size(scriptFile)
    for index in range(1, size + 1):
        run_job(scriptFile, workDir, index if size > 1 else None, env=env)

def require_tools(*tools):
    '''
    Skips the calling test if any tool is not on PATH.
    '''
    missing = [ t for t in tools if shutil.which(t) == None ]
    if missing != []:
        pytest.skip(f"requires {', '.join(missing)} on PATH")

def tool_dir(tool):
    return os.path.dirname(shutil.which(tool))

@pytest.fixture
def fake_qsub(tmp_path):
    '''
    Puts a fake qsub on PATH which logs its arguments and reports sequential job IDs.

    Returns:
        (env, logFile) -- a dict of environment variables for run_psst(), and a string of
                          the file which logs each qsub call as one line.
    '''
    binDir = tmp_path / "fakebin"
    binDir.mkdir()
    logFile = tmp_path / "qsub.log"
    qsub = binDir / "qsub"
    qsub.write_text("#!/bin/bash\n" +
                    f'echo "$*" >> "{logFile}"\n' +
                    f'echo "$(( $(wc -l < "{logFile}") + 100 ))[].pbs"\n')
    qsub.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{binDir}{os.pathsep}{env['PATH']}"
    return env, str(logFile)

@pytest.fixture
def stub_bcftools(tmp_path):
    '''
    Puts a stub bcftools on PATH (see stubs/bcftools) which records the regions each
    calling subjob was given and writes placeholder outputs.

    Returns:
        (env, stubDir) -- a dict of environment variables for run_job(), and a string of
                          the directory the stub records regions in.
    '''
    stubDir = tmp_path / "stub_records"
    stubDir.mkdir()
    env = dict(os.environ)
    env["PATH"] = f"{STUBS_DIR}{os.pathsep}{env['PATH']}"
    env["STUB_DIR"] = str(stubDir)
    return env, str(stubDir)
