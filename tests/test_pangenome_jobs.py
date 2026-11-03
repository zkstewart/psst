# Tests of the "pangenome create" job script's own logic (i.e., writing the seqfile, and
# resuming or refusing to repeat a cactus run), run against stub conda and cactus-pangenome
# (see stubs/cactus).

import os

import pytest

from conftest import STUBS_DIR, run_psst, run_job

@pytest.fixture
def setup(tmp_path):
    '''
    Creates a reference and genome FASTAs, a work directory,
    and an environment with the stubs on PATH.

    Returns:
        (workDir, genomes, env, stubDir) -- the work directory, a dict of FASTA locations
            by file name, a dict of environment variables for run_job(), and the directory
            the stub records its calls in.
    '''
    genomes = {}
    for name in ["ref.fa", "A.hap1.fa", "A.hap2.fa", "B1.fa"]:
        genomes[name] = str(tmp_path / name)
        open(genomes[name], "w").write(">c\nACGT\n")
    workDir = tmp_path / "work"
    workDir.mkdir()
    stubDir = tmp_path / "stub_records"
    stubDir.mkdir()
    env = dict(os.environ)
    env["PATH"] = f"{os.path.join(STUBS_DIR, 'cactus')}{os.pathsep}{env['PATH']}"
    env["STUB_DIR"] = str(stubDir)
    return str(workDir), genomes, env, str(stubDir)

def generate(setup, args=[]):
    workDir, genomes = setup[:2]
    run_psst(["pangenome", "create", "-i", genomes["ref.fa"], genomes["A.hap1.fa"], genomes["A.hap2.fa"],
              genomes["B1.fa"], *args], workDir)

def run(setup, env=None, check=True):
    workDir, _, stubEnv, stubDir = setup
    result = run_job(os.path.join(workDir, "pangenome.sh"), workDir, env=env or stubEnv, check=check)
    logFile = os.path.join(stubDir, "cactus.log")
    calls = [ line.split("|") for line in open(logFile).read().splitlines() ] if os.path.exists(logFile) else []
    return result, calls

def test_create_writes_seqfile_and_runs_cactus(setup):
    workDir, genomes, _, stubDir = setup
    generate(setup)
    result, calls = run(setup)
    outDir = os.path.join(workDir, "cactus_output")
    expected = f"reference\t{genomes['ref.fa']}\nA.1\t{genomes['A.hap1.fa']}\nA.2\t{genomes['A.hap2.fa']}\nB1.0\t{genomes['B1.fa']}\n"
    assert open(os.path.join(outDir, "cactus_samples.txt")).read() == expected
    assert open(os.path.join(stubDir, "seqfile.txt")).read() == expected
    assert len(calls) == 1
    call = calls[0]
    assert call[:2] == ["cactus_output/jobstore", "cactus_output/cactus_samples.txt"]
    for flag, value in [("--outDir", "cactus_output/graph"), ("--outName", "pangenome"), ("--reference", "reference"),
                        ("--giraffe", "filter"), ("--maxCores", "16"), ("--maxMemory", "120G"),
                        ("--binariesMode", "singularity")]:
        assert call[call.index(flag)+1] == value
    assert "--vcf" in call and "--vcfwave" in call and "--restart" not in call
    assert sorted(os.listdir(outDir)) == ["cactus.ok", "cactus_samples.txt", "graph"]
    assert os.listdir(os.path.join(outDir, "graph")) == ["pangenome.d2.gbz"]

def test_create_options(setup):
    generate(setup, ["--giraffe", "clip", "--no-vcfwave", "--binaries", "local", "--mem", "64G", "--cpus", "4",
                     "-p", "pg", "--refname", "myref"])
    result, calls = run(setup)
    call = calls[0]
    for flag, value in [("--outName", "pg"), ("--reference", "myref"), ("--giraffe", "clip"), ("--maxCores", "4"),
                        ("--maxMemory", "64G"), ("--binariesMode", "local")]:
        assert call[call.index(flag)+1] == value
    assert "--vcfwave" not in call
    assert open(os.path.join(setup[0], "cactus_output", "cactus_samples.txt")).read().startswith("myref\t")

def test_create_resumes_and_refuses_to_repeat(setup):
    workDir, _, env, _ = setup
    generate(setup)
    failEnv = dict(env)
    failEnv["STUB_FAIL"] = "1"
    result, calls = run(setup, env=failEnv, check=False)
    assert result.returncode != 0
    assert os.path.isdir(os.path.join(workDir, "cactus_output", "jobstore"))
    assert not os.path.exists(os.path.join(workDir, "cactus_output", "cactus.ok"))

    result, calls = run(setup) # resumes with --restart, reusing the identical seqfile
    assert "--restart" in calls[1]
    assert os.path.exists(os.path.join(workDir, "cactus_output", "cactus.ok"))

    result, calls = run(setup, check=False)
    assert result.returncode != 0
    assert "already holds a completed pangenome" in result.stdout
    assert len(calls) == 2

def test_create_refuses_changed_seqfile(setup):
    workDir, _, env, _ = setup
    generate(setup)
    failEnv = dict(env)
    failEnv["STUB_FAIL"] = "1"
    run(setup, env=failEnv, check=False)
    with open(os.path.join(workDir, "cactus_output", "cactus_samples.txt"), "a") as fileOut:
        fileOut.write("C\t/elsewhere/C.fa\n")
    result, calls = run(setup, check=False)
    assert result.returncode != 0
    assert "differs from the genomes of this job" in result.stdout
    assert len(calls) == 1
