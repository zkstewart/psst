# Tests of the "qc library trim" job script's own logic (i.e., finding each subjob's
# reads, naming outputs, skipping trimmed samples and never leaving partial output under
# a final name), run against stub conda and bbduk.sh (see stubs/bbduk).

import os

import pytest

from conftest import STUBS_DIR, run_psst, run_job, run_array_job

@pytest.fixture
def setup(tmp_path):
    '''
    Creates paired and single-end reads directories, a work directory, and an
    environment with the stubs on PATH.

    Returns:
        (workDir, env, stubDir) -- the work directory, a dict of environment variables
            for run_job(), and the directory the stub records its calls in.
    '''
    for name in ["A_R1.fq.gz", "A_R2.fq.gz", "B.1.fq.gz", "B.2.fq.gz"]:
        (tmp_path / "paired").mkdir(exist_ok=True)
        (tmp_path / "paired" / name).write_text(name)
    for name in ["C.fq.gz", "D.fq.gz"]:
        (tmp_path / "single").mkdir(exist_ok=True)
        (tmp_path / "single" / name).write_text(name)
    workDir = tmp_path / "work"
    workDir.mkdir()
    stubDir = tmp_path / "stub_records"
    stubDir.mkdir()
    env = dict(os.environ)
    env["PATH"] = f"{os.path.join(STUBS_DIR, 'bbduk')}{os.pathsep}{env['PATH']}"
    env["STUB_DIR"] = str(stubDir)
    return str(workDir), env, str(stubDir)

def generate(setup, args):
    workDir = setup[0]
    run_psst(["qc", "library", "trim", *args], workDir)
    return os.path.join(workDir, "trim_reads.sh")

def calls(stubDir):
    logFile = os.path.join(stubDir, "bbduk.log")
    return [ line.split("|") for line in open(logFile).read().splitlines() ] if os.path.exists(logFile) else []

def test_trim_paired(setup, tmp_path):
    workDir, env, stubDir = setup
    script = generate(setup, ["-d", str(tmp_path / "paired"), "--mem", "8G", "--trimq", "10", "--minlen", "30"])
    run_array_job(script, workDir, env=env)
    outDir = os.path.join(workDir, "trimmed")
    assert sorted(os.listdir(outDir)) == ["A.bbduk_stats.txt", "A_1.fq.gz", "A_2.fq.gz",
                                          "B.bbduk_stats.txt", "B_1.fq.gz", "B_2.fq.gz"]
    assert open(os.path.join(outDir, "A_2.fq.gz")).read() == "A_R2.fq.gz" # pairs stay in order
    assert open(os.path.join(outDir, "B_1.fq.gz")).read() == "B.1.fq.gz"
    a = calls(stubDir)[0]
    assert a[0] == "-Xmx6144m"
    for arg in [f"in1={tmp_path}/paired/A_R1.fq.gz", f"in2={tmp_path}/paired/A_R2.fq.gz",
                "out1=trimmed/A.incomplete_1.fq.gz", "out2=trimmed/A.incomplete_2.fq.gz",
                "tpe", "tbo", "ref=adapters", "threads=4", "ktrim=r", "k=23", "mink=11", "hdist=1",
                "qtrim=rl", "trimq=10", "minlength=30", "stats=trimmed/A.bbduk_stats.txt"]:
        assert arg in a

def test_trim_single_end(setup, tmp_path):
    workDir, env, stubDir = setup
    adapters = tmp_path / "my adapters.fa" # a space tests quoting
    adapters.write_text(">a\nACGT\n")
    script = generate(setup, ["-d", str(tmp_path / "single"), "--singleEnd", "--adapters", str(adapters)])
    run_array_job(script, workDir, env=env)
    assert sorted(f for f in os.listdir(os.path.join(workDir, "trimmed")) if f.endswith(".fq.gz")) == ["C.fq.gz", "D.fq.gz"]
    c = calls(stubDir)[0]
    for arg in [f"in={tmp_path}/single/C.fq.gz", "out=trimmed/C.incomplete.fq.gz", f"ref={adapters}"]:
        assert arg in c
    assert "tpe" not in c and "tbo" not in c and not any(a.startswith("in2=") for a in c)

def test_trim_skips_trimmed_samples(setup, tmp_path):
    workDir, env, stubDir = setup
    script = generate(setup, ["-d", str(tmp_path / "paired")])
    run_job(script, workDir, index=1, env=env)
    result = run_job(script, workDir, index=1, env=env)
    assert "skipping A" in result.stdout
    assert len(calls(stubDir)) == 1

def test_trim_refuses_partial_output(setup, tmp_path):
    workDir, env, stubDir = setup
    script = generate(setup, ["-d", str(tmp_path / "paired")])
    os.makedirs(os.path.join(workDir, "trimmed"))
    open(os.path.join(workDir, "trimmed", "A_1.fq.gz"), "w").close()
    result = run_job(script, workDir, index=1, env=env, check=False)
    assert result.returncode != 0
    assert "Only some of" in result.stdout
    assert calls(stubDir) == []

def test_trim_failure_leaves_no_final_output(setup, tmp_path):
    workDir, env, stubDir = setup
    script = generate(setup, ["-d", str(tmp_path / "paired")])
    result = run_job(script, workDir, index=1, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert sorted(f for f in os.listdir(os.path.join(workDir, "trimmed")) if f.endswith(".fq.gz")) == \
        ["A.incomplete_1.fq.gz", "A.incomplete_2.fq.gz"]
    run_job(script, workDir, index=1, env=env) # a rerun replaces the incomplete output
    assert sorted(f for f in os.listdir(os.path.join(workDir, "trimmed")) if f.endswith(".fq.gz")) == \
        ["A_1.fq.gz", "A_2.fq.gz"]

def test_trim_refuses_outdir_as_readsdir(setup, tmp_path):
    workDir, env, stubDir = setup
    script = generate(setup, ["-d", str(tmp_path / "single"), "--singleEnd"])
    text = open(script).read().replace("OUTDIR=trimmed\n", f"OUTDIR={tmp_path}/single/../single\n")
    open(script, "w").write(text) # as if edited by hand after validation
    os.remove(os.path.join(tmp_path, "single", "D.fq.gz")) # so C's output name is free
    result = run_job(script, workDir, env=env, check=False)
    assert result.returncode != 0
    assert "must differ from the reads directory" in result.stdout
    assert calls(stubDir) == []
