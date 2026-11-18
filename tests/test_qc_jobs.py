# Tests of the "qc library" job scripts' own logic (i.e., finding each subjob's reads, naming
# outputs, skipping finished samples and never leaving partial output under a final name),
# run against stubs (see stubs/bbduk and stubs/fastqc, and a stand-in tabulate_fastqc.py).

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

# A stand-in for Various_scripts' tabulate_fastqc.py, which records its arguments and
# writes the table, failing if $STUB_FAIL is set
TABULATE_STUB = """import os, sys
open(os.path.join(os.environ["STUB_DIR"], "tabulate.log"), "a").write("|".join(sys.argv[1:]) + "\\n")
open(sys.argv[sys.argv.index("-o") + 1], "w").write("table")
sys.exit(1 if os.environ.get("STUB_FAIL") else 0)
"""

@pytest.fixture
def fastqc_setup(setup, tmp_path):
    '''
    Extends setup with the fastqc stubs on PATH and a Various_scripts directory holding a
    stand-in tabulate_fastqc.py.

    Returns:
        (workDir, env, stubDir, varScriptDir) -- as for setup, plus the Various_scripts
            directory to give --varscripts.
    '''
    workDir, env, stubDir = setup
    varScriptDir = tmp_path / "Various_scripts"
    (varScriptDir / "QC" / "reads").mkdir(parents=True)
    (varScriptDir / "QC" / "reads" / "tabulate_fastqc.py").write_text(TABULATE_STUB)
    env = dict(env)
    env["PATH"] = f"{os.path.join(STUBS_DIR, 'fastqc')}{os.pathsep}{env['PATH']}"
    return workDir, env, stubDir, str(varScriptDir)

def generate_fastqc(fastqc_setup, args):
    workDir, _, _, varScriptDir = fastqc_setup
    run_psst(["qc", "library", "fastqc", "--varscripts", varScriptDir, *args], workDir)
    return os.path.join(workDir, "fastqc.sh"), os.path.join(workDir, "fastqc_tabulate.sh")

def fastqc_calls(stubDir, tool="fastqc"):
    logFile = os.path.join(stubDir, f"{tool}.log")
    return [ line.split("|") for line in open(logFile).read().splitlines() ] if os.path.exists(logFile) else []

def test_fastqc_then_tabulate(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired"), "--fastqc-cpus", "2"])
    run_array_job(fastqcScript, workDir, env=env)
    outDir = os.path.join(workDir, "fastqc")
    assert sorted(os.listdir(outDir)) == ["A", "B"]
    assert sorted(os.listdir(os.path.join(outDir, "A"))) == ["A_R1_fastqc.html", "A_R1_fastqc.zip",
                                                             "A_R2_fastqc.html", "A_R2_fastqc.zip"]
    assert fastqc_calls(stubDir)[0] == ["--outdir", "fastqc/A.incomplete", "--threads", "2",
                                        f"{tmp_path}/paired/A_R1.fq.gz", f"{tmp_path}/paired/A_R2.fq.gz"]
    run_job(tabulateScript, workDir, env=env)
    assert fastqc_calls(stubDir, "tabulate") == [["-i", "fastqc", "-o", "fastqc/fastqc_summary.incomplete.xlsx"]]
    assert sorted(os.listdir(outDir)) == ["A", "B", "fastqc_summary.xlsx"]
    result = run_job(tabulateScript, workDir, env=env, check=False) # refuses to overwrite the table
    assert result.returncode != 0 and "already exists" in result.stdout

def test_fastqc_single_end(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "single"), "--singleEnd",
                                                                  "-o", "summary.xlsx"])
    run_array_job(fastqcScript, workDir, env=env)
    assert sorted(os.listdir(os.path.join(workDir, "fastqc", "C"))) == ["C_fastqc.html", "C_fastqc.zip"]
    run_job(tabulateScript, workDir, env=env)
    assert os.path.exists(os.path.join(workDir, "summary.xlsx"))

def test_fastqc_skips_checked_samples(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, _ = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    run_job(fastqcScript, workDir, index=1, env=env)
    result = run_job(fastqcScript, workDir, index=1, env=env)
    assert "skipping A" in result.stdout
    assert len(fastqc_calls(stubDir)) == 1

def test_fastqc_failure_leaves_no_final_output(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, _ = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    result = run_job(fastqcScript, workDir, index=1, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert os.listdir(os.path.join(workDir, "fastqc")) == ["A.incomplete"]
    run_job(fastqcScript, workDir, index=1, env=env) # a rerun replaces the incomplete output
    assert os.listdir(os.path.join(workDir, "fastqc")) == ["A"]
    assert len(os.listdir(os.path.join(workDir, "fastqc", "A"))) == 4

@pytest.mark.parametrize("prepare,message", [
    (lambda outDir: os.makedirs(os.path.join(outDir, "A")), "fastqc/B was not found"),
    (lambda outDir: [ os.makedirs(os.path.join(outDir, d)) for d in ["A", "B", "C.incomplete"] ], "contains incomplete FastQC output"),
])
def test_tabulate_checks_reports(prepare, message, fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    _, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    prepare(os.path.join(workDir, "fastqc"))
    result = run_job(tabulateScript, workDir, env=env, check=False)
    assert result.returncode != 0
    assert message in result.stdout
    assert fastqc_calls(stubDir, "tabulate") == []

def test_tabulate_failure_leaves_no_table(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    run_array_job(fastqcScript, workDir, env=env)
    result = run_job(tabulateScript, workDir, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert not os.path.exists(os.path.join(workDir, "fastqc", "fastqc_summary.xlsx"))
    run_job(tabulateScript, workDir, env=env) # a rerun replaces the incomplete table
    assert sorted(os.listdir(os.path.join(workDir, "fastqc"))) == ["A", "B", "fastqc_summary.xlsx"]
