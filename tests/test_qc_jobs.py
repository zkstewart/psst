# Tests of the "qc library" job scripts' own logic (i.e., finding each subjob's reads, naming
# outputs, skipping finished samples and never leaving partial output under a final name),
# run against stubs (see stubs/bbduk, stubs/fastqc, stubs/kraken2 and stubs/smudgeplot, and stand-in
# Various_scripts).

import gzip, os

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

# A stand-in for Various_scripts' kreport_reformat.py, which records its arguments in
# $STUB_DIR/tabulate.log and writes its -o file, failing if $STUB_FAIL is set
TABULATE_STUB = """import os, sys
open(os.path.join(os.environ["STUB_DIR"], "tabulate.log"), "a").write("|".join(sys.argv[1:]) + "\\n")
open(sys.argv[sys.argv.index("-o") + 1], "w").write("table")
sys.exit(1 if os.environ.get("STUB_FAIL") else 0)
"""

@pytest.fixture
def fastqc_setup(setup, tmp_path):
    '''
    Extends setup with the fastqc stubs on PATH, with the stub annotarium (see
    stubs/annotarium) to tabulate reports.

    Returns:
        (workDir, env, stubDir, annotariumDir) -- as for setup, plus the stub annotarium
            directory to give --annotarium.
    '''
    workDir, env, stubDir = setup
    env = dict(env)
    env["PATH"] = f"{os.path.join(STUBS_DIR, 'fastqc')}{os.pathsep}{env['PATH']}"
    return workDir, env, stubDir, os.path.join(STUBS_DIR, "annotarium")

def generate_fastqc(fastqc_setup, args):
    workDir, _, _, annotariumDir = fastqc_setup
    run_psst(["qc", "library", "fastqc", "--annotarium", annotariumDir, *args], workDir)
    return os.path.join(workDir, "fastqc.sh"), os.path.join(workDir, "fastqc_tabulate.sh")

def stub_calls(stubDir, tool="fastqc"):
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
    assert stub_calls(stubDir)[0] == ["--outdir", "fastqc/A.incomplete", "--threads", "2",
                                        f"{tmp_path}/paired/A_R1.fq.gz", f"{tmp_path}/paired/A_R2.fq.gz"]
    run_job(tabulateScript, workDir, env=env)
    assert stub_calls(stubDir, "annotarium") == [["qc", "fastqc", "tabulate", "-i", "fastqc/A", "fastqc/B",
                                                  "-o", "fastqc/fastqc_summary.incomplete.xlsx"]]
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
    assert len(stub_calls(stubDir)) == 1

def test_fastqc_failure_leaves_no_final_output(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, _ = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    result = run_job(fastqcScript, workDir, index=1, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert os.listdir(os.path.join(workDir, "fastqc")) == ["A.incomplete"]
    run_job(fastqcScript, workDir, index=1, env=env) # a rerun replaces the incomplete output
    assert os.listdir(os.path.join(workDir, "fastqc")) == ["A"]
    assert len(os.listdir(os.path.join(workDir, "fastqc", "A"))) == 4

def test_tabulate_checks_reports(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    _, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    os.makedirs(os.path.join(workDir, "fastqc", "A")) # as if only sample A were checked
    result = run_job(tabulateScript, workDir, env=env, check=False)
    assert result.returncode != 0
    assert "fastqc/B was not found" in result.stdout
    assert stub_calls(stubDir, "annotarium") == []

def test_tabulate_only_given_samples(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    _, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    for d in ["A", "B", "C", "D.incomplete"]: # C and D aren't samples of -d
        os.makedirs(os.path.join(workDir, "fastqc", d))
    run_job(tabulateScript, workDir, env=env)
    call = stub_calls(stubDir, "annotarium")[0]
    assert call[call.index("-i") + 1:call.index("-o")] == ["fastqc/A", "fastqc/B"]

def test_tabulate_failure_leaves_no_table(fastqc_setup, tmp_path):
    workDir, env, stubDir, _ = fastqc_setup
    fastqcScript, tabulateScript = generate_fastqc(fastqc_setup, ["-d", str(tmp_path / "paired")])
    run_array_job(fastqcScript, workDir, env=env)
    result = run_job(tabulateScript, workDir, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert not os.path.exists(os.path.join(workDir, "fastqc", "fastqc_summary.xlsx"))
    run_job(tabulateScript, workDir, env=env) # a rerun replaces the incomplete table
    assert sorted(os.listdir(os.path.join(workDir, "fastqc"))) == ["A", "B", "fastqc_summary.xlsx"]

@pytest.fixture
def kraken2_setup(setup, tmp_path):
    '''
    Extends setup with the kraken2 stubs on PATH, a database directory, and a
    Various_scripts directory holding a stand-in kreport_reformat.py.

    Returns:
        (workDir, env, stubDir, args) -- as for setup, plus the psst arguments giving the
            database and Various_scripts locations.
    '''
    workDir, env, stubDir = setup
    varScriptDir = tmp_path / "Various_scripts"
    (varScriptDir / "kraken").mkdir(parents=True)
    (varScriptDir / "kraken" / "kreport_reformat.py").write_text(TABULATE_STUB)
    db = tmp_path / "k2db"
    db.mkdir()
    (db / "hash.k2d").write_text("")
    env = dict(env)
    env["PATH"] = f"{os.path.join(STUBS_DIR, 'kraken2')}{os.pathsep}{env['PATH']}"
    return workDir, env, stubDir, ["--db", str(db), "--varscripts", str(varScriptDir)]

def generate_kraken2(kraken2_setup, args):
    workDir, _, _, setupArgs = kraken2_setup
    run_psst(["qc", "library", "kraken2", *setupArgs, *args], workDir)
    return os.path.join(workDir, "kraken2.sh"), os.path.join(workDir, "kraken2_tabulate.sh")

def test_kraken2_then_tabulate(kraken2_setup, tmp_path):
    workDir, env, stubDir, _ = kraken2_setup
    classifyScript, tabulateScript = generate_kraken2(kraken2_setup, ["-d", str(tmp_path / "paired"), "--classify-cpus", "2"])
    assert "#PBS -J" not in open(classifyScript).read() # samples run serially, not as an array
    result = run_job(classifyScript, workDir, env=env)
    assert result.stdout.count("Reading the kraken2 database") == 1
    outDir = os.path.join(workDir, "kraken2")
    assert sorted(os.listdir(outDir)) == ["A", "B"]
    assert sorted(os.listdir(os.path.join(outDir, "A"))) == [
        "A_1_classified.fastq", "A_1_unclassified.fastq", "A_2_classified.fastq", "A_2_unclassified.fastq",
        "A_output.txt", "A_report.html", "A_report.krona", "A_report.txt"]
    assert stub_calls(stubDir, "k2")[0] == [
        "classify", "--threads", "2", "--db", f"{tmp_path}/k2db", "--memory-mapping",
        "--output", "kraken2/A.incomplete/A_output.txt",
        "--unclassified-out", "kraken2/A.incomplete/A#_unclassified.fastq",
        "--classified-out", "kraken2/A.incomplete/A#_classified.fastq",
        "--report", "kraken2/A.incomplete/A_report.txt",
        "--paired", f"{tmp_path}/paired/A_R1.fq.gz", f"{tmp_path}/paired/A_R2.fq.gz"]
    assert [ c[-1] for c in stub_calls(stubDir, "k2") ] == [f"{tmp_path}/paired/A_R2.fq.gz", f"{tmp_path}/paired/B.2.fq.gz"]
    assert stub_calls(stubDir, "ktImportText")[0] == [
        "kraken2/A.incomplete/A_report.krona", "-o", "kraken2/A.incomplete/A_report.html"]
    run_job(tabulateScript, workDir, env=env)
    assert stub_calls(stubDir, "combine_kreports.py") == [[
        "-r", "kraken2/A/A_report.txt", "kraken2/B/B_report.txt", "-o", "kraken2/kraken_summary.incomplete.tsv"]]
    assert stub_calls(stubDir, "tabulate") == [[
        "-i", "kraken2/kraken_summary.incomplete.tsv", "-o", "kraken2/kraken_summary.reformat.incomplete.tsv"]]
    assert sorted(os.listdir(outDir)) == ["A", "B", "kraken_summary.reformat.tsv", "kraken_summary.tsv"]
    result = run_job(tabulateScript, workDir, env=env, check=False) # refuses to overwrite the tables
    assert result.returncode != 0 and "already exists" in result.stdout

def test_kraken2_single_end(kraken2_setup, tmp_path):
    workDir, env, stubDir, _ = kraken2_setup
    classifyScript, tabulateScript = generate_kraken2(kraken2_setup, ["-d", str(tmp_path / "single"), "--singleEnd",
                                                                      "-p", "run1"])
    run_array_job(classifyScript, workDir, env=env)
    assert sorted(os.listdir(os.path.join(workDir, "kraken2", "C"))) == [
        "C_classified.fastq", "C_output.txt", "C_report.html", "C_report.krona", "C_report.txt", "C_unclassified.fastq"]
    call = stub_calls(stubDir, "k2")[0]
    assert "--paired" not in call and call[-1] == f"{tmp_path}/single/C.fq.gz"
    run_job(tabulateScript, workDir, env=env)
    assert os.path.exists(os.path.join(workDir, "kraken2", "run1.reformat.tsv"))

def test_kraken2_skips_classified_samples(kraken2_setup, tmp_path):
    workDir, env, stubDir, _ = kraken2_setup
    classifyScript, _ = generate_kraken2(kraken2_setup, ["-d", str(tmp_path / "paired")])
    os.makedirs(os.path.join(workDir, "kraken2", "A"))
    result = run_job(classifyScript, workDir, env=env)
    assert "skipping A" in result.stdout
    assert [ c[-1] for c in stub_calls(stubDir, "k2") ] == [f"{tmp_path}/paired/B.2.fq.gz"] # only B is classified
    result = run_job(classifyScript, workDir, env=env) # with every sample done, the database isn't read
    assert "skipping B" in result.stdout and "Reading the kraken2 database" not in result.stdout
    assert len(stub_calls(stubDir, "k2")) == 1

def test_kraken2_failure_leaves_no_final_output(kraken2_setup, tmp_path):
    workDir, env, stubDir, _ = kraken2_setup
    classifyScript, _ = generate_kraken2(kraken2_setup, ["-d", str(tmp_path / "paired")])
    result = run_job(classifyScript, workDir, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert os.listdir(os.path.join(workDir, "kraken2")) == ["A.incomplete"] # the job stops at the first failure
    assert stub_calls(stubDir, "kreport2krona.py") == [] # nothing runs after a failed classification
    run_job(classifyScript, workDir, env=env) # a rerun replaces the incomplete output
    assert sorted(os.listdir(os.path.join(workDir, "kraken2"))) == ["A", "B"]
    assert len(os.listdir(os.path.join(workDir, "kraken2", "A"))) == 8

def test_kraken2_tabulate_checks_reports(kraken2_setup, tmp_path):
    workDir, env, stubDir, _ = kraken2_setup
    classifyScript, tabulateScript = generate_kraken2(kraken2_setup, ["-d", str(tmp_path / "paired")])
    os.makedirs(os.path.join(workDir, "kraken2", "A")) # as if only sample A were classified
    open(os.path.join(workDir, "kraken2", "A", "A_report.txt"), "w").close()
    result = run_job(tabulateScript, workDir, env=env, check=False)
    assert result.returncode != 0
    assert "kraken2/B/B_report.txt was not found" in result.stdout
    assert stub_calls(stubDir, "combine_kreports.py") == []

def test_kraken2_tabulate_failure_leaves_no_tables(kraken2_setup, tmp_path):
    workDir, env, stubDir, _ = kraken2_setup
    classifyScript, tabulateScript = generate_kraken2(kraken2_setup, ["-d", str(tmp_path / "paired")])
    run_array_job(classifyScript, workDir, env=env)
    result = run_job(tabulateScript, workDir, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert sorted(f for f in os.listdir(os.path.join(workDir, "kraken2")) if not f.endswith("incomplete.tsv")) == ["A", "B"]
    run_job(tabulateScript, workDir, env=env) # a rerun replaces the incomplete tables
    assert sorted(os.listdir(os.path.join(workDir, "kraken2"))) == ["A", "B", "kraken_summary.reformat.tsv", "kraken_summary.tsv"]

@pytest.fixture
def smudgeplot_setup(setup, tmp_path):
    '''
    Extends setup with the FastK and smudgeplot stubs on PATH, and FASTQ reads: paired
    samples A and B of 50 100-bp reads per file, gzipped.

    Returns:
        (workDir, env, stubDir, readsDir) -- as for setup, plus the reads directory.
    '''
    workDir, env, stubDir = setup
    readsDir = tmp_path / "fastq"
    readsDir.mkdir()
    for name in ["A_1.fq.gz", "A_2.fq.gz", "B_1.fq.gz", "B_2.fq.gz"]:
        with gzip.open(readsDir / name, "wt") as fileOut:
            fileOut.write("".join(f"@r{i}\n{'A' * 100}\n+\n{'I' * 100}\n" for i in range(50)))
    env = dict(env)
    stubs = os.path.join(STUBS_DIR, "smudgeplot")
    env["PATH"] = os.pathsep.join([stubs, os.path.join(stubs, "v0.5"), env["PATH"]])
    env["STUB_SMUDGES"] = "AB=10 AAB=30 4A1B=5"
    return workDir, env, stubDir, str(readsDir)

def generate_smudgeplot(smudgeplot_setup, args):
    workDir = smudgeplot_setup[0]
    run_psst(["qc", "library", "smudgeplot", *args], workDir)
    return os.path.join(workDir, "smudgeplot.sh"), os.path.join(workDir, "smudgeplot_tabulate.sh")

def read_tsv(fileName):
    return [ line.split("\t") for line in open(fileName).read().splitlines() ]

PLOIDY_HEADER = ["sample", "bases", "depth", "error_threshold", "smudge", "smudge_fraction", "ploidy"]

def test_smudgeplot_then_tabulate(smudgeplot_setup, tmp_path):
    workDir, env, stubDir, readsDir = smudgeplot_setup
    script, tabulateScript = generate_smudgeplot(smudgeplot_setup, ["-d", readsDir, "-g", "100", "--smudgeplot-cpus", "2",
                                                                    "--smudgeplot-mem", "8G"])
    result = run_array_job(script, workDir, env=env)
    # 10,000 bases give a depth of 100 and k-mer depth of 70 (as 70 of each read's 100 positions
    # start a 31-mer), which is 35 per genome copy, half of which is 17.5
    assert stub_calls(stubDir, "FastK")[0] == [
        "-t4", "-k31", "-M6", "-T2", f"{readsDir}/A_1.fq.gz", f"{readsDir}/A_2.fq.gz",
        "-Nsmudgeplot/A.incomplete/A_table", "-Psmudgeplot/A.incomplete/tmp"]
    assert stub_calls(stubDir, "smudgeplot")[:2] == [
        ["hetmers", "-L", "18", "-t", "2", "-o", "smudgeplot/A.incomplete/A_kmerpairs", "-tmp", "smudgeplot/A.incomplete/tmp",
         "smudgeplot/A.incomplete/A_table.ktab"],
        ["all", "-t", "A", "-o", "smudgeplot/A.incomplete/A", "smudgeplot/A.incomplete/A_kmerpairs.smu"]]
    outDir = os.path.join(workDir, "smudgeplot")
    assert sorted(os.listdir(outDir)) == ["A", "B"]
    assert sorted(os.listdir(os.path.join(outDir, "A"))) == [ # the k-mer table is removed
        "A.ploidy.tsv", "A.smudge_report.tsv", "A_kmerpairs.smu", "A_smudgeplot.png", "A_table.hist"]
    assert read_tsv(os.path.join(outDir, "A", "A.ploidy.tsv")) == [
        PLOIDY_HEADER, ["A", "10000", "100.0", "18", "AAB", "0.6667", "3"]]
    run_job(tabulateScript, workDir, env=env)
    assert read_tsv(os.path.join(outDir, "ploidy_numbers.tsv")) == [
        PLOIDY_HEADER, ["A", "10000", "100.0", "18", "AAB", "0.6667", "3"], ["B", "10000", "100.0", "18", "AAB", "0.6667", "3"]]
    result = run_job(tabulateScript, workDir, env=env, check=False) # refuses to overwrite the table
    assert result.returncode != 0 and "already exists" in result.stdout

def test_smudgeplot_estimates_large_files(smudgeplot_setup, tmp_path):
    workDir, env, stubDir, _ = smudgeplot_setup
    readsDir = tmp_path / "large"
    readsDir.mkdir()
    with open(readsDir / "C.fq", "w") as fileOut: # ~25 MB, so the bases are scaled from the first 10 MB
        fileOut.write("".join("@r\n" + "A" * 100 + "\n+\n" + "I" * 100 + "\n" for _ in range(120000)))
    script, _ = generate_smudgeplot(smudgeplot_setup, ["-d", str(readsDir), "--singleEnd", "-rs", ".fq", "-g", "100k",
                                                       "--ploidy", "3", "--keep-tables"])
    result = run_job(script, workDir, env=env)
    sample = read_tsv(os.path.join(workDir, "smudgeplot", "C", "C.ploidy.tsv"))[1]
    assert abs(int(sample[1]) - 12000000) < 12000 # within 0.1% of the 12 Mb sequenced
    # A depth of 120 gives a k-mer depth of 84, which is 28 per genome copy at ploidy 3, half of which is 14
    assert sample[2:] == ["120.0", "14", "AAB", "0.6667", "3"]
    assert "using -L 14" in result.stdout
    assert sorted(f for f in os.listdir(os.path.join(workDir, "smudgeplot", "C")) if "table" in f) == [
        ".C_table.ktab.1", "C_table.hist", "C_table.ktab"] # --keep-tables keeps the k-mer table

def test_smudgeplot_given_threshold_and_ploidy_names(smudgeplot_setup, tmp_path):
    workDir, env, stubDir, readsDir = smudgeplot_setup
    script, _ = generate_smudgeplot(smudgeplot_setup, ["-d", readsDir, "--L", "7"])
    os.makedirs(os.path.join(workDir, "smudgeplot", "A"))
    assert "skipping A" in run_job(script, workDir, index=1, env=env).stdout
    run_job(script, workDir, index=2, env=dict(env, STUB_SMUDGES="AB=1 AAB=1 4A1B=3")) # 4A1B is 5 copies
    assert stub_calls(stubDir, "smudgeplot")[0][:3] == ["hetmers", "-L", "7"]
    assert read_tsv(os.path.join(workDir, "smudgeplot", "B", "B.ploidy.tsv"))[1] == ["B", ".", ".", "7", "4A1B", "0.6000", "5"]
    assert len(stub_calls(stubDir, "FastK")) == 1 # A was skipped

def test_smudgeplot_earlier_versions(smudgeplot_setup, tmp_path):
    workDir, env, stubDir, readsDir = smudgeplot_setup
    stubs = os.path.join(STUBS_DIR, "smudgeplot")
    path = [ d for d in env["PATH"].split(os.pathsep) if not os.path.exists(os.path.join(d, "smudgeplot")) ]
    env = dict(env, PATH=os.pathsep.join([stubs, os.path.join(stubs, "v0.4"), *path]))
    script, _ = generate_smudgeplot(smudgeplot_setup, ["-d", readsDir, "--L", "7"])
    run_job(script, workDir, index=1, env=dict(env, STUB_SMUDGES="AB=0.2 AABB=0.7 AAB=0.1"))
    assert stub_calls(stubDir, "smudgeplot")[1][-1] == "smudgeplot/A.incomplete/A_kmerpairs_text.smu"
    assert read_tsv(os.path.join(workDir, "smudgeplot", "A", "A.ploidy.tsv"))[1] == ["A", ".", ".", "7", "AABB", "0.7000", "4"]
    run_job(script, workDir, index=2, env=dict(env, STUB_SMUDGES="")) # no smudges were found
    assert read_tsv(os.path.join(workDir, "smudgeplot", "B", "B.ploidy.tsv"))[1] == ["B", ".", ".", "7", ".", ".", "undetermined"]

def test_smudgeplot_failure_leaves_no_final_output(smudgeplot_setup, tmp_path):
    workDir, env, stubDir, readsDir = smudgeplot_setup
    script, _ = generate_smudgeplot(smudgeplot_setup, ["-d", readsDir, "-g", "100"])
    result = run_job(script, workDir, index=1, env=dict(env, STUB_FAIL="1"), check=False)
    assert result.returncode != 0
    assert os.listdir(os.path.join(workDir, "smudgeplot")) == ["A.incomplete"]
    assert stub_calls(stubDir, "smudgeplot") == [] # nothing runs after a failed k-mer count
    run_job(script, workDir, index=1, env=env) # a rerun replaces the incomplete output
    assert os.listdir(os.path.join(workDir, "smudgeplot")) == ["A"]

def test_smudgeplot_tabulate_checks_summaries(smudgeplot_setup, tmp_path):
    workDir, env, stubDir, readsDir = smudgeplot_setup
    script, tabulateScript = generate_smudgeplot(smudgeplot_setup, ["-d", readsDir, "-g", "100", "-o", "ploidy.tsv"])
    run_job(script, workDir, index=1, env=env) # as if only sample A were analysed
    result = run_job(tabulateScript, workDir, env=env, check=False)
    assert result.returncode != 0
    assert "smudgeplot/B/B.ploidy.tsv was not found" in result.stdout
    assert not os.path.exists(os.path.join(workDir, "ploidy.tsv"))
