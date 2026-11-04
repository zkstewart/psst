# Tests of the "assembly braker3" job script's own logic (i.e., linking local reads under
# the names BRAKER3 looks for, combining SRA accession files, checking RNAseq set IDs, and
# refusing to repeat a run), run against a stub singularity (see stubs/braker3).

import os

import pytest

from conftest import STUBS_DIR, run_psst, run_job

@pytest.fixture
def setup(tmp_path):
    '''
    Creates a genome, proteins, a reads directory whose files are links to another directory
    (which must also be bound into the container), two SRA accession files, a placeholder
    image, a work directory, and an environment with the stub on PATH.

    Returns:
        (workDir, files, env, stubDir) -- the work directory, a dict of input locations, a dict
            of environment variables for run_job(), and the directory the stub records in.
    '''
    (tmp_path / "inputs").mkdir() # bound into the container, but not holding the reads
    files = { name: str(tmp_path / "inputs" / name) for name in ["genome.fasta.masked", "proteins.fa", "braker3.sif",
                                                                  "reticulata.txt", "glauca.txt"] }
    for name in ["genome.fasta.masked", "proteins.fa", "braker3.sif"]:
        open(files[name], "w").write(">x\nACGT\n")
    open(files["reticulata.txt"], "w").write("SRR100\nSRR101\n\n")
    open(files["glauca.txt"], "w").write("SRR200\r\nSRR101\r\n") # a repeat, with Windows line endings
    store = tmp_path / "store"
    store.mkdir()
    files["reads"] = str(tmp_path / "reads")
    os.mkdir(files["reads"])
    for sample in ["leaf", "root"]:
        for r in [1, 2]:
            (store / f"{sample}.trimmed_{r}P.fq.gz").write_text("")
            os.symlink(store / f"{sample}.trimmed_{r}P.fq.gz", os.path.join(files["reads"], f"{sample}.trimmed_{r}P.fq.gz"))
    workDir = tmp_path / "work"
    workDir.mkdir()
    stubDir = tmp_path / "stub_records"
    stubDir.mkdir()
    env = dict(os.environ)
    env["PATH"] = f"{os.path.join(STUBS_DIR, 'braker3')}{os.pathsep}{env['PATH']}"
    env["STUB_DIR"] = str(stubDir)
    return str(workDir), files, env, str(stubDir)

def generate(setup, args):
    workDir, files = setup[:2]
    run_psst(["assembly", "braker3", "-g", files["genome.fasta.masked"], "--proteins", files["proteins.fa"],
              "--sif", files["braker3.sif"], *args], workDir)

def run(setup, env=None, check=True):
    workDir, _, stubEnv, stubDir = setup
    result = run_job(os.path.join(workDir, "braker3.sh"), workDir, env=env or stubEnv, check=check)
    logFile = os.path.join(stubDir, "singularity.log")
    calls = [ line.split("|") for line in open(logFile).read().splitlines() ] if os.path.exists(logFile) else []
    return result, calls

def rnaseq_sets(stubDir):
    return dict( line.split("\t") for line in open(os.path.join(stubDir, "rnaseq.tsv")).read().splitlines() )

def test_braker3_local_and_sra(setup):
    workDir, files, _, stubDir = setup
    generate(setup, ["-d", files["reads"], "-rs", "P.fq.gz", "--sra", files["reticulata.txt"], files["glauca.txt"]])
    result, calls = run(setup)
    outDir = os.path.join(workDir, "braker3")
    assert open(os.path.join(outDir, "sra_accessions.txt")).read() == "SRR100\nSRR101\nSRR200\n"
    assert sorted(os.listdir(os.path.join(outDir, "rnaseq_reads"))) == ["leaf.trimmed_1.fq.gz", "leaf.trimmed_2.fq.gz",
                                                                         "root.trimmed_1.fq.gz", "root.trimmed_2.fq.gz"]
    assert rnaseq_sets(stubDir) == { "leaf.trimmed": "leaf.trimmed_1.fq.gz,leaf.trimmed_2.fq.gz",
                                     "root.trimmed": "root.trimmed_1.fq.gz,root.trimmed_2.fq.gz",
                                     "SRR100": "SRA", "SRR101": "SRA", "SRR200": "SRA" }
    assert [ c[4] for c in calls ] == ["cp", "braker.pl"]
    braker = calls[1]
    for arg in [f"--species=genome", f"--busco_lineage=eudicots", "--threads=24", "--gff3", "--nocleanup",
                f"--genome={os.path.realpath(files['genome.fasta.masked'])}",
                f"--rnaseq_sets_dirs={os.path.realpath(outDir)}/rnaseq_reads"]:
        assert arg in braker
    for ext in ["gff3", "gtf", "aa", "codingseq"]:
        assert open(os.path.join(outDir, f"genome.braker.{ext}")).read() == f"{ext}\n"
    assert not [ f for f in os.listdir(outDir) if f.endswith(".incomplete") ]

def test_braker3_sra_only(setup):
    workDir, files, _, stubDir = setup
    generate(setup, ["--sra", files["glauca.txt"], "-p", "glauca"])
    result, calls = run(setup)
    assert rnaseq_sets(stubDir) == { "SRR200": "SRA", "SRR101": "SRA" }
    assert not any(a.startswith("--rnaseq_sets_dirs") for a in calls[1])
    assert os.listdir(os.path.join(workDir, "braker3", "rnaseq_reads")) == []
    assert os.path.exists(os.path.join(workDir, "braker3", "glauca.braker.gtf"))

def test_braker3_refuses_to_resume_or_repeat(setup):
    workDir, files, env, _ = setup
    generate(setup, ["-d", files["reads"], "-rs", "P.fq.gz"])
    failEnv = dict(env)
    failEnv["STUB_FAIL"] = "1"
    result, calls = run(setup, env=failEnv, check=False)
    assert result.returncode != 0
    result, calls = run(setup, check=False)
    assert result.returncode != 0 and "which BRAKER3 cannot resume" in result.stdout
    assert len(calls) == 2

    for name in ["braker", "augustus_config"]:
        os.system(f"rm -rf {os.path.join(workDir, 'braker3', name)}")
    run(setup)
    result, calls = run(setup, check=False)
    assert result.returncode != 0 and "already holds BRAKER3 results for genome" in result.stdout
    assert len(calls) == 4

@pytest.mark.parametrize("change, message", [
    (lambda f: open(os.path.join(f["reads"], "groot.trimmed_1P.fq.gz"), "w").close() or
               open(os.path.join(f["reads"], "groot.trimmed_2P.fq.gz"), "w").close(),
     "'root.trimmed' matches the end of 'groot.trimmed'"),
    (lambda f: open(f["glauca.txt"], "a").write("leaf.trimmed\n"), "is not an SRA run accession"),
    (lambda f: os.remove(os.path.join(f["reads"], "leaf.trimmed_2P.fq.gz")), "has no matching reverse read file"),
])
def test_braker3_checks_inputs_at_runtime(setup, change, message):
    workDir, files, _, stubDir = setup
    generate(setup, ["-d", files["reads"], "-rs", "P.fq.gz", "--sra", files["glauca.txt"]])
    change(files) # after psst validated them
    result, calls = run(setup, check=False)
    assert result.returncode != 0 and message in result.stdout
    assert calls == []
