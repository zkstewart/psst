# Tests of psst's command line: every mode generates valid job scripts with only --rpid,
# validation errors are raised before anything is written, and --qsub submits (or
# refuses to submit) as expected.

import glob, gzip, os, re, subprocess, sys

import pytest

from conftest import PSST, RPID, run_psst

def discover_modes(path=()):
    '''
    Finds every leaf mode (e.g., ('map', 'index', 'bwa')) by walking psst's -h output,
    so that new modes are tested without updating this file.
    '''
    result = subprocess.run([sys.executable, PSST, *path, "-h"], capture_output=True, text=True)
    match = re.search(r"positional arguments:\n\s+\{([a-z0-9_,]+)\}", result.stdout)
    if match == None:
        return [ path ]
    return [ mode for choice in match.group(1).split(",") for mode in discover_modes(path + (choice,)) ]

# Modes which run immediately rather than writing a job script
NON_JOB_MODES = [ ("basespace", "prepare") ]
MODES = [ m for m in discover_modes() if m not in NON_JOB_MODES ]

def written_scripts(directory):
    return sorted(glob.glob(os.path.join(directory, "**", "*.sh"), recursive=True))

def reported_placeholders(stdout):
    match = re.search(r"placeholders to fill in before submission: (.+)", stdout)
    return set(match.group(1).split(", ")) if match else set()

def test_modes_are_discovered():
    assert ("blast", "run", "mmseqs2") in MODES
    assert ("call", "short") in MODES
    assert ("map", "index", "giraffe") in MODES
    assert ("hic", "contacts") in MODES

@pytest.mark.parametrize("mode", MODES, ids=" ".join)
def test_mode_generates_with_only_rpid(mode, tmp_path):
    result = run_psst(list(mode), tmp_path)
    scripts = written_scripts(str(tmp_path))
    assert scripts != [], "no job script was written"
    placeholders = set()
    for script in scripts:
        text = open(script).read()
        assert text.startswith("#!/bin/bash -l\n")
        assert f"#PBS -P {RPID}" in text
        subprocess.run(["bash", "-n", script], check=True) # valid bash syntax
        placeholders.update(re.findall(r"@@([A-Za-z_][A-Za-z0-9_]*)@@", text))
    assert placeholders == reported_placeholders(result.stdout)

@pytest.mark.parametrize("mode", MODES, ids=" ".join)
def test_mode_refuses_existing_scripts(mode, tmp_path):
    run_psst(list(mode), tmp_path)
    result = run_psst(list(mode), tmp_path, check=False)
    assert result.returncode != 0
    assert "already exists" in result.stderr

@pytest.mark.parametrize("mode", MODES, ids=" ".join)
def test_qsub_refused_with_placeholders(mode, tmp_path, fake_qsub):
    env, logFile = fake_qsub
    previewDir, workDir = tmp_path / "preview", tmp_path / "work"
    previewDir.mkdir()
    workDir.mkdir()
    if mode == ("basecalling", "dorado"):
        (workDir / "pod5").mkdir() # dorado refuses --qsub without its input directory
    if reported_placeholders(run_psst(list(mode), previewDir).stdout) == set():
        pytest.skip("mode has no placeholders by default")
    result = run_psst(list(mode) + ["--qsub"], workDir, env=env, check=False)
    assert result.returncode != 0
    assert "--qsub cannot be used while job scripts have unfilled placeholders" in result.stderr
    assert written_scripts(str(workDir)) == [] # nothing is written
    assert not os.path.exists(logFile)

def test_invalid_rpid_rejected(tmp_path):
    result = subprocess.run([sys.executable, PSST, "call", "short", "--rpid", "BAD"],
                            cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode != 0
    assert "must be 6 letters followed by 4 numbers" in result.stderr

@pytest.fixture
def inputs(tmp_path):
    '''
    Creates small input files for validation tests, returning a dict of their locations.
    '''
    files = {}
    files["genome"] = tmp_path / "genome.fa"
    files["genome"].write_text(">chr1\nACGT\n>chr2\nACGT\n")
    (tmp_path / "genome.fa.fai").write_text("chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\n")
    files["nofai"] = tmp_path / "nofai.fa"
    files["nofai"].write_text(">chr1\nACGT\n")
    files["bams"] = tmp_path / "bams"
    files["bams"].mkdir()
    for sample in ["s1", "s2"]:
        (files["bams"] / f"{sample}.sorted.md.bam").write_text("not a real BAM")
        (files["bams"] / f"{sample}.sorted.md.bam.bai").write_text("")
    files["unindexed"] = tmp_path / "unindexed"
    files["unindexed"].mkdir()
    (files["unindexed"] / "s1.sorted.md.bam").write_text("")
    files["reads"] = tmp_path / "reads"
    files["reads"].mkdir()
    for sample in ["s1", "s2"]:
        for r in [1, 2]:
            (files["reads"] / f"{sample}_{r}.fq.gz").write_text("")
    files["unpaired"] = tmp_path / "unpaired"
    files["unpaired"].mkdir()
    (files["unpaired"] / "s1_1.fq.gz").write_text("")
    files["gbz"] = tmp_path / "graph.gbz"
    files["gbz"].write_text("")
    files["giraffe"] = tmp_path / "pg_giraffe_index"
    files["giraffe"].mkdir()
    for name in ["pangenome.gbz", "pangenome.dist", "pangenome.min", "pangenome.zipcodes",
                 "reference.fa", "targets.tsv.gz"]:
        (files["giraffe"] / name).write_text("")
    (files["giraffe"] / "reference.fa.fai").write_text("reference#0#chr1\t4\t6\t4\t5\n")
    with gzip.open(files["giraffe"] / "sites.vcf.gz", "wt") as fileOut:
        fileOut.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\ts1\n")
    files["vgdir"] = tmp_path / "vgdir"
    files["vgdir"].mkdir()
    (files["vgdir"] / "vg").write_text("")
    return { k: str(v) for k, v in files.items() }

VALIDATION_ERRORS = [
    # call modes
    (["call", "short", "-r", "{nofai}"], "is not indexed; run 'samtools faidx'"),
    (["call", "short", "-r", "{genome}", "-b", "{reads}"], "No BAM files ending with '.sorted.md.bam'"),
    (["call", "short", "-b", "{unindexed}"], "BAM file(s) are not indexed"),
    (["call", "short", "--chunks", "0"], "--chunks must be >= 1"),
    (["call", "short", "-q", "-1"], "-q must be >= 0"),
    (["call", "short", "-Q", "-1"], "-Q must be >= 0"),
    (["call", "short", "--max-depth", "0"], "--max-depth must be >= 1"),
    (["call", "short", "-p", "a/b"], "must not contain whitespace or '/'"),
    (["call", "short", "-bs", ""], "-bs must be given a non-empty value"),
    (["call", "short", "--call-mem", "8X"], "is not a valid PBS memory value"),
    (["call", "short", "--concat-walltime", "49"], "exceeds the HPC limit of 48 hours"),
    (["call", "short", "--ploidy", "4"], "invalid choice"),
    (["call", "nanopore", "--hic"], "unrecognized arguments: --hic"),
    (["call", "giraffe", "--hic"], "unrecognized arguments: --hic"),
    (["call", "giraffe", "-i", "{unpaired}"], "is missing reference.fa"),
    # map modes
    (["map", "short", "--hic", "--singleEnd"], "--hic and --singleEnd are mutually exclusive"),
    (["map", "short", "-r", "{nofai}"], "run 'bwa index' on it first"),
    (["map", "short", "-d", "{unpaired}"], "has no matching reverse read file"),
    (["map", "giraffe", "-i", "{unpaired}"], "is missing pangenome.gbz"),
    (["map", "giraffe", "--vg", "{unpaired}"], "vg not found in --vg"),
    (["map", "index", "giraffe", "-g", "{genome}.missing"], "Pangenome GBZ file"),
    (["map", "index", "giraffe", "--refsample", "a b"], "must not contain whitespace or '#'"),
    (["map", "index", "giraffe", "-g", "{gbz}", "-o", "{giraffe}"], "already exists"),
]

@pytest.mark.parametrize("args,message", VALIDATION_ERRORS, ids=lambda x: " ".join(x) if isinstance(x, list) else None)
def test_validation_error(args, message, inputs, tmp_path):
    workDir = tmp_path / "work"
    workDir.mkdir()
    result = run_psst([ a.format(**inputs) for a in args ], workDir, check=False)
    assert result.returncode != 0
    assert message in result.stderr
    assert written_scripts(str(workDir)) == [] # nothing is written on error

VALIDATION_WARNINGS = [
    (["call", "short", "-b", "{unindexed}", "--afterok", "1.pbs"], "are not indexed yet"),
    (["call", "short", "-b", "{tmp}/nowhere", "--afterok", "1.pbs"], "the --afterok job(s) must create them"),
    (["call", "giraffe", "-i", "{tmp}/nowhere", "--afterok", "1.pbs"], "giraffe index"),
    (["map", "giraffe", "-i", "{tmp}/nowhere", "--afterok", "1.pbs"], "giraffe index"),
    (["call", "giraffe", "-i", "{giraffe}", "-b", "{bams}"], "sample(s) s1 are also samples in"),
]

@pytest.mark.parametrize("args,message", VALIDATION_WARNINGS, ids=lambda x: " ".join(x) if isinstance(x, list) else None)
def test_validation_warning(args, message, inputs, tmp_path, monkeypatch):
    # Sample names are read from BAM headers, so give s1's BAM a real header for the clash test
    with gzip.open(os.path.join(inputs["bams"], "s1.sorted.md.bam"), "wb") as fileOut:
        text = b"@RG\tID:s1\tSM:s1\n"
        fileOut.write(b"BAM\x01" + len(text).to_bytes(4, "little") + text)
    workDir = tmp_path / "work"
    workDir.mkdir()
    result = run_psst([ a.format(tmp=tmp_path, **inputs) for a in args ], workDir)
    assert message in result.stdout

def test_call_chunks_capped_to_contigs(inputs, tmp_path):
    result = run_psst(["call", "short", "-r", inputs["genome"], "--chunks", "10"], tmp_path)
    assert "Reducing --chunks from 10 to the genome's 2 contig(s)" in result.stdout
    assert "#PBS -J 1-2" in open(tmp_path / "call_short.sh").read()

def test_call_derives_prefixes(inputs, tmp_path):
    run_psst(["call", "short", "-r", inputs["genome"]], tmp_path)
    assert "PREFIX=genome\n" in open(tmp_path / "call_short.sh").read()
    run_psst(["call", "giraffe", "-i", inputs["giraffe"]], tmp_path)
    text = open(tmp_path / "call_giraffe.sh").read()
    assert "PREFIX=pg\n" in text
    assert f"TARGETS={inputs['giraffe']}/targets.tsv.gz\n" in text
    assert f"SITES={inputs['giraffe']}/sites.vcf.gz\n" in open(tmp_path / "call_giraffe_concat.sh").read()

def test_call_short_hic_sets_orphans(tmp_path):
    run_psst(["call", "short", "--hic"], tmp_path)
    assert "ORPHANS=true" in open(tmp_path / "call_short.sh").read()

def test_map_short_hic_sets_flag(tmp_path):
    run_psst(["map", "short", "--hic"], tmp_path)
    assert "HIC=true" in open(tmp_path / "map_short.sh").read()

def test_map_index_giraffe_derives_index(inputs, tmp_path):
    run_psst(["map", "index", "giraffe", "-g", inputs["gbz"]], tmp_path)
    assert f"INDEX={inputs['gbz'][:-4]}_giraffe_index\n" in open(tmp_path / "giraffe_index.sh").read()

def test_qsub_chains_call_jobs(inputs, tmp_path, fake_qsub):
    env, logFile = fake_qsub
    result = run_psst(["call", "short", "-r", inputs["genome"], "-b", inputs["bams"], "--qsub",
                       "--afterok", "99[].pbs"], tmp_path, env=env)
    calls = open(logFile).read().splitlines()
    assert len(calls) == 2
    assert calls[0].endswith("call_short.sh")
    assert calls[1] == f"-W depend=afterok:101[].pbs {tmp_path}/call_short_concat.sh"
    assert "#PBS -W depend=afterok:99[].pbs" in open(tmp_path / "call_short.sh").read()
    assert "depend" not in open(tmp_path / "call_short_concat.sh").read()

def test_qsub_single_job(inputs, tmp_path, fake_qsub):
    env, logFile = fake_qsub
    run_psst(["map", "short", "-r", inputs["genome"], "-d", inputs["reads"], "--qsub"], tmp_path, env=env,
             check=False) # the genome has no bwa index
    assert not os.path.exists(logFile)
    run_psst(["map", "giraffe", "-i", inputs["giraffe"], "-d", inputs["reads"], "--qsub"], tmp_path, env=env)
    assert open(logFile).read() == f"{tmp_path}/map_giraffe.sh\n"
    assert "#PBS -J 1-2" in open(tmp_path / "map_giraffe.sh").read()
