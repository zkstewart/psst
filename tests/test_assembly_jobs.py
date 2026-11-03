# Tests of the "assembly revcomp" and "assembly scaffolding" job scripts' own logic (i.e.,
# phasing or reusing a phasing, then processing each partition), run against a stub
# annotarium (see stubs/annotarium).

import os

import pytest

from conftest import STUBS_DIR, run_psst, run_job

STUB_ANNOTARIUM = os.path.join(STUBS_DIR, "annotarium")

@pytest.fixture
def setup(tmp_path):
    '''
    Creates contig and reference FASTAs (one with a space in its name), a directory of
    placeholder 'minimap2' and 'paf2dotplot.r' files to pass validation, and a work directory.

    Returns:
        (workDir, contigs, reference, env, logFile, toolDir) -- the work directory, a list
            of the contig FASTA locations, the reference FASTA location, a dict of
            environment variables for run_job(), the stub's log file location, and the
            placeholder tool directory.
    '''
    contigs = [ str(tmp_path / "hap1.fa"), str(tmp_path / "hap 2.fa") ]
    reference = str(tmp_path / "ref.fa")
    for fileName in contigs + [reference]:
        open(fileName, "w").write(">x\nACGT\n")
    toolDir = tmp_path / "tools"
    toolDir.mkdir()
    for tool in ["minimap2", "paf2dotplot.r"]:
        (toolDir / tool).write_text("")
    workDir = tmp_path / "work"
    workDir.mkdir()
    stubDir = tmp_path / "stub_records"
    stubDir.mkdir()
    env = dict(os.environ)
    env["STUB_DIR"] = str(stubDir)
    return str(workDir), contigs, reference, env, str(stubDir / "annotarium.log"), str(toolDir)

def run_mode(mode, setup, args=[], check=True, env=None):
    workDir, contigs, reference, stubEnv, logFile, toolDir = setup
    if not os.path.exists(os.path.join(workDir, f"{mode}.sh")):
        extra = ["--paf2dotplot", toolDir] if mode == "scaffolding" else []
        run_psst(["assembly", mode, "-i", *contigs, "-r", reference, "--annotarium", STUB_ANNOTARIUM,
                  "--minimap2", toolDir, *extra, *args], workDir)
    result = run_job(os.path.join(workDir, f"{mode}.sh"), workDir, env=env or stubEnv, check=check)
    calls = [ line.split("|") for line in open(logFile).read().splitlines() ] if os.path.exists(logFile) else []
    return result, calls

def test_revcomp_phases_then_pairs_each_partition(setup):
    workDir, contigs, reference = setup[:3]
    result, calls = run_mode("revcomp", setup, ["--ploidy", "3", "--preset", "asm10"])
    assert [ c[:2] for c in calls ] == [["assembly", "phase"]] + [["assembly", "homology_pair"]] * 3
    phase = calls[0]
    assert phase[phase.index("-i")+1:phase.index("-i")+3] == contigs # file names with spaces stay intact
    assert phase[phase.index("-o")+1] == "output/phase"
    assert phase[phase.index("--ploidy")+1] == "3" and phase[phase.index("--preset")+1] == "asm10"
    assert phase[phase.index("--minimap2")+1] == os.path.join(setup[5], "minimap2")
    for partition, call in zip(["1", "2", "3"], calls[1:]):
        assert call[call.index("-i")+1] == "output/phase/phasing.tsv"
        assert call[call.index("-f")+1:call.index("-f")+3] == contigs
        assert call[call.index("-r")+1] == reference
        assert call[call.index("--partition")+1] == partition
        assert open(os.path.join(workDir, "output", f"partition{partition}.fasta")).read() == f">partition{partition}\n"
    assert sorted(os.listdir(os.path.join(workDir, "output"))) == ["partition1.fasta", "partition2.fasta",
                                                                   "partition3.fasta", "phase"]

def test_revcomp_reuses_given_phasing(setup, tmp_path):
    phaseDir = tmp_path / "earlier"
    phaseDir.mkdir()
    (phaseDir / "phasing.tsv").write_text("#contig\tlength\tpartition\tstatus\n" +
                                          "a\t1\t2\tassigned\nb\t1\t2\tassigned\nc\t1\t.\tredundant\n")
    result, calls = run_mode("revcomp", setup, ["--phase", str(phaseDir)])
    assert [ c[1] for c in calls ] == ["homology_pair"] # only partition 2 has contigs
    assert calls[0][calls[0].index("-i")+1] == str(phaseDir / "phasing.tsv")
    assert sorted(os.listdir(os.path.join(setup[0], "output"))) == ["partition2.fasta"]

def test_revcomp_resumes_after_failure(setup):
    workDir, env = setup[0], dict(setup[3])
    env["STUB_FAIL_PARTITION"] = "2"
    result, calls = run_mode("revcomp", setup, env=env, check=False)
    assert result.returncode != 0
    assert os.path.exists(os.path.join(workDir, "output", "partition1.fasta"))
    assert not os.path.exists(os.path.join(workDir, "output", "partition2.fasta"))
    open(os.path.join(workDir, "output", "partition2.incomplete.fasta"), "w").close() # as a killed job might leave

    result, calls = run_mode("revcomp", setup)
    assert "Reusing the completed phasing in output/phase" in result.stdout
    assert "Skipping partition 1" in result.stdout
    assert [ (c[1], c[c.index("--partition")+1]) for c in calls if c[1] == "homology_pair" ] == \
        [("homology_pair", "1"), ("homology_pair", "2"), ("homology_pair", "2")]
    assert [ c[1] for c in calls ].count("phase") == 1
    assert sorted(os.listdir(os.path.join(workDir, "output"))) == ["partition1.fasta", "partition2.fasta", "phase"]

def test_revcomp_fails_without_assigned_contigs(setup, tmp_path):
    phaseDir = tmp_path / "earlier"
    phaseDir.mkdir()
    (phaseDir / "phasing.tsv").write_text("#contig\tlength\tpartition\tstatus\na\t1\t.\tunaligned\n")
    result, calls = run_mode("revcomp", setup, ["--phase", str(phaseDir)], check=False)
    assert result.returncode != 0
    assert "No contigs were assigned to a partition" in result.stdout
    assert calls == []

def test_scaffolding_phases_then_scaffolds(setup):
    workDir, contigs, reference = setup[:3]
    result, calls = run_mode("scaffolding", setup, ["--preset", "asm5"])
    assert [ c[1] for c in calls ] == ["phase", "scaffold"]
    scaffold = calls[1]
    assert scaffold[scaffold.index("-i")+1] == "output/phase/phasing.tsv"
    assert scaffold[scaffold.index("-f")+1:scaffold.index("-f")+3] == contigs
    assert scaffold[scaffold.index("-r")+1] == reference
    assert scaffold[scaffold.index("-o")+1] == "output"
    assert scaffold[scaffold.index("--paf2dotplot")+1] == os.path.join(setup[5], "paf2dotplot.r")
    assert scaffold[scaffold.index("--preset")+1] == "asm5"
    assert sorted(os.listdir(os.path.join(workDir, "output"))) == ["partition1.fasta", "partition2.fasta", "phase"]

def test_scaffolding_reuses_given_phasing(setup, tmp_path):
    phaseDir = tmp_path / "earlier"
    phaseDir.mkdir()
    (phaseDir / "phasing.tsv").write_text("#contig\tlength\tpartition\tstatus\na\t1\t1\tassigned\n")
    result, calls = run_mode("scaffolding", setup, ["--phase", str(phaseDir)])
    assert [ c[1] for c in calls ] == ["scaffold"]
    assert calls[0][calls[0].index("-i")+1] == str(phaseDir / "phasing.tsv")
