# Integration tests of job steps run with annotarium's 'fasta' modes: "basecalling dorado"'s
# concatenation of corrected read blocks ('fasta concat', then 'fasta stats'), and
# "assembly hifiasm"'s statistics of its contig FASTAs ('fasta stats'), with hifiasm and
# compleasm stubbed out. Requires annotarium at $PSST_TEST_ANNOTARIUM, or else beside the
# psst repository.

import os, random

import pytest

import sim
from conftest import TESTS_DIR, run_psst, run_job

pytestmark = pytest.mark.integration

ANNOTARIUM = os.environ.get("PSST_TEST_ANNOTARIUM",
                            os.path.join(os.path.dirname(os.path.dirname(TESTS_DIR)), "annotarium"))

# A stand-in for hifiasm, writing primary and haplotype contig graphs (GFA) whose contigs
# are $STUB_CONTIG_LENGTH long; 'hifiasm --version' prints a version
HIFIASM_STUB = """#!/bin/bash
if [[ $1 == "--version" ]]; then echo "0.0-stub"; exit 0; fi
PREFIX=$2
for graph in bp.hap1.p_ctg bp.hap2.p_ctg bp.p_ctg; do
    printf "S\\tctg1\\t%s\\nS\\tctg2\\t%s\\n" $(head -c ${STUB_CONTIG_LENGTH} /dev/zero | tr '\\0' A) \\
        $(head -c ${STUB_CONTIG_LENGTH} /dev/zero | tr '\\0' C) > ${PREFIX}.${graph}.gfa
done
"""

@pytest.fixture
def annotarium():
    if not os.path.isfile(os.path.join(ANNOTARIUM, "annotarium.py")):
        pytest.skip(f"requires annotarium at {ANNOTARIUM} (or set $PSST_TEST_ANNOTARIUM)")
    return ANNOTARIUM

def stats_totals(statsFile):
    '''
    Returns:
        (numRows, total) -- the number of sequence rows in an annotarium 'fasta stats' file,
                            and the total length it reports.
    '''
    lines = [ line.rstrip("\n").split("\t") for line in open(statsFile) ]
    rows = [ l for l in lines[1:] if not l[0].startswith("#") ]
    total = [ l for l in lines if l[0] == "#Total" ][0]
    return len(rows), int(total[1])

def test_dorado_concat(annotarium, tmp_path):
    run_psst(["basecalling", "dorado", "-p", "sp", "--annotarium", annotarium, "--concat-minimum", "100"], tmp_path)
    correctedDir = tmp_path / "fastq_corrected"
    (correctedDir / "step2").mkdir()
    rng = random.Random(4)
    reads = {}
    for block in range(3):
        blockReads = { f"r{block}_{i}": "".join(rng.choice(sim.BASES) for _ in range(rng.randint(50, 400))) for i in range(20) }
        sim.write_fasta(str(correctedDir / "step2" / f"reads.fastq.block_{block}.fasta"), blockReads)
        reads.update(blockReads)
    kept = { name: seq for name, seq in reads.items() if len(seq) >= 100 }

    for _ in range(2): # a rerun replaces the statistics
        run_job(str(correctedDir / "concat_corrected_fastq.sh"), str(correctedDir))
        outputDir = correctedDir / "output"
        assert sim.read_fasta(str(outputDir / "sp_corrected_reads.fasta")) == kept
        assert stats_totals(outputDir / "sp_corrected_reads.stats") == (len(kept), sum(len(s) for s in kept.values()))
        assert sorted(os.listdir(outputDir)) == ["sp_corrected_reads.fasta", "sp_corrected_reads.stats"]

def test_hifiasm_stats(annotarium, tmp_path):
    (tmp_path / "fastq_corrected" / "output").mkdir(parents=True)
    (tmp_path / "fastq_corrected" / "output" / "sp_corrected_reads.fasta").write_text(">r\nACGT\n")
    workDir = tmp_path / "assembly"
    binDir = tmp_path / "bin"
    compleasmDir = tmp_path / "compleasm"
    for d in [workDir, binDir, compleasmDir]:
        d.mkdir()
    (binDir / "hifiasm").write_text(HIFIASM_STUB)
    (binDir / "hifiasm").chmod(0o755)
    (compleasmDir / "compleasm.py").write_text("") # does nothing
    run_psst(["assembly", "hifiasm", "--annotarium", annotarium, "--compleasm", str(compleasmDir), "--cpus", "2"], workDir)

    for contigLength in [100, 250]: # a rerun replaces the statistics of the regenerated contigs
        env = dict(os.environ, PATH=f"{binDir}{os.pathsep}{os.environ['PATH']}", STUB_CONTIG_LENGTH=str(contigLength))
        run_job(str(workDir / "hifiasm.sh"), str(workDir), env=env)
        for graph in ["bp.hap1.p_ctg", "bp.hap2.p_ctg", "bp.p_ctg"]:
            assert stats_totals(workDir / f"sp.asm.{graph}.stats") == (2, 2 * contigLength)
        assert [ f for f in os.listdir(workDir) if f.endswith(".incomplete") ] == []
