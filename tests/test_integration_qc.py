# Integration tests of "qc library" modes on simulated reads, checking their outputs against
# the simulation's truth: "trim" runs bbduk on reads which read through into adapters or end
# in low-quality bases (requiring bbduk.sh on PATH), and "fastqc" runs FastQC then Various_scripts'
# tabulate_fastqc.py (requiring fastqc on PATH, pandas and xlsxwriter, and Various_scripts at
# $PSST_TEST_VARSCRIPTS, or else beside the psst repository).

import gzip, os, random, shutil

import pytest

import sim
from conftest import STUBS_DIR, TESTS_DIR, require_tools, run_psst, run_job, run_array_job

pytestmark = pytest.mark.integration

SAMPLES = ["S1", "S2"]
MINLEN = 25
VARSCRIPTS = os.environ.get("PSST_TEST_VARSCRIPTS",
                            os.path.join(os.path.dirname(os.path.dirname(TESTS_DIR)), "Various_scripts"))

def read_fastq(fileName):
    with gzip.open(fileName, "rt") as fileIn:
        lines = fileIn.read().splitlines()
    return { lines[i][1:].split()[0]: lines[i+1] for i in range(0, len(lines), 4) }

@pytest.fixture(scope="module")
def trimmed(tmp_path_factory):
    '''
    Simulates reads for two samples and trims them with "qc library trim", with conda
    stubbed out.

    Returns:
        (base, truth) -- the work directory, and a dict of samples to the simulation's
            truth as returned by sim.simulate_trimming_reads().
    '''
    require_tools("bbduk.sh")
    base = tmp_path_factory.mktemp("qc")
    rng = random.Random(3)
    genome = sim.random_genome(rng, { "chr1": 50000, "chr2": 30000 })
    genome = { c: "".join(s) for c, s in genome.items() }
    (base / "reads").mkdir()
    truth = { s: sim.simulate_trimming_reads(rng, genome, str(base / "reads" / s), 3000) for s in SAMPLES }

    binDir = base / "bin" # only the conda stub, since stubs/bbduk also holds a stub bbduk.sh
    binDir.mkdir()
    os.symlink(os.path.join(STUBS_DIR, "bbduk", "conda"), binDir / "conda")
    env = dict(os.environ)
    env["PATH"] = f"{binDir}{os.pathsep}{env['PATH']}"
    run_psst(["qc", "library", "trim", "-d", str(base / "reads"), "--cpus", "2", "--mem", "4G",
              "--minlen", str(MINLEN)], base)
    run_array_job(str(base / "trim_reads.sh"), str(base), env=env)
    return base, truth

@pytest.mark.parametrize("sample", SAMPLES)
def test_trim(trimmed, sample):
    base, truth = trimmed
    truth = truth[sample]
    reads1 = read_fastq(base / "trimmed" / f"{sample}_1.fq.gz")
    reads2 = read_fastq(base / "trimmed" / f"{sample}_2.fq.gz")
    assert list(reads1) == list(reads2) # pairs are kept together, in order
    assert os.path.getsize(base / "trimmed" / f"{sample}.bbduk_stats.txt") > 0

    readThrough = [ n for n in truth if int(n[1:]) % 3 == 0 ] # see sim.simulate_trimming_reads()
    others = [ n for n in truth if int(n[1:]) % 3 != 0 ]
    for name, (t1, t2) in truth.items():
        if min(len(t1), len(t2)) < MINLEN: # too short once trimmed, so discarded
            assert name not in reads1
        elif name in reads1: # no adapter or low-quality base remains
            assert t1.startswith(reads1[name]) and t2.startswith(reads2[name])

    # Untouched and low-quality pairs are trimmed exactly; read-through pairs are trimmed
    # to their insert, though bbduk's k-mer matching (hdist=1) often removes a few more bases
    exact = [ n for n in others if n in reads1 and (reads1[n], reads2[n]) == truth[n] ]
    assert len(exact) >= 0.99 * len(others)
    kept = [ n for n in readThrough if len(truth[n][0]) >= MINLEN + 10 ]
    close = [ n for n in kept if n in reads1 and len(truth[n][0]) - len(reads1[n]) <= 5 ]
    assert len(close) >= 0.95 * len(kept)

@pytest.fixture(scope="module")
def fastqc_reports(tmp_path_factory):
    '''
    Simulates paired reads for a sample with adapter read-through and a clean sample (plus the
    former's forward reads as a single-end sample), then checks and tabulates them with
    "qc library fastqc", with the module command stubbed out.

    Returns:
        (base, numPairs) -- the work directory, and a dict of samples to their number of pairs.
    '''
    require_tools("fastqc")
    pytest.importorskip("pandas")
    pytest.importorskip("xlsxwriter")
    if not os.path.isfile(os.path.join(VARSCRIPTS, "QC", "reads", "tabulate_fastqc.py")):
        pytest.skip(f"requires Various_scripts at {VARSCRIPTS} (or set $PSST_TEST_VARSCRIPTS)")
    base = tmp_path_factory.mktemp("fastqc")
    rng = random.Random(5)
    genome = { c: "".join(s) for c, s in sim.random_genome(rng, { "chr1": 40000 }).items() }
    (base / "reads").mkdir()
    (base / "single").mkdir()
    sim.simulate_trimming_reads(rng, genome, str(base / "reads" / "adapters"), 1500)
    sim.simulate_paired_reads(rng, [genome], str(base / "reads" / "clean"), coverage=5)
    shutil.copy(base / "reads" / "adapters_1.fq.gz", base / "single" / "adapters.fq.gz")
    numPairs = { s: len(read_fastq(base / "reads" / f"{s}_1.fq.gz")) for s in ["adapters", "clean"] }

    binDir = base / "bin"
    binDir.mkdir()
    os.symlink(os.path.join(STUBS_DIR, "fastqc", "module"), binDir / "module")
    env = dict(os.environ)
    env["PATH"] = f"{binDir}{os.pathsep}{env['PATH']}"
    for readsDir, extraArgs in [("reads", []), ("single", ["--singleEnd", "--outdir", "fastqc_se", "--name", "se"])]:
        run_psst(["qc", "library", "fastqc", "-d", str(base / readsDir), "--varscripts", VARSCRIPTS,
                  "--fastqc-cpus", "2", *extraArgs], base)
        name = "se" if extraArgs else "fastqc"
        run_array_job(str(base / f"{name}.sh"), str(base), env=env)
        run_job(str(base / f"{name}_tabulate.sh"), str(base), env=env)
    return base, numPairs

def test_fastqc(fastqc_reports):
    import pandas as pd
    base, numPairs = fastqc_reports
    tables = pd.read_excel(base / "fastqc" / "fastqc_summary.xlsx", sheet_name=None, index_col=0)
    assert list(tables) == ["Forward Reads", "Reverse Reads"]
    for table in tables.values():
        assert sorted(table.index) == ["adapters", "clean"]
        for sample, n in numPairs.items():
            assert table.loc[sample, "total_reads"] == n
            assert table.loc[sample, "approx_bp"] == n * 150
        # A third of the adapters sample's pairs read through into adapters
        assert table.loc["adapters", "adapter_content"] == "FAIL"
        assert table.loc["clean", "adapter_content"] == "PASS"
    single = pd.read_excel(base / "fastqc_se" / "fastqc_summary.xlsx", sheet_name=None, index_col=0)
    assert list(single) == ["Forward Reads"] # single-end reads have no reverse reads
    assert single["Forward Reads"].loc["adapters", "total_reads"] == numPairs["adapters"]
