# Integration tests of "qc library" modes on simulated reads, checking their outputs against
# the simulation's truth: "trim" runs bbduk on reads which read through into adapters or end
# in low-quality bases (requiring bbduk.sh on PATH), "fastqc" runs FastQC then Various_scripts'
# tabulate_fastqc.py (requiring fastqc on PATH, pandas, xlsxwriter and openpyxl), and "kraken2" classifies
# reads from known taxa against a small database built here, then combines and reformats the
# reports (requiring kraken2, KrakenTools and Krona on PATH, and pandas). Various_scripts is
# found at $PSST_TEST_VARSCRIPTS, or else beside the psst repository.

import gzip, os, random, shutil, subprocess

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
    pytest.importorskip("openpyxl") # for pandas to read the table
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

# A small taxonomy for a kraken2 database: (taxid, parent, rank, name)
TAXA = [ (1, 1, "no rank", "root"), (2, 1, "superkingdom", "Bacteria"), (3, 2, "genus", "Simulobacter"),
         (101, 3, "species", "Simulobacter alpha"), (102, 3, "species", "Simulobacter beta"),
         (4, 1, "superkingdom", "Eukaryota"), (201, 4, "species", "Simulata host") ]

# The number of read pairs simulated from each species' genome (None being random reads,
# which kraken2 leaves unclassified)
KRAKEN2_SAMPLES = { "S1": { 101: 600, 201: 300, None: 100 }, "S2": { 101: 500, 102: 500 } }

def build_kraken2_db(rng, dbDir):
    '''
    Builds a kraken2 database of a random genome per species in TAXA, with a taxonomy in
    the format of NCBI's nodes.dmp and names.dmp.

    Returns:
        genomes -- a dict of taxids to genome sequence strings.
    '''
    os.makedirs(os.path.join(dbDir, "taxonomy"))
    with open(os.path.join(dbDir, "taxonomy", "nodes.dmp"), "w") as nodes, \
         open(os.path.join(dbDir, "taxonomy", "names.dmp"), "w") as names:
        for taxid, parent, rank, name in TAXA:
            nodes.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|\n")
            names.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    genomes = { taxid: "".join(sim.random_genome(rng, { "c": 30000 })["c"]) for taxid in [101, 102, 201] }
    libraryFasta = os.path.join(dbDir, "genomes.fa")
    sim.write_fasta(libraryFasta, { f"g{taxid}|kraken:taxid|{taxid}": seq for taxid, seq in genomes.items() })
    for args in [["--add-to-library", libraryFasta], ["--build", "--threads", "2"]]:
        subprocess.run(["kraken2-build", *args, "--db", dbDir], check=True, capture_output=True)
    return genomes

@pytest.fixture(scope="module")
def kraken2_results(tmp_path_factory):
    '''
    Builds a small kraken2 database, simulates reads for samples mixing its species (and
    random reads), then classifies and tabulates them with "qc library kraken2", with conda
    stubbed out.

    Returns:
        base -- the work directory.
    '''
    require_tools("k2", "kraken2-build", "combine_kreports.py", "kreport2krona.py", "ktImportText")
    pytest.importorskip("pandas")
    if not os.path.isfile(os.path.join(VARSCRIPTS, "kraken", "kreport_reformat.py")):
        pytest.skip(f"requires Various_scripts at {VARSCRIPTS} (or set $PSST_TEST_VARSCRIPTS)")
    base = tmp_path_factory.mktemp("kraken2")
    rng = random.Random(2)
    genomes = build_kraken2_db(rng, str(base / "db"))
    (base / "reads").mkdir()
    for sample, mix in KRAKEN2_SAMPLES.items():
        pairs = []
        for taxid, n in mix.items():
            genome = genomes[taxid] if taxid != None else "".join(rng.choice(sim.BASES) for _ in range(30000))
            for _ in range(n):
                start = rng.randint(0, len(genome) - 400)
                fragment = genome[start:start+400]
                pairs.append((fragment[:150], sim.revcomp(fragment)[:150]))
        rng.shuffle(pairs)
        with gzip.open(base / "reads" / f"{sample}_1.fq.gz", "wt") as out1, gzip.open(base / "reads" / f"{sample}_2.fq.gz", "wt") as out2:
            for i, (r1, r2) in enumerate(pairs):
                out1.write(f"@p{i}\n{r1}\n+\n{'I' * len(r1)}\n")
                out2.write(f"@p{i}\n{r2}\n+\n{'I' * len(r2)}\n")

    binDir = base / "bin"
    binDir.mkdir()
    os.symlink(os.path.join(STUBS_DIR, "kraken2", "conda"), binDir / "conda")
    env = dict(os.environ)
    env["PATH"] = f"{binDir}{os.pathsep}{env['PATH']}"
    run_psst(["qc", "library", "kraken2", "-d", str(base / "reads"), "--db", str(base / "db"),
              "--varscripts", VARSCRIPTS, "--classify-cpus", "2"], base)
    run_array_job(str(base / "kraken2.sh"), str(base), env=env)
    run_job(str(base / "kraken2_tabulate.sh"), str(base), env=env)
    return base

def test_kraken2(kraken2_results):
    import pandas as pd
    base = kraken2_results
    for sample in KRAKEN2_SAMPLES:
        sampleDir = base / "kraken2" / sample
        assert sorted(os.listdir(sampleDir)) == sorted([f"{sample}_{r}_{c}.fastq" for r in [1, 2] for c in ["classified", "unclassified"]] +
                                                       [f"{sample}_output.txt", f"{sample}_report.html", f"{sample}_report.krona", f"{sample}_report.txt"])
        assert len(read_fastq_plain(sampleDir / f"{sample}_1_unclassified.fastq")) == KRAKEN2_SAMPLES[sample].get(None, 0)
    
    # Each sample's fraction of reads per taxon matches the simulated mix
    table = pd.read_csv(base / "kraken2" / "kraken_summary.reformat.tsv", sep="\t")
    table["taxon"] = table["taxon"].str.strip()
    table = table.set_index("taxon")
    names = { taxid: name for taxid, _, _, name in TAXA }
    for sample, mix in KRAKEN2_SAMPLES.items():
        total = sum(mix.values())
        assert table.loc["unclassified", sample] == mix.get(None, 0) / total
        for taxid in [101, 102, 201]:
            assert table.loc[names[taxid], sample] == mix.get(taxid, 0) / total
    combined = open(base / "kraken2" / "kraken_summary.tsv").read()
    assert "#Number of Samples: 2" in combined

def read_fastq_plain(fileName):
    with open(fileName) as fileIn:
        return fileIn.read().splitlines()[::4]
