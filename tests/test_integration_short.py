# Integration tests which run the "map index bwa", "map short", "map nanopore" and "call"
# job scripts with real tools on simulated data. Run with 'pytest -m integration'.

import os, random, subprocess

import pytest

import sim
from conftest import run_psst, run_job, run_array_job, require_tools, tool_dir

pytestmark = pytest.mark.integration

SAMPLES = ["sA", "sB", "sC"]
HIC_DUPLICATE_FRACTION = 0.05

@pytest.fixture(scope="module")
def data(tmp_path_factory):
    '''
    Simulates a genome and three diploid samples sequenced with Illumina WGS, Hi-C and
    nanopore reads, then indexes the genome with 'map index bwa'.
    '''
    require_tools("bwa", "samtools", "bcftools", "bgzip", "tabix")
    base = tmp_path_factory.mktemp("short")
    rng = random.Random(7)
    genome = sim.random_genome(rng, { "chr1": 150000, "chr2": 100000, "chr3": 60000,
                                      "scaf4": 20000, "scaf5": 8000, "scaf6": 3000 }, homopolymers=True)
    sim.write_fasta(base / "ref.fa", { c: "".join(s) for c, s in genome.items() })
    population = sim.simulate_population(rng, genome, SAMPLES)
    for readType in ["wgs", "hic", "ont"]:
        (base / readType).mkdir()
    for sample, haplotypes in population.items():
        sim.simulate_paired_reads(rng, haplotypes, str(base / "wgs" / sample), coverage=15)
        sim.simulate_paired_reads(rng, haplotypes, str(base / "hic" / sample), coverage=15, hic=True,
                                  duplicateFraction=HIC_DUPLICATE_FRACTION)
        sim.simulate_long_reads(rng, haplotypes, str(base / "ont" / f"{sample}.fq.gz"), coverage=12)

    run_psst(["map", "index", "bwa", "-r", str(base / "ref.fa")], base)
    run_job(str(base / "bwa_index.sh"), str(base))
    subprocess.run(["samtools", "faidx", str(base / "ref.fa")], check=True)
    return base

@pytest.fixture(scope="module")
def mapped(data):
    '''
    Maps the WGS and Hi-C reads with "map short", and the nanopore reads with "map nanopore"
    if minimap2 is available.
    '''
    for readType, extraArgs in [("wgs", []), ("hic", ["--hic"])]:
        run_psst(["map", "short", "-r", str(data / "ref.fa"), "-d", str(data / readType), "--outdir",
                  f"bam_{readType}", "--cpus", "2", "--sort-mem", "500M", "--name", readType,
                  "--script", f"map_{readType}.sh", *extraArgs], data)
        run_array_job(str(data / f"map_{readType}.sh"), str(data))
    return data

def call(workDir, mode, args, name):
    '''
    Runs a "call" mode's jobs, returning the output VCF's location.
    '''
    run_psst(["call", mode, *args, "-p", name, "--name", name], workDir)
    run_array_job(os.path.join(workDir, f"{name}.sh"), str(workDir))
    run_job(os.path.join(workDir, f"{name}_concat.sh"), str(workDir))
    return os.path.join(workDir, "variants", f"{name}.vcf.gz")

def records(vcfFile, extra=""):
    result = subprocess.run(["bcftools", "query", "-f", f"%CHROM\\t%POS\\t%REF\\t%ALT{extra}\\n", vcfFile],
                            check=True, capture_output=True, text=True)
    return result.stdout.splitlines()

def flagstat(bamFile):
    return subprocess.run(["samtools", "flagstat", bamFile], check=True, capture_output=True, text=True).stdout

@pytest.fixture(scope="module")
def wgs_calls(mapped):
    return call(mapped, "short", ["-r", str(mapped / "ref.fa"), "-b", str(mapped / "bam_wgs"), "--chunks", "4"], "wgs")

def test_map_short(mapped):
    for sample in SAMPLES:
        bamFile = str(mapped / "bam_wgs" / f"{sample}.sorted.md.bam")
        assert os.path.exists(bamFile + ".bai")
        header = subprocess.run(["samtools", "view", "-H", bamFile], check=True, capture_output=True, text=True).stdout
        assert f"SM:{sample}" in header
        stats = flagstat(bamFile)
        properlyPaired = float(stats.split(" properly paired (")[1].split("%")[0])
        assert properlyPaired > 95
    result = run_job(str(mapped / "map_wgs.sh"), str(mapped), 1) # completed samples are skipped
    assert "already exists; skipping sA" in result.stdout

def test_map_short_hic(mapped):
    bamFile = str(mapped / "bam_hic" / "sA.sorted.md.bam")
    header = subprocess.run(["samtools", "view", "-H", bamFile], check=True, capture_output=True, text=True).stdout
    assert "bwa mem -t 2 -5SP" in header
    stats = flagstat(bamFile)
    assert "0 + 0 properly paired" in stats
    total = int(stats.split(" + ")[0])
    duplicates = int(stats.split(" + 0 duplicates")[0].split("\n")[-1])
    expected = round(total * HIC_DUPLICATE_FRACTION / (1 + HIC_DUPLICATE_FRACTION))
    assert 0.95 * expected <= duplicates <= 1.05 * expected # the duplicated pairs, and few others

def test_call_short_matches_whole_genome_calling(mapped, wgs_calls):
    # Chunked calling must give exactly what one whole-genome pipeline gives
    bams = sorted(str(mapped / "bam_wgs" / f) for f in os.listdir(mapped / "bam_wgs") if f.endswith(".bam"))
    direct = str(mapped / "direct.vcf.gz")
    subprocess.run(f"bcftools mpileup -f {mapped}/ref.fa -q 20 -Q 20 -d 250 -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD " +
                   f"-Ou {' '.join(bams)} | bcftools call -m -v -a GQ,GP -Ou | bcftools norm -f {mapped}/ref.fa " +
                   f"-Oz -o {direct}", shell=True, check=True, stderr=subprocess.DEVNULL)
    extra = "\\t%QUAL[\\t%GT:%PL:%AD]"
    assert records(wgs_calls, extra) == records(direct, extra)
    assert len(records(wgs_calls)) > 150 # ~1 variant per 1.5 kb was simulated
    assert subprocess.run(["bcftools", "query", "-l", wgs_calls], check=True, capture_output=True,
                          text=True).stdout.split() == SAMPLES

def test_call_short_output_is_normalised(mapped, wgs_calls):
    result = subprocess.run(["bcftools", "norm", "-f", str(mapped / "ref.fa"), "-Ou", "-o", "/dev/null", wgs_calls],
                            check=True, capture_output=True, text=True)
    assert "/0/0/0/0/0/0\n" in result.stderr.replace("\t", "/") # nothing left to split, join or realign

def test_call_short_hic(mapped, wgs_calls):
    expected = set(records(wgs_calls))
    hic = set(records(call(mapped, "short", ["-r", str(mapped / "ref.fa"), "-b", str(mapped / "bam_hic"),
                                             "--hic", "--chunks", "2"], "hic")))
    noOrphans = set(records(call(mapped, "short", ["-r", str(mapped / "ref.fa"), "-b", str(mapped / "bam_hic"),
                                                   "--chunks", "2"], "hic_no_orphans")))
    assert len(hic & expected) >= 0.95 * len(expected)
    assert len(noOrphans & expected) < 0.5 * len(expected) # without -A, Hi-C pairs are skipped

def test_call_nanopore(mapped, wgs_calls):
    require_tools("minimap2")
    run_psst(["map", "nanopore", "-r", str(mapped / "ref.fa"), "-d", str(mapped / "ont"), "--outdir", "bam_ont",
              "--minimap2", tool_dir("minimap2"), "--cpus", "2", "--sort-mem", "500M"], mapped)
    run_array_job(str(mapped / "map_nanopore.sh"), str(mapped))
    expected = set(r for r in records(wgs_calls) if len(r.split("\t")[2]) == len(r.split("\t")[3]) == 1)
    ont = set(records(call(mapped, "nanopore", ["-r", str(mapped / "ref.fa"), "-b", str(mapped / "bam_ont"),
                                                "--chunks", "2"], "ont")))
    assert len(ont & expected) >= 0.9 * len(expected) # SNPs, as nanopore indels are less reliable
