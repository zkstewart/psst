# Integration tests which run the "pangenome create" job script with real cactus on a
# simulated pangenome of a reference (with a second haplotype of its own), two diploid samples
# and a haploid sample, check its
# VCF genotypes against the simulation, then index the graph with "map index giraffe".
# Run with 'pytest -m integration'; cactus runs from a Singularity image given by
# $PSST_TEST_CACTUS_SIF (e.g., from 'singularity pull docker://quay.io/comparative-genomics-toolkit/cactus:v3.3.0'),
# and takes ~10 minutes.

import os, random, subprocess

import pytest

import sim
from conftest import STUBS_DIR, run_psst, run_job, require_tools, tool_dir

pytestmark = pytest.mark.integration

CACTUS_SIF = os.environ.get("PSST_TEST_CACTUS_SIF")

@pytest.fixture(scope="module")
def pangenome(tmp_path_factory):
    if CACTUS_SIF == None or not os.path.isfile(CACTUS_SIF):
        pytest.skip("requires a cactus Singularity image at $PSST_TEST_CACTUS_SIF")
    require_tools("singularity", "bcftools")
    base = tmp_path_factory.mktemp("pangenome")
    rng = random.Random(13)
    genome = sim.random_genome(rng, { "chr1": 60000, "chr2": 40000 })
    records = sim.simulate_pangenome(rng, genome)
    haplotypes = sim.pangenome_haplotypes(genome, records)
    genomeFiles = [str(base / "R.hap1.fa"), str(base / "R.hap2.fa")] # the reference, and its other haplotype
    sim.write_fasta(genomeFiles[0], { c: "".join(s) for c, s in genome.items() })
    sim.write_fasta(genomeFiles[1], haplotypes[0]) # sharing P1's first haplotype
    for h, haplotype in enumerate(haplotypes):
        genomeFiles.append(str(base / f"P{h//2+1}.hap{h%2+1}.fa"))
        sim.write_fasta(genomeFiles[-1], haplotype)
    genomeFiles.append(str(base / "Iso.fa"))
    sim.write_fasta(genomeFiles[-1], haplotypes[3]) # a haploid sample, sharing P2's second haplotype

    # Run cactus from the image (/tmp, holding the test data, is bound by default), with
    # conda stubbed out
    binDir = base / "bin"
    binDir.mkdir()
    (binDir / "cactus-pangenome").write_text(f'#!/bin/bash\nexec singularity exec "{CACTUS_SIF}" cactus-pangenome "$@"\n')
    (binDir / "cactus-pangenome").chmod(0o755)
    os.symlink(os.path.join(STUBS_DIR, "cactus", "conda"), binDir / "conda")
    env = dict(os.environ)
    env["PATH"] = f"{binDir}{os.pathsep}{env['PATH']}"

    run_psst(["pangenome", "create", "-i", *genomeFiles, "-p", "pg",
              "--binaries", "local", "--cpus", "4", "--mem", "8G"], base)
    run_job(str(base / "pangenome.sh"), str(base), env=env)
    return { "base": base, "graph": base / "cactus_output" / "graph", "records": records }

def bcftools(*args):
    return subprocess.run(["bcftools", *args], check=True, capture_output=True, text=True).stdout

def test_outputs(pangenome):
    outDir = pangenome["base"] / "cactus_output"
    assert sorted(os.listdir(outDir)) == ["cactus.ok", "cactus_samples.txt", "graph"] # the job store is removed
    assert open(outDir / "cactus_samples.txt").read().split()[::2] == ["reference", "R.0", "P1.1", "P1.2", "P2.1", "P2.2", "Iso.0"]
    for name in ["pg.d2.gbz", "pg.d2.dist", "pg.full.hal", "pg.vcf.gz", "pg.wave.vcf.gz", "cactus_samples.txt"]:
        assert os.path.exists(pangenome["graph"] / name)

def test_snp_genotypes(pangenome):
    vcf = str(pangenome["graph"] / "pg.wave.vcf.gz")
    assert sorted(bcftools("query", "-l", vcf).split()) == ["Iso", "P1", "P2", "R"]
    calls = {}
    for line in bcftools("query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", vcf).splitlines():
        chrom, pos, ref, alt, *gts = line.split("\t")
        calls[(chrom.split("#")[-1], int(pos), ref, alt)] = gts
    order = bcftools("query", "-l", vcf).split()

    # Compare biallelic SNPs away from other variants, whose representation is unambiguous
    spans = [ (c, p, p + max(len(ref), 1)) for c, p, ref, _, _ in pangenome["records"] ]
    correct, wrong = 0, []
    for c, p, ref, alts, carriers in pangenome["records"]:
        if len(ref) != 1 or len(alts) != 1 or len(alts[0]) != 1:
            continue
        if any(c2 == c and (c2, p2, e2) != (c, p, p + 1) and p2 - 20 < p < e2 + 20 for c2, p2, e2 in spans):
            continue
        truth = { "P1": f"{carriers.get(0, 0)}|{carriers.get(1, 0)}", "P2": f"{carriers.get(2, 0)}|{carriers.get(3, 0)}",
                  "Iso": str(carriers.get(3, 0)), "R": str(carriers.get(0, 0)) }
        gts = calls.get((c, p + 1, ref, alts[0]))
        if gts != None and dict(zip(order, gts)) == truth:
            correct += 1
        else:
            wrong.append((c, p + 1, ref, alts[0], truth, gts))
    assert correct > 50
    assert wrong == []

def test_giraffe_index_of_graph(pangenome):
    require_tools("vg", "tabix", "bgzip")
    base, graph = pangenome["base"], pangenome["graph"]
    run_psst(["map", "index", "giraffe", "-g", str(graph / "pg.d2.gbz"), "--vcf", str(graph / "pg.vcf.gz"),
              "-o", "index", "--vg", tool_dir("vg"), "--cpus", "2", "--sort-mem", "500M"], base)
    run_job(str(base / "giraffe_index.sh"), str(base))
    assert "Reference sample: reference" in open(base / "index" / "index_info.txt").read()
    assert [ l.split()[0] for l in open(base / "index" / "reference.fa.fai") ] == ["reference#0#chr1", "reference#0#chr2"]
