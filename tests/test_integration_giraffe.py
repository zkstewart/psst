# Integration tests which run the "map index giraffe", "map giraffe" and "call giraffe"
# job scripts with real tools on a simulated pangenome, and check that progeny are
# genotyped correctly at every pangenome variant site. Run with 'pytest -m integration'.

import gzip, json, os, random, subprocess

import pytest

import sim
from conftest import run_psst, run_job, run_array_job, require_tools, tool_dir

pytestmark = pytest.mark.integration

PROGENY = ["F0", "F1", "F2", "F3"]

@pytest.fixture(scope="module")
def pangenome(tmp_path_factory):
    '''
    Simulates a pangenome of a reference and two diploid parents (P1 and P2), progeny
    each inheriting one haplotype from each parent per contig, and progeny reads; then
    runs the index, mapping and calling jobs.
    '''
    require_tools("vg", "samtools", "bcftools", "bgzip", "tabix")
    base = tmp_path_factory.mktemp("giraffe")
    rng = random.Random(11)
    genome = sim.random_genome(rng, { "chr1": 120000, "chr2": 70000, "chr3": 25000 })
    sim.write_fasta(base / "ref.fa", { c: "".join(s) for c, s in genome.items() })
    records = sim.simulate_pangenome(rng, genome)
    sim.write_pangenome_vcf(str(base / "pangenome.vcf.gz"), genome, records)
    sim.build_pangenome_gbz(str(base), str(base / "ref.fa"), str(base / "pangenome.vcf.gz"), "pangenome.gbz")

    haplotypes = sim.pangenome_haplotypes(genome, records)
    (base / "reads").mkdir()
    inheritance = {}
    for sample in PROGENY:
        inheritance[sample] = { c: (rng.randint(0, 1), rng.randint(2, 3)) for c in genome }
        progenyHaplotypes = [ { c: haplotypes[inheritance[sample][c][i]][c] for c in genome } for i in (0, 1) ]
        sim.simulate_paired_reads(rng, progenyHaplotypes, str(base / "reads" / sample), coverage=15)

    vgDir = tool_dir("vg")
    run_psst(["map", "index", "giraffe", "-g", str(base / "pangenome.gbz"), "--vcf", str(base / "pangenome.vcf.gz"),
              "-o", "index", "--vg", vgDir, "--cpus", "2", "--sort-mem", "500M"], base)
    run_job(str(base / "giraffe_index.sh"), str(base))
    run_psst(["map", "giraffe", "-i", "index", "-d", str(base / "reads"), "--vg", vgDir, "--cpus", "2",
              "--sort-mem", "500M"], base)
    run_array_job(str(base / "map_giraffe.sh"), str(base))
    run_psst(["call", "giraffe", "-i", "index", "-b", "mapped", "--chunks", "2"], base)
    run_array_job(str(base / "call_giraffe.sh"), str(base))
    run_job(str(base / "call_giraffe_concat.sh"), str(base))
    return { "base": base, "inheritance": inheritance }

def bcftools(*args, stdin=None):
    return subprocess.run(["bcftools", *args], check=True, capture_output=True, text=True, input=stdin).stdout

def haplotype_atoms(base):
    '''
    Returns:
        atoms -- a list, for each pangenome haplotype, of the set of (CHROM, POS, REF, ALT)
                 variants it carries, named, split and atomised as in the giraffe index.
    '''
    renames = "".join(f"{line.split()[0].split('#')[-1]}\t{line.split()[0]}\n"
                      for line in open(base / "index" / "reference.fa.fai"))
    (base / "rename.txt").write_text(renames)
    vcf = bcftools("view", str(base / "pangenome.vcf.gz"))
    atoms = []
    for h in range(4):
        lines = []
        for line in vcf.splitlines():
            if line.startswith("##"):
                lines.append(line)
                continue
            f = line.split("\t")
            if line.startswith("#"):
                lines.append("\t".join(f[:9] + ["H"]))
                continue
            allele = (f[9].split("|") + f[10].split("|"))[h]
            if allele != "0":
                lines.append("\t".join(f[:4] + [f[4].split(",")[int(allele)-1]] + f[5:9] + ["1"]))
        normalised = subprocess.run(
            f"bcftools annotate --rename-chrs {base}/rename.txt -Ou | bcftools norm -f {base}/index/reference.fa " +
            "-a --atom-overlaps . -m -any -Ou | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'",
            shell=True, check=True, capture_output=True, text=True, input="\n".join(lines) + "\n").stdout
        atoms.append(set(tuple(l.split("\t")) for l in normalised.splitlines()))
    return atoms

def test_index_contents(pangenome):
    index = pangenome["base"] / "index"
    assert sorted(os.listdir(index)) == sorted(["index_info.txt", "pangenome.dist", "pangenome.gbz", "pangenome.min",
        "pangenome.zipcodes", "reference.fa", "reference.fa.fai", "sites.vcf.gz", "sites.vcf.gz.tbi",
        "targets.tsv.gz", "targets.tsv.gz.tbi"])
    assert "Reference sample: reference" in open(index / "index_info.txt").read() # detected from the graph
    assert [ l.split()[0] for l in open(index / "reference.fa.fai") ] == ["reference#0#chr1", "reference#0#chr2",
                                                                           "reference#0#chr3"]
    sites = str(index / "sites.vcf.gz")
    assert bcftools("view", "-H", "-v", "mnps", sites) == "" # MNPs are atomised into SNPs
    assert bcftools("view", "-H", "-m3", sites) == "" # multiallelic records are split
    assert bcftools("view", "-H", "-i", 'TYPE="other"', sites) != "" # complex alleles are kept
    targets = [ l.split("\t") for l in gzip.open(index / "targets.tsv.gz", "rt").read().splitlines() ]
    complexPositions = set(tuple(l.split("\t")[:2]) for l in bcftools("view", "-H", "-i", 'TYPE="other"', sites).splitlines())
    assert not any((chrom, pos) in complexPositions and len(set(len(a) for a in alleles.split(","))) > 1
                   and alleles.split(",")[0][0] != alleles.split(",")[1][0] for chrom, pos, alleles in targets)
    snpTargets = [ (chrom, pos) for chrom, pos, alleles in targets if set(len(a) for a in alleles.split(",")) == {1} ]
    assert len(snpTargets) == len(set(snpTargets)) # SNP alleles at one position are joined into one target

def test_bams_use_reference_paths(pangenome):
    base = pangenome["base"]
    header = subprocess.run(["samtools", "view", "-H", str(base / "mapped" / "F0.sorted.md.bam")],
                            check=True, capture_output=True, text=True).stdout
    contigs = [ l.split("\t")[1][3:] for l in header.splitlines() if l.startswith("@SQ") ]
    assert contigs == [ l.split()[0] for l in open(base / "index" / "reference.fa.fai") ]

def test_genotypes(pangenome):
    base, inheritance = pangenome["base"], pangenome["inheritance"]
    atoms = haplotype_atoms(base)
    output = str(base / "variants" / "index.vcf.gz")
    sites = bcftools("query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(base / "index" / "sites.vcf.gz")).splitlines()
    assert bcftools("query", "-l", output).split() == PROGENY + ["P1", "P2"]

    rows = bcftools("query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%GT]\n", output).splitlines()
    assert sorted("\t".join(r.split("\t")[:4]) for r in rows) == sorted(sites) # every site, once
    def in_own_deletion(key, haps):
        # Whether a site lies within a deletion carried by one of a sample's haplotypes
        position = int(key[1])
        return any(chrom == key[0] and len(ref) > len(alt) and int(pos) < position < int(pos) + len(ref)
                   for h in haps for chrom, pos, ref, alt in atoms[h])
    correct, wrong, invariant = 0, [], 0
    for row in rows:
        f = row.split("\t")
        key, variantType, genotypes = tuple(f[:4]), f[4], dict(zip(PROGENY + ["P1", "P2"], f[5:]))
        contig = key[0].split("#")[-1]
        if variantType == "OTHER": # complex alleles aren't genotyped
            assert all(genotypes[s] == "./." for s in PROGENY)
            continue
        for parent, haps in [("P1", (0, 1)), ("P2", (2, 3))]:
            assert genotypes[parent].replace("|", "/").count("1") == sum(key in atoms[h] for h in haps)
        allHomRef = True
        for sample in PROGENY:
            truth = sum(key in atoms[h] for h in inheritance[sample][contig])
            called = genotypes[sample].count("1")
            allHomRef = allHomRef and genotypes[sample] == "0/0"
            if called == truth:
                correct += 1
            else:
                wrong.append((key, sample, truth, genotypes[sample], in_own_deletion(key, inheritance[sample][contig])))
        invariant += allHomRef
    # The only errors expected are at sites within one of a sample's own deletions, where
    # it is hemizygous but a diploid caller must report two alleles
    assert [ w for w in wrong if not w[4] ] == []
    assert correct > 2500
    assert invariant > 0 # sites where all progeny are homozygous for the reference are kept
