# Integration tests which run the "assembly revcomp" and "assembly scaffolding" job scripts
# with annotarium, minimap2 and RagTag on two simulated haplotypes of a reference, and check
# that each partition holds one haplotype of every reference contig, named after and
# oriented to match it. Run with 'pytest -m integration'; annotarium is found at
# $PSST_TEST_ANNOTARIUM, or else beside the psst repository.

import os, random, re

import pytest

import sim
from conftest import TESTS_DIR, run_psst, run_job, require_tools, tool_dir

pytestmark = pytest.mark.integration

ANNOTARIUM = os.environ.get("PSST_TEST_ANNOTARIUM",
                            os.path.join(os.path.dirname(os.path.dirname(TESTS_DIR)), "annotarium"))

@pytest.fixture(scope="module")
def haplotypes(tmp_path_factory):
    '''
    Simulates a reference genome and two haplotypes, each with 1% SNP divergence from it.
    '''
    require_tools("minimap2")
    if not os.path.isfile(os.path.join(ANNOTARIUM, "annotarium.py")):
        pytest.skip(f"requires annotarium at {ANNOTARIUM} (or set $PSST_TEST_ANNOTARIUM)")
    base = tmp_path_factory.mktemp("assembly")
    rng = random.Random(5)
    genome = sim.random_genome(rng, { "chr1": 300000, "chr2": 200000, "chr3": 120000 })
    sim.write_fasta(base / "ref.fa", { c: "".join(s) for c, s in genome.items() })
    haps = [ sim.diverged_haplotype(rng, genome, 0.01) for _ in range(2) ]
    return { "base": base, "rng": rng, "haps": haps }

def write_contigs(fileName, contigs, rng):
    '''
    Writes contigs in a shuffled order, alternately reverse complemented.

    Parameters:
        contigs -- a dict of contig IDs to sequence strings.
    '''
    order = list(contigs)
    rng.shuffle(order)
    sim.write_fasta(fileName, { c: sim.revcomp(contigs[c]) if i % 2 == 0 else contigs[c]
                                for i, c in enumerate(order) })

def run_mode(mode, haplotypes, contigFiles, args=[]):
    base = haplotypes["base"]
    workDir = base / mode
    workDir.mkdir()
    run_psst(["assembly", mode, "-i", *contigFiles, "-r", str(base / "ref.fa"), "--annotarium", ANNOTARIUM,
              "--minimap2", tool_dir("minimap2"), "--cpus", "2", *args], workDir)
    run_job(str(workDir / f"{mode}.sh"), str(workDir))
    return workDir / "output"

def test_revcomp(haplotypes):
    base, rng, haps = haplotypes["base"], haplotypes["rng"], haplotypes["haps"]
    contigFiles = []
    for h, hap in enumerate(haps):
        contigFiles.append(str(base / f"hap{h+1}.fa"))
        write_contigs(contigFiles[-1], { f"h{h+1}_tig{i}": hap[c] for i, c in enumerate(hap) }, rng)
    outputDir = run_mode("revcomp", haplotypes, contigFiles)

    assert sorted(os.listdir(outputDir)) == ["partition1.fasta", "partition2.fasta", "phase"]
    used = { c: [] for c in haps[0] }
    for partition in [1, 2]:
        seqs = sim.read_fasta(outputDir / f"partition{partition}.fasta")
        assert list(seqs) == list(haps[0]) # named after the reference contigs, in their order
        for c, seq in seqs.items():
            matches = [ h for h, hap in enumerate(haps) if hap[c] == seq ] # oriented as the reference
            assert len(matches) == 1
            used[c] += matches
    assert all(sorted(hs) == [0, 1] for hs in used.values()) # each haplotype is in one partition

def test_scaffolding(haplotypes, tmp_path):
    require_tools("ragtag.py")
    base, rng, haps = haplotypes["base"], haplotypes["rng"], haplotypes["haps"]
    # Break each haplotype at the same points per contig, so that either haplotype's piece
    # of a region can be phased with either haplotype's pieces of other regions; pieces are
    # longer than 'assembly phase' requires alignments to be (10 kb)
    cuts = { c: [rng.randrange(20000, len(seq) // 2 - 10000), rng.randrange(len(seq) // 2 + 10000, len(seq) - 20000)]
             for c, seq in haps[0].items() }
    pieces = {}
    contigFiles = []
    for h, hap in enumerate(haps):
        hapPieces = {}
        for c, seq in hap.items():
            bounds = [0] + cuts[c] + [len(seq)]
            for i in range(3):
                hapPieces[f"h{h+1}_{c}_{i}"] = seq[bounds[i]:bounds[i+1]]
        pieces.update(hapPieces)
        contigFiles.append(str(base / f"hap{h+1}_pieces.fa"))
        write_contigs(contigFiles[-1], hapPieces, rng)
    p2dpDir = tmp_path / "paf2dotplot" # a stub which writes the plots paf2dotplot.r would
    p2dpDir.mkdir()
    (p2dpDir / "paf2dotplot.r").write_text('#!/bin/bash\nwhile [[ $1 != -o ]]; do shift; done\n' +
                                           'touch "$2.pdf" "$2.png"\n')
    (p2dpDir / "paf2dotplot.r").chmod(0o755)
    outputDir = run_mode("scaffolding", haplotypes, contigFiles, ["--paf2dotplot", str(p2dpDir)])

    for partition in [1, 2]:
        for suffix in [".fasta", ".dotplot.paf", ".dotplot.pdf", ".dotplot.png"]:
            assert os.path.exists(outputDir / f"partition{partition}{suffix}")
    used = []
    for partition in [1, 2]:
        seqs = sim.read_fasta(outputDir / f"partition{partition}.fasta")
        assert list(seqs) == list(haps[0]) # every piece was placed; scaffolds are named after the reference
        for c, seq in seqs.items():
            parts = re.split("N+", seq)
            assert len(parts) == 3
            for i, part in enumerate(parts): # in reference order and orientation
                matches = [ name for name in (f"h1_{c}_{i}", f"h2_{c}_{i}") if pieces[name] == part ]
                assert len(matches) == 1
                used += matches
    assert sorted(used) == sorted(pieces) # each piece is in one partition
