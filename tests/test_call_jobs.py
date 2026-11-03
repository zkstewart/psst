# Tests of the "call" job scripts' own logic (chunking, resuming, concatenation and the
# arguments given to bcftools), using a stub bcftools so that no real data are needed.

import gzip, os, random

import pytest

from conftest import run_psst, run_job, array_size

@pytest.fixture
def genome(tmp_path):
    '''
    A fake genome .fai like a chromosome-level assembly with many small scaffolds, plus
    a directory of (fake) indexed BAMs.
    '''
    rng = random.Random(1)
    contigs = [ (f"chr{i}", rng.randint(30, 60) * 10**6) for i in range(1, 13) ]
    contigs += [ (f"scaffold_{i}", rng.randint(2000, 200000)) for i in range(300) ]
    genomeFile = tmp_path / "genome.fa"
    genomeFile.write_text(">placeholder\nA\n")
    offset = 0
    with open(tmp_path / "genome.fa.fai", "w") as fileOut:
        for name, length in contigs:
            fileOut.write(f"{name}\t{length}\t{offset}\t60\t61\n")
            offset += length
    bamDir = tmp_path / "bams"
    bamDir.mkdir()
    for sample in ["s1", "s2"]:
        (bamDir / f"{sample}.sorted.md.bam").write_text("")
        (bamDir / f"{sample}.sorted.md.bam.bai").write_text("")
    workDir = tmp_path / "work"
    workDir.mkdir()
    return { "genome": str(genomeFile), "bams": str(bamDir), "contigs": contigs, "work": str(workDir) }

def run_calling(genome, env, chunks, extraArgs=[]):
    run_psst(["call", "short", "-r", genome["genome"], "-b", genome["bams"], "-p", "out",
              "--chunks", str(chunks), *extraArgs], genome["work"])
    script = os.path.join(genome["work"], "call_short.sh")
    for index in range(1, array_size(script) + 1):
        run_job(script, genome["work"], index, env=env)
    return script

def test_chunks_cover_genome_in_order(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    run_calling(genome, env, 50)
    chunkDir = os.path.join(genome["work"], "variants", "out_chunks")
    regions, emptyChunks, sizes = [], [], []
    for chunk in range(1, 51):
        regionsFile = os.path.join(stubDir, f"regions_{chunk}.txt")
        if os.path.exists(regionsFile):
            lines = [ l.split("\t") for l in open(regionsFile).read().splitlines() ]
            regions += [ (name, int(end)) for name, start, end in lines ]
            sizes.append(sum(int(end) for _, _, end in lines))
            assert os.path.exists(os.path.join(chunkDir, f"chunk_{chunk}.vcf.gz"))
        else:
            emptyChunks.append(chunk)
            assert os.path.exists(os.path.join(chunkDir, f"chunk_{chunk}.empty"))
    assert regions == genome["contigs"] # every contig once, in genome order
    assert emptyChunks != [] # 12 large chromosomes can't fill 50 chunks
    assert max(sizes) == max(length for _, length in genome["contigs"]) # chunks are no larger than needed
    assert not any(f.endswith(".regions.txt") or ".incomplete." in f for f in os.listdir(chunkDir))

def test_single_chunk_has_no_array(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    script = run_calling(genome, env, 1)
    assert array_size(script) == 1
    assert len(open(os.path.join(stubDir, "regions_1.txt")).read().splitlines()) == len(genome["contigs"])

def test_completed_chunks_are_skipped(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    script = run_calling(genome, env, 5)
    result = run_job(script, genome["work"], 1, env=env)
    assert "Chunk 1 has already been called; skipping it" in result.stdout

def test_concat_joins_chunks_in_order(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    run_calling(genome, env, 50)
    run_job(os.path.join(genome["work"], "call_short_concat.sh"), genome["work"], env=env)
    variants = os.path.join(genome["work"], "variants")
    chunks = [ int(line.split()[1]) for line in open(os.path.join(variants, "out.vcf.gz")).read().splitlines() ]
    assert chunks == sorted(chunks) and len(chunks) == len(set(chunks))
    assert sorted(os.listdir(variants)) == ["out.vcf.gz", "out.vcf.gz.csi", "out.vcf.gz.tbi"] # chunks removed
    assert "merge" not in open(os.path.join(stubDir, "calls.log")).read()

def test_concat_refuses_missing_chunk(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    run_calling(genome, env, 5)
    os.remove(os.path.join(genome["work"], "variants", "out_chunks", "chunk_3.vcf.gz"))
    result = run_job(os.path.join(genome["work"], "call_short_concat.sh"), genome["work"], env=env, check=False)
    assert result.returncode != 0
    assert "Chunk 3 has not been called" in result.stdout
    assert not os.path.exists(os.path.join(genome["work"], "variants", "out.vcf.gz"))

def test_concat_refuses_existing_output(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    run_calling(genome, env, 5)
    open(os.path.join(genome["work"], "variants", "out.vcf.gz"), "w").close()
    result = run_job(os.path.join(genome["work"], "call_short_concat.sh"), genome["work"], env=env, check=False)
    assert result.returncode != 0
    assert "already exists" in result.stdout

@pytest.mark.parametrize("extraArgs,expected,unexpected", [
    ([], ["-q 20", "-Q 20", "-d 250", "INFO/AD", "call -m -v"], [" -A ", "--indels-cns", "-X ", "-T "]),
    (["--hic", "--indels-cns", "-q", "30"], [" -A ", "--indels-cns", "-q 30"], []),
    (["--ploidy", "1"], ["--ploidy 1"], []),
])
def test_call_short_arguments(extraArgs, expected, unexpected, genome, stub_bcftools):
    env, stubDir = stub_bcftools
    run_calling(genome, env, 1, extraArgs)
    log = open(os.path.join(stubDir, "calls.log")).read()
    for text in expected:
        assert text in log
    for text in unexpected:
        assert text not in log

def test_call_nanopore_arguments(genome, stub_bcftools):
    env, stubDir = stub_bcftools
    run_psst(["call", "nanopore", "-r", genome["genome"], "-b", genome["bams"], "-bs", ".sorted.md.bam",
              "--chunks", "1"], genome["work"])
    run_job(os.path.join(genome["work"], "call_nanopore.sh"), genome["work"], env=env)
    log = open(os.path.join(stubDir, "calls.log")).read().splitlines()
    mpileup = next(l for l in log if l.startswith("mpileup")) # pipeline processes log in any order
    assert " -X ont-sup -q 20 " in mpileup # the profile comes first, so later options override it
    assert " -Q " not in mpileup

def test_call_giraffe_genotypes_targets_and_merges_sites(genome, stub_bcftools, tmp_path):
    env, stubDir = stub_bcftools
    index = tmp_path / "pg_giraffe_index"
    index.mkdir()
    os.symlink(genome["genome"], index / "reference.fa")
    os.symlink(genome["genome"] + ".fai", index / "reference.fa.fai")
    (index / "targets.tsv.gz").write_text("")
    with gzip.open(index / "sites.vcf.gz", "wt") as fileOut:
        fileOut.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\n")
    run_psst(["call", "giraffe", "-i", str(index), "-b", genome["bams"], "--chunks", "1"], genome["work"])
    run_job(os.path.join(genome["work"], "call_giraffe.sh"), genome["work"], env=env)
    run_job(os.path.join(genome["work"], "call_giraffe_concat.sh"), genome["work"], env=env)
    log = open(os.path.join(stubDir, "calls.log")).read().splitlines()
    targets, sites = str(index / "targets.tsv.gz"), str(index / "sites.vcf.gz")
    mpileup = next(l for l in log if l.startswith("mpileup"))
    assert f"-T {targets}" in mpileup and "INFO/AD" not in mpileup
    assert " -A " in mpileup # giraffe can leave pairs spanning deletions unflagged as proper
    assert any(l.startswith(f"call -m -C alleles -T {targets} -A ") for l in log)
    assert any(l.startswith("norm -f ") and " -m -any " in l for l in log)
    merge = next(l for l in log if l.startswith("merge"))
    assert merge.startswith("merge -m none --force-samples") and merge.endswith(sites)
    variants = os.path.join(genome["work"], "variants")
    assert sorted(os.listdir(variants)) == ["pg.vcf.gz", "pg.vcf.gz.csi", "pg.vcf.gz.tbi"] # intermediate removed
