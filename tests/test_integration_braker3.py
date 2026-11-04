# Integration test which runs the "assembly braker3" job script with the real BRAKER3 image on
# the example data inside it (a 1 Mb genome, proteins and paired RNAseq, as local read files),
# checking its genes against the image's reference results for the same data. Gene prediction
# needs realistic genes for GeneMark-ETP to train on, so this uses BRAKER3's example rather
# than simulated data.
# Run with 'pytest -m integration'; BRAKER3 runs from a Singularity image given by
# $PSST_TEST_BRAKER3_SIF (e.g., from 'singularity pull docker://teambraker/braker3:latest'),
# needs internet access (compleasm downloads a BUSCO lineage), and takes ~6 minutes.

import os, re, subprocess

import pytest

from conftest import run_psst, run_job, require_tools

pytestmark = pytest.mark.integration

BRAKER3_SIF = os.environ.get("PSST_TEST_BRAKER3_SIF")
EXAMPLE_DIR = "/opt/BRAKER/example"

def genes(gtfFile):
    '''
    Returns:
        genes -- a set of (contig, strand, start, end) tuples of each gene's CDS span.
    '''
    spans = {}
    for line in open(gtfFile):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        gene = re.search(r'gene_id "([^"]+)"', fields[8]).group(1)
        key = (fields[0], fields[6], gene)
        start, end = int(fields[3]), int(fields[4])
        spans[key] = (min(start, spans[key][0]), max(end, spans[key][1])) if key in spans else (start, end)
    return { (c, s, start, end) for (c, s, _), (start, end) in spans.items() }

@pytest.fixture(scope="module")
def braker3(tmp_path_factory):
    if BRAKER3_SIF == None or not os.path.isfile(BRAKER3_SIF):
        pytest.skip("requires a BRAKER3 Singularity image at $PSST_TEST_BRAKER3_SIF")
    require_tools("singularity")
    base = tmp_path_factory.mktemp("braker3")
    (base / "reads").mkdir()

    # Copy the example data out of the image, as a softmasked genome and paired read files
    subprocess.run(["singularity", "exec", BRAKER3_SIF, "bash", "-c",
                    f"cp {EXAMPLE_DIR}/genome.fa genome.fasta.masked && cp {EXAMPLE_DIR}/proteins.fa {EXAMPLE_DIR}/results/test3/braker.gtf . && " +
                    f"samtools view -b -f 1 -F 0x900 {EXAMPLE_DIR}/RNAseq.bam | samtools collate -O -u - | " +
                    "samtools fastq -n -1 reads/DRR095543_1.fq.gz -2 reads/DRR095543_2.fq.gz -0 /dev/null -s /dev/null -"],
                   cwd=base, check=True, capture_output=True)

    run_psst(["assembly", "braker3", "-g", str(base / "genome.fasta.masked"), "--proteins", str(base / "proteins.fa"),
              "-d", str(base / "reads"), "--sif", BRAKER3_SIF, "--lineage", "eukaryota_odb10", "--cpus", "8"], base)
    # Speed up the run with the options of BRAKER3's own test (which must not be used on real data)
    script = base / "braker3.sh"
    text = script.read_text()
    assert "--gff3 --nocleanup )" in text
    script.write_text(text.replace("--gff3 --nocleanup )", "--gff3 --nocleanup --skipOptimize --gm_max_intergenic 10000 )"))
    run_job(str(script), str(base))
    return base

def test_outputs(braker3):
    outDir = braker3 / "braker3"
    for ext in ["gff3", "gtf", "aa", "codingseq"]:
        assert os.path.getsize(outDir / f"genome.braker.{ext}") > 0
    assert sorted(os.listdir(outDir / "rnaseq_reads")) == ["DRR095543_1.fq.gz", "DRR095543_2.fq.gz"]
    assert "Found paired RNA-Seq reads" in open(outDir / "braker" / "braker.log").read() # used the local reads

def test_genes_match_reference_results(braker3):
    predicted = genes(braker3 / "braker3" / "genome.braker.gtf")
    reference = genes(braker3 / "braker.gtf")
    def overlaps(gene, others):
        c, s, start, end = gene
        return any(c == c2 and s == s2 and start <= end2 and start2 <= end for c2, s2, start2, end2 in others)
    # Locally, 256 of 265 reference genes were overlapped, and 248 had identical CDS spans
    assert sum(overlaps(g, predicted) for g in reference) > 0.9 * len(reference)
    assert len(predicted & reference) > 0.85 * len(reference)
    assert len(predicted) < 1.1 * len(reference)
