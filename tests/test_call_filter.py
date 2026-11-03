# Tests of "call filter": the bcftools commands its job builds from each option (using a
# stub bcftools), and the sites real bcftools keeps from a VCF of known sites.

import os, subprocess

import pytest

from conftest import run_psst, run_job, require_tools

# Sites designed to fail particular filters, with genotypes for samples a-d plus a haploid
# sample h (as pangenome haplotypes may be)
VCF_HEADER = """##fileformat=VCFv4.2
##contig=<ID=c,length=1000>
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\th
"""
SITES = [ # position, REF, ALT, QUAL, FILTER, genotypes, what the site tests
    (1, "A", "C", "50", ".", "0/0 0/1 0/0 0/0 0", "passes defaults; minor allele frequency 1/9"),
    (2, "A", "C", "50", "PASS", "0/1 0/1 1/1 0/0 1", "passes everything"),
    (3, "A", "C", ".", ".", "0/1 0/1 0/0 0/0 1", "no QUAL"),
    (4, "A", "C", "20", ".", "0/1 0/1 0/0 0/0 1", "low QUAL"),
    (5, "A", "C", "50", "LowQual", "0/1 0/1 0/0 0/0 1", "failed FILTER"),
    (6, "A", "C", "50", ".", "./. ./. 0/1 0/0 .", "3 of 5 samples missing"),
    (7, "A", "C", "50", ".", "1/1 1/1 1/1 1/1 1", "no minor allele"),
    (8, "A", "C,G", "50", ".", "0/1 0/2 0/1 1/2 0", "multiallelic record"),
    (9, "A", "C", "50", ".", "0/1 0/1 0/0 0/0 1", "split multiallelic site (with 9 below)"),
    (9, "A", "G", "50", ".", "0/0 0/0 0/1 0/1 0", "split multiallelic site (with 9 above)"),
    (10, "AT", "A", "50", ".", "0/1 0/1 0/0 0/0 1", "indel"),
    (11, "ATC", "GA", "50", ".", "0/1 0/1 0/0 0/0 1", "complex variant"),
]

@pytest.fixture
def vcf(tmp_path):
    lines = [ f"c\t{p}\t.\t{ref}\t{alt}\t{qual}\t{flt}\t.\tGT\t" + "\t".join(gts.split()) for p, ref, alt, qual, flt, gts, _ in SITES ]
    (tmp_path / "sites.vcf").write_text(VCF_HEADER + "\n".join(lines) + "\n")
    return str(tmp_path / "sites.vcf")

def run_filter(vcf, tmp_path, args, env=None):
    workDir = tmp_path / "work"
    workDir.mkdir()
    run_psst(["call", "filter", "-i", vcf, *args], workDir)
    run_job(str(workDir / "call_filter.sh"), str(workDir), env=env)
    return vcf[:-len(".vcf")] + ".filtered.vcf.gz"

@pytest.mark.parametrize("args,expected,unexpected", [
    ([], ["view -i F_MISSING<=0.5 && QUAL>=30 && MAC>=1 -f PASS,. -Ou"], ["norm", "-V ", "MAF"]),
    (["--max-missing", "1", "--minq", "0", "--mac", "0", "--keep-filtered"], ["view -Ou "], ["-i ", "-f PASS"]),
    (["--maf", "0.05", "--max-missing", "0.2", "--minq", "20.5", "--mac", "3"],
     ["-i F_MISSING<=0.2 && QUAL>=20.5 && MAC>=3 && MAF>=0.05 -f PASS,."], []),
    (["--biallelic"], ["norm -m +any -Ou", "view -M 2 -Ou"], ["-V "]),
    (["--no-indels"], ["view -V indels,other -Ou"], ["norm"]),
])
def test_filter_commands(args, expected, unexpected, vcf, tmp_path, stub_bcftools):
    env, stubDir = stub_bcftools
    output = run_filter(vcf, tmp_path, args, env=env)
    log = open(os.path.join(stubDir, "calls.log")).read()
    for text in expected:
        assert text in log
    for text in unexpected:
        assert text not in log
    for suffix in [".vcf.gz", ".vcf.gz.csi", ".vcf.gz.tbi", ".stats.txt"]:
        assert os.path.exists(output.replace(".vcf.gz", suffix))
    assert not os.path.exists(output.replace(".vcf.gz", ".incomplete.vcf.gz"))

def test_filter_refuses_existing_output(vcf, tmp_path, stub_bcftools):
    env, stubDir = stub_bcftools
    workDir = tmp_path / "work"
    workDir.mkdir()
    run_psst(["call", "filter", "-i", vcf], workDir)
    open(vcf[:-len(".vcf")] + ".filtered.vcf.gz", "w").close()
    result = run_job(str(workDir / "call_filter.sh"), str(workDir), env=env, check=False)
    assert result.returncode != 0
    assert "already exists" in result.stdout

@pytest.mark.integration
@pytest.mark.parametrize("args,kept", [
    ([], [1, 2, 8, 9, 9, 10, 11]),
    (["--maf", "0.25"], [2, 9, 10, 11]), # the rarest allele of site 8 has frequency 2/9
    (["--biallelic"], [1, 2, 10, 11]),
    (["--no-indels"], [1, 2, 8, 9, 9]),
    (["--biallelic", "--no-indels", "--maf", "0.15"], [2]),
    (["--keep-filtered"], [1, 2, 5, 8, 9, 9, 10, 11]),
    (["--max-missing", "1", "--minq", "0", "--mac", "0", "--keep-filtered"], [1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 10, 11]),
    (["--max-missing", "0.6"], [1, 2, 6, 8, 9, 9, 10, 11]),
    (["--minq", "10"], [1, 2, 4, 8, 9, 9, 10, 11]),
], ids=lambda x: " ".join(x) if isinstance(x, list) and x and isinstance(x[0], str) else None)
def test_filter_sites(args, kept, vcf, tmp_path):
    require_tools("bcftools")
    output = run_filter(vcf, tmp_path, args)
    positions = subprocess.run(["bcftools", "query", "-f", "%POS\\n", output], check=True,
                               capture_output=True, text=True).stdout.split()
    assert [ int(p) for p in positions ] == kept
    assert "number of records:" in open(output.replace(".vcf.gz", ".stats.txt")).read()
