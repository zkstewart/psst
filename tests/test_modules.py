# Tests of psst's helper functions.

import gzip, os, sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.calling import get_bam_samples, get_vcf_samples
from modules.mapping import get_array_size
from modules.pbs import fill_template, format_pbs_header
from modules.reads import get_read_files
from modules.validation import parse_mem_gb

def touch(directory, *names):
    for name in names:
        open(os.path.join(directory, name), "w").close()

def test_get_read_files_paired(tmp_path):
    touch(tmp_path, "a_R1.fq.gz", "a_R2.fq.gz", "b.1.fq.gz", "b.2.fq.gz", "c-1.fq.gz", "c-2.fq.gz")
    samples = get_read_files(str(tmp_path), ".fq.gz", False)
    assert [ s[0] for s in samples ] == ["a", "b", "c"] # separators (and 'R') are stripped
    assert samples[0][1:] == (str(tmp_path / "a_R1.fq.gz"), str(tmp_path / "a_R2.fq.gz"))

def test_get_read_files_single_end(tmp_path):
    touch(tmp_path, "a.fq.gz", "b.fq.gz", "notes.txt")
    assert get_read_files(str(tmp_path), ".fq.gz", True) == [
        ("a", str(tmp_path / "a.fq.gz"), None), ("b", str(tmp_path / "b.fq.gz"), None) ]

@pytest.mark.parametrize("files,message", [
    (["a_1.fq.gz"], "has no matching reverse read file"),
    (["a_1.fq.gz", "a_2.fq.gz", "b_2.fq.gz"], "have no matching forward read file"),
    (["a_1.fq.gz", "a_2.fq.gz", "a.fq.gz"], "is not preceded by a 1 or 2"),
    (["a_1.fq.gz", "a_2.fq.gz", "a.1.fq.gz", "a.2.fq.gz"], "is shared by multiple read files"),
    (["notes.txt"], "No files ending with"),
])
def test_get_read_files_errors(files, message, tmp_path):
    touch(tmp_path, *files)
    with pytest.raises(Exception, match=message):
        get_read_files(str(tmp_path), ".fq.gz", False)

def write_bam_header(fileName, text):
    with gzip.open(fileName, "wb") as fileOut: # BAM is BGZF, which gzip can read
        fileOut.write(b"BAM\x01" + len(text).to_bytes(4, "little") + text)

def test_get_bam_samples(tmp_path):
    bamFile = str(tmp_path / "a.bam")
    write_bam_header(bamFile, b"@HD\tVN:1.6\n@RG\tID:x\tSM:s2\n@RG\tID:y\tSM:s1\tLB:l\n@RG\tID:z\tSM:s1\n")
    assert get_bam_samples(bamFile) == ["s1", "s2"]
    write_bam_header(bamFile, b"@HD\tVN:1.6\n")
    assert get_bam_samples(bamFile) == []
    with open(bamFile, "w") as fileOut:
        fileOut.write("not a BAM")
    assert get_bam_samples(bamFile) == None

def test_get_vcf_samples(tmp_path):
    vcfFile = str(tmp_path / "a.vcf.gz")
    with gzip.open(vcfFile, "wt") as fileOut:
        fileOut.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\n")
    assert get_vcf_samples(vcfFile) == ["P1", "P2"]

def test_get_array_size():
    assert get_array_size(None) == "@@numSamples@@"
    assert get_array_size(1) == None
    assert get_array_size(5) == 5

def test_parse_mem_gb():
    assert parse_mem_gb("70G") == 70
    assert parse_mem_gb("70gb") == 70
    assert parse_mem_gb("2048M") == 2

def test_fill_template_quotes_and_keeps_placeholders():
    filled, missing = fill_template("call/run_bcftools_concat.sh", {
        "outputDir": "my dir", "prefix": None, "sites": "", "chunks": 3, "cpus": 2 })
    assert "OUTDIR='my dir'\n" in filled
    assert "PREFIX=@@prefix@@\n" in filled
    assert "SITES=''\n" in filled
    assert missing == ["prefix"]

def test_fill_template_rejects_mismatched_fields():
    with pytest.raises(KeyError, match="fields without values = \\['sites'\\], values without fields = \\['typo'\\]"):
        fill_template("call/run_bcftools_concat.sh", {
            "outputDir": "x", "prefix": "x", "chunks": 3, "cpus": 2, "typo": 1 })

def test_format_pbs_header():
    header = format_pbs_header("job", "01:00:00", "4G", 2, "ABCDEF1234", afterok=["1.pbs", "2.pbs"], arraySize=7)
    assert "#PBS -J 1-7" in header
    assert "#PBS -W depend=afterok:1.pbs:2.pbs" in header
    assert "#PBS -P ABCDEF1234" in header
