# Tests of psst's helper functions.

import gzip, os, sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.assembly import braker_prefix, check_rnaseq_ids, read_sra_accessions
from modules.calling import get_bam_samples, get_vcf_samples
from modules.mapping import get_array_size
from modules.pangenome import group_haplotypes, name_genomes, parse_seqfile, seqfile_name
from modules.pbs import fill_template, format_pbs_header, q
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

def test_q_quotes_lists_item_by_item():
    assert q("plain") == "plain"
    assert q(["a.fa", "my file.fa", "it's.fa"]) == """a.fa 'my file.fa' 'it'"'"'s.fa'"""
    assert q(2) == "2"

@pytest.mark.parametrize("files,expected", [
    (["Busk.hap1.fasta", "Busk.hap2.fasta"], [("Busk", [1, 2])]),
    (["x/A_1.fa.gz", "x/A_2.fa.gz", "x/A_3.fa.gz", "x/A_4.fa.gz"], [("A", [1, 2, 3, 4])]), # any ploidy
    (["14Q021.hap1.fasta", "14Q021.hap2.fasta", "18Q053.fasta"], [("14Q021", [1, 2]), ("18Q053", [0])]),
    (["s.asm.bp.hap1.p_ctg.fa", "s.asm.bp.hap2.p_ctg.fa"], [("s.asm.bp", [1, 2])]), # text may follow the number
    (["P1.h1.fa", "P1.h2.fa", "P2.h1.fa", "P2.h2.fa"], [("P1", [1, 2]), ("P2", [1, 2])]), # only one number may differ
    (["A.hap2.fa", "A.hap1.fa"], [("A", [2, 1])]), # numbered by name, not order
    (["A.hap1.fa", "B.fa", "A.hap2.fa"], [("A.hap1", [0]), ("B", [0]), ("A.hap2", [0])]), # must be consecutive
    (["A.hap1.fa", "A.hap01.fa"], [("A.hap1", [0]), ("A.hap01", [0])]), # the same number isn't a haplotype
    (["iso1.fa", "iso2.fa"], [("iso", [1, 2])]), # unavoidably grouped; see --haploid
])
def test_group_haplotypes(files, expected):
    samples = group_haplotypes(files)
    assert [ (sample, [ h for h, _ in haplotypes ]) for sample, haplotypes in samples ] == expected
    assert [ f for _, haplotypes in samples for _, f in haplotypes ] == files

def test_group_haplotypes_haploid():
    assert group_haplotypes(["iso1.fa", "iso2.fa"], haploid=True) == [("iso1", [(0, "iso1.fa")]), ("iso2", [(0, "iso2.fa")])]

@pytest.mark.parametrize("files,message", [
    (["1.fa", "2.fa"], "have no sample name"),
    (["hap1.fa", "hap2.fa"], "have no sample name"),
    (["A.hap0.fa", "A.hap1.fa"], "distinct haplotype numbers >= 1"),
])
def test_group_haplotypes_errors(files, message):
    with pytest.raises(ValueError, match=message):
        group_haplotypes(files)

def test_seqfile_names():
    assert seqfile_name("reference", None) == "reference"
    assert seqfile_name("18Q049", 0) == "18Q049.0"
    assert seqfile_name("13Q019", 2) == "13Q019.2"

@pytest.mark.parametrize("files,expected", [
    (["s.hap1.fa", "s.hap2.fa", "o.hap1.fa", "o.hap2.fa"], ["reference", "s.0", "o.1", "o.2"]), # diploid reference
    (["s.hap1.fa", "s.hap2.fa", "s.hap3.fa", "o.hap1.fa", "o.hap2.fa", "o.hap3.fa"],
     ["reference", "s.1", "s.2", "o.1", "o.2", "o.3"]), # triploid reference; the rest are renumbered
    (["ref.fa", "o.hap1.fa", "o.hap2.fa", "iso.fa"], ["reference", "o.1", "o.2", "iso.0"]), # haploid reference
    (["s.hap2.fa", "s.hap1.fa", "o.fa"], ["reference", "s.0", "o.0"]), # the first file listed is the reference
])
def test_name_genomes(files, expected):
    entries = name_genomes(group_haplotypes(files), "reference")
    assert [ seqfile_name(sample, haplotype) for sample, haplotype, _ in entries ] == expected
    assert [ f for _, _, f in entries ] == files

def test_parse_seqfile(tmp_path):
    (tmp_path / "seqfile.txt").write_text("reference       /g/18Q049.hap2.fasta\n18Q049.0        /g/18Q049.hap1.fasta\n" +
                                          "13Q019.1\t/g/13Q019.hap1.fasta\n\n13Q019#2\trel.fa\nv1.5.2 /g/v.fa\n")
    assert parse_seqfile(str(tmp_path / "seqfile.txt")) == [
        ("reference", 0, "/g/18Q049.hap2.fasta"), ("18Q049", 0, "/g/18Q049.hap1.fasta"),
        ("13Q019", 1, "/g/13Q019.hap1.fasta"), ("13Q019", 2, os.path.abspath("rel.fa")), ("v1.5", 2, "/g/v.fa") ]
    (tmp_path / "bad.txt").write_text("a b c\n")
    with pytest.raises(ValueError, match="must have 2 columns"):
        parse_seqfile(str(tmp_path / "bad.txt"))

@pytest.mark.parametrize("fileName, expected", [
    ("/x/Fingerlime_genome.fna.masked", "Fingerlime_genome"),
    ("genome.softmasked.fasta", "genome"),
    ("a.b.fa.gz", "a.b"),
    ("plain", "plain"),
])
def test_braker_prefix(fileName, expected):
    assert braker_prefix(fileName) == expected

def test_read_sra_accessions(tmp_path):
    (tmp_path / "ok.txt").write_text("SRR1\n\n  ERR22 \r\nDRR3\n")
    assert read_sra_accessions(str(tmp_path / "ok.txt")) == ["SRR1", "ERR22", "DRR3"]
    (tmp_path / "empty.txt").write_text("\n")
    with pytest.raises(ValueError, match="lists no SRA accessions"):
        read_sra_accessions(str(tmp_path / "empty.txt"))
    (tmp_path / "bad.txt").write_text("SRR1\nSRX5\n")
    with pytest.raises(ValueError, match="Line 2 .* \\('SRX5'\\) is not an SRA run accession"):
        read_sra_accessions(str(tmp_path / "bad.txt"))

@pytest.mark.parametrize("ids, message", [
    (["SRR1", "SRR11", "leaf"], None),
    (["leaf", "SRR1", "leaf"], "'leaf' is given more than once"),
    (["A1", "BA1"], "'A1' matches the end of 'BA1'"),
    (["SRR11", "SRR1"], None), # neither ends the other
    (["x.1", "xy1"], "'x.1' matches the end of 'xy1'"), # IDs are regexes to BRAKER3
])
def test_check_rnaseq_ids(ids, message):
    if message == None:
        check_rnaseq_ids(ids)
    else:
        with pytest.raises(ValueError, match=message):
            check_rnaseq_ids(ids)
