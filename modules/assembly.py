# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, re

from .pangenome import genome_name
from .pbs import finalise_job

# SRA run accessions (e.g., 'SRR1234567') which BRAKER3 can download
SRA_ACCESSION_REGEX = re.compile(r"^[SED]RR\d+$")

# Characters allowed in the IDs of BRAKER3's RNAseq sets (i.e., local sample names), which
# it uses within regexes to find their files
RNASEQ_ID_REGEX = re.compile(r"^[A-Za-z0-9._-]+$")

# Read file suffixes whose files can be linked under names BRAKER3 recognises as FASTQ
BRAKER_READS_SUFFIX_REGEX = re.compile(r"\.(fq|fastq)(\.gz)?$")

# Text left after a softmasked genome's FASTA extension, or preceding it
MASKED_REGEX = re.compile(r"\.(soft)?masked$", re.IGNORECASE)

def assembly_hifiasm(args):
    '''
    Generates (and optionally submits) a PBS job script which runs a hifiasm genome
    assembly, followed by FASTA conversion, statistics and completeness scoring. The
    Hi-C template is used if --hic was given (with or without files).

    Parameters:
        args -- an argparse.Namespace containing validated "assembly hifiasm" arguments.
    '''
    values = {
        "genScriptDir": args.genScriptDir,
        "compleasmDir": args.compleasmDir,
        "cpus": args.cpus,
        "lineage": args.lineage
    }
    if args.hicFiles == None:
        finalise_job(args, "assembly/run_hifiasm_assembly.sh", values)
    else:
        hic1, hic2 = args.hicFiles if args.hicFiles != [] else (None, None)
        finalise_job(args, "assembly/run_hifiasm_assembly.hic.sh", {
            **values,
            "ploidy": args.ploidy,
            "hic1": hic1,
            "hic2": hic2
        })

def phase_values(args):
    '''
    Parameters:
        args -- an argparse.Namespace containing validated arguments of an "assembly" mode
                which phases contigs with annotarium (see validation.validate_phasing()).
    Returns:
        values -- a dict of the template values shared by these modes.
    '''
    return {
        "annotariumDir": args.annotariumDir,
        "minimap2Dir": args.minimap2Dir,
        "contigFiles": args.contigFiles,
        "refGenome": args.refGenome,
        "phaseDir": args.phaseDir if args.phaseDir != None else "", # empty phases the contigs
        "ploidy": args.ploidy,
        "preset": args.preset,
        "outputDir": args.outputDir,
        "cpus": args.cpus
    }

def assembly_revcomp(args):
    '''
    Generates (and optionally submits) a PBS job script which phases contigs into haplotype
    partitions against a reference genome with annotarium (unless given an earlier phasing),
    then writes each partition's contigs paired 1:1 with the reference contigs, renamed and
    reverse complemented to match them.

    Parameters:
        args -- an argparse.Namespace containing validated "assembly revcomp" arguments.
    '''
    finalise_job(args, "assembly/run_revcomp.sh", phase_values(args))

def assembly_scaffolding(args):
    '''
    Generates (and optionally submits) a PBS job script which phases contigs into haplotype
    partitions against a reference genome with annotarium (unless given an earlier phasing),
    then scaffolds each partition against the reference with RagTag and dotplots it.

    Parameters:
        args -- an argparse.Namespace containing validated "assembly scaffolding" arguments.
    '''
    finalise_job(args, "assembly/run_scaffolding.sh", {
        **phase_values(args),
        "paf2dotplotDir": args.paf2dotplotDir
    })

def braker_prefix(genomeFile):
    '''
    Parameters:
        genomeFile -- a string of a (softmasked) genome FASTA file location.
    Returns:
        prefix -- a string of the file name without its directory, FASTA extension, or
                  '.masked'/'.softmasked' before or after it.
    '''
    name = MASKED_REGEX.sub("", os.path.basename(genomeFile))
    return MASKED_REGEX.sub("", genome_name(name))

def read_sra_accessions(sraFile):
    '''
    Parameters:
        sraFile -- a string of a text file listing SRA accessions, one per line; blank lines
                   are ignored.
    Returns:
        accessions -- a list of the accession strings, in file order.
    '''
    accessions = []
    with open(sraFile, "r") as fileIn:
        for lineNum, line in enumerate(fileIn, start=1):
            line = line.strip()
            if line == "":
                continue
            if SRA_ACCESSION_REGEX.match(line) == None:
                raise ValueError(f"Line {lineNum} of --sra file {sraFile} ('{line}') is not an SRA run " +
                                 "accession (e.g., SRR1234567)")
            accessions.append(line)
    if accessions == []:
        raise ValueError(f"--sra file {sraFile} lists no SRA accessions")
    return accessions

def check_rnaseq_ids(rnaseqIDs):
    '''
    Checks that BRAKER3 can tell RNAseq sets apart. It finds the files of each set ID by
    matching the ID, as a regex, to the end of file names (before '_1.fq.gz' etc.), so an ID
    matching the end of another ID would be given that set's files.

    Parameters:
        rnaseqIDs -- a list of strings of RNAseq set IDs (local sample names and SRA accessions).
    '''
    for i, id1 in enumerate(rnaseqIDs):
        for j, id2 in enumerate(rnaseqIDs):
            if i == j:
                continue
            if id1 == id2:
                raise ValueError(f"RNAseq set '{id1}' is given more than once (as a local sample name or SRA accession)")
            if re.search(f"{id1}$", id2) != None:
                raise ValueError(f"RNAseq set '{id1}' matches the end of '{id2}', so BRAKER3 would confuse " +
                                 "their files; rename the local read files")

def assembly_braker3(args):
    '''
    Generates (and optionally submits) a PBS job script which predicts genes in a softmasked
    genome with BRAKER3 from its Singularity image, using protein evidence and RNAseq from
    local paired read files and/or SRA accessions (combined from one or more files).

    Parameters:
        args -- an argparse.Namespace containing validated "assembly braker3" arguments.
    '''
    finalise_job(args, "assembly/run_braker3.sh", {
        "brakerSif": args.brakerSif,
        "genome": args.genome,
        "proteins": args.proteins,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "sraFiles": args.sraFiles,
        "prefix": args.prefix,
        "lineage": args.lineage,
        "outputDir": args.outputDir,
        "cpus": args.cpus
    })
