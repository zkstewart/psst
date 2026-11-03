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

from .pbs import finalise_job

# Genome FASTA extensions removed to name a genome after its file
GENOME_EXTENSION_REGEX = re.compile(r"\.(fasta|fa|fna|fas)(\.gz)?$", re.IGNORECASE)

# Text between a sample name and its haplotype number (e.g., the '.hap' of 'Busk.hap1')
HAPLOTYPE_LABEL_REGEX = re.compile(r"[._-]*(haplotype|hap|h)?[._-]*$", re.IGNORECASE)

def genome_name(fileName):
    '''
    Parameters:
        fileName -- a string of a genome FASTA file location.
    Returns:
        name -- a string of the file name without its directory or FASTA extension.
    '''
    return GENOME_EXTENSION_REGEX.sub("", os.path.basename(fileName))

def haplotype_position(name1, name2):
    '''
    Parameters:
        name1, name2 -- strings of genome names (see genome_name()).
    Returns:
        position -- an int index, within re.split(r"(\\d+)", name), of the only number at
                    which the names differ, or None if they differ in any other way.
    '''
    tokens1, tokens2 = re.split(r"(\d+)", name1), re.split(r"(\d+)", name2)
    if len(tokens1) != len(tokens2):
        return None
    differences = [ i for i, (a, b) in enumerate(zip(tokens1, tokens2)) if a != b ]
    if len(differences) != 1 or differences[0] % 2 == 0: # numbers are at odd indices
        return None
    position = differences[0]
    return position if int(tokens1[position]) != int(tokens2[position]) else None

def group_haplotypes(genomeFiles, haploid=False):
    '''
    Groups genome FASTA files into samples. Consecutive files whose names are identical
    except for one number are haplotypes of one sample, numbered by it (e.g.,
    'Busk.hap1.fasta' and 'Busk.hap2.fasta' are haplotypes 1 and 2 of sample 'Busk');
    other files are haploid samples (haplotype 0) named after their file.

    Parameters:
        genomeFiles -- a list of strings of genome FASTA file locations.
        haploid -- OPTIONAL; a boolean indicating whether to treat every file as a
                   haploid sample, without grouping; default == False.
    Returns:
        samples -- a list of (sample name, [(haplotype number, genome file), ...]) tuples
                   in input order.
    '''
    names = [ genome_name(f) for f in genomeFiles ]
    samples = []
    i = 0
    while i < len(names):
        # Extend the group while files differ from its first file in the same number
        j = i + 1
        position = haplotype_position(names[i], names[j]) if j < len(names) and not haploid else None
        if position != None:
            while j < len(names) and haplotype_position(names[i], names[j]) == position:
                j += 1
        if position == None:
            samples.append((names[i], [(0, genomeFiles[i])]))
        else:
            sample = HAPLOTYPE_LABEL_REGEX.sub("", "".join(re.split(r"(\d+)", names[i])[:position]))
            haplotypes = [ (int(re.split(r"(\d+)", names[k])[position]), genomeFiles[k]) for k in range(i, j) ]
            numbers = [ h for h, _ in haplotypes ]
            files = ", ".join([ os.path.basename(f) for _, f in haplotypes ])
            if sample == "":
                raise ValueError(f"Haplotype files ({files}) have no sample name before their haplotype " +
                                 "number; rename them (e.g., 'sample.hap1.fasta') or use --seqfile")
            if 0 in numbers or len(set(numbers)) != len(numbers):
                raise ValueError(f"Haplotype files ({files}) must have distinct haplotype numbers >= 1")
            samples.append((sample, haplotypes))
        i = j
    return samples

def name_genomes(samples, refName):
    '''
    Names grouped genomes for a cactus seqfile, taking the first genome as the reference.
    Any other haplotypes of the reference's sample form a sample of their own: haploid
    (haplotype 0) if one remains, or else haplotypes 1, 2, ... in order. Other samples keep
    their haplotype numbers.

    Parameters:
        samples -- a list of samples as produced by group_haplotypes().
        refName -- a string of the name to give the reference genome.
    Returns:
        entries -- a list of (sample name, haplotype number, genome file) tuples in seqfile
                   order, where the reference's haplotype is None.
    '''
    (firstSample, firstHaplotypes), otherSamples = samples[0], samples[1:]
    entries = [(refName, None, firstHaplotypes[0][1])]
    remaining = [ genomeFile for _, genomeFile in firstHaplotypes[1:] ]
    if len(remaining) == 1:
        entries.append((firstSample, 0, remaining[0]))
    else:
        entries += [ (firstSample, haplotype, genomeFile) for haplotype, genomeFile in enumerate(remaining, start=1) ]
    entries += [ (sample, haplotype, genomeFile) for sample, haplotypes in otherSamples
                 for haplotype, genomeFile in haplotypes ]
    return entries

def seqfile_name(sample, haplotype):
    '''
    Names a genome in a cactus seqfile as 'SAMPLE.HAPLOTYPE' (haplotype 0 being a haploid
    sample), or 'SAMPLE' for the reference, which cactus requires to be named as given to
    --reference. Cactus also documents PanSN-style 'SAMPLE#HAPLOTYPE' names, but v3.3.0
    fails to match a reference named 'SAMPLE#0', and fails to align 'SAMPLE#1' contigs, so
    these are not used.

    Parameters:
        sample -- a string of the sample name.
        haplotype -- an int of the haplotype number, or None for the reference.
    Returns:
        name -- a string of the seqfile name.
    '''
    return sample if haplotype == None else f"{sample}.{haplotype}"

def parse_seqfile(seqFile):
    '''
    Parses a cactus seqfile of '<name> <FASTA>' lines, where names are 'SAMPLE.HAPLOTYPE',
    or PanSN-style 'SAMPLE#HAPLOTYPE', or 'SAMPLE' for a haploid sample (haplotype 0).

    Parameters:
        seqFile -- a string of the seqfile location; relative FASTA locations within it
                   are resolved against the current directory, as cactus would.
    Returns:
        entries -- a list of (sample name, haplotype number, absolute FASTA location)
                   tuples in file order.
    '''
    entries = []
    with open(seqFile) as fileIn:
        for lineNumber, line in enumerate(fileIn, start=1):
            fields = line.split()
            if fields == []:
                continue
            if len(fields) != 2:
                raise ValueError(f"Line {lineNumber} of --seqfile {seqFile} must have 2 columns " +
                                 f"(a name and a FASTA file), not {len(fields)}")
            name, fastaFile = fields
            match = re.match(r"^(.+?)(?:#(\d+)|\.(\d+))?$", name) # as cactus reads names
            sample, haplotype = match.group(1), match.group(2) or match.group(3) or "0"
            if "#" in sample:
                raise ValueError(f"Name '{name}' on line {lineNumber} of --seqfile {seqFile} must be " +
                                 "'SAMPLE' or 'SAMPLE#HAPLOTYPE' with a numeric haplotype")
            entries.append((sample, int(haplotype), os.path.abspath(fastaFile)))
    return entries

def pangenome_create(args):
    '''
    Generates (and optionally submits) a PBS job script which writes a cactus seqfile of
    the genomes and builds a pangenome graph from them with cactus-pangenome, producing
    a VCF against the reference and giraffe indexes.

    Parameters:
        args -- an argparse.Namespace containing validated "pangenome create" arguments.
    '''
    finalise_job(args, "pangenome/run_cactus_pangenome.sh", {
        "condaEnv": args.condaEnv,
        "binariesMode": args.binariesMode,
        "seqNames": args.seqNames,
        "genomeFiles": args.genomeFiles,
        "refName": args.refName,
        "outputDir": args.outputDir,
        "prefix": args.prefix,
        "giraffe": args.giraffe,
        "vcfwave": "false" if args.noVcfwave else "true",
        "cpus": args.cpus,
        "maxMemory": args.maxMemory
    })
