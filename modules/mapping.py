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

from .pbs import finalise_job

# Files written alongside a FASTA file by 'bwa index'
BWA_INDEX_SUFFIXES = [".amb", ".ann", ".bwt", ".pac", ".sa"]

# Files within a 'map index giraffe' directory; the templates use these names too
GIRAFFE_MAP_FILES = ["pangenome.gbz", "pangenome.dist", "pangenome.min", "pangenome.zipcodes"]
GIRAFFE_CALL_FILES = ["reference.fa", "reference.fa.fai", "sites.vcf.gz", "targets.tsv.gz"]

def get_array_size(numSamples):
    '''
    Parameters:
        numSamples -- an int of the number of samples to process with one subjob each,
                      or None if unknown.
    Returns:
        arraySize -- an int of the PBS array size, None if only one sample needs no
                     array, or an '@@numSamples@@' placeholder if numSamples is unknown.
    '''
    if numSamples == None:
        return "@@numSamples@@" # left for the user to fill in
    return numSamples if numSamples > 1 else None

def map_short(args):
    '''
    Generates (and optionally submits) a PBS array job script which maps short DNA reads
    (e.g., Illumina WGS or Hi-C) with bwa mem, piping its output through samtools to sort
    and mark duplicates. Each subjob maps one sample; the array size is a placeholder if
    no reads directory was given, and no array is requested if there is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "map short" arguments.
    '''
    finalise_job(args, "map/run_bwa_short.sh", {
        "genome": args.genome,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "hic": "true" if args.hic else "false",
        "outputDir": args.outputDir,
        "library": args.library,
        "cpus": args.cpus,
        "sortMem": args.sortMem
    }, arraySize=get_array_size(args.numSamples))

def map_index_bwa(args):
    '''
    Generates (and optionally submits) a PBS job script which indexes a reference
    genome FASTA file with bwa index, writing index files alongside it.

    Parameters:
        args -- an argparse.Namespace containing validated "map index bwa" arguments.
    '''
    finalise_job(args, "map/run_bwa_index.sh", {
        "genome": args.genome
    })

def map_index_salmon(args):
    '''
    Generates (and optionally submits) a PBS job script which indexes a transcripts
    FASTA file with salmon index, optionally with the genome as decoys.

    Parameters:
        args -- an argparse.Namespace containing validated "map index salmon" arguments.
    '''
    finalise_job(args, "map/run_salmon_index.sh", {
        "salmonDir": args.salmonDir,
        "transcripts": args.transcripts,
        "genome": args.genome if args.genome != None else "", # empty means no decoys
        "indexDir": args.indexDir,
        "kmer": args.kmer,
        "ramLimit": args.ramLimit,
        "cpus": args.cpus
    })

def map_rna_salmon(args):
    '''
    Generates (and optionally submits) a PBS array job script which quantifies RNAseq
    reads against a salmon index. Each subjob quantifies one sample; the array size is
    a placeholder if no reads directory was given, and no array is requested if there
    is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "map rna salmon" arguments.
    '''
    finalise_job(args, "map/run_salmon_quant.sh", {
        "salmonDir": args.salmonDir,
        "indexDir": args.indexDir,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir,
        "libType": args.libType,
        "geneMap": args.geneMap if args.geneMap != None else "", # empty means no gene map
        "gcBias": "true" if args.bias in ["gc", "both"] else "false",
        "seqBias": "true" if args.bias in ["seq", "both"] else "false",
        "bootstraps": args.bootstraps,
        "cpus": args.cpus
    }, arraySize=get_array_size(args.numSamples))

def map_index_star(args):
    '''
    Generates (and optionally submits) a PBS job script which indexes a reference
    genome FASTA file and its GTF annotation with STAR genomeGenerate.

    Parameters:
        args -- an argparse.Namespace containing validated "map index star" arguments.
    '''
    finalise_job(args, "map/run_star_index.sh", {
        "starModule": args.starModule,
        "genome": args.genome,
        "gtf": args.gtf,
        "indexDir": args.indexDir,
        "overhang": args.overhang,
        "cpus": args.cpus,
        "ramLimit": args.ramLimit
    })

def map_rna_star(args):
    '''
    Generates (and optionally submits) a PBS array job script which maps RNAseq reads
    to a genome with STAR, producing a coordinate-sorted BAM, a transcriptome BAM and
    per-gene read counts. Each subjob maps one sample; the array size is a placeholder
    if no reads directory was given, and no array is requested if there is only one
    sample.

    Parameters:
        args -- an argparse.Namespace containing validated "map rna star" arguments.
    '''
    finalise_job(args, "map/run_star_align.sh", {
        "starModule": args.starModule,
        "indexDir": args.indexDir,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir,
        "twopass": "true" if args.twopass else "false",
        "cpus": args.cpus
    }, arraySize=get_array_size(args.numSamples))

def map_nanopore(args):
    '''
    Generates (and optionally submits) a PBS array job script which maps nanopore reads
    with minimap2, optionally filtering them with chopper first, and pipes alignments
    into a sorted BAM. Each subjob maps one sample; the array size is a placeholder if
    no reads directory was given, and no array is requested if there is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "map nanopore" arguments.
    '''
    finalise_job(args, "map/run_minimap2_nanopore.sh", {
        "minimap2Dir": args.minimap2Dir,
        "chopperDir": args.chopperDir,
        "genome": args.genome,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "outputDir": args.outputDir,
        "preset": args.preset,
        "library": args.library,
        "minLength": args.minLength,
        "maxLength": args.maxLength,
        "minQual": args.minQual,
        "cpus": args.cpus,
        "sortMem": args.sortMem
    }, arraySize=get_array_size(args.numSamples))

def map_index_giraffe(args):
    '''
    Generates (and optionally submits) a PBS job script which prepares a pangenome graph
    for read mapping with vg giraffe and variant calling with "call giraffe", creating
    an index directory containing the graph, its giraffe indexes, the reference
    sample's paths as a linear reference, and the pangenome's variant sites.

    Parameters:
        args -- an argparse.Namespace containing validated "map index giraffe" arguments.
    '''
    finalise_job(args, "map/run_giraffe_index.sh", {
        "vgDir": args.vgDir,
        "gbz": args.gbz,
        "vcf": args.vcf,
        "refSample": args.refSample,
        "indexDir": args.indexDir,
        "cpus": args.cpus,
        "sortMem": args.sortMem
    })

def map_giraffe(args):
    '''
    Generates (and optionally submits) a PBS array job script which maps short reads to
    a pangenome graph with vg giraffe, projecting alignments onto the reference paths and
    piping them through samtools to sort and mark duplicates. Each subjob maps one sample;
    the array size is a placeholder if no reads directory was given, and no array is
    requested if there is only one sample.

    Parameters:
        args -- an argparse.Namespace containing validated "map giraffe" arguments.
    '''
    finalise_job(args, "map/run_giraffe_map.sh", {
        "vgDir": args.vgDir,
        "indexDir": args.indexDir,
        "readsDir": args.readsDir,
        "readsSuffix": args.readsSuffix,
        "paired": "false" if args.singleEnd else "true",
        "outputDir": args.outputDir,
        "cpus": args.cpus,
        "sortMem": args.sortMem
    }, arraySize=get_array_size(args.numSamples))
