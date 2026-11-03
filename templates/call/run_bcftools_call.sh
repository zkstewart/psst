####

# Specify the reference genome FASTA, which must be indexed with 'samtools faidx'
GENOME=@@genome@@

# Specify the directory of BAM files to call variants from jointly, and the suffix
# identifying them; BAMs must be indexed, and their read groups (SM) name the samples
BAMDIR=@@bamDir@@
SUFFIX=@@bamSuffix@@

# Specify the output directory and the prefix for output files
OUTDIR=@@outputDir@@
PREFIX=@@prefix@@

# Specify the number of chunks the genome is divided into (one per subjob)
NCHUNKS=@@chunks@@

# Specify pileup parameters
PROFILE=@@profile@@ # bcftools mpileup -X platform profile; empty for none
MINMAPQ=@@minMapQ@@
MINBASEQ=@@minBaseQ@@ # empty to use the profile's (or bcftools') default
MAXDEPTH=@@maxDepth@@ # per BAM file
ORPHANS=@@countOrphans@@ # true or false; include pairs not flagged as properly paired (e.g., Hi-C)
INDELSCNS=@@indelsCns@@ # true or false; use the --indels-cns indel calling model

# Specify calling parameters
PLOIDY=@@ploidy@@ # 1 or 2

# Specify computational parameters
CPUS=@@cpus@@

####

set -eo pipefail
shopt -s nullglob

# STEP 1: Skip chunks which have already been called
CHUNK=${PBS_ARRAY_INDEX:-1} # PBS_ARRAY_INDEX is unset if there is only one chunk
CHUNKDIR=${OUTDIR}/${PREFIX}_chunks
CHUNKVCF=${CHUNKDIR}/chunk_${CHUNK}.vcf.gz
EMPTYFLAG=${CHUNKDIR}/chunk_${CHUNK}.empty
if [[ -f ${CHUNKVCF} || -f ${EMPTYFLAG} ]]; then
    echo "Chunk ${CHUNK} has already been called; skipping it"
    exit 0
fi
mkdir -p "${CHUNKDIR}"

# STEP 2: Find the contigs of this chunk; contigs are assigned to chunks of similar total
# length in genome order, so that chunk outputs concatenate in genome order
if [[ ! -f ${GENOME}.fai ]]; then
    echo "${GENOME}.fai does not exist; index the genome with 'samtools faidx' first"
    exit 1
fi
REGIONS=${CHUNKDIR}/chunk_${CHUNK}.regions.txt
awk -v nchunks=${NCHUNKS} -v chunk=${CHUNK} 'BEGIN{OFS="\t"}
    NR==FNR{total+=$2; next}
    {
        c = int((sum + $2 / 2) * nchunks / total) + 1; # the chunk holding this contig'"'"'s midpoint
        if (c > nchunks) c = nchunks;
        sum += $2;
        if (c == chunk) print $1, 1, $2;
    }' "${GENOME}.fai" "${GENOME}.fai" > "${REGIONS}"
if [[ ! -s ${REGIONS} ]]; then
    echo "Chunk ${CHUNK} has no contigs (the genome has fewer large contigs than chunks)"
    rm "${REGIONS}"
    touch "${EMPTYFLAG}"
    exit 0
fi

# STEP 3: Find the BAM files
BAMS=( "${BAMDIR}"/*"${SUFFIX}" )
if (( ${#BAMS[@]} == 0 )); then
    echo "No BAM files ending with ${SUFFIX} were found in ${BAMDIR}"
    exit 1
fi

# STEP 4: Pileup, call and normalise variants without intermediate files, giving the
# output its final name once complete; norm left-aligns indels and trims redundant bases
MPILEUPARGS=()
[[ -n ${PROFILE} ]] && MPILEUPARGS+=( -X ${PROFILE} ) # first, so later options override it
MPILEUPARGS+=( -q ${MINMAPQ} -d ${MAXDEPTH} )
[[ -n ${MINBASEQ} ]] && MPILEUPARGS+=( -Q ${MINBASEQ} )
[[ ${ORPHANS} == "true" ]] && MPILEUPARGS+=( -A )
[[ ${INDELSCNS} == "true" ]] && MPILEUPARGS+=( --indels-cns )
TMPVCF=${CHUNKDIR}/chunk_${CHUNK}.incomplete.vcf.gz
bcftools mpileup -f "${GENOME}" -R "${REGIONS}" "${MPILEUPARGS[@]}" \
    -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD -Ou "${BAMS[@]}" \
    | bcftools call -m -v --ploidy ${PLOIDY} -a GQ,GP -Ou \
    | bcftools norm -f "${GENOME}" --threads $(( CPUS - 1 )) -Oz -o "${TMPVCF}"
mv "${TMPVCF}" "${CHUNKVCF}"
rm "${REGIONS}"
