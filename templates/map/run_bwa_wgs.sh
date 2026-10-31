####

# Specify the reference genome FASTA, which must be indexed with 'bwa index'
GENOME=@@genome@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory
OUTDIR=@@outputDir@@

# Specify the library name for read groups
LIBRARY=@@library@@

# Specify computational parameters
CPUS=@@cpus@@
SORTMEM=@@sortMem@@ # memory per samtools sort thread

####

set -eo pipefail
shopt -s nullglob

# STEP 1: Find the read files for this array index (one sample per subjob)
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
INDEX=$(( ${PBS_ARRAY_INDEX:-1} - 1 )) # PBS_ARRAY_INDEX is unset if there is only one sample
if (( INDEX >= ${#READFILES[@]} )); then
    echo "Array index ${PBS_ARRAY_INDEX} exceeds the ${#READFILES[@]} read files found in ${READSDIR}"
    exit 1
fi
R1=${READFILES[${INDEX}]}
BASE=$(basename "${R1}")

# STEP 2: Locate the read pair (if paired) and name the sample
if [[ ${PAIRED} == "true" ]]; then
    PREFIX=${BASE%1"${SUFFIX}"}
    R2=${READSDIR}/${PREFIX}2${SUFFIX}
    if [[ ! -f ${R2} ]]; then
        echo "${R1} has no matching reverse read file (${R2})"
        exit 1
    fi
    READS=( "${R1}" "${R2}" )
    SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
else
    READS=( "${R1}" )
    SAMPLE=${BASE%"${SUFFIX}"}
fi

# STEP 3: Skip samples which have already been mapped
OUTBAM=${OUTDIR}/${SAMPLE}.sorted.md.bam
if [[ -f ${OUTBAM} ]]; then
    echo "${OUTBAM} already exists; skipping ${SAMPLE}"
    exit 0
fi
mkdir -p "${OUTDIR}"

# STEP 4: Map, sort and mark duplicates without intermediate files
TMPBAM=${OUTDIR}/${SAMPLE}.incomplete.bam
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${LIBRARY}\tPU:${SAMPLE}"
bwa mem -t ${CPUS} -R "${RG}" "${GENOME}" "${READS[@]}" \
    | samtools fixmate -m -u - - \
    | samtools sort -u -@ ${CPUS} -m ${SORTMEM} -T "${OUTDIR}/${SAMPLE}.sort_tmp" - \
    | samtools markdup -@ ${CPUS} -f "${OUTDIR}/${SAMPLE}.markdup.stats" - "${TMPBAM}"

# STEP 5: Index and summarise the BAM, then give it its final name once complete
samtools index -@ ${CPUS} "${TMPBAM}" "${OUTBAM}.bai"
samtools flagstat -@ ${CPUS} "${TMPBAM}" > "${OUTDIR}/${SAMPLE}.sorted.md.flagstat"
mv "${TMPBAM}" "${OUTBAM}"
