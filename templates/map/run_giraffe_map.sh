####

# Specify the vg executable location
VGDIR=@@vgDir@@

# Specify the index directory created by 'psst map index giraffe'
INDEX=@@indexDir@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory
OUTDIR=@@outputDir@@

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
INDEX_NUM=$(( ${PBS_ARRAY_INDEX:-1} - 1 )) # PBS_ARRAY_INDEX is unset if there is only one sample
if (( INDEX_NUM >= ${#READFILES[@]} )); then
    echo "Array index ${PBS_ARRAY_INDEX} exceeds the ${#READFILES[@]} read files found in ${READSDIR}"
    exit 1
fi
R1=${READFILES[${INDEX_NUM}]}
BASE=$(basename "${R1}")

# STEP 2: Locate the read pair (if paired) and name the sample
if [[ ${PAIRED} == "true" ]]; then
    PREFIX=${BASE%1"${SUFFIX}"}
    R2=${READSDIR}/${PREFIX}2${SUFFIX}
    if [[ ! -f ${R2} ]]; then
        echo "${R1} has no matching reverse read file (${R2})"
        exit 1
    fi
    READARGS=( -f "${R1}" -f "${R2}" )
    SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
else
    READARGS=( -f "${R1}" )
    SAMPLE=${BASE%"${SUFFIX}"}
fi

# STEP 3: Skip samples which have already been mapped
OUTBAM=${OUTDIR}/${SAMPLE}.sorted.md.bam
if [[ -f ${OUTBAM} ]]; then
    echo "${OUTBAM} already exists; skipping ${SAMPLE}"
    exit 0
fi
mkdir -p "${OUTDIR}"

# STEP 4: Map reads to the graph, outputting alignments projected onto the reference paths,
# then sort and mark duplicates without intermediate files; indels are left-aligned so
# that they sit where bcftools expects them (e.g., at pangenome variant sites); giraffe
# outputs mates together, so they need no collation before fixmate
TMPBAM=${OUTDIR}/${SAMPLE}.incomplete.bam
${VGDIR}/vg giraffe -t ${CPUS} --left-align \
    -Z "${INDEX}/pangenome.gbz" -d "${INDEX}/pangenome.dist" \
    -m "${INDEX}/pangenome.min" -z "${INDEX}/pangenome.zipcodes" \
    "${READARGS[@]}" -N "${SAMPLE}" -R "${SAMPLE}" -o BAM \
    | samtools fixmate -m -u - - \
    | samtools sort -u -@ ${CPUS} -m ${SORTMEM} -T "${OUTDIR}/${SAMPLE}.sort_tmp" - \
    | samtools markdup -@ ${CPUS} -f "${OUTDIR}/${SAMPLE}.markdup.stats" - "${TMPBAM}"

# STEP 5: Index and summarise the BAM, then give it its final name once complete
samtools index -@ ${CPUS} "${TMPBAM}" "${OUTBAM}.bai"
samtools flagstat -@ ${CPUS} "${TMPBAM}" > "${OUTDIR}/${SAMPLE}.sorted.md.flagstat"
mv "${TMPBAM}" "${OUTBAM}"
