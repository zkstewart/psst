####

# Specify program locations
MINIMAP2DIR=@@minimap2Dir@@
CHOPPERDIR=@@chopperDir@@

# Specify the reference genome FASTA
GENOME=@@genome@@

# Specify the reads directory and the suffix identifying read files
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@

# Specify the output directory
OUTDIR=@@outputDir@@

# Specify mapping parameters
PRESET=@@preset@@ # minimap2 preset
LIBRARY=@@library@@

# Specify read filtering parameters for chopper (0 to disable each filter)
MINLENGTH=@@minLength@@
MAXLENGTH=@@maxLength@@
MINQUAL=@@minQual@@

# Specify computational parameters
CPUS=@@cpus@@
SORTMEM=@@sortMem@@ # memory per samtools sort thread

####

set -eo pipefail
shopt -s nullglob

# STEP 1: Find the read file for this array index (one sample per subjob)
READFILES=( "${READSDIR}"/*"${SUFFIX}" )
INDEX_NUM=$(( ${PBS_ARRAY_INDEX:-1} - 1 )) # PBS_ARRAY_INDEX is unset if there is only one sample
if (( INDEX_NUM >= ${#READFILES[@]} )); then
    echo "Array index ${PBS_ARRAY_INDEX} exceeds the ${#READFILES[@]} read files found in ${READSDIR}"
    exit 1
fi
READS=${READFILES[${INDEX_NUM}]}
SAMPLE=$(basename "${READS}" "${SUFFIX}")

# STEP 2: Skip samples which have already been mapped
OUTBAM=${OUTDIR}/${SAMPLE}.sorted.bam
if [[ -f ${OUTBAM} ]]; then
    echo "${OUTBAM} already exists; skipping ${SAMPLE}"
    exit 0
fi
mkdir -p "${OUTDIR}"

# STEP 3: Set up read filtering with chopper, if any filter is enabled
CHOPPERARGS=()
(( MINLENGTH > 0 )) && CHOPPERARGS+=( --minlength ${MINLENGTH} )
(( MAXLENGTH > 0 )) && CHOPPERARGS+=( --maxlength ${MAXLENGTH} )
(( MINQUAL > 0 )) && CHOPPERARGS+=( --quality ${MINQUAL} )
filter_reads() {
    if (( ${#CHOPPERARGS[@]} > 0 )); then
        ${CHOPPERDIR}/chopper --threads ${CPUS} "${CHOPPERARGS[@]}" 2> "${OUTDIR}/${SAMPLE}.chopper.log"
    else
        cat
    fi
}

# STEP 4: Filter, map and sort without intermediate files
TMPBAM=${OUTDIR}/${SAMPLE}.incomplete.bam
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ONT\tLB:${LIBRARY}"
zcat -f "${READS}" \
    | filter_reads \
    | ${MINIMAP2DIR}/minimap2 -t ${CPUS} -ax ${PRESET} -R "${RG}" "${GENOME}" - \
    | samtools sort -@ ${CPUS} -m ${SORTMEM} -T "${OUTDIR}/${SAMPLE}.sort_tmp" -o "${TMPBAM}" -

# STEP 5: Index and summarise the BAM, then give it its final name once complete
samtools index -@ ${CPUS} "${TMPBAM}" "${OUTBAM}.bai"
samtools flagstat -@ ${CPUS} "${TMPBAM}" > "${OUTDIR}/${SAMPLE}.sorted.flagstat"
mv "${TMPBAM}" "${OUTBAM}"
