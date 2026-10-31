####

# Specify the STAR module to load
STARMODULE=@@starModule@@

# Specify the STAR index directory
INDEX=@@indexDir@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory, which will contain a subdirectory for each sample
OUTDIR=@@outputDir@@

# Specify mapping parameters
TWOPASS=@@twopass@@ # true or false; re-map each sample using its novel splice junctions

# Specify computational parameters
CPUS=@@cpus@@

####

module load ${STARMODULE}
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
    READS=( "${R1}" "${R2}" )
    SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
else
    READS=( "${R1}" )
    SAMPLE=${BASE%"${SUFFIX}"}
fi

# STEP 3: Skip samples which have already been mapped
OUTSAMPLE=${OUTDIR}/${SAMPLE}
if [[ -d ${OUTSAMPLE} ]]; then
    echo "${OUTSAMPLE} already exists; skipping ${SAMPLE}"
    exit 0
fi
mkdir -p "${OUTDIR}"

# STEP 4: Map reads, giving the output its final name once complete
OPTARGS=()
[[ ${READS[0]} == *.gz ]] && OPTARGS+=( --readFilesCommand zcat )
[[ ${TWOPASS} == "true" ]] && OPTARGS+=( --twopassMode Basic )
TMPSAMPLE=${OUTSAMPLE}.incomplete
rm -rf "${TMPSAMPLE}" # remove any output left by a previous failed run
mkdir -p "${TMPSAMPLE}"
STAR --runThreadN ${CPUS} \
     --genomeDir "${INDEX}" \
     --readFilesIn "${READS[@]}" \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattrRGline ID:${SAMPLE} SM:${SAMPLE} PL:ILLUMINA \
     --quantMode TranscriptomeSAM GeneCounts \
     "${OPTARGS[@]}" \
     --outFileNamePrefix "${TMPSAMPLE}/"
mv "${TMPSAMPLE}" "${OUTSAMPLE}"
