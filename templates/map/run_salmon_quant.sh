####

# Specify the salmon executable location
SALMONDIR=@@salmonDir@@

# Specify the salmon index directory
INDEX=@@indexDir@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory, which will contain a subdirectory for each sample
OUTDIR=@@outputDir@@

# Specify quantification parameters
LIBTYPE=@@libType@@ # A to detect automatically
GENEMAP=@@geneMap@@ # GTF/GFF or transcript-to-gene TSV for gene-level estimates; empty for none
GCBIAS=@@gcBias@@ # true or false; fragment GC bias correction (paired reads only)
SEQBIAS=@@seqBias@@ # true or false; sequence-specific bias correction
BOOTSTRAPS=@@bootstraps@@

# Specify computational parameters
CPUS=@@cpus@@

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
    READARGS=( -1 "${R1}" -2 "${R2}" )
    SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
else
    READARGS=( -r "${R1}" )
    SAMPLE=${BASE%"${SUFFIX}"}
fi

# STEP 3: Skip samples which have already been quantified
OUTSAMPLE=${OUTDIR}/${SAMPLE}
if [[ -d ${OUTSAMPLE} ]]; then
    echo "${OUTSAMPLE} already exists; skipping ${SAMPLE}"
    exit 0
fi
mkdir -p "${OUTDIR}"

# STEP 4: Quantify, giving the output its final name once complete; salmon runs one
# thread for reading input beyond --threads, so leave it one of the requested CPUs
OPTARGS=()
[[ -n ${GENEMAP} ]] && OPTARGS+=( --geneMap "${GENEMAP}" )
[[ ${GCBIAS} == "true" ]] && OPTARGS+=( --gcBias )
[[ ${SEQBIAS} == "true" ]] && OPTARGS+=( --seqBias )
TMPSAMPLE=${OUTSAMPLE}.incomplete
rm -rf "${TMPSAMPLE}" # remove any output left by a previous failed run
${SALMONDIR}/salmon quant --threads $(( CPUS > 1 ? CPUS - 1 : 1 )) \
    --index "${INDEX}" --libType ${LIBTYPE} "${READARGS[@]}" \
    "${OPTARGS[@]}" --numBootstraps ${BOOTSTRAPS} \
    --output "${TMPSAMPLE}"
mv "${TMPSAMPLE}" "${OUTSAMPLE}"
