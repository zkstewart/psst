####

# Specify the FastQC module to load
FASTQCMODULE=@@fastqcModule@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory; each sample's reports are written to <outdir>/<sample>/
OUTDIR=@@outputDir@@

# Specify computational parameters
CPUS=@@cpus@@

####

module load ${FASTQCMODULE}
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

# STEP 3: Skip samples which have already been checked
SAMPLEDIR=${OUTDIR}/${SAMPLE}
if [[ -d ${SAMPLEDIR} ]]; then
    echo "${SAMPLEDIR} already exists; skipping ${SAMPLE}"
    exit 0
fi

# STEP 4: Run FastQC, giving the reports their final location once complete
INCOMPLETE=${OUTDIR}/${SAMPLE}.incomplete
rm -rf "${INCOMPLETE}" # left by a failed run
mkdir -p "${INCOMPLETE}"
fastqc --outdir "${INCOMPLETE}" --threads ${CPUS} "${READS[@]}" # long options work in FastQC 0.12 and 0.13
mv "${INCOMPLETE}" "${SAMPLEDIR}"
