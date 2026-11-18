####

# Specify the conda environment containing kraken2, KrakenTools and Krona
CONDAENV=@@condaEnv@@

# Specify the kraken2 database directory
DBDIR=@@db@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory; each sample's results are written to <outdir>/<sample>/
OUTDIR=@@outputDir@@

# Specify computational parameters
CPUS=@@cpus@@

####

conda activate ${CONDAENV}
set -eo pipefail
shopt -s nullglob

# STEP 1: Find the read files of every sample
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
if (( ${#READFILES[@]} == 0 )); then
    echo "No read files ending with ${SUFFIX} were found in ${READSDIR}"
    exit 1
fi

# STEP 2: Read the database once into the node's file cache, so that each sample's
# classification maps it from memory (--memory-mapping) rather than loading its own copy;
# this is only done once a sample needs classifying
DBCACHED=false
cache_database() {
    if [[ ${DBCACHED} == "false" ]]; then
        echo "Reading the kraken2 database (${DBDIR}) into the file cache"
        cat "${DBDIR}"/*.k2d > /dev/null
        DBCACHED=true
    fi
}

# STEP 3: Classify each sample in turn
for R1 in "${READFILES[@]}"; do
    # Locate the read pair (if paired) and name the sample; paired classified and unclassified
    # reads are written as <sample>_1_classified.fastq, <sample>_2_classified.fastq, ...
    BASE=$(basename "${R1}")
    if [[ ${PAIRED} == "true" ]]; then
        PREFIX=${BASE%1"${SUFFIX}"}
        R2=${READSDIR}/${PREFIX}2${SUFFIX}
        if [[ ! -f ${R2} ]]; then
            echo "${R1} has no matching reverse read file (${R2})"
            exit 1
        fi
        SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
        READARGS=( --paired "${R1}" "${R2}" )
        PAIRTAG="#" # replaced by kraken2 with _1 and _2
    else
        SAMPLE=${BASE%"${SUFFIX}"}
        READARGS=( "${R1}" )
        PAIRTAG=""
    fi
    
    # Skip samples which have already been classified
    SAMPLEDIR=${OUTDIR}/${SAMPLE}
    if [[ -d ${SAMPLEDIR} ]]; then
        echo "${SAMPLEDIR} already exists; skipping ${SAMPLE}"
        continue
    fi
    INCOMPLETE=${OUTDIR}/${SAMPLE}.incomplete
    rm -rf "${INCOMPLETE}" # left by a failed run
    mkdir -p "${INCOMPLETE}"
    
    # Classify the reads
    cache_database
    echo "Classifying ${SAMPLE} ($(date))"
    k2 classify --threads ${CPUS} --db "${DBDIR}" --memory-mapping \
        --output "${INCOMPLETE}/${SAMPLE}_output.txt" \
        --unclassified-out "${INCOMPLETE}/${SAMPLE}${PAIRTAG}_unclassified.fastq" \
        --classified-out "${INCOMPLETE}/${SAMPLE}${PAIRTAG}_classified.fastq" \
        --report "${INCOMPLETE}/${SAMPLE}_report.txt" \
        "${READARGS[@]}"
    
    # Chart the classifications with Krona
    kreport2krona.py -r "${INCOMPLETE}/${SAMPLE}_report.txt" -o "${INCOMPLETE}/${SAMPLE}_report.krona"
    ktImportText "${INCOMPLETE}/${SAMPLE}_report.krona" -o "${INCOMPLETE}/${SAMPLE}_report.html"
    
    # Give the results their final location once complete
    mv "${INCOMPLETE}" "${SAMPLEDIR}"
done
