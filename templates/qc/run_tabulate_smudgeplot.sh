####

# Specify the reads directory and the suffix identifying read files, as given to smudgeplot
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the directory of smudgeplot results (one subdirectory per sample) and the table
# of every sample's most likely ploidy to write
OUTDIR=@@outputDir@@
TABLE=@@outputFile@@

####

set -eo pipefail
shopt -s nullglob

# STEP 1: Refuse to overwrite an existing table
if [[ -f ${TABLE} ]]; then
    echo "${TABLE} already exists; remove it to tabulate again"
    exit 1
fi

# STEP 2: Locate every sample's ploidy summary
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
if (( ${#READFILES[@]} == 0 )); then
    echo "No read files ending with ${SUFFIX} were found in ${READSDIR}"
    exit 1
fi
SUMMARIES=()
for R1 in "${READFILES[@]}"; do
    BASE=$(basename "${R1}")
    if [[ ${PAIRED} == "true" ]]; then
        SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${BASE%1"${SUFFIX}"}")
    else
        SAMPLE=${BASE%"${SUFFIX}"}
    fi
    SUMMARY=${OUTDIR}/${SAMPLE}/${SAMPLE}.ploidy.tsv
    if [[ ! -f ${SUMMARY} ]]; then
        echo "${SUMMARY} was not found; has smudgeplot finished for every sample?"
        exit 1
    fi
    SUMMARIES+=( "${SUMMARY}" )
done

# STEP 3: Combine the summaries under one header, giving the table its final name once complete
TMPTABLE=${TABLE}.incomplete
awk 'FNR > 1 || NR == 1' "${SUMMARIES[@]}" > "${TMPTABLE}"
mv "${TMPTABLE}" "${TABLE}"
