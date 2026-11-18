####

# Specify the Various_scripts location
VARSCRIPTDIR=@@varScriptDir@@

# Specify the reads directory and the suffix identifying read files, as given to FastQC
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the directory of FastQC reports (one subdirectory per sample) and the Excel
# summary table to write
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

# STEP 2: Check that every sample has FastQC reports, since tabulate_fastqc.py summarises
# every subdirectory of the reports directory as a sample
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
if (( ${#READFILES[@]} == 0 )); then
    echo "No read files ending with ${SUFFIX} were found in ${READSDIR}"
    exit 1
fi
for R1 in "${READFILES[@]}"; do
    BASE=$(basename "${R1}")
    if [[ ${PAIRED} == "true" ]]; then
        SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${BASE%1"${SUFFIX}"}")
    else
        SAMPLE=${BASE%"${SUFFIX}"}
    fi
    if [[ ! -d ${OUTDIR}/${SAMPLE} ]]; then
        echo "${OUTDIR}/${SAMPLE} was not found; has FastQC finished for every sample?"
        exit 1
    fi
done
INCOMPLETE=( "${OUTDIR}"/*.incomplete/ )
if (( ${#INCOMPLETE[@]} > 0 )); then
    echo "${OUTDIR} contains incomplete FastQC output (${INCOMPLETE[*]}); remove it to tabulate"
    exit 1
fi

# STEP 3: Tabulate the reports, giving the table its final name once complete
TMPTABLE=${TABLE%.xlsx}.incomplete.xlsx
rm -f "${TMPTABLE}" # left by a failed run
python "${VARSCRIPTDIR}/QC/reads/tabulate_fastqc.py" -i "${OUTDIR}" -o "${TMPTABLE}"
mv "${TMPTABLE}" "${TABLE}"
