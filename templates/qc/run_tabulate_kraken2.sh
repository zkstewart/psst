####

# Specify the conda environment containing KrakenTools (and pandas), and the
# Various_scripts location
CONDAENV=@@condaEnv@@
VARSCRIPTDIR=@@varScriptDir@@

# Specify the reads directory and the suffix identifying read files, as given to kraken2
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the directory of kraken2 results (one subdirectory per sample), and the prefix
# of the combined tables written there (<prefix>.tsv and <prefix>.reformat.tsv)
OUTDIR=@@outputDir@@
PREFIX=@@prefix@@

####

conda activate ${CONDAENV}
set -eo pipefail
shopt -s nullglob

# STEP 1: Refuse to overwrite existing tables
COMBINED=${OUTDIR}/${PREFIX}.tsv
REFORMAT=${OUTDIR}/${PREFIX}.reformat.tsv
for f in "${COMBINED}" "${REFORMAT}"; do
    if [[ -f ${f} ]]; then
        echo "${f} already exists; remove it to tabulate again"
        exit 1
    fi
done

# STEP 2: Locate every sample's kraken2 report
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
if (( ${#READFILES[@]} == 0 )); then
    echo "No read files ending with ${SUFFIX} were found in ${READSDIR}"
    exit 1
fi
REPORTS=()
for R1 in "${READFILES[@]}"; do
    BASE=$(basename "${R1}")
    if [[ ${PAIRED} == "true" ]]; then
        SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${BASE%1"${SUFFIX}"}")
    else
        SAMPLE=${BASE%"${SUFFIX}"}
    fi
    REPORT=${OUTDIR}/${SAMPLE}/${SAMPLE}_report.txt
    if [[ ! -f ${REPORT} ]]; then
        echo "${REPORT} was not found; has kraken2 finished for every sample?"
        exit 1
    fi
    REPORTS+=( "${REPORT}" )
done

# STEP 3: Combine the reports, then reformat them as the fraction of each sample's reads
# assigned to each taxon, giving the tables their final names once both are complete
TMPCOMBINED=${OUTDIR}/${PREFIX}.incomplete.tsv
TMPREFORMAT=${OUTDIR}/${PREFIX}.reformat.incomplete.tsv
rm -f "${TMPCOMBINED}" "${TMPREFORMAT}" # left by a failed run
combine_kreports.py -r "${REPORTS[@]}" -o "${TMPCOMBINED}"
python "${VARSCRIPTDIR}/kraken/kreport_reformat.py" -i "${TMPCOMBINED}" -o "${TMPREFORMAT}"
mv "${TMPCOMBINED}" "${COMBINED}"
mv "${TMPREFORMAT}" "${REFORMAT}"
