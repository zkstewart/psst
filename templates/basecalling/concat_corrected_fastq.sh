####

# Specify the annotarium location
ANNOTARIUMDIR=@@annotariumDir@@

# Specify the prefix of the corrected reads file
PREFIX=@@prefix@@

# Specify script behaviours
MINIMUM=@@concatMinimum@@
MULTILINE=@@concatMultiline@@

####

# STEP 0: Set hard-coded variables
INDIR=step2
OUTDIR=output
mkdir -p ${OUTDIR}

# STEP 1: Concatenate individual FASTA files
if [[ ! -n "$(ls -A ${OUTDIR} 2>/dev/null)" ]]; then
    python ${ANNOTARIUMDIR}/annotarium.py fasta concat \
        -i ${INDIR} \
        -o ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
        --minimum ${MINIMUM} --multiline ${MULTILINE};
fi

# STEP 2: Obtain statistics for the reads (a row per read, then summary rows), replacing any
# from a previous run; annotarium refuses to overwrite, so they're written under a temporary
# name first
STATS=${OUTDIR}/${PREFIX}_corrected_reads.stats
rm -f "${STATS}.incomplete"
python ${ANNOTARIUMDIR}/annotarium.py fasta stats \
    -i ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
    -o "${STATS}.incomplete" && mv -f "${STATS}.incomplete" "${STATS}"
