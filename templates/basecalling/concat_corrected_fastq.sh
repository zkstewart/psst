####

# Specify the location of the Various_scripts and Genome_analysis_scripts folders
VARSCRIPTDIR=@@varScriptDir@@
GENSCRIPTDIR=@@genScriptDir@@

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
    python ${VARSCRIPTDIR}/fasta_concat.py \
        -i ${INDIR} \
        -o ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
        --minimum ${MINIMUM} --multiline ${MULTILINE};
fi

# STEP 2: Obtain statistics for the reads
python ${GENSCRIPTDIR}/genome_stats.py \
    -i ${OUTDIR}/${PREFIX}_corrected_reads.fasta \
    -o ${OUTDIR}/${PREFIX}_corrected_reads.stats
