####

# Specify program locations
MMSEQDIR=@@mmseqsDir@@
VARSCRIPTDIR=@@varScriptDir@@

# Specify query FASTA location
QUERY=@@query@@

# Specify target FASTA location
DATABASE=@@database@@

# Specify output details
OUTPREFIX=@@outputPrefix@@

# Specify tmp directory
MMS2TMP=@@tmpDir@@

# Specify program parameters
CPUS=@@cpus@@
NUM_ITERS=@@iterations@@
SENS=@@sensitivity@@
SEARCHTYPE=@@searchType@@

####

# Step 1: Add MMSEQS to path
export PATH=${MMSEQDIR}/:${PATH}

# Step 2: Run MMseqs2
python ${VARSCRIPTDIR}/run_mmseqs2.py -st ${SEARCHTYPE} \
    -q ${QUERY} -t ${DATABASE} \
    -o ${OUTPREFIX} \
    -m ${MMSEQDIR} \
    -c ${CPUS} \
    -n ${NUM_ITERS} -s ${SENS} \
    --resume \
    --blast_sort \
    --tmpDir ${MMS2TMP}
