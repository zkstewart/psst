####

# Specify program locations
MMSEQDIR=@@mmseqsDir@@
VARSCRIPTDIR=@@varScriptDir@@

# Specify query FASTA location
QUERY=@@query@@

# Specify reference FASTA location
DB=@@databaseTag@@
DATABASE=@@database@@

# Specify database annotation file locations
IM=@@idmappingFile@@
IO=@@oboFile@@

# Specify program parameters
CPUS=@@cpus@@
NUMHITS=@@numHits@@
NUM_ITERS=@@iterations@@

# Specify prefix for output files
OUTPREFIX=@@outputPrefix@@

# Specify tmp directory
MMS2TMP=@@tmpDir@@

####

# STEP 0: Set PATH variable
export PATH=${MMSEQDIR}/:${PATH}

# STEP 1: Run MMseqs2
python ${VARSCRIPTDIR}/run_mmseqs2.py -q ${QUERY} \
    -t ${DATABASE} \
    -o ${OUTPREFIX} \
    -m ${MMSEQDIR} -c ${CPUS} \
    -n ${NUM_ITERS} \
    --resume --blast_sort \
    --tmpDir ${MMS2TMP} \
    -st auto

# STEP 2: Generate the annotation table
python ${VARSCRIPTDIR}/annotation_table/generate_annotation_table.py -rf ${QUERY} \
    -bf ${DATABASE} -bo ${OUTPREFIX}_mms2SEARCH_sorted.m8 \
    -id ${IM} -io ${IO} \
    -o ${OUTPREFIX}.annotation.tsv \
    --numhits ${NUMHITS} \
    --db ${DB}
