module load zlib/1.3.1
module load SAMtools/1.18
export LD_LIBRARY_PATH=/mnt/weka/pkg/rhel94/AuthenticAMD-25/software/zlib/1.3.1/lib:${LD_LIBRARY_PATH}

####

# Specify the program locations
DORADODIR=@@doradoDir@@
VARSCRIPTDIR=@@varScriptDir@@

# Specify the prefix for output files
OUTPREFIX=@@prefix@@

# Specify program behavioural values
MODEL=@@model@@
MINREADLEN=@@minReadLen@@

####

# STEP 0: Export and set all hard-coded variables
POD5DIR=pod5
BAMDIR=bam_sup
FQDIR=fastq
MODELSUFFIX=$(echo ${MODEL} | tr , - | tr @ _)
export LD_LIBRARY_PATH=${DORADODIR}/lib:${LD_LIBRARY_PATH}

# STEP 2: Run samtools to convert BAM to FASTQ
mkdir -p ${FQDIR}
BAMRESULTFILE=${BAMDIR}/${OUTPREFIX}.${MODELSUFFIX}.bam
FQRESULTFILE=${FQDIR}/${OUTPREFIX}.fastq
if [[ ! -f ${FQRESULTFILE} ]]; then
    samtools fastq -T "*" ${BAMRESULTFILE} > ${FQRESULTFILE};
fi

# STEP 3: Filter FASTQ reads to minimum size cutoff
FQFILTERFILE=${FQDIR}/${OUTPREFIX}.min${MINREADLEN}.fastq
if [[ ! -f ${FQDIR}/${OUTPREFIX}.min${MINREADLEN}.fastq ]]; then
    python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f cullbelow -i ${FQRESULTFILE} -n ${MINREADLEN} -o ${FQFILTERFILE};
fi

# STEP 4: Index the FASTQ file
FQINDEXFILE=${FQFILTERFILE}.fai
if [[ ! -f ${FQINDEXFILE} ]]; then
    samtools faidx ${FQFILTERFILE};
fi

# STEP 5: Figure out how many blocks we can chunk the FASTQ into for correction
NUMBLOCKS=$(${DORADODIR}/bin/dorado correct ${FQFILTERFILE} --compute-num-blocks)
echo "NUMBLOCKS = ${NUMBLOCKS}"
echo "${NUMBLOCKS}" > ${FQDIR}/numblocks.txt
