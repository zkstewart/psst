module load zlib/1.3.1
export LD_LIBRARY_PATH=/mnt/weka/pkg/rhel94/AuthenticAMD-25/software/zlib/1.3.1/lib:${LD_LIBRARY_PATH}
module load CUDA/12.8.0

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

# STEP 1: Run dorado
mkdir -p ${BAMDIR}
BAMRESULTFILE=${BAMDIR}/${OUTPREFIX}.${MODELSUFFIX}.bam
if [[ ! -f ${BAMRESULTFILE} ]]; then
    ${DORADODIR}/bin/dorado basecaller ${MODEL} ${POD5DIR}/ > ${BAMRESULTFILE};
fi
