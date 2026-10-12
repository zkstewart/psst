module load zlib/1.3.1
export LD_LIBRARY_PATH=/mnt/weka/pkg/rhel94/AuthenticAMD-25/software/zlib/1.3.1/lib:${LD_LIBRARY_PATH}

####

# Specify parameters for this sample / basecalling run
PREFIX=@@prefix@@
MINREADLEN=@@minReadLen@@

# Specify the location of the dorado bin and lib folders
DORADODIR=@@doradoDir@@

# Specify the location of the Various_scripts folder
VARSCRIPTDIR=@@varScriptDir@@

# Specify input and output read file locations
FQDIR=fastq
OUTDIR=fastq_corrected

# Specify the job name prefix and resources of the step 1 (CPU) and step 2 (GPU) correction jobs
JOBPREFIX=@@jobPrefix@@
S1CPUS=@@s1Cpus@@
S1MEM=@@s1Mem@@
S1WALLTIME=@@s1Walltime@@
S2CPUS=@@s2Cpus@@
S2MEM=@@s2Mem@@
S2WALLTIME=@@s2Walltime@@

####

# STEP 0: Derive file locations (following on from basecaller_part2.sh)
FQFILE=${PREFIX}.min${MINREADLEN}.fastq
NUMBLOCKS=$(cat ${FQDIR}/numblocks.txt)
mkdir -p ${OUTDIR}

# STEP 1: Run the pipeliner script which submits the step 1 and step 2 correction jobs
python ${VARSCRIPTDIR}/basecalling/dorado_correct_pipeline.py -f ${FQDIR}/${FQFILE} \
    -n ${NUMBLOCKS} -o ${OUTDIR} \
    --dorado ${DORADODIR}/bin/dorado \
    --jobPrefix ${JOBPREFIX} \
    --cpuCPU ${S1CPUS} --memCPU ${S1MEM} --walltimeCPU ${S1WALLTIME} \
    --cpuGPU ${S2CPUS} --memGPU ${S2MEM} --walltimeGPU ${S2WALLTIME} | tee ${OUTDIR}/dorado_correct_pipeline.log

# STEP 2: Submit the concatenation job to run after the step 2 correction job succeeds
STEP2JOBID=$(grep -oP "step2=\K\S+" ${OUTDIR}/dorado_correct_pipeline.log)
if [[ -z ${STEP2JOBID} ]]; then
    echo "Unable to find the step 2 job ID; submit ${OUTDIR}/concat_corrected_fastq.sh manually"
    exit 1
fi
cd ${OUTDIR}
qsub -W depend=afterok:${STEP2JOBID} concat_corrected_fastq.sh
