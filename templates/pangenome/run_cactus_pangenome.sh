####

# Specify the conda environment containing cactus-pangenome, and how cactus runs its tools
# (singularity, docker or local)
CONDAENV=@@condaEnv@@
BINARIESMODE=@@binariesMode@@

# Specify the genomes as seqfile names ('SAMPLE.N' for haplotype N of a sample, where 0 is
# a haploid sample, or just the reference's name) and their FASTA files, in matching order
NAMES=( @@seqNames@@ )
GENOMES=( @@genomeFiles@@ )

# Specify the name of the reference genome, which must be haploid and chromosome scale
REFERENCE=@@refName@@

# Specify the output directory and prefix of output files
OUTDIR=@@outputDir@@
PREFIX=@@prefix@@

# Specify the graph to make giraffe indexes of (filter, clip or full), and whether to
# realign VCF alleles with vcfwave
GIRAFFE=@@giraffe@@
VCFWAVE=@@vcfwave@@

# Specify computational parameters
CPUS=@@cpus@@
MAXMEMORY=@@maxMemory@@

####

conda activate ${CONDAENV}
set -eo pipefail

# STEP 1: Refuse to rebuild a completed pangenome
SEQFILE=${OUTDIR}/cactus_samples.txt
JOBSTORE=${OUTDIR}/jobstore
GRAPHDIR=${OUTDIR}/graph
if [[ -f ${OUTDIR}/cactus.ok ]]; then
    echo "${OUTDIR} already holds a completed pangenome; remove it to rebuild"
    exit 1
fi

# STEP 2: Write the seqfile; a seqfile from a previous run of this job must match it, as
# cactus resumes that run
if (( ${#NAMES[@]} != ${#GENOMES[@]} )); then
    echo "NAMES and GENOMES must have the same number of values"
    exit 1
fi
mkdir -p "${OUTDIR}"
SEQLINES=$(for i in "${!NAMES[@]}"; do printf '%s\t%s\n' "${NAMES[$i]}" "${GENOMES[$i]}"; done)
if [[ -f ${SEQFILE} ]]; then
    if [[ $(cat "${SEQFILE}") != "${SEQLINES}" ]]; then
        echo "${SEQFILE} differs from the genomes of this job; remove ${OUTDIR} to start again"
        exit 1
    fi
else
    echo "${SEQLINES}" > "${SEQFILE}"
fi

# STEP 3: Build the pangenome, resuming a previous run of this job if its job store
# remains (cactus removes it on success); cactus writes outputs to a separate directory
# since it exports a modified copy of the seqfile there
CACTUSARGS=( --outDir "${GRAPHDIR}" --outName ${PREFIX} --reference ${REFERENCE} \
    --vcf --giraffe ${GIRAFFE} --maxCores ${CPUS} --maxMemory ${MAXMEMORY} \
    --binariesMode ${BINARIESMODE} )
[[ ${VCFWAVE} == "true" ]] && CACTUSARGS+=( --vcfwave )
[[ -d ${JOBSTORE} ]] && CACTUSARGS+=( --restart )
cactus-pangenome "${JOBSTORE}" "${SEQFILE}" "${CACTUSARGS[@]}"
touch "${OUTDIR}/cactus.ok"
