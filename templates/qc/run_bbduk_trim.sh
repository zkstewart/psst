####

# Specify the conda environment containing BBTools (e.g., conda create -n readqc bioconda::bbmap)
CONDAENV=@@condaEnv@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory; trimmed reads are written as <sample>_1.fq.gz and
# <sample>_2.fq.gz (paired) or <sample>.fq.gz (single-end)
OUTDIR=@@outputDir@@

# Specify the adapter sequences to trim; 'adapters' uses those bundled with BBTools
ADAPTERS=@@adapters@@

# Specify the quality to trim read ends to, and the minimum length of reads to keep
TRIMQ=@@trimq@@
MINLEN=@@minLength@@

# Specify computational parameters
CPUS=@@cpus@@
XMX=@@xmx@@ # Java heap for bbduk; it otherwise sizes this from the whole node's memory

####

conda activate ${CONDAENV}
set -eo pipefail
shopt -s nullglob

# STEP 1: Find the read files for this array index (one sample per subjob)
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
INDEX=$(( ${PBS_ARRAY_INDEX:-1} - 1 )) # PBS_ARRAY_INDEX is unset if there is only one sample
if (( INDEX >= ${#READFILES[@]} )); then
    echo "Array index ${PBS_ARRAY_INDEX} exceeds the ${#READFILES[@]} read files found in ${READSDIR}"
    exit 1
fi
R1=${READFILES[${INDEX}]}
BASE=$(basename "${R1}")

# STEP 2: Locate the read pair (if paired), name the sample and its outputs
if [[ ${PAIRED} == "true" ]]; then
    PREFIX=${BASE%1"${SUFFIX}"}
    R2=${READSDIR}/${PREFIX}2${SUFFIX}
    if [[ ! -f ${R2} ]]; then
        echo "${R1} has no matching reverse read file (${R2})"
        exit 1
    fi
    SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
    OUTFILES=( "${OUTDIR}/${SAMPLE}_1.fq.gz" "${OUTDIR}/${SAMPLE}_2.fq.gz" )
    TMPFILES=( "${OUTDIR}/${SAMPLE}.incomplete_1.fq.gz" "${OUTDIR}/${SAMPLE}.incomplete_2.fq.gz" )
    IOARGS=( in1="${R1}" in2="${R2}" out1="${TMPFILES[0]}" out2="${TMPFILES[1]}" tpe tbo )
else
    SAMPLE=${BASE%"${SUFFIX}"}
    OUTFILES=( "${OUTDIR}/${SAMPLE}.fq.gz" )
    TMPFILES=( "${OUTDIR}/${SAMPLE}.incomplete.fq.gz" )
    IOARGS=( in="${R1}" out="${TMPFILES[0]}" )
fi

# STEP 3: Skip samples which have already been trimmed, refusing to overwrite partial output
# (or the input reads, which single-end output would share names with)
if [[ $(realpath -m "${OUTDIR}") == $(realpath -m "${READSDIR}") ]]; then
    echo "The output directory (${OUTDIR}) must differ from the reads directory"
    exit 1
fi
EXISTING=0
for f in "${OUTFILES[@]}"; do
    [[ -f ${f} ]] && EXISTING=$(( EXISTING + 1 ))
done
if (( EXISTING == ${#OUTFILES[@]} )); then
    echo "${OUTFILES[*]} already exist; skipping ${SAMPLE}"
    exit 0
elif (( EXISTING > 0 )); then
    echo "Only some of ${OUTFILES[*]} exist; remove them to trim ${SAMPLE} again"
    exit 1
fi
mkdir -p "${OUTDIR}"

# STEP 4: Trim adapters from the right (ktrim=r; with tbo also by pair overlap, and tpe
# trimming both reads of a pair to the same length), then trim both ends to quality TRIMQ,
# discarding reads (and their mates) left shorter than MINLEN
bbduk.sh -Xmx${XMX} "${IOARGS[@]}" ref="${ADAPTERS}" threads=${CPUS} \
    ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=${TRIMQ} minlength=${MINLEN} \
    stats="${OUTDIR}/${SAMPLE}.bbduk_stats.txt" overwrite=t

# STEP 5: Give the trimmed reads their final names once complete
for i in "${!TMPFILES[@]}"; do
    mv "${TMPFILES[${i}]}" "${OUTFILES[${i}]}"
done
