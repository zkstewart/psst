####

# Specify the salmon executable location
SALMONDIR=@@salmonDir@@

# Specify the transcripts FASTA to index
TRANSCRIPTS=@@transcripts@@

# Optionally, specify the genome FASTA to use as decoys (leave empty for no decoys)
GENOME=@@genome@@

# Specify the index directory to create
INDEX=@@indexDir@@

# Specify index parameters
KMER=@@kmer@@

# Specify computational parameters
CPUS=@@cpus@@
RAMLIMIT=@@ramLimit@@ # in GiB; salmon otherwise sizes this from the whole node's memory

####

set -eo pipefail

# STEP 1: Refuse to overwrite an existing index
if [[ -e ${INDEX} ]]; then
    echo "${INDEX} already exists; remove it to rebuild the index"
    exit 1
fi

# STEP 2: Prepare decoys (if any) by appending the genome to the transcripts
TMPPREFIX=${INDEX}.incomplete
if [[ -n ${GENOME} ]]; then
    grep '^>' <(zcat -f "${GENOME}") | cut -d ' ' -f 1 | sed 's/^>//' > ${TMPPREFIX}.decoys.txt
    zcat -f "${TRANSCRIPTS}" "${GENOME}" > ${TMPPREFIX}.gentrome.fasta
    DECOYARGS=( --decoys ${TMPPREFIX}.decoys.txt )
    TARGETS=${TMPPREFIX}.gentrome.fasta
else
    DECOYARGS=()
    TARGETS=${TRANSCRIPTS}
fi

# STEP 3: Build the index, giving it its final name once complete
${SALMONDIR}/salmon index --threads ${CPUS} -k ${KMER} --ramLimit ${RAMLIMIT} \
    --transcripts "${TARGETS}" "${DECOYARGS[@]}" \
    --index ${TMPPREFIX}
mv ${TMPPREFIX} ${INDEX}
rm -f ${TMPPREFIX}.decoys.txt ${TMPPREFIX}.gentrome.fasta
