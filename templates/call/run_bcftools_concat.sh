####

# Specify the output directory and the prefix for output files, as given to the calling job
OUTDIR=@@outputDir@@
PREFIX=@@prefix@@

# Specify the number of chunks the genome was divided into by the calling job
NCHUNKS=@@chunks@@

# Specify computational parameters
CPUS=@@cpus@@

####

set -eo pipefail

# STEP 1: Refuse to overwrite a previous output
OUTVCF=${OUTDIR}/${PREFIX}.vcf.gz
if [[ -f ${OUTVCF} ]]; then
    echo "${OUTVCF} already exists; remove it to rebuild it"
    exit 1
fi

# STEP 2: Check that every chunk was called, collecting them in genome order
CHUNKDIR=${OUTDIR}/${PREFIX}_chunks
CHUNKVCFS=()
for (( CHUNK=1; CHUNK<=NCHUNKS; CHUNK++ )); do
    if [[ -f ${CHUNKDIR}/chunk_${CHUNK}.vcf.gz ]]; then
        CHUNKVCFS+=( "${CHUNKDIR}/chunk_${CHUNK}.vcf.gz" )
    elif [[ ! -f ${CHUNKDIR}/chunk_${CHUNK}.empty ]]; then
        echo "Chunk ${CHUNK} has not been called; rerun the calling job (completed chunks are skipped)"
        exit 1
    fi
done

# STEP 3: Concatenate chunks without recompressing them, then index the output; a .tbi
# index cannot be made if a contig is longer than 2^29 bp, but the .csi index can
TMPVCF=${OUTDIR}/${PREFIX}.incomplete.vcf.gz
bcftools concat --naive -o "${TMPVCF}" "${CHUNKVCFS[@]}"
bcftools index --threads ${CPUS} -o "${OUTVCF}.csi" "${TMPVCF}"
bcftools index --threads ${CPUS} -t -o "${OUTVCF}.tbi" "${TMPVCF}" \
    || echo "WARNING: a .tbi index could not be made; use the .csi index instead"
mv "${TMPVCF}" "${OUTVCF}"

# STEP 4: Clean up the chunk outputs, which are now contained in the output
rm -r "${CHUNKDIR}"
