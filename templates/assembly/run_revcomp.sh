####

# Specify the locations of the annotarium repository and minimap2
ANNOTARIUMDIR=@@annotariumDir@@
MINIMAP2DIR=@@minimap2Dir@@

# Specify the contig FASTA file(s) to phase, and the reference genome to match
CONTIGS=( @@contigFiles@@ )
REFGENOME=@@refGenome@@

# Specify an existing 'assembly phase' output directory to reuse; empty to phase the
# contigs into ${OUTDIR}/phase
PHASEDIR=@@phaseDir@@

# Specify phasing parameters
PLOIDY=@@ploidy@@
PRESET=@@preset@@

# Specify the output directory
OUTDIR=@@outputDir@@

# Specify computational parameters
CPUS=@@cpus@@

####

set -eo pipefail

# STEP 1: Phase contigs into haplotype partitions, unless an earlier phasing was given or
# a previous run of this job completed it ('assembly phase' writes coverage.tsv last)
mkdir -p "${OUTDIR}"
if [[ -z ${PHASEDIR} ]]; then
    PHASEDIR=${OUTDIR}/phase
    if [[ -f ${PHASEDIR}/coverage.tsv ]]; then
        echo "Reusing the completed phasing in ${PHASEDIR}"
    else
        python "${ANNOTARIUMDIR}/annotarium.py" assembly phase -i "${CONTIGS[@]}" -r "${REFGENOME}" \
            -o "${PHASEDIR}" --ploidy ${PLOIDY} --preset ${PRESET} \
            --threads ${CPUS} --minimap2 "${MINIMAP2DIR}/minimap2"
    fi
fi

# STEP 2: Find the partitions that contigs were assigned to
PARTITIONS=$(awk -F '\t' '$4 == "assigned" { print $3 }' "${PHASEDIR}/phasing.tsv" | sort -nu)
if [[ -z ${PARTITIONS} ]]; then
    echo "No contigs were assigned to a partition in ${PHASEDIR}/phasing.tsv"
    exit 1
fi

# STEP 3: Pair each partition's contigs 1:1 with the reference contigs they best match,
# renamed and reverse complemented to match them, in reference order; partitions finished
# by a previous run of this job are skipped
for PARTITION in ${PARTITIONS}; do
    OUTFILE=${OUTDIR}/partition${PARTITION}.fasta
    if [[ -f ${OUTFILE} ]]; then
        echo "Skipping partition ${PARTITION}; ${OUTFILE} already exists"
        continue
    fi
    TMPFILE=${OUTFILE%.fasta}.incomplete.fasta
    rm -f "${TMPFILE}"
    python "${ANNOTARIUMDIR}/annotarium.py" assembly homology_pair -i "${PHASEDIR}/phasing.tsv" \
        -f "${CONTIGS[@]}" -r "${REFGENOME}" --partition ${PARTITION} -o "${TMPFILE}"
    mv "${TMPFILE}" "${OUTFILE}"
done
