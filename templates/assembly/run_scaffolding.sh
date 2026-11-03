####

# Specify the locations of the annotarium repository, minimap2 and paf2dotplot; RagTag
# (ragtag.py) must be in the PATH
ANNOTARIUMDIR=@@annotariumDir@@
MINIMAP2DIR=@@minimap2Dir@@
P2DPDIR=@@paf2dotplotDir@@

# Specify the contig FASTA file(s) to phase, and the reference genome to scaffold against
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

# STEP 2: Scaffold each partition against the reference with RagTag, writing
# ${OUTDIR}/partition<N>.fasta and a dotplot of it against the reference; RagTag runs
# finished by a previous run of this job are reused
python "${ANNOTARIUMDIR}/annotarium.py" assembly scaffold -i "${PHASEDIR}/phasing.tsv" \
    -f "${CONTIGS[@]}" -r "${REFGENOME}" -o "${OUTDIR}" \
    --threads ${CPUS} --minimap2 "${MINIMAP2DIR}/minimap2" \
    --paf2dotplot "${P2DPDIR}/paf2dotplot.r" --preset ${PRESET}
