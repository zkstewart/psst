####

# Specify the annotarium location
ANNOTARIUMDIR=@@annotariumDir@@

# Specify compleasm directory
COMPLEASMDIR=@@compleasmDir@@
COMPLEASMDB=${COMPLEASMDIR}/db

# Specify computational parameters
CPUS=@@cpus@@

# Specify lineage to search
LINEAGE=@@lineage@@

####

# STEP 0: Get the file name and prefix
INDIR=$(dirname ${PBS_O_WORKDIR})/fastq_corrected/output
FILEIN=$(ls ${INDIR}/*.fasta)
PREFIX=$(basename ${FILEIN} _corrected_reads.fasta)

# STEP 1: Run hifiasm
hifiasm --version > hifiasm.version # leave a record of the version used to create this assembly
hifiasm -o ${PREFIX}.asm -t ${CPUS} ${FILEIN}

# STEP 2: Generate the FASTA files
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.bp.hap1.p_ctg.gfa > ${PREFIX}.asm.bp.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.bp.hap2.p_ctg.gfa > ${PREFIX}.asm.bp.hap2.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.bp.p_ctg.gfa > ${PREFIX}.asm.bp.p_ctg.fasta

# STEP 3: Get the assembly statistics, replacing any from a previous run (the FASTA files
# are regenerated above); annotarium refuses to overwrite, so each is written under a
# temporary name first
fasta_stats() {
    rm -f "$2.incomplete"
    python "${ANNOTARIUMDIR}/annotarium.py" fasta stats -i "$1" -o "$2.incomplete" && mv -f "$2.incomplete" "$2"
}
fasta_stats ${PREFIX}.asm.bp.hap1.p_ctg.fasta ${PREFIX}.asm.bp.hap1.p_ctg.stats
fasta_stats ${PREFIX}.asm.bp.hap2.p_ctg.fasta ${PREFIX}.asm.bp.hap2.p_ctg.stats
fasta_stats ${PREFIX}.asm.bp.p_ctg.fasta ${PREFIX}.asm.bp.p_ctg.stats

# STEP 4: Run completeness scoring of assemblies
for assembly in ${PREFIX}.asm.bp.hap1.p_ctg.fasta ${PREFIX}.asm.bp.hap2.p_ctg.fasta ${PREFIX}.asm.bp.p_ctg.fasta; do
    mkdir -p ${assembly}.compleasm
    python ${COMPLEASMDIR}/compleasm.py run -a ${assembly} \
        -L ${COMPLEASMDB} -l ${LINEAGE} -t ${CPUS} \
        -o ${assembly}.compleasm;
done;
