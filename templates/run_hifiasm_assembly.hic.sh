####

# Specify the Genome_analysis_scripts location
GENSCRIPTDIR=@@genScriptDir@@

# Specify compleasm directory
COMPLEASMDIR=@@compleasmDir@@
COMPLEASMDB=${COMPLEASMDIR}/db

# Specify computational and behavioural parameters
CPUS=@@cpus@@
PLOIDY=@@ploidy@@ # if ploidy > 2, you need to manually add some changes below

# Specify lineage to search
LINEAGE=@@lineage@@

# Specify the Hi-C read files
HIC1=@@hic1@@
HIC2=@@hic2@@

####

# STEP 0: Get the file name and prefix
INDIR=$(dirname ${PBS_O_WORKDIR})/fastq_corrected/output
FILEIN=$(ls ${INDIR}/*.fasta)
PREFIX=$(basename ${FILEIN} _corrected_reads.fasta)

# STEP 1: Run hifiasm
hifiasm --version > hifiasm.version # leave a record of the version used to create this assembly
hifiasm -o ${PREFIX}.asm \
    --n-hap ${PLOIDY} \
    -t ${CPUS} \
    --h1 ${HIC1} --h2 ${HIC2} \
    ${FILEIN}

# STEP 2: Generate the FASTA files
## Add hap3... if needed
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.hic.hap1.p_ctg.gfa > ${PREFIX}.asm.hic.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.hic.hap2.p_ctg.gfa > ${PREFIX}.asm.hic.hap2.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' ${PREFIX}.asm.hic.p_ctg.gfa > ${PREFIX}.asm.hic.p_ctg.fasta

# STEP 3: Get the assembly statistics
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.asm.hic.hap1.p_ctg.fasta -o ${PREFIX}.asm.hic.hap1.p_ctg.stats
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.asm.hic.hap2.p_ctg.fasta -o ${PREFIX}.asm.hic.hap2.p_ctg.stats
python ${GENSCRIPTDIR}/genome_stats.py -i ${PREFIX}.asm.hic.p_ctg.fasta -o ${PREFIX}.asm.hic.p_ctg.stats

# STEP 4: Run completeness scoring of assemblies
## Add hap3... if needed
for assembly in ${PREFIX}.asm.hic.hap1.p_ctg.fasta ${PREFIX}.asm.hic.hap2.p_ctg.fasta ${PREFIX}.asm.hic.p_ctg.fasta; do
    mkdir -p ${assembly}.compleasm
    python ${COMPLEASMDIR}/compleasm.py run -a ${assembly} \
        -L ${COMPLEASMDB} -l ${LINEAGE} -t ${CPUS} \
        -o ${assembly}.compleasm;
done;
