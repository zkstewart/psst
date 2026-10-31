####

# Specify the location of the Various_scripts repository
VARSCRIPTDIR=@@varScriptDir@@

# Specify the location of the reference genome
REFGENOME=@@refGenome@@

# Specify how many CPUs to use
CPUS=@@cpus@@

# Specify behavioural params
OUTDIR=@@outputDir@@
PRESET=@@preset@@

####

# STEP 0: Specify the location of the input files
ASSEMBLYDIR=$(realpath ../scaffolding/output)
HAP1=${ASSEMBLYDIR}/hap1.fasta
HAP2=${ASSEMBLYDIR}/hap2.fasta

# STEP 1: Run the pipeline
python ${VARSCRIPTDIR}/assembly/revcomp_to_match_reference.py -i ${HAP1} ${HAP2} \
    -r ${REFGENOME} \
    -o ${OUTDIR} \
    --threads ${CPUS} --preset ${PRESET}
