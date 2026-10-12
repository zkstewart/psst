####

# Specify the location of the Various_scripts repository
VARSCRIPTDIR=@@varScriptDir@@

# Specify the location of the diploid reference genome haplotypes
REFHAP1=@@refHap1@@
REFHAP2=@@refHap2@@

# Specify how many CPUs to use
CPUS=@@cpus@@

# Specify behavioural params
OUTDIR=@@outputDir@@
PRESET=@@preset@@

####

# STEP 0: Specify the location of the input files
ASSEMBLYDIR=$(realpath ../assembly)
HAP1=$(ls ${ASSEMBLYDIR}/*.hap1.p_ctg.fasta)
HAP2=$(ls ${ASSEMBLYDIR}/*.hap2.p_ctg.fasta)

# STEP 1: Run the pipeline
python ${VARSCRIPTDIR}/assembly/haplotype_scaffolding_pipeline.py -i ${HAP1} ${HAP2} \
    -r ${REFHAP1} ${REFHAP2} \
    -o ${OUTDIR} --threads ${CPUS} --preset ${PRESET}
