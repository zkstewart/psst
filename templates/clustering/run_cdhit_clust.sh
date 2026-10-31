####

# Specify the program locations
VARSCRIPTDIR=@@varScriptDir@@
CDHITDIR=@@cdhitDir@@

# Specify the location of the input FASTA file and its molecule type
FASTA=@@fasta@@
MOLECULE=@@molecule@@ # protein or nucleotide

# Specify the output file name
OUTPUTFILE=@@outputFile@@

# Specify program parameters
CPUS=@@cpus@@
IDENTITY=@@identity@@
MEM=@@cdhitMem@@ # in MB
ALIGN=@@align@@ # local or global
SHORTERCOV=@@shorterCov@@
LONGERCOV=@@longerCov@@

####

# Run wrapper script
python ${VARSCRIPTDIR}/run_clustering.py -i ${FASTA} -m ${MOLECULE} -o ${OUTPUTFILE} -t ${CPUS} \
    -p cd-hit --cdhit ${CDHITDIR} \
    --identity ${IDENTITY} --mem ${MEM} --align ${ALIGN} \
    --shorter_coverage ${SHORTERCOV} --longer_coverage ${LONGERCOV}
