####

# Specify the program locations
VARSCRIPTDIR=@@varScriptDir@@
MMS2DIR=@@mmseqsDir@@

# Specify the location of the input FASTA file and its molecule type
FASTA=@@fasta@@
MOLECULE=@@molecule@@ # protein or nucleotide

# Specify the output file name
OUTPUTFILE=@@outputFile@@

# Specify the MMseqs2 tmp directory
MMS2TMP=@@tmpDir@@

# Specify clustering type
PROGRAM=@@program@@ # mmseqs-cascade or mmseqs-linclust

# Specify program parameters
CPUS=@@cpus@@
IDENTITY=@@identity@@
COVERAGE=@@coverage@@
EVALUE=@@evalue@@
MODE=@@mode@@ # set-cover, connected-component, or greedy
SENSITIVITY=@@sensitivity@@ # only relevant for cascade
STEPS=@@steps@@ # only relevant for cascade

####

# Run wrapper script
python ${VARSCRIPTDIR}/run_clustering.py -i ${FASTA} -m ${MOLECULE} -o ${OUTPUTFILE} -t ${CPUS} \
    -p ${PROGRAM} --mmseqs ${MMS2DIR} --tmpDir ${MMS2TMP} \
    --identity ${IDENTITY} --coverage ${COVERAGE} --evalue ${EVALUE} \
    --mode ${MODE} --sensitivity ${SENSITIVITY} --steps ${STEPS}
