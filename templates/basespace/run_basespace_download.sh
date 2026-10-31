####

# Specify basespace executable location
BSDIR=@@bsDir@@

# Specify project ID
## Find this with: ${BSDIR}/bs projects list
PROJECTID=@@projectID@@

# Specify output directory
OUTDIR=@@outputDir@@

####

# Run basespace download
${BSDIR}/bs download project -i ${PROJECTID} -o ${OUTDIR}
