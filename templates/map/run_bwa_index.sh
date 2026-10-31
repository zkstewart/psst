####

# Specify the reference genome FASTA to index; index files are written alongside it
GENOME=@@genome@@

####

set -e

# STEP 1: Refuse to overwrite an existing index
for SUFFIX in .amb .ann .bwt .pac .sa; do
    if [[ -f ${GENOME}${SUFFIX} ]]; then
        echo "${GENOME}${SUFFIX} already exists; remove the existing index files to rebuild it"
        exit 1
    fi
done

# STEP 2: Index the genome
bwa index "${GENOME}"
