####

# Specify the STAR module to load
STARMODULE=@@starModule@@

# Specify the reference genome FASTA (uncompressed)
GENOME=@@genome@@

# Specify the gene annotation GTF (uncompressed)
## Convert GFF3 to GTF like: agat_convert_sp_gff2gtf.pl --gff input.gff3 -o output.gtf
GTF=@@gtf@@

# Specify the index directory to create
INDEX=@@indexDir@@

# Specify index parameters
OVERHANG=@@overhang@@ # ideally read length - 1

# Specify computational parameters
CPUS=@@cpus@@
RAMLIMIT=@@ramLimit@@ # in bytes; STAR's default limit is 31 GB

####

module load ${STARMODULE}
set -eo pipefail

# STEP 1: Refuse to overwrite an existing index
if [[ -e ${INDEX} ]]; then
    echo "${INDEX} already exists; remove it to rebuild the index"
    exit 1
fi

# STEP 2: Scale the suffix array index to the genome size, as recommended by the
# STAR manual: min(14, log2(genome length)/2 - 1)
GENOMELEN=$(grep -v '^>' "${GENOME}" | tr -d '\n\r' | wc -c)
SAINDEX=$(awk -v len=${GENOMELEN} 'BEGIN{v=int(log(len)/log(2)/2 - 1); print (v < 14) ? v : 14}')
echo "Genome length is ${GENOMELEN} bp; using --genomeSAindexNbases ${SAINDEX}"

# STEP 3: Build the index, giving it its final name once complete
TMPINDEX=${INDEX}.incomplete
rm -rf "${TMPINDEX}" # remove any output left by a previous failed run
mkdir -p "${TMPINDEX}"
STAR --runThreadN ${CPUS} \
     --runMode genomeGenerate \
     --genomeDir "${TMPINDEX}" \
     --genomeFastaFiles "${GENOME}" \
     --sjdbGTFfile "${GTF}" \
     --sjdbOverhang ${OVERHANG} \
     --genomeSAindexNbases ${SAINDEX} \
     --limitGenomeGenerateRAM ${RAMLIMIT} \
     --outFileNamePrefix "${TMPINDEX}/"
mv "${TMPINDEX}" "${INDEX}"
