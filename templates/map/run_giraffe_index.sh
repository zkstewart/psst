####

# Specify the vg executable location
VGDIR=@@vgDir@@

# Specify the pangenome graph to map to (e.g., cactus-pangenome's .d2.gbz) and the VCF of
# its variants (e.g., cactus-pangenome's vcfwave-normalised .wave.vcf.gz)
GBZ=@@gbz@@
VCF=@@vcf@@

# Specify the sample whose paths form the linear reference; empty to detect the graph's
# only reference sample
REFSAMPLE=@@refSample@@

# Specify the index directory to create
INDEX=@@indexDir@@

# Specify computational parameters
CPUS=@@cpus@@
SORTMEM=@@sortMem@@ # memory for bcftools sort

####

set -eo pipefail

# STEP 1: Refuse to overwrite an existing index
if [[ -e ${INDEX} ]]; then
    echo "${INDEX} already exists; remove it to rebuild the index"
    exit 1
fi
TMPINDEX=${INDEX}.incomplete
rm -rf "${TMPINDEX}" # remove any output left by a previous failed run
mkdir -p "${TMPINDEX}/tmp"

# STEP 2: Copy the graph into the index, recording its origin and the vg version
cp "${GBZ}" "${TMPINDEX}/pangenome.gbz"
{
    echo "Pangenome graph: $(realpath "${GBZ}")"
    echo "Pangenome variants: $(realpath "${VCF}")"
    ${VGDIR}/vg version | head -n 1
} > "${TMPINDEX}/index_info.txt"

# STEP 3: Find the reference sample if not given
if [[ -z ${REFSAMPLE} ]]; then
    REFSAMPLES=( $(${VGDIR}/vg paths -x "${TMPINDEX}/pangenome.gbz" -M -R | awk 'NR>1{print $3}' | sort -u) )
    if (( ${#REFSAMPLES[@]} != 1 )); then
        echo "The graph has ${#REFSAMPLES[@]} reference samples (${REFSAMPLES[*]}); specify one as REFSAMPLE"
        exit 1
    fi
    REFSAMPLE=${REFSAMPLES[0]}
fi
echo "Reference sample: ${REFSAMPLE}" >> "${TMPINDEX}/index_info.txt"

# STEP 4: Write the reference sample's paths as the linear reference, which mapped reads
# are projected onto and variants are called against
${VGDIR}/vg paths -x "${TMPINDEX}/pangenome.gbz" -F -S "${REFSAMPLE}" > "${TMPINDEX}/reference.fa"
if [[ ! -s ${TMPINDEX}/reference.fa ]]; then
    echo "The graph has no paths for sample '${REFSAMPLE}'"
    exit 1
fi
samtools faidx "${TMPINDEX}/reference.fa"

# STEP 5: Build the distance and minimizer indexes for vg giraffe
${VGDIR}/vg index -t ${CPUS} -j "${TMPINDEX}/pangenome.dist" "${TMPINDEX}/pangenome.gbz"
${VGDIR}/vg minimizer -t ${CPUS} -d "${TMPINDEX}/pangenome.dist" -z "${TMPINDEX}/pangenome.zipcodes" \
    -o "${TMPINDEX}/pangenome.min" "${TMPINDEX}/pangenome.gbz"

# STEP 6: Name the variants' contigs like the reference paths (e.g., 'chr1' becomes
# 'reference#0#chr1'), then split multiallelic records and left-align them
cut -f1 "${TMPINDEX}/reference.fa.fai" | awk '{n=$1; sub(/^[^#]*#[^#]*#/, "", n); print n"\t"$1}' \
    > "${TMPINDEX}/tmp/rename.txt"
bcftools annotate --rename-chrs "${TMPINDEX}/tmp/rename.txt" -Ou "${VCF}" \
    | bcftools norm -f "${TMPINDEX}/reference.fa" -c e -m -any -Ou \
    | bcftools sort -T "${TMPINDEX}/tmp/sort" -m ${SORTMEM} -Ob -o "${TMPINDEX}/tmp/split.bcf"

# STEP 7: Decompose MNPs into SNPs, matching how mapped reads show them, so that samples
# can be genotyped at every site; complex alleles (e.g., ATTCC>TCG) are kept whole, since
# reads represent them inconsistently and they can't be genotyped reliably
bcftools view -i 'TYPE="other"' -Ob -o "${TMPINDEX}/tmp/complex.bcf" --write-index "${TMPINDEX}/tmp/split.bcf"
bcftools view -e 'TYPE="other"' -Ou "${TMPINDEX}/tmp/split.bcf" \
    | bcftools norm -f "${TMPINDEX}/reference.fa" -a --atom-overlaps . -Ou \
    | bcftools view -e 'ALT="*"' -Ou \
    | bcftools sort -T "${TMPINDEX}/tmp/sort" -m ${SORTMEM} -Ob -o "${TMPINDEX}/tmp/simple.bcf" --write-index
bcftools concat -a -Oz -o "${TMPINDEX}/sites.vcf.gz" --write-index=tbi \
    "${TMPINDEX}/tmp/simple.bcf" "${TMPINDEX}/tmp/complex.bcf"

# STEP 8: List the alleles to genotype samples at, one line per position and variant type
# (as bcftools mpileup reports SNPs and indels separately)
bcftools view -e 'TYPE="other"' -Ou "${TMPINDEX}/sites.vcf.gz" \
    | bcftools norm -m +both -Ou \
    | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\n' \
    | bgzip > "${TMPINDEX}/targets.tsv.gz"
tabix -s1 -b2 -e2 "${TMPINDEX}/targets.tsv.gz"

# STEP 9: Give the index its final name once complete
rm -r "${TMPINDEX}/tmp"
mv "${TMPINDEX}" "${INDEX}"
