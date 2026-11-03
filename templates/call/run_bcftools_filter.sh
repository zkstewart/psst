####

# Specify the VCF to filter and the filtered VCF to create
VCF=@@vcf@@
OUTVCF=@@outputFile@@

# Specify site filters, each of which can be disabled as noted
MAXMISSING=@@maxMissing@@ # maximum fraction of samples with a missing genotype; 1 to disable
MINQ=@@minQ@@ # minimum QUAL; 0 to disable
MAC=@@mac@@ # minimum minor allele count; 0 to disable
MAF=@@maf@@ # minimum minor allele frequency; empty to disable
PASSONLY=@@passOnly@@ # true or false; remove sites whose FILTER is not PASS (or '.')
BIALLELIC=@@biallelic@@ # true or false; remove multiallelic sites
NOINDELS=@@noIndels@@ # true or false; remove indels and complex variants

# Specify computational parameters
CPUS=@@cpus@@

####

set -eo pipefail

# STEP 1: Refuse to overwrite a previous output
if [[ -f ${OUTVCF} ]]; then
    echo "${OUTVCF} already exists; remove it to refilter"
    exit 1
fi

# STEP 2: Combine the enabled site filters; bcftools calculates F_MISSING, MAC and MAF from
# the genotypes, and a site without QUAL fails a QUAL filter
CLAUSES=()
awk -v x=${MAXMISSING} 'BEGIN{exit !(x < 1)}' && CLAUSES+=( "F_MISSING<=${MAXMISSING}" )
awk -v x=${MINQ} 'BEGIN{exit !(x > 0)}' && CLAUSES+=( "QUAL>=${MINQ}" )
(( MAC > 0 )) && CLAUSES+=( "MAC>=${MAC}" )
[[ -n ${MAF} ]] && CLAUSES+=( "MAF>=${MAF}" )
VIEWARGS=()
if (( ${#CLAUSES[@]} > 0 )); then
    printf -v EXPRESSION '%s && ' "${CLAUSES[@]}"
    VIEWARGS+=( -i "${EXPRESSION% && }" )
fi
[[ ${PASSONLY} == "true" ]] && VIEWARGS+=( -f PASS,. )

# STEP 3: Define the optional filters as pipeline stages; multiallelic sites split across
# records (e.g., by 'call giraffe') are joined so that they are removed too
keep_biallelic() {
    if [[ ${BIALLELIC} == "true" ]]; then
        bcftools norm -m +any -Ou | bcftools view -M 2 -Ou
    else
        cat
    fi
}
remove_indels() {
    if [[ ${NOINDELS} == "true" ]]; then
        bcftools view -V indels,other -Ou
    else
        cat
    fi
}

# STEP 4: Filter without intermediate files, then index the output and give it its final
# name once complete; a .tbi index cannot be made if a contig is longer than 2^29 bp, but
# the .csi index can
TMPVCF=${OUTVCF%.vcf.gz}.incomplete.vcf.gz
bcftools view "${VIEWARGS[@]}" -Ou "${VCF}" \
    | keep_biallelic \
    | remove_indels \
    | bcftools view --threads ${CPUS} -Oz -o "${TMPVCF}"
bcftools index --threads ${CPUS} -o "${OUTVCF}.csi" "${TMPVCF}"
bcftools index --threads ${CPUS} -t -o "${OUTVCF}.tbi" "${TMPVCF}" \
    || echo "WARNING: a .tbi index could not be made; use the .csi index instead"
mv "${TMPVCF}" "${OUTVCF}"

# STEP 5: Summarise the filtered variants
bcftools stats -s - "${OUTVCF}" > "${OUTVCF%.vcf.gz}.stats.txt"
