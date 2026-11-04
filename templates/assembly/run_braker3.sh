####

# Specify the BRAKER3 Singularity image
BRAKERSIF=@@brakerSif@@

# Specify the softmasked genome FASTA file, and the protein evidence FASTA file
GENOME=@@genome@@
PROTEINS=@@proteins@@

# Specify RNAseq evidence as a directory of paired read files and the suffix identifying
# them (the character before it must be a 1 or 2, e.g., sample_1.fq.gz and sample_2.fq.gz),
# and/or text files listing SRA accessions (one per line) for BRAKER3 to download; leave
# READSDIR empty or SRAFILES=( ) to use only the other
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
SRAFILES=( @@sraFiles@@ )

# Specify the species name to train AUGUSTUS as, which also prefixes output files
PREFIX=@@prefix@@

# Specify the BUSCO lineage BRAKER3 uses to assess and improve its predictions
LINEAGE=@@lineage@@

# Specify the output directory
OUTDIR=@@outputDir@@

# Specify computational parameters
CPUS=@@cpus@@

####

set -eo pipefail
shopt -s nullglob

# STEP 1: Refuse to repeat a BRAKER3 run, or to resume one, which BRAKER3 cannot do
mkdir -p "${OUTDIR}"
OUTDIR=$(cd "${OUTDIR}" && pwd -P) # absolute, as paths are given to the container
BRAKERDIR=${OUTDIR}/braker
CONFIGDIR=${OUTDIR}/augustus_config
if [[ -f ${OUTDIR}/${PREFIX}.braker.gtf ]]; then
    echo "${OUTDIR} already holds BRAKER3 results for ${PREFIX}"
    exit 1
fi
for DIR in "${BRAKERDIR}" "${CONFIGDIR}"; do
    if [[ -e ${DIR} ]]; then
        echo "${DIR} remains from a previous run, which BRAKER3 cannot resume; remove ${BRAKERDIR} and ${CONFIGDIR} to start again"
        exit 1
    fi
done
GENOME=$(readlink -f "${GENOME}")
PROTEINS=$(readlink -f "${PROTEINS}")
BINDS=( "${OUTDIR}" "$(dirname "${GENOME}")" "$(dirname "${PROTEINS}")" )

# STEP 2: Link each pair of local read files into one directory as <ID>_1.fq.gz and
# <ID>_2.fq.gz (or .fq), the names BRAKER3 looks for, binding the files' real locations
# into the container
RNAIDS=()
LINKDIR=${OUTDIR}/rnaseq_reads
rm -rf "${LINKDIR}" # holds only links made by a previous run of this job
mkdir "${LINKDIR}"
if [[ -n ${READSDIR} ]]; then
    [[ ${SUFFIX} == *.gz ]] && EXT=fq.gz || EXT=fq
    R1FILES=( "${READSDIR}"/*1"${SUFFIX}" )
    if (( ${#R1FILES[@]} == 0 )); then
        echo "No read files ending with 1${SUFFIX} were found in ${READSDIR}"
        exit 1
    fi
    for R1 in "${R1FILES[@]}"; do
        PAIRPREFIX=$(basename "${R1%1"${SUFFIX}"}")
        R2=${READSDIR}/${PAIRPREFIX}2${SUFFIX}
        if [[ ! -f ${R2} ]]; then
            echo "${R1} has no matching reverse read file (${R2})"
            exit 1
        fi
        ID=$(sed -E 's/[._-]R?$//' <<< "${PAIRPREFIX}")
        if [[ ! ${ID} =~ ^[A-Za-z0-9._-]+$ ]]; then
            echo "Sample name '${ID}' of ${R1} may only contain letters, numbers, '.', '_' and '-'"
            exit 1
        fi
        R1=$(readlink -f "${R1}")
        R2=$(readlink -f "${R2}")
        ln -s "${R1}" "${LINKDIR}/${ID}_1.${EXT}"
        ln -s "${R2}" "${LINKDIR}/${ID}_2.${EXT}"
        RNAIDS+=( "${ID}" )
        BINDS+=( "$(dirname "${R1}")" "$(dirname "${R2}")" )
    done
fi

# STEP 3: Combine the SRA accessions of all files, without blank lines or repeats
SRALIST=${OUTDIR}/sra_accessions.txt
rm -f "${SRALIST}"
if (( ${#SRAFILES[@]} > 0 )); then
    cat "${SRAFILES[@]}" | tr -d '\r' | awk 'NF && !seen[$1]++ { print $1 }' > "${SRALIST}"
    readarray -t SRAIDS < "${SRALIST}"
    for ID in "${SRAIDS[@]}"; do
        if [[ ! ${ID} =~ ^[SED]RR[0-9]+$ ]]; then
            echo "'${ID}' in the SRA accessions files is not an SRA run accession (e.g., SRR1234567)"
            exit 1
        fi
    done
    RNAIDS+=( "${SRAIDS[@]}" )
fi

# STEP 4: Check that BRAKER3 can tell the RNAseq sets apart; it finds each set's files by
# matching its ID to the end of file names, so an ID matching the end of another ID
# would be given that set's files
if (( ${#RNAIDS[@]} == 0 )); then
    echo "No RNAseq sets were given; specify READSDIR and/or SRAFILES"
    exit 1
fi
for i in "${!RNAIDS[@]}"; do
    for j in "${!RNAIDS[@]}"; do
        (( i == j )) && continue
        if [[ ${RNAIDS[$i]} == "${RNAIDS[$j]}" ]]; then
            echo "RNAseq set '${RNAIDS[$i]}' is given more than once (as a local sample name or SRA accession)"
            exit 1
        fi
        if [[ ${RNAIDS[$j]} =~ ${RNAIDS[$i]}$ ]]; then
            echo "RNAseq set '${RNAIDS[$i]}' matches the end of '${RNAIDS[$j]}', so BRAKER3 would confuse their files"
            exit 1
        fi
    done
done
RNAIDLIST=$(IFS=,; echo "${RNAIDS[*]}")
echo "RNAseq sets: ${RNAIDLIST}"

# STEP 5: Copy AUGUSTUS's configuration out of the container, where BRAKER3 can write the
# parameters it trains
BINDLIST=$(printf '%s\n' "${BINDS[@]}" | sort -u | paste -sd, -)
singularity exec --bind "${BINDLIST}" "${BRAKERSIF}" cp -R /opt/Augustus/config "${CONFIGDIR}"

# STEP 6: Give compleasm (run by BRAKER3 in its working directory) BUSCO's list of files
# without placement files, which it fails to parse in newer BUSCO releases (e.g.,
# 'odb12.2' names, since May 2026); they are only used to choose a lineage automatically,
# not with --busco_lineage. compleasm still checks the lineage it downloads against this list
mkdir -p "${BRAKERDIR}/mb_downloads"
curl -fsSL https://busco-data.ezlab.org/v5/data/file_versions.tsv \
    | awk -F '\t' '$5 != "placement_files"' > "${BRAKERDIR}/mb_downloads/file_versions.tsv"
touch "${BRAKERDIR}/mb_downloads/file_versions.tsv.done"

# STEP 7: Run BRAKER3
BRAKERARGS=( --genome="${GENOME}" --prot_seq="${PROTEINS}" --species="${PREFIX}" \
    --rnaseq_sets_ids="${RNAIDLIST}" --busco_lineage="${LINEAGE}" \
    --workingdir="${BRAKERDIR}" --AUGUSTUS_CONFIG_PATH="${CONFIGDIR}" \
    --threads=${CPUS} --gff3 --nocleanup )
[[ -n ${READSDIR} ]] && BRAKERARGS+=( --rnaseq_sets_dirs="${LINKDIR}" )
singularity exec --bind "${BINDLIST}" "${BRAKERSIF}" braker.pl "${BRAKERARGS[@]}"

# STEP 8: Copy the predictions out of BRAKER3's working directory under the species name;
# the GTF is moved into place last, as it marks the job as complete
for TYPE in gff3 aa codingseq gtf; do
    cp "${BRAKERDIR}/braker.${TYPE}" "${OUTDIR}/${PREFIX}.braker.${TYPE}.incomplete"
done
for TYPE in gff3 aa codingseq gtf; do
    mv "${OUTDIR}/${PREFIX}.braker.${TYPE}.incomplete" "${OUTDIR}/${PREFIX}.braker.${TYPE}"
done
