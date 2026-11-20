####

# Specify the conda environment containing smudgeplot and FastK
CONDAENV=@@condaEnv@@

# Specify the reads directory and the suffix identifying read files; for paired reads,
# the character before the suffix must be a 1 or 2 (e.g., sample_1.fq.gz and sample_2.fq.gz)
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@
PAIRED=@@paired@@ # true or false

# Specify the output directory; each sample's results are written to <outdir>/<sample>/
OUTDIR=@@outputDir@@

# Specify the count below which k-mers are considered erroneous (smudgeplot hetmers -L);
# if blank, it is estimated for each sample as half its k-mer coverage per genome copy,
# from the bases sequenced, the approximate (haploid) genome size in bp, and the highest
# expected ploidy (over-estimating the ploidy only lowers the threshold)
ERRORTHRESHOLD=@@errorThreshold@@
GENOMESIZE=@@genomeSize@@
PLOIDY=@@ploidy@@

# Specify whether to keep each sample's FastK k-mer table (needed to rerun smudgeplot
# hetmers with another threshold, but large); its k-mer histogram is always kept
KEEPTABLE=@@keepTable@@ # true or false

# Specify computational parameters
CPUS=@@cpus@@
FASTKMEM=@@fastkMem@@ # GB of memory for FastK's sorting

####

conda activate ${CONDAENV}
set -eo pipefail
shopt -s nullglob

KMER=31 # k-mer size counted by FastK
MINTHRESHOLD=4 # FastK's table only holds k-mers seen at least this many times
SAMPLEBYTES=10000000 # bytes read from the start of each read file to estimate its bases

# STEP 1: Find the read files for this array index (one sample per subjob)
if [[ ${PAIRED} == "true" ]]; then
    READFILES=( "${READSDIR}"/*1"${SUFFIX}" )
else
    READFILES=( "${READSDIR}"/*"${SUFFIX}" )
fi
INDEX=$(( ${PBS_ARRAY_INDEX:-1} - 1 )) # PBS_ARRAY_INDEX is unset if there is only one sample
if (( INDEX >= ${#READFILES[@]} )); then
    echo "Array index ${PBS_ARRAY_INDEX} exceeds the ${#READFILES[@]} read files found in ${READSDIR}"
    exit 1
fi
R1=${READFILES[${INDEX}]}
BASE=$(basename "${R1}")

# STEP 2: Locate the read pair (if paired) and name the sample
if [[ ${PAIRED} == "true" ]]; then
    PREFIX=${BASE%1"${SUFFIX}"}
    R2=${READSDIR}/${PREFIX}2${SUFFIX}
    if [[ ! -f ${R2} ]]; then
        echo "${R1} has no matching reverse read file (${R2})"
        exit 1
    fi
    READS=( "${R1}" "${R2}" )
    SAMPLE=$(sed -E 's/[._-]R?$//' <<< "${PREFIX}")
else
    READS=( "${R1}" )
    SAMPLE=${BASE%"${SUFFIX}"}
fi

# STEP 3: Skip samples which have already been analysed
SAMPLEDIR=${OUTDIR}/${SAMPLE}
if [[ -d ${SAMPLEDIR} ]]; then
    echo "${SAMPLEDIR} already exists; skipping ${SAMPLE}"
    exit 0
fi
INCOMPLETE=${OUTDIR}/${SAMPLE}.incomplete
rm -rf "${INCOMPLETE}" # left by a failed run
mkdir -p "${INCOMPLETE}/tmp"

# STEP 4: Estimate the bases sequenced, unless the error threshold was given. The reads in
# the first ${SAMPLEBYTES} bytes of each file give the bases per (compressed) byte, which
# is scaled up to the file's size; their mean length gives the fraction of each read's
# positions that start a k-mer
estimate_bases() {
    # Prints the number of reads and bases in the first SAMPLEBYTES bytes of a FASTQ file,
    # and the file's size
    local SIZE=$(stat -L -c %s "$1")
    if [[ $1 == *.gz ]]; then
        { head -c ${SAMPLEBYTES} "$1" | gzip -dc 2>/dev/null || true; } # a truncated stream fails
    else
        head -c ${SAMPLEBYTES} "$1"
    fi | awk -v size=${SIZE} 'NR % 4 == 2 { len = length($0) } NR % 4 == 0 { reads++; bases += len }
        END { printf "%d\t%.0f\t%.0f\n", reads, bases, size }'
}
if [[ -z ${ERRORTHRESHOLD} ]]; then
    read BASES DEPTH ERRORTHRESHOLD < <(
        for READFILE in "${READS[@]}"; do
            estimate_bases "${READFILE}"
        done | awk -v bytes=${SAMPLEBYTES} -v genome=${GENOMESIZE} -v ploidy=${PLOIDY} \
                   -v k=${KMER} -v minimum=${MINTHRESHOLD} '
            { reads += $1; sampled += $2; bases += $2 * ($3 > bytes ? $3 / bytes : 1) }
            END {
                if (reads == 0) { print "0 0 0"; exit }
                readLength = sampled / reads
                depth = bases / genome
                kmerDepth = depth * (readLength - k + 1) / readLength
                threshold = int(kmerDepth / ploidy / 2 + 0.5)
                printf "%.0f %.1f %d\n", bases, depth, (threshold < minimum ? minimum : threshold)
            }'
    )
    if (( BASES == 0 )); then
        echo "No reads were found at the start of ${READS[*]}"
        exit 1
    fi
    echo "${SAMPLE}: estimated ${BASES} bases (depth ${DEPTH}x of a ${GENOMESIZE} bp genome); using -L ${ERRORTHRESHOLD}"
else
    BASES=.
    DEPTH=.
fi

# STEP 5: Count k-mers with FastK
TABLE=${INCOMPLETE}/${SAMPLE}_table
FastK -t${MINTHRESHOLD} -k${KMER} -M${FASTKMEM} -T${CPUS} "${READS[@]}" -N"${TABLE}" -P"${INCOMPLETE}/tmp"

# STEP 6: Find heterozygous k-mer pairs, then infer the ploidy and plot the smudgeplot;
# smudgeplot >= 0.5 is run as 'smudgeplot', with earlier versions as 'smudgeplot.py'
if command -v smudgeplot > /dev/null; then
    SMUDGEPLOT=smudgeplot
else
    SMUDGEPLOT=smudgeplot.py
fi
PAIRS=${INCOMPLETE}/${SAMPLE}_kmerpairs
${SMUDGEPLOT} hetmers -L ${ERRORTHRESHOLD} -t ${CPUS} -o "${PAIRS}" -tmp "${INCOMPLETE}/tmp" "${TABLE}.ktab"
if [[ -f ${PAIRS}.smu ]]; then # smudgeplot >= 0.5
    SMU=${PAIRS}.smu
elif [[ -f ${PAIRS}_text.smu ]]; then # earlier versions
    SMU=${PAIRS}_text.smu
else
    echo "smudgeplot hetmers did not write ${PAIRS}.smu (or ${PAIRS}_text.smu)"
    exit 1
fi
${SMUDGEPLOT} all -t "${SAMPLE}" -o "${INCOMPLETE}/${SAMPLE}" "${SMU}"

# STEP 7: Summarise the largest smudge, whose number of k-mer copies (e.g., 3 for AAB or
# 2A1B) is the most likely ploidy; smudgeplot >= 0.5 reports each smudge's k-mer pairs,
# and earlier versions each smudge's fraction of them
REPORT=${INCOMPLETE}/${SAMPLE}.smudge_report.tsv
SIZES=${INCOMPLETE}/${SAMPLE}_smudge_sizes.txt
if [[ -f ${REPORT} ]]; then
    SMUDGES=$(awk -F '\t' 'NR == 1 { for (i = 4; i <= NF; i++) { name[i] = $i } }
        NR == 2 { for (i = 4; i <= NF; i++) { if ($i > 0) { print name[i] "\t" $i } } }' "${REPORT}")
elif [[ -f ${SIZES} ]]; then
    SMUDGES=$(awk -F '\t' '$1 != "" && $2 > 0 { print $1 "\t" $2 }' "${SIZES}")
else
    echo "smudgeplot all did not write ${REPORT} (or ${SIZES})"
    exit 1
fi
awk -F '\t' -v sample="${SAMPLE}" -v bases=${BASES} -v depth=${DEPTH} -v threshold=${ERRORTHRESHOLD} '
    { total += $2; if ($2 > best) { best = $2; smudge = $1 } }
    END {
        if (best == 0) { smudge = "."; fraction = "."; ploidy = "undetermined" }
        else {
            fraction = sprintf("%.4f", best / total)
            if (smudge ~ /^[AB]+$/) { ploidy = length(smudge) }
            else { n = split(smudge, counts, /[AB]/); ploidy = 0; for (i = 1; i <= n; i++) { ploidy += counts[i] } }
        }
        print "sample\tbases\tdepth\terror_threshold\tsmudge\tsmudge_fraction\tploidy"
        print sample "\t" bases "\t" depth "\t" threshold "\t" smudge "\t" fraction "\t" ploidy
    }' <<< "${SMUDGES}" > "${INCOMPLETE}/${SAMPLE}.ploidy.tsv"

# STEP 8: Tidy up, giving the results their final location once complete
rm -rf "${INCOMPLETE}/tmp"
if [[ ${KEEPTABLE} == "false" ]]; then
    rm -f "${TABLE}.ktab" "${INCOMPLETE}/.${SAMPLE}_table.ktab."*
fi
mv "${INCOMPLETE}" "${SAMPLEDIR}"
