####

# Specify the conda environment containing pairtools, and the 3D-DNA location
CONDAENV=@@condaEnv@@
THREEDDNADIR=@@threeddnaDir@@

# Specify the reference genome FASTA, which must be indexed with 'bwa index'
GENOME=@@genome@@

# Specify the Hi-C reads directory and the suffix identifying read files; the character
# before the suffix must be a 1 or 2 (e.g., hic_1.fq.gz and hic_2.fq.gz); all read pairs
# in the directory are mapped together as one library
READSDIR=@@readsDir@@
SUFFIX=@@readsSuffix@@

# Specify the output directory and the prefix for output files
OUTDIR=@@outputDir@@
PREFIX=@@prefix@@

# Specify the minimum MAPQ of both reads for a contact to be shown in the .hic file
MAPQ=@@mapq@@

# Specify computational parameters
CPUS=@@cpus@@
SORTMEM=@@sortMem@@ # memory for pairtools sort

####

conda activate ${CONDAENV}
set -eo pipefail
shopt -s nullglob

# STEP 1: Find all read pairs
R1FILES=( "${READSDIR}"/*1"${SUFFIX}" )
if (( ${#R1FILES[@]} == 0 )); then
    echo "No read files ending with 1${SUFFIX} were found in ${READSDIR}"
    exit 1
fi
R2FILES=()
for R1 in "${R1FILES[@]}"; do
    BASE=$(basename "${R1}")
    R2=${READSDIR}/${BASE%1"${SUFFIX}"}2${SUFFIX}
    if [[ ! -f ${R2} ]]; then
        echo "${R1} has no matching reverse read file (${R2})"
        exit 1
    fi
    R2FILES+=( "${R2}" )
done

# STEP 2: Refuse to overwrite a previous .hic file
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"
[[ ${MAPQ} -eq 1 ]] && HICFILE=${PREFIX}.hic || HICFILE=${PREFIX}_${MAPQ}.hic # 3D-DNA's naming
if [[ -f ${HICFILE} ]]; then
    echo "${OUTDIR}/${HICFILE} already exists; remove it to rebuild it"
    exit 1
fi

# STEP 3: Describe the reference contigs for pairtools and 3D-DNA
samtools faidx "${GENOME}" --fai-idx ${PREFIX}.fai
cut -f1,2 ${PREFIX}.fai > ${PREFIX}.chrom.sizes
{ awk '{print ">"$1" "NR" "$2}' ${PREFIX}.fai; awk '{print NR}' ${PREFIX}.fai; } > ${PREFIX}.assembly

# STEP 4: Map reads, then parse, sort and deduplicate contacts; MAPQ 0 contacts are kept
# with their MAPQ values so that the .hic file can apply (and later change) the threshold
if [[ ! -f ${PREFIX}.pairs.gz ]]; then
    mkdir -p ${PREFIX}.sort_tmp
    bwa mem -5SP -t ${CPUS} "${GENOME}" <(cat "${R1FILES[@]}") <(cat "${R2FILES[@]}") \
        | pairtools parse --chroms-path ${PREFIX}.chrom.sizes --min-mapq 0 --walks-policy 5unique \
            --add-columns mapq --drop-sam --drop-seq --nproc-in ${CPUS} --nproc-out ${CPUS} \
        | pairtools sort --nproc ${CPUS} --memory ${SORTMEM} --tmpdir ${PREFIX}.sort_tmp \
        | pairtools dedup --nproc-in ${CPUS} --nproc-out ${CPUS} --output-stats ${PREFIX}.dedup.stats \
            -o ${PREFIX}.incomplete.pairs.gz
    mv ${PREFIX}.incomplete.pairs.gz ${PREFIX}.pairs.gz
    rmdir ${PREFIX}.sort_tmp
fi

# STEP 5: Convert mapped contacts to Juicer's merged_nodups (mnd) format, keeping MAPQ values
if [[ ! -f ${PREFIX}.mnd.txt ]]; then
    zcat ${PREFIX}.pairs.gz | awk '
    BEGIN{OFS="\t"}
    /^#columns:/{for(i=2;i<=NF;i++){col[$i]=i-1}; next}
    /^#/{next}
    $(col["chrom1"])!="!" && $(col["chrom2"])!="!" {
        print ($(col["strand1"])=="+")?0:16, $(col["chrom1"]), $(col["pos1"]), 0,
              ($(col["strand2"])=="+")?0:16, $(col["chrom2"]), $(col["pos2"]), 1,
              $(col["mapq1"]), "-", "-", $(col["mapq2"]), "-", "-", $(col["readID"]), $(col["readID"]);
    }' > ${PREFIX}.incomplete.mnd.txt
    mv ${PREFIX}.incomplete.mnd.txt ${PREFIX}.mnd.txt
fi

# STEP 6: Build the .hic file, suitable for editing with Juicebox Assembly Tools
bash ${THREEDDNADIR}/visualize/run-assembly-visualizer.sh -p false -c -q ${MAPQ} \
    ${PREFIX}.assembly ${PREFIX}.mnd.txt
