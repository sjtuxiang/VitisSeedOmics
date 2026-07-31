#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

trap 'printf "[ERROR] Line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"

if [[ ! -d "$PROJECT_DIR" ]]; then
    printf "[ERROR] Project directory does not exist: %s\n" "$PROJECT_DIR" >&2
    exit 1
fi

cd "$PROJECT_DIR"

GENOME="Vv.fa"
GFF3="Vv.gff3"
RAW_DIR="rawdata"

ANNOTATION_DIR="00.annotation"
INDEX_DIR="00.index"
CLEAN_DIR="01.clean_data"
BAM_DIR="03.bam"
COUNT_DIR="04.counts"
LOG_DIR="05.logs"
EXPRESSION_DIR="06.expression"
QC_DIR="07.qc"

GTF="${ANNOTATION_DIR}/Vv.gtf"
COUNT_GTF="${ANNOTATION_DIR}/Vv.featureCounts.gtf"
MISSING_GENE_ID_EXONS="${ANNOTATION_DIR}/exons_missing_gene_id.gtf"
SPLICE_SITES="${ANNOTATION_DIR}/Vv.splice_sites.txt"
INDEX_PREFIX="${INDEX_DIR}/Vv_hisat2"
INDEX_MARKER="${INDEX_DIR}/Vv_hisat2.build.ok"

TOTAL_THREADS="${SLURM_CPUS_PER_TASK:-40}"
FASTP_THREADS=16
HISAT2_THREADS=24
SAMTOOLS_THREADS=12
FEATURECOUNTS_THREADS=32
SAMTOOLS_SORT_MEM="2G"

FC_STRAND="${FC_STRAND:-0}"
MIN_MAPQ=10

ALL_COUNT="${COUNT_DIR}/gene_counts.featureCounts.txt"
COUNT_SETTINGS="${COUNT_DIR}/gene_counts.settings.txt"
SAMPLE_INFO="${LOG_DIR}/sample_info.tsv"

JOB_TAG="${SLURM_JOB_ID:-$$}"

###############################################################################
# 1. 输入、参数和软件检查
###############################################################################

printf "[INFO] Project directory: %s\n" "$PROJECT_DIR"

for file in "$GENOME" "$GFF3"; do
    if [[ ! -s "$file" ]]; then
        printf "[ERROR] Input file not found or empty: %s\n" "$file" >&2
        exit 1
    fi
done

if [[ ! -d "$RAW_DIR" ]]; then
    printf "[ERROR] Raw-data directory not found: %s\n" "$RAW_DIR" >&2
    exit 1
fi

if [[ "$FC_STRAND" != "0" && "$FC_STRAND" != "1" && "$FC_STRAND" != "2" ]]; then
    printf "[ERROR] FC_STRAND must be 0, 1, or 2; current value: %s\n" \
        "$FC_STRAND" >&2
    exit 1
fi

mkdir -p \
    "$ANNOTATION_DIR" \
    "$INDEX_DIR" \
    "$CLEAN_DIR" \
    "$BAM_DIR" \
    "$COUNT_DIR" \
    "$LOG_DIR" \
    "$EXPRESSION_DIR" \
    "$QC_DIR"

required_commands=(
    fastp
    gzip
    hisat2
    hisat2-build
    hisat2-inspect
    hisat2_extract_splice_sites.py
    samtools
    featureCounts
    gffread
    Rscript
)

for cmd in "${required_commands[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        printf "[ERROR] Command not found in current environment: %s\n" "$cmd" >&2
        exit 1
    fi
done

{
    printf "date\t%s\n" "$(date '+%F %T')"
    printf "conda_environment\t%s\n" "${CONDA_DEFAULT_ENV:-unknown}"
    printf "FC_STRAND\t%s\n" "$FC_STRAND"
    fastp --version 2>&1
    hisat2 --version 2>&1
    samtools --version 2>&1
    featureCounts -v 2>&1
    gffread --version 2>&1
    Rscript --version 2>&1
} > "${LOG_DIR}/software_versions.txt"

###############################################################################
# 2. GFF3 转换成 GTF，并验证注释
###############################################################################

if [[ ! -s "$GTF" || "$GFF3" -nt "$GTF" ]]; then
    printf "[INFO] Converting GFF3 to GTF with gffread...\n"
    TMP_GTF="${GTF}.tmp.${JOB_TAG}"
    gffread "$GFF3" -T -o "$TMP_GTF"

    if [[ ! -s "$TMP_GTF" ]]; then
        printf "[ERROR] gffread produced an empty GTF file.\n" >&2
        exit 1
    fi

    mv -f "$TMP_GTF" "$GTF"
else
    printf "[INFO] Existing GTF is current; skipping conversion.\n"
fi

EXON_COUNT="$(
    awk -F '\t' '$0 !~ /^#/ && $3 == "exon" {n++} END {print n+0}' "$GTF"
)"
EXON_WITH_GENE_ID="$(
    awk -F '\t' \
        '$0 !~ /^#/ && $3 == "exon" && $9 ~ /gene_id "[^"]+"/ {n++}
         END {print n+0}' \
        "$GTF"
)"

if (( EXON_COUNT == 0 )); then
    printf "[ERROR] No exon records were found in converted GTF: %s\n" "$GTF" >&2
    exit 1
fi

MISSING_GENE_ID_COUNT=$(( EXON_COUNT - EXON_WITH_GENE_ID ))

# 少量没有 gene_id 的孤立/异常 exon 不适合随意伪造基因编号：
# 1) 完整 GTF 继续用于 HISAT2 提取剪接位点；
# 2) 缺少 gene_id 的 exon 单独保存；
# 3) featureCounts 使用仅含有效 gene_id exon 的专用 GTF。
awk -F '\t' \
    '$0 !~ /^#/ && $3 == "exon" && $9 !~ /gene_id "[^"]+"/' \
    "$GTF" \
    > "$MISSING_GENE_ID_EXONS"

if (( MISSING_GENE_ID_COUNT > 0 )); then
    printf "[WARN] %s of %s exons lack gene_id and will be excluded only from gene-level counting.\n" \
        "$MISSING_GENE_ID_COUNT" "$EXON_COUNT" >&2
    printf "[WARN] These records were saved to: %s\n" \
        "$MISSING_GENE_ID_EXONS" >&2

    # 若超过全部 exon 的 1%，说明注释结构可能存在系统性问题，不再自动继续。
    if (( MISSING_GENE_ID_COUNT * 100 > EXON_COUNT )); then
        printf "[ERROR] More than 1%% of exons lack gene_id; inspect the annotation before continuing.\n" >&2
        exit 1
    fi
fi

if [[ ! -s "$COUNT_GTF" || "$GTF" -nt "$COUNT_GTF" ]]; then
    TMP_COUNT_GTF="${COUNT_GTF}.tmp.${JOB_TAG}"
    awk -F '\t' '
        BEGIN {OFS = "\t"}
        /^#/ {print; next}
        $3 == "exon" && $9 ~ /gene_id "[^"]+"/ {print}
    ' "$GTF" > "$TMP_COUNT_GTF"

    COUNT_EXON_COUNT="$(
        awk -F '\t' '$0 !~ /^#/ && $3 == "exon" {n++} END {print n+0}' \
            "$TMP_COUNT_GTF"
    )"

    if (( COUNT_EXON_COUNT != EXON_WITH_GENE_ID )); then
        printf "[ERROR] featureCounts GTF validation failed: expected=%s, observed=%s\n" \
            "$EXON_WITH_GENE_ID" "$COUNT_EXON_COUNT" >&2
        exit 1
    fi

    mv -f "$TMP_COUNT_GTF" "$COUNT_GTF"
fi

printf "[INFO] GTF exon records: total=%s, usable_for_gene_counting=%s, excluded=%s\n" \
    "$EXON_COUNT" "$EXON_WITH_GENE_ID" "$MISSING_GENE_ID_COUNT"

if [[ ! -s "${GENOME}.fai" || "$GENOME" -nt "${GENOME}.fai" ]]; then
    samtools faidx "$GENOME"
fi

awk '{print $1}' "${GENOME}.fai" \
    | LC_ALL=C sort -u \
    > "${ANNOTATION_DIR}/genome.seqids.txt"

awk -F '\t' '$0 !~ /^#/ {print $1}' "$GTF" \
    | LC_ALL=C sort -u \
    > "${ANNOTATION_DIR}/annotation.seqids.txt"

comm -23 \
    "${ANNOTATION_DIR}/annotation.seqids.txt" \
    "${ANNOTATION_DIR}/genome.seqids.txt" \
    > "${ANNOTATION_DIR}/annotation_seqids_missing_from_genome.txt"

MISSING_SEQID_COUNT="$(
    wc -l < "${ANNOTATION_DIR}/annotation_seqids_missing_from_genome.txt"
)"

if (( MISSING_SEQID_COUNT > 0 )); then
    printf "[ERROR] %s annotation sequence IDs are absent from Vv.fa.\n" \
        "$MISSING_SEQID_COUNT" >&2
    printf "[ERROR] First missing IDs:\n" >&2
    sed -n '1,20p' \
        "${ANNOTATION_DIR}/annotation_seqids_missing_from_genome.txt" >&2
    exit 1
fi

###############################################################################
# 3. 提取已知剪接位点，建立 HISAT2 线性索引
###############################################################################

if [[ ! -s "$SPLICE_SITES" || "$GTF" -nt "$SPLICE_SITES" ]]; then
    printf "[INFO] Extracting known splice sites from GTF...\n"
    TMP_SPLICE="${SPLICE_SITES}.tmp.${JOB_TAG}"
    hisat2_extract_splice_sites.py "$GTF" > "$TMP_SPLICE"

    if [[ ! -s "$TMP_SPLICE" ]]; then
        printf "[ERROR] No splice sites were extracted from GTF.\n" >&2
        exit 1
    fi

    mv -f "$TMP_SPLICE" "$SPLICE_SITES"
fi

INDEX_OK=false
if [[ -s "$INDEX_MARKER" && ! "$GENOME" -nt "$INDEX_MARKER" ]]; then
    if hisat2-inspect -n "$INDEX_PREFIX" >/dev/null 2>&1; then
        INDEX_OK=true
    fi
fi

if [[ "$INDEX_OK" == false ]]; then
    printf "[INFO] Building HISAT2 genome index...\n"
    TMP_INDEX_PREFIX="${INDEX_PREFIX}.tmp.${JOB_TAG}"

    hisat2-build \
        -p "$TOTAL_THREADS" \
        "$GENOME" \
        "$TMP_INDEX_PREFIX" \
        > "${LOG_DIR}/hisat2-build.log" 2>&1

    built_index_files=(
        "${TMP_INDEX_PREFIX}".*.ht2
        "${TMP_INDEX_PREFIX}".*.ht2l
    )

    if (( ${#built_index_files[@]} == 0 )); then
        printf "[ERROR] HISAT2 index files were not generated.\n" >&2
        exit 1
    fi

    if (( ${#built_index_files[@]} != 8 )); then
        printf "[ERROR] Expected 8 HISAT2 index files, but found %s temporary files.\n" \
            "${#built_index_files[@]}" >&2
        exit 1
    fi

    # 只删除精确命名的旧正式索引。不能使用
    # "${INDEX_PREFIX}".*.ht2，因为该模式也会匹配并删除刚生成的
    # "${INDEX_PREFIX}.tmp.${JOB_TAG}".*.ht2 临时索引。
    for index_part in {1..8}; do
        rm -f \
            "${INDEX_PREFIX}.${index_part}.ht2" \
            "${INDEX_PREFIX}.${index_part}.ht2l"
    done

    for index_file in "${built_index_files[@]}"; do
        suffix="${index_file#${TMP_INDEX_PREFIX}}"
        mv -f "$index_file" "${INDEX_PREFIX}${suffix}"
    done

    touch "$INDEX_MARKER"

    if ! hisat2-inspect -n "$INDEX_PREFIX" >/dev/null 2>&1; then
        printf "[ERROR] The completed HISAT2 index could not be inspected.\n" >&2
        exit 1
    fi
else
    printf "[INFO] HISAT2 index already exists and is valid; skipping build.\n"
fi

###############################################################################
# 4. 识别双端样本并生成样本信息表
###############################################################################

r1_files=("${RAW_DIR}"/*_R1.fq.gz)

if (( ${#r1_files[@]} == 0 )); then
    printf "[ERROR] No files matching %s/*_R1.fq.gz were found.\n" \
        "$RAW_DIR" >&2
    exit 1
fi

samples_unsorted=()
pair_error=0

for r1 in "${r1_files[@]}"; do
    r1_name="$(basename "$r1")"
    sample="${r1_name%_R1.fq.gz}"
    r2="${RAW_DIR}/${sample}_R2.fq.gz"

    if [[ ! -s "$r2" ]]; then
        printf "[ERROR] Missing R2 mate for sample %s: %s\n" "$sample" "$r2" >&2
        pair_error=1
    else
        samples_unsorted+=("$sample")
    fi
done

for r2 in "${RAW_DIR}"/*_R2.fq.gz; do
    r2_name="$(basename "$r2")"
    sample="${r2_name%_R2.fq.gz}"
    r1="${RAW_DIR}/${sample}_R1.fq.gz"

    if [[ ! -s "$r1" ]]; then
        printf "[ERROR] Missing R1 mate for sample %s: %s\n" "$sample" "$r1" >&2
        pair_error=1
    fi
done

if (( pair_error != 0 )); then
    exit 1
fi

mapfile -t samples < <(
    printf '%s\n' "${samples_unsorted[@]}" | LC_ALL=C sort -V
)

if (( ${#samples[@]} == 0 )); then
    printf "[ERROR] No complete paired-end samples were found.\n" >&2
    exit 1
fi

TMP_SAMPLE_INFO="${SAMPLE_INFO}.tmp.${JOB_TAG}"
{
    printf "sample\tgroup\treplicate\n"
    for sample in "${samples[@]}"; do
        group="${sample%-*}"
        replicate="${sample##*-}"
        printf "%s\t%s\t%s\n" "$sample" "$group" "$replicate"
    done
} > "$TMP_SAMPLE_INFO"
mv -f "$TMP_SAMPLE_INFO" "$SAMPLE_INFO"

printf "[INFO] Total paired-end samples: %s\n" "${#samples[@]}"
printf '%s\n' "${samples[@]}" > "${LOG_DIR}/sample_list.txt"

###############################################################################
# 5. fastp 质控、HISAT2 比对、BAM 排序与检查
###############################################################################

bam_list=()

for sample in "${samples[@]}"; do
    printf "[INFO] Processing sample: %s\n" "$sample"

    RAW_R1="${RAW_DIR}/${sample}_R1.fq.gz"
    RAW_R2="${RAW_DIR}/${sample}_R2.fq.gz"

    CLEAN_R1="${CLEAN_DIR}/${sample}_clean_R1.fq.gz"
    CLEAN_R2="${CLEAN_DIR}/${sample}_clean_R2.fq.gz"
    FASTP_HTML="${LOG_DIR}/${sample}.fastp.html"
    FASTP_JSON="${LOG_DIR}/${sample}.fastp.json"

    BAM="${BAM_DIR}/${sample}.sorted.bam"
    BAI="${BAM}.bai"
    bam_list+=("$BAM")

    CLEAN_OK=false
    if [[ -s "$CLEAN_R1" && -s "$CLEAN_R2" && \
          -s "$FASTP_HTML" && -s "$FASTP_JSON" ]] && \
       gzip -t "$CLEAN_R1" "$CLEAN_R2" >/dev/null 2>&1; then
        CLEAN_OK=true
    fi

    if [[ "$CLEAN_OK" == true ]]; then
        printf "[INFO] fastp outputs exist; skipping sample %s.\n" "$sample"
    else
        if [[ -s "$CLEAN_R1" || -s "$CLEAN_R2" ]]; then
            printf "[WARN] Existing clean FASTQ files for %s are incomplete or not valid gzip; regenerating them.\n" \
                "$sample" >&2
        fi

        # 临时文件名必须仍以 .fq.gz 结尾，否则 fastp 会写出未压缩文本；
        # 之后若直接重命名为 .fq.gz，HISAT2 会将其误当成 gzip 并读到0条序列。
        TMP_CLEAN_R1="${CLEAN_DIR}/.${sample}_clean_R1.tmp.${JOB_TAG}.fq.gz"
        TMP_CLEAN_R2="${CLEAN_DIR}/.${sample}_clean_R2.tmp.${JOB_TAG}.fq.gz"

        fastp \
            --thread "$FASTP_THREADS" \
            --in1 "$RAW_R1" \
            --in2 "$RAW_R2" \
            --out1 "$TMP_CLEAN_R1" \
            --out2 "$TMP_CLEAN_R2" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --n_base_limit 5 \
            --length_required 50 \
            --compression 6 \
            --html "$FASTP_HTML" \
            --json "$FASTP_JSON"

        if ! gzip -t "$TMP_CLEAN_R1" "$TMP_CLEAN_R2"; then
            printf "[ERROR] fastp output failed gzip validation for sample %s.\n" \
                "$sample" >&2
            exit 1
        fi

        mv -f "$TMP_CLEAN_R1" "$CLEAN_R1"
        mv -f "$TMP_CLEAN_R2" "$CLEAN_R2"
    fi

    BAM_OK=false
    if [[ -s "$BAM" ]] && samtools quickcheck -q "$BAM"; then
        PAIRED_ALIGNMENT_RECORDS="$(samtools view -c -f 1 "$BAM")"
        if (( PAIRED_ALIGNMENT_RECORDS > 0 )); then
            BAM_OK=true
        else
            printf "[WARN] Existing BAM contains no paired alignment records and will be regenerated: %s\n" \
                "$BAM" >&2
        fi
    fi

    if [[ "$BAM_OK" == false ]]; then
        printf "[INFO] Running HISAT2 for sample %s...\n" "$sample"
        TMP_BAM="${BAM}.tmp.${JOB_TAG}"

        hisat2 \
            -p "$HISAT2_THREADS" \
            --known-splicesite-infile "$SPLICE_SITES" \
            --no-unal \
            --rg-id "$sample" \
            --rg "SM:${sample}" \
            --rg "LB:${sample}" \
            --rg "PL:ILLUMINA" \
            --new-summary \
            --summary-file "${LOG_DIR}/${sample}.hisat2.summary.txt" \
            -x "$INDEX_PREFIX" \
            -1 "$CLEAN_R1" \
            -2 "$CLEAN_R2" \
            2> "${LOG_DIR}/${sample}.hisat2.stderr.log" \
        | samtools sort \
            -@ "$SAMTOOLS_THREADS" \
            -m "$SAMTOOLS_SORT_MEM" \
            -o "$TMP_BAM" \
            -

        samtools quickcheck -v "$TMP_BAM"

        TMP_PAIRED_ALIGNMENT_RECORDS="$(samtools view -c -f 1 "$TMP_BAM")"
        if (( TMP_PAIRED_ALIGNMENT_RECORDS == 0 )); then
            printf "[ERROR] HISAT2 produced no paired alignment records for sample %s.\n" \
                "$sample" >&2
            printf "[ERROR] Inspect %s/%s.hisat2.summary.txt and the clean FASTQ files.\n" \
                "$LOG_DIR" "$sample" >&2
            exit 1
        fi

        mv -f "$TMP_BAM" "$BAM"
    else
        printf "[INFO] Valid BAM exists; skipping alignment for %s.\n" "$sample"
    fi

    if [[ ! -s "$BAI" || "$BAM" -nt "$BAI" ]]; then
        TMP_BAI="${BAI}.tmp.${JOB_TAG}"
        samtools index -@ "$SAMTOOLS_THREADS" "$BAM" "$TMP_BAI"
        mv -f "$TMP_BAI" "$BAI"
    fi

    samtools flagstat \
        -@ "$SAMTOOLS_THREADS" \
        "$BAM" \
        > "${LOG_DIR}/${sample}.samtools.flagstat.txt"

    samtools idxstats \
        "$BAM" \
        > "${LOG_DIR}/${sample}.samtools.idxstats.txt"
done

###############################################################################
# 6. featureCounts：按基因统计 paired-end fragments
###############################################################################

CURRENT_COUNT_SETTINGS="$(
    printf 'FC_STRAND=%s;MIN_MAPQ=%s;GTF=%s' \
        "$FC_STRAND" "$MIN_MAPQ" "$COUNT_GTF"
)"

COUNT_OK=false
if [[ -s "$ALL_COUNT" && -s "${ALL_COUNT}.summary" && -s "$COUNT_SETTINGS" ]]; then
    SAVED_COUNT_SETTINGS="$(<"$COUNT_SETTINGS")"
    if [[ "$SAVED_COUNT_SETTINGS" == "$CURRENT_COUNT_SETTINGS" && \
          ! "$COUNT_GTF" -nt "$ALL_COUNT" ]]; then
        COUNT_OK=true
        for bam in "${bam_list[@]}"; do
            if [[ "$bam" -nt "$ALL_COUNT" ]]; then
                COUNT_OK=false
                break
            fi
        done
    fi
fi

if [[ "$COUNT_OK" == false ]]; then
    printf "[INFO] Running featureCounts for all samples...\n"
    TMP_COUNT="${ALL_COUNT}.tmp.${JOB_TAG}"

    featureCounts \
        -T "$FEATURECOUNTS_THREADS" \
        -F GTF \
        -a "$COUNT_GTF" \
        -t exon \
        -g gene_id \
        -p \
        --countReadPairs \
        -B \
        -C \
        -Q "$MIN_MAPQ" \
        -s "$FC_STRAND" \
        -o "$TMP_COUNT" \
        "${bam_list[@]}" \
        > "${LOG_DIR}/featureCounts.log" 2>&1

    if [[ ! -s "$TMP_COUNT" || ! -s "${TMP_COUNT}.summary" ]]; then
        printf "[ERROR] featureCounts did not generate complete output.\n" >&2
        exit 1
    fi

    mv -f "$TMP_COUNT" "$ALL_COUNT"
    mv -f "${TMP_COUNT}.summary" "${ALL_COUNT}.summary"
    printf '%s\n' "$CURRENT_COUNT_SETTINGS" > "$COUNT_SETTINGS"
else
    printf "[INFO] featureCounts output is current; skipping recount.\n"
fi

###############################################################################
# 7. 由 featureCounts 的 gene counts 和 exon-union Length 计算 TPM/FPKM
#
# FPKM = count / (gene_length_kb * assigned_fragments_million)
# TPM  = RPK / sum(RPK) * 1e6
###############################################################################

Rscript --vanilla - \
    "$ALL_COUNT" \
    "$SAMPLE_INFO" \
    "$EXPRESSION_DIR" \
    "$QC_DIR" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
sample_info_file <- args[2]
expression_dir <- args[3]
qc_dir <- args[4]

dir.create(expression_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

normalize_sample_name <- function(x) {
    x <- basename(x)
    sub("\\.sorted\\.bam$", "", x)
}

fc <- read.delim(
    count_file,
    header = TRUE,
    comment.char = "#",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
)

required_columns <- c("Geneid", "Length")
missing_columns <- setdiff(required_columns, colnames(fc))
if (length(missing_columns) > 0L) {
    stop(
        "Missing required featureCounts columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

length_column <- match("Length", colnames(fc))
if (length_column >= ncol(fc)) {
    stop("No sample count columns were found after the Length column.")
}

count_columns <- seq.int(length_column + 1L, ncol(fc))
counts <- as.matrix(fc[, count_columns, drop = FALSE])
storage.mode(counts) <- "numeric"

sample_names <- normalize_sample_name(colnames(counts))
if (anyDuplicated(sample_names)) {
    stop("Duplicated sample names after removing BAM path/suffix.")
}
colnames(counts) <- sample_names

gene_ids <- as.character(fc$Geneid)
gene_length_bp <- as.numeric(fc$Length)

if (anyDuplicated(gene_ids)) {
    stop("Duplicated Geneid values were found in featureCounts output.")
}
if (any(!is.finite(gene_length_bp)) || any(gene_length_bp <= 0)) {
    stop("Invalid non-positive or non-numeric gene lengths were found.")
}
if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Invalid counts were found.")
}

rownames(counts) <- gene_ids
length_kb <- gene_length_bp / 1000
assigned_library_size <- colSums(counts)

if (any(assigned_library_size <= 0)) {
    stop(
        "At least one sample has zero assigned gene counts: ",
        paste(names(assigned_library_size)[assigned_library_size <= 0],
              collapse = ", ")
    )
}

rpk <- sweep(counts, 1, length_kb, "/")
rpk_sum <- colSums(rpk)

if (any(rpk_sum <= 0)) {
    stop("At least one sample has a zero RPK sum.")
}

tpm <- sweep(rpk, 2, rpk_sum / 1e6, "/")
fpkm <- sweep(rpk, 2, assigned_library_size / 1e6, "/")

write_expression_matrix <- function(matrix_data, output_file, digits = NULL) {
    if (!is.null(digits)) {
        matrix_data <- round(matrix_data, digits)
    }
    output <- data.frame(
        GeneID = gene_ids,
        Length_bp = gene_length_bp,
        matrix_data,
        check.names = FALSE
    )
    write.table(
        output,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
    )
}

write_expression_matrix(
    counts,
    file.path(expression_dir, "gene_counts_matrix.tsv")
)
write_expression_matrix(
    tpm,
    file.path(expression_dir, "gene_TPM_matrix.tsv"),
    digits = 6
)
write_expression_matrix(
    fpkm,
    file.path(expression_dir, "gene_FPKM_matrix.tsv"),
    digits = 6
)

sample_info <- read.delim(
    sample_info_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

required_info <- c("sample", "group", "replicate")
if (!all(required_info %in% colnames(sample_info))) {
    stop("sample_info.tsv lacks sample/group/replicate columns.")
}

sample_info <- sample_info[
    match(colnames(counts), sample_info$sample),
    ,
    drop = FALSE
]

if (any(is.na(sample_info$sample))) {
    stop("Some count-matrix samples were not found in sample_info.tsv.")
}

sample_stats <- data.frame(
    sample = colnames(counts),
    group = sample_info$group,
    replicate = sample_info$replicate,
    assigned_gene_fragments = as.numeric(assigned_library_size),
    detected_genes_count_gt_0 = colSums(counts > 0),
    detected_genes_count_ge_10 = colSums(counts >= 10),
    TPM_sum = colSums(tpm),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

summary_file <- paste0(count_file, ".summary")
if (file.exists(summary_file)) {
    fc_summary <- read.delim(
        summary_file,
        header = TRUE,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    status <- fc_summary[[1]]
    summary_matrix <- as.matrix(fc_summary[, -1, drop = FALSE])
    storage.mode(summary_matrix) <- "numeric"
    rownames(summary_matrix) <- status
    colnames(summary_matrix) <- normalize_sample_name(colnames(summary_matrix))

    summary_matrix <- summary_matrix[
        ,
        match(colnames(counts), colnames(summary_matrix)),
        drop = FALSE
    ]

    if ("Assigned" %in% rownames(summary_matrix)) {
        total_fragments <- colSums(summary_matrix)
        sample_stats$featureCounts_total_fragments <- total_fragments
        sample_stats$featureCounts_assigned_fragments <-
            summary_matrix["Assigned", ]
        sample_stats$featureCounts_assignment_rate <-
            summary_matrix["Assigned", ] / total_fragments
    }
}

write.table(
    sample_stats,
    file = file.path(expression_dir, "sample_expression_statistics.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

log_tpm <- log2(tpm + 1)
sample_correlation <- cor(log_tpm, method = "pearson")

write.table(
    data.frame(
        Sample = rownames(sample_correlation),
        sample_correlation,
        check.names = FALSE
    ),
    file = file.path(qc_dir, "sample_correlation_log2TPM.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

gene_variance <- apply(log_tpm, 1, var)
valid_genes <- which(is.finite(gene_variance) & gene_variance > 0)

if (length(valid_genes) >= 2L &&
    requireNamespace("ggplot2", quietly = TRUE)) {
    top_n <- min(5000L, length(valid_genes))
    top_genes <- valid_genes[
        order(gene_variance[valid_genes], decreasing = TRUE)[seq_len(top_n)]
    ]

    pca <- prcomp(t(log_tpm[top_genes, , drop = FALSE]), scale. = FALSE)
    variance_percent <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

    pca_data <- data.frame(
        sample = rownames(pca$x),
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        group = sample_info$group,
        stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(
        pca_data,
        ggplot2::aes(x = PC1, y = PC2, color = group, label = sample)
    ) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_text(vjust = -0.8, size = 2.8, show.legend = FALSE) +
        ggplot2::labs(
            x = sprintf("PC1 (%.2f%%)", variance_percent[1]),
            y = sprintf("PC2 (%.2f%%)", variance_percent[2]),
            color = "Group",
            title = "PCA based on log2(TPM + 1)"
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(hjust = 0.5)
        )

    ggplot2::ggsave(
        filename = file.path(qc_dir, "PCA_log2TPM.pdf"),
        plot = p,
        width = 8,
        height = 6
    )
    ggplot2::ggsave(
        filename = file.path(qc_dir, "PCA_log2TPM.png"),
        plot = p,
        width = 8,
        height = 6,
        dpi = 300
    )
}

if (requireNamespace("pheatmap", quietly = TRUE)) {
    annotation_col <- data.frame(
        Group = sample_info$group,
        row.names = sample_info$sample,
        check.names = FALSE
    )

    pheatmap::pheatmap(
        sample_correlation,
        annotation_col = annotation_col,
        annotation_row = annotation_col,
        border_color = NA,
        main = "Sample correlation: log2(TPM + 1)",
        filename = file.path(qc_dir, "sample_correlation_log2TPM.pdf"),
        width = 9,
        height = 8
    )

    pheatmap::pheatmap(
        sample_correlation,
        annotation_col = annotation_col,
        annotation_row = annotation_col,
        border_color = NA,
        main = "Sample correlation: log2(TPM + 1)",
        filename = file.path(qc_dir, "sample_correlation_log2TPM.png"),
        width = 9,
        height = 8
    )
}

message("Counts, TPM, FPKM and expression QC outputs were generated successfully.")
RSCRIPT

###############################################################################
# 8. 汇总质控报告
###############################################################################

if command -v multiqc >/dev/null 2>&1; then
    if ! multiqc \
        "$LOG_DIR" \
        "$COUNT_DIR" \
        --outdir "${QC_DIR}/multiqc" \
        --force \
        > "${LOG_DIR}/multiqc.log" 2>&1; then
        printf "[WARN] MultiQC failed; see %s/multiqc.log\n" "$LOG_DIR" >&2
    fi
else
    printf "[WARN] multiqc is not installed; skipping the combined QC report.\n" >&2
fi

printf "[INFO] RNA-seq workflow completed successfully.\n"
printf "[INFO] Raw count matrix: %s/gene_counts_matrix.tsv\n" "$EXPRESSION_DIR"
printf "[INFO] TPM matrix: %s/gene_TPM_matrix.tsv\n" "$EXPRESSION_DIR"
printf "[INFO] FPKM matrix: %s/gene_FPKM_matrix.tsv\n" "$EXPRESSION_DIR"
printf "[INFO] Sample metadata: %s\n" "$SAMPLE_INFO"
