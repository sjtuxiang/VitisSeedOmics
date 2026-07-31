#!/usr/bin/env bash

set -Eeuo pipefail
shopt -s nullglob

trap 'printf "[ERROR] Line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR

###############################################################################
# Configuration
###############################################################################

PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
THREADS="${SLURM_CPUS_PER_TASK:-20}"
SAMPLE_GLOB="${SAMPLE_GLOB:-V*/}"
HAPLOTYPES=(hap1 hap2)

###############################################################################
# Helper functions
###############################################################################

extract_anchors() {
    local collinearity_file="$1"
    local output_prefix="$2"

    if [[ ! -s "$collinearity_file" ]]; then
        printf "[ERROR] Collinearity file not found or empty: %s\n" \
            "$collinearity_file" >&2
        return 1
    fi

    python3 - "$collinearity_file" "$output_prefix" <<'PYTHON_ANCHORS'
import sys
from pathlib import Path


def extract_anchors(collinearity_file: Path, output_prefix: str) -> None:
    anchors_simple = Path(f"{output_prefix}.anchors.simple")
    anchors_with_block = Path(f"{output_prefix}.anchors")

    print(
        f"   Extracting anchors: {collinearity_file} -> "
        f"{anchors_simple}, {anchors_with_block}"
    )

    block_id = -1
    simple_count = 0

    with (
        collinearity_file.open("r", encoding="utf-8") as source,
        anchors_simple.open("w", encoding="utf-8") as simple_out,
        anchors_with_block.open("w", encoding="utf-8") as block_out,
    ):
        for raw_line in source:
            line = raw_line.strip()

            if line.startswith("## Alignment"):
                block_id += 1
                continue

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            # A standard MCScanX anchor record contains at least four fields:
            # alignment position, pair index, gene A, and gene B.
            if len(fields) < 4:
                continue

            gene_a = fields[2]
            gene_b = fields[3]

            simple_out.write(f"{gene_a}\t{gene_b}\n")
            block_out.write(f"{gene_a}\t{gene_b}\t{block_id}\n")
            simple_count += 1

    print(
        f"   Anchor extraction completed: {simple_count} gene pairs "
        f"across {block_id + 1} collinearity blocks."
    )


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: embedded-anchor-parser "
            "<collinearity_file> <output_prefix>"
        )

    collinearity_file = Path(sys.argv[1])
    output_prefix = sys.argv[2]

    if not collinearity_file.is_file():
        raise FileNotFoundError(
            f"Collinearity file does not exist: {collinearity_file}"
        )

    extract_anchors(collinearity_file, output_prefix)


if __name__ == "__main__":
    main()
PYTHON_ANCHORS
}

extract_mcscanx_gff() {
    local input_gff3="$1"
    local output_gff="$2"

    awk -F '\t' '
        BEGIN {
            OFS = "\t"
        }
        $0 !~ /^#/ && $3 == "gene" {
            gene_id = ""
            attribute_count = split($9, attributes, ";")

            for (i = 1; i <= attribute_count; i++) {
                if (attributes[i] ~ /^ID=/) {
                    gene_id = attributes[i]
                    sub(/^ID=/, "", gene_id)
                    break
                }
            }

            if (gene_id != "") {
                print $1, gene_id, $4, $5
            }
        }
    ' "$input_gff3" > "$output_gff"

    if [[ ! -s "$output_gff" ]]; then
        printf "[ERROR] No gene records with ID attributes were extracted from %s\n" \
            "$input_gff3" >&2
        return 1
    fi
}

###############################################################################
# Environment and input checks
###############################################################################

if [[ ! -d "$PROJECT_DIR" ]]; then
    printf "[ERROR] Project directory does not exist: %s\n" "$PROJECT_DIR" >&2
    exit 1
fi

required_commands=(
    awk
    diamond
    MCScanX
    duplicate_gene_classifier
    python3
)

for command_name in "${required_commands[@]}"; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        printf "[ERROR] Required command not found in PATH: %s\n" \
            "$command_name" >&2
        exit 1
    fi
done

cd "$PROJECT_DIR"

sample_directories=($SAMPLE_GLOB)

if (( ${#sample_directories[@]} == 0 )); then
    printf "[ERROR] No sample directories matched '%s' under %s\n" \
        "$SAMPLE_GLOB" "$PROJECT_DIR" >&2
    exit 1
fi

###############################################################################
# Main workflow
###############################################################################

printf "===== MCScanX batch workflow started =====\n"
printf "Project directory: %s\n" "$PROJECT_DIR"
printf "Sample directory pattern: %s\n" "$SAMPLE_GLOB"
printf "Threads: %s\n" "$THREADS"
printf "Samples detected: %s\n" "${#sample_directories[@]}"

completed_runs=0
skipped_runs=0
failed_runs=0

for sample_directory in "${sample_directories[@]}"; do
    sample_directory="${sample_directory%/}"
    sample="$(basename "$sample_directory")"

    printf "\n====== Processing sample: %s ======\n" "$sample"

    for haplotype in "${HAPLOTYPES[@]}"; do
        printf -- "--- %s %s ---\n" "$sample" "$haplotype"

        gff3="${sample_directory}/${sample}.${haplotype}.gff3"
        protein_fasta="${sample_directory}/${sample}.${haplotype}.pep.fa"
        output_prefix="${sample}_${haplotype}"
        output_base="${sample_directory}/${output_prefix}"

        if [[ ! -s "$gff3" ]]; then
            printf "[WARN] Missing or empty GFF3 file; skipping: %s\n" \
                "$gff3" >&2
            ((skipped_runs += 1))
            continue
        fi

        if [[ ! -s "$protein_fasta" ]]; then
            printf "[WARN] Missing or empty protein FASTA file; skipping: %s\n" \
                "$protein_fasta" >&2
            ((skipped_runs += 1))
            continue
        fi

        if ! (
            cd "$sample_directory"

            printf "   Creating MCScanX GFF: %s.gff\n" "$output_prefix"
            extract_mcscanx_gff \
                "${sample}.${haplotype}.gff3" \
                "${output_prefix}.gff" \
                || exit 1

            printf "   Building DIAMOND database...\n"
            diamond makedb \
                --in "${sample}.${haplotype}.pep.fa" \
                --db "$output_prefix" \
                --quiet \
                || exit 1

            printf "   Running DIAMOND self-alignment...\n"
            diamond blastp \
                --query "${sample}.${haplotype}.pep.fa" \
                --db "$output_prefix" \
                --evalue 1e-5 \
                --max-target-seqs 5 \
                --threads "$THREADS" \
                --outfmt 6 \
                --out "${output_prefix}.blast" \
                --quiet \
                || exit 1

            if [[ ! -s "${output_prefix}.blast" ]]; then
                printf "[ERROR] DIAMOND produced an empty BLAST file: %s.blast\n" \
                    "$output_prefix" >&2
                exit 1
            fi

            printf "   Running MCScanX...\n"
            MCScanX "$output_prefix" || exit 1

            if [[ ! -s "${output_prefix}.collinearity" ]]; then
                printf "[ERROR] MCScanX did not generate a valid collinearity file.\n" \
                    >&2
                exit 1
            fi

            printf "   Classifying gene-duplication types...\n"
            duplicate_gene_classifier "$output_prefix" || exit 1

            extract_anchors \
                "${output_prefix}.collinearity" \
                "$output_prefix" \
                || exit 1
        ); then
            printf "[ERROR] Processing failed for %s %s\n" \
                "$sample" "$haplotype" >&2
            ((failed_runs += 1))
            continue
        fi

        if [[ ! -s "${output_base}.anchors.simple" ||
              ! -s "${output_base}.anchors" ]]; then
            printf "[ERROR] Anchor outputs are missing for %s %s\n" \
                "$sample" "$haplotype" >&2
            ((failed_runs += 1))
            continue
        fi

        printf "[OK] Completed %s %s\n" "$sample" "$haplotype"
        ((completed_runs += 1))
    done
done

###############################################################################
# Final summary
###############################################################################

printf "\n===== MCScanX batch workflow finished =====\n"
printf "Completed: %s\n" "$completed_runs"
printf "Skipped:   %s\n" "$skipped_runs"
printf "Failed:    %s\n" "$failed_runs"

if (( failed_runs > 0 )); then
    exit 1
fi
