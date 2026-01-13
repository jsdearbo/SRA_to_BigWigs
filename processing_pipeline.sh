#!/bin/bash
set -euo pipefail

# ============================================================
# Configuration
# ============================================================

METADATA_FILE="/mnt/netfiles/jsdearbo/public_data/pipeline_v2/metadata.csv"
BASE_DIR="/mnt/netfiles/jsdearbo/public_data/pipeline_v2"
SRA_SOURCE_DIR="${BASE_DIR}/sra_files"

STAR_EXEC="/mnt/netfiles/conda_envs/rna-seq/bin/STAR"
STAR_INDEX="/mnt/netfiles/genomes/STAR/mm10/star_index"
TRIMMOMATIC_EXEC="/mnt/netfiles/conda_envs/rna-seq/bin/trimmomatic"
ADAPTERS="/mnt/netfiles/conda_envs/rna-seq/share/trimmomatic/adapters/TruSeq3-PE.fa"

ALIGN_THREADS=6
IO_THREADS=4
STAR_OVERHANG=100
export _JAVA_OPTIONS="-Xmx8G"

# ============================================================
# Processing Function
# ============================================================

process_sample() {
    local cell_type=$1
    local srr=$2
    local description=$3

    local out_root="${BASE_DIR}/${cell_type}"
    mkdir -p "$out_root"/{fastq,trimmed,bam,bw,logs,tmp}

    # ------------------------------------------------
    # Per-sample temp directories
    # ------------------------------------------------
    local sample_tmp="${out_root}/tmp/${srr}"
    local fq_tmp="${sample_tmp}/fasterq"
    local star_tmp="${sample_tmp}/star"
    mkdir -p "$fq_tmp" "$star_tmp"

    echo "=================================================="
    echo "Processing: $srr ($cell_type)"
    echo "Desc: $description"
    echo "=================================================="

    # ------------------------------------------------
    # Step 1: SRA -> FASTQ
    # ------------------------------------------------
    local fq1="${out_root}/fastq/${srr}_1.fastq"
    local fq2="${out_root}/fastq/${srr}_2.fastq"

    if [[ ! -f "$fq1" || ! -f "$fq2" ]]; then
        echo "[FASTQ] Dumping ${srr}"
        
        # Handle cases where prefetch puts files in ./SRR/SRR.sra or just ./SRR.sra
        local sra_path="${SRA_SOURCE_DIR}/${srr}/${srr}.sra"
        [[ -f "$sra_path" ]] || sra_path="${SRA_SOURCE_DIR}/${srr}.sra"

        if [[ ! -f "$sra_path" ]]; then
            echo "ERROR: SRA file not found for $srr. Did you run download_data.sh?"
            return 1
        fi

        fasterq-dump --split-3 --threads "$IO_THREADS" --outdir "${out_root}/fastq" --temp "$fq_tmp" "$sra_path"
        sync
    else
        echo "[FASTQ] Exists — skipping"
    fi

    # ------------------------------------------------
    # Step 2: Trimming
    # ------------------------------------------------
    local trim_fq1="${out_root}/trimmed/${srr}_1_trimmed.fastq"
    local trim_fq2="${out_root}/trimmed/${srr}_2_trimmed.fastq"

    if [[ ! -f "$trim_fq1" || ! -f "$trim_fq2" ]]; then
        echo "[TRIM] Trimming reads"
        "$TRIMMOMATIC_EXEC" PE -threads "$ALIGN_THREADS" \
            "$fq1" "$fq2" \
            "$trim_fq1" "${sample_tmp}/1_unpaired.fq" \
            "$trim_fq2" "${sample_tmp}/2_unpaired.fq" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        
        rm -f "${sample_tmp}/"*_unpaired.fq
        sync
    else
        echo "[TRIM] Exists — skipping"
    fi

    # ------------------------------------------------
    # Step 3: STAR alignment
    # ------------------------------------------------
    local final_bam="${out_root}/bam/${srr}.bam"

    if [[ ! -f "$final_bam" ]]; then
        echo "[STAR] Aligning ${srr}"
        "$STAR_EXEC" \
            --runThreadN "$ALIGN_THREADS" \
            --genomeDir "$STAR_INDEX" \
            --readFilesIn "$trim_fq1" "$trim_fq2" \
            --outFileNamePrefix "${star_tmp}/" \
            --outTmpDir "${star_tmp}/tmp" \
            --outSAMtype BAM SortedByCoordinate \
            --sjdbOverhang "$STAR_OVERHANG" \
            --outBAMcompression 1 \
            --outSAMattributes All

        mv "${star_tmp}/Aligned.sortedByCoord.out.bam" "$final_bam"
        # Move logs if they exist
        mv "${star_tmp}/Log."* "${out_root}/logs/" 2>/dev/null || true
        mv "${star_tmp}/SJ.out.tab" "${out_root}/logs/" 2>/dev/null || true

        samtools index "$final_bam"
        sync
    else
        echo "[STAR] BAM exists — skipping"
    fi

    # ------------------------------------------------
    # Step 4: BigWig
    # ------------------------------------------------
    local bw_out="${out_root}/bw/${srr}.bw"

    if [[ ! -f "$bw_out" ]]; then
        echo "[BW] Generating BigWig"
        bamCoverage -b "$final_bam" -o "$bw_out" \
            --binSize 10 \
            --normalizeUsing CPM \
            --effectiveGenomeSize 2652783500 \
            --numberOfProcessors "$IO_THREADS"
        sync
    else
        echo "[BW] Exists — skipping"
    fi

    # Cleanup
    rm -rf "$sample_tmp"
    echo "[DONE] ${srr}"
}

# ============================================================
# Main Execution
# ============================================================

# Read CSV, skip header, parse line by line
while IFS=, read -r cell_type run_id study_id description; do
    # clean up any carriage returns if edited on windows
    run_id=$(echo "$run_id" | tr -d '\r')
    cell_type=$(echo "$cell_type" | tr -d '\r')
    
    if [[ "$cell_type" != "CellType" ]]; then
        process_sample "$cell_type" "$run_id" "$description"
    fi
done < "$METADATA_FILE"

echo "All pipeline tasks complete."