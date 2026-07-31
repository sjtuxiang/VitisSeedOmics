#!/bin/bash

# ===== 1. 核心目录规划 =====
JOB_STORE="$(pwd)/03_JobStore"

# 临时目录在计算节点的本地固态硬盘
LOCAL_TMP="/tmp/${USER}_cactus_${SLURM_JOB_ID}"
mkdir -p "$LOCAL_TMP"

export TOIL_WORKDIR="$LOCAL_TMP"
export TMPDIR="$LOCAL_TMP"
export SINGULARITY_TMPDIR="$LOCAL_TMP"
export APPTAINER_TMPDIR="$LOCAL_TMP"

# ===== 2. 自动断点续传检测 =====
RESTART_FLAG=""
# 注意：if 和 [ 之间必须有空格！
if [ -d "$JOB_STORE" ]; then
    echo "=========================================================="
    echo "检测到已存在的 JobStore: 03_JobStore"
    echo "自动启用 --restart 模式，从上次中断的地方继续..."
    echo "=========================================================="
    RESTART_FLAG="--restart"
fi

# ===== 3. 运行 Cactus =====
singularity run \
    --bind ${LOCAL_TMP}:${LOCAL_TMP} \
    --bind $(pwd):$(pwd) \
    cactus_v3.2.0.sif cactus-pangenome \
    $JOB_STORE \
    $(pwd)/02_Config/vitis_seqfile.txt \
    $RESTART_FLAG \
    --outDir $(pwd)/Output \
    --outName VitisPan \
    --maxCores 64 \
    --maxMemory 480G \
    --mgSplit \
    --mgCores 16 \
    --mapCores 16 \
    --consCores 16 \
    --indexCores 16 \
    --vcfwaveCores 16 \
    --vcfwave \
    --gbz clip \
    --odgi clip \
    --chrom-vg clip \
    --chrom-og clip \
    --reference V148 \
    --giraffe clip \
    --vcf clip \
    --logFile $(pwd)/all_log.txt

# ===== 4. 任务正常跑完时的清理 =====
echo "Cactus finished before the 7-day limit. Cleaning up..."
rm -rf "$LOCAL_TMP"
echo "All done!"
