#!/bin/bash
set -euo pipefail

THREADS=${SLURM_CPUS_PER_TASK:-32}

# ===== 基础路径 =====
WORKDIR=/dssg/home/xiang/Output
cd ${WORKDIR}

GRAPH=VitisPan.gbz
GFF=PN40024.graphnames.gff3
OUTDIR=${WORKDIR}/RNAseq_graph
mkdir -p ${OUTDIR}

# ===== 环境检查 =====
command -v vg >/dev/null 2>&1 || { echo "ERROR: vg 未安装或不在 PATH"; exit 1; }
vg version || true

# ===== 记录图中参考路径 =====
vg paths --list --reference-paths --generic-paths -x ${GRAPH} > ${OUTDIR}/graph.refpaths.txt

# ===== 简单检查 GFF3 是否已经改成图路径命名 =====
head -n 50 ${GFF} | awk 'BEGIN{ok=1} !/^#/ && $1 !~ /^V148#0#Chr/ {ok=0} END{if(ok==0) exit 1}'
echo "[INFO] GFF3 第一列命名看起来与图一致"

# ===== 1. 构建 pantranscriptome =====
vg rna --progress \
  --proj-embed-paths \
  --gbz-format ${GRAPH} \
  --transcripts ${GFF} \
  --transcript-tag Parent \
  --write-info ${OUTDIR}/pantranscript.info.tsv \
  --gbwt-bidirectional \
  --write-gbwt ${OUTDIR}/pantranscript.transcripts.gbwt \
  > ${OUTDIR}/pantranscript.pg

# ===== 2. 构建供 mpmap 使用的索引 =====
vg index -t ${THREADS} \
  -x ${OUTDIR}/pantranscript.xg \
  ${OUTDIR}/pantranscript.pg

vg index -t ${THREADS} \
  -g ${OUTDIR}/pantranscript.gcsa \
  -k 16 \
  ${OUTDIR}/pantranscript.pg

vg index -t ${THREADS} \
  -j ${OUTDIR}/pantranscript.dist \
  ${OUTDIR}/pantranscript.pg

# ===== 3. 生成一个图文件清单，便于后面样本脚本调用 =====
cat > ${OUTDIR}/pantranscript.files.txt <<EOF
PANTRAN_PG=${OUTDIR}/pantranscript.pg
PANTRAN_XG=${OUTDIR}/pantranscript.xg
PANTRAN_GCSA=${OUTDIR}/pantranscript.gcsa
PANTRAN_DIST=${OUTDIR}/pantranscript.dist
PANTRAN_GBWT=${OUTDIR}/pantranscript.transcripts.gbwt
PANTRAN_INFO=${OUTDIR}/pantranscript.info.tsv
EOF

echo "[INFO] pantranscriptome 构建完成"
echo "[INFO] 输出目录: ${OUTDIR}"
