#!/bin/bash
#SBATCH --job-name=sracha-big
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=sracha-big-%j.out

set -euo pipefail

OUTDIR=/tmp/sracha-big-${SLURM_JOB_ID}

echo "Starting sracha get SRR17778108 at $(date)"
echo "Cores: $SLURM_CPUS_PER_TASK"

time ./target/release/sracha get SRR17778108 \
    -O ${OUTDIR} \
    --force -t ${SLURM_CPUS_PER_TASK} --gzip-level 1

echo "Done at $(date)"
ls -lh ${OUTDIR}/
zcat ${OUTDIR}/*_1.fastq.gz 2>/dev/null | head -4 || zcat ${OUTDIR}/*_0.fastq.gz 2>/dev/null | head -4
