#!/bin/bash
#SBATCH --job-name=sracha-big
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=sracha-big-%j.out

set -euo pipefail

export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK}
OUTDIR=/tmp/sracha-big-${SLURM_JOB_ID}

echo "Starting sracha get SRR17778108 at $(date)"
echo "Cores: $SLURM_CPUS_PER_TASK"

time ./target/release/sracha get SRR17778108 \
    -O ${OUTDIR} \
    --force -t ${SLURM_CPUS_PER_TASK}

echo "Done at $(date)"
ls -lh ${OUTDIR}/
head -4 ${OUTDIR}/*_1.fastq.gz 2>/dev/null || zcat ${OUTDIR}/*_1.fastq.gz 2>/dev/null | head -4 || head -4 ${OUTDIR}/*_0.fastq* 2>/dev/null
