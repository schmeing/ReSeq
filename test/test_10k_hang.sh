#!/bin/bash
# Integration test for the 10k read hang fix (GitHub issue #24).
#
# The seqToIllumina/illuminaPE simulation would hang after generating exactly
# 10,000 reads because written_blocks_ was not incremented after flushing,
# causing condition_variable waits to never be satisfied for subsequent blocks.
#
# This test runs illuminaPE with --numReads 20000 and verifies:
#   1. The process completes (does not hang)
#   2. Output files contain approximately 20,000 read pairs
#
# Usage: test_10k_hang.sh <reseq_binary> <reference.fa> <stats.reseq>
#
# Example:
#   ./test/test_10k_hang.sh build/bin/reseq \
#     test/ecoli-GCF_000005845.2_ASM584v2_genomic.fa \
#     /path/to/stats.reseq

set -euo pipefail

RESEQ="${1:?Usage: $0 <reseq_binary> <reference.fa> <stats.reseq>}"
REF="${2:?Usage: $0 <reseq_binary> <reference.fa> <stats.reseq>}"
STATS="${3:?Usage: $0 <reseq_binary> <reference.fa> <stats.reseq>}"

NUM_READS=20000
TIMEOUT_SECONDS=300
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

R1="$TMPDIR/test_R1.fq"
R2="$TMPDIR/test_R2.fq"

echo "=== 10k hang regression test ==="
echo "Binary:    $RESEQ"
echo "Reference: $REF"
echo "Stats:     $STATS"
echo "Requesting $NUM_READS read pairs..."

# Run with a timeout to catch hangs
if ! timeout "$TIMEOUT_SECONDS" "$RESEQ" illuminaPE \
    -r "$REF" \
    -s "$STATS" \
    -1 "$R1" \
    -2 "$R2" \
    --numReads "$NUM_READS" \
    --threads 2 \
    --seed 42 \
    2>&1; then
    echo "FAIL: reseq illuminaPE did not complete within ${TIMEOUT_SECONDS}s (likely hung at 10k)"
    exit 1
fi

# Count output reads
R1_COUNT=$(grep -c "^@" "$R1" || true)
R2_COUNT=$(grep -c "^@" "$R2" || true)

echo "Generated $R1_COUNT R1 reads and $R2_COUNT R2 reads"

# Verify we got more than 10,000 (the hang point)
if [ "$R1_COUNT" -le 10000 ]; then
    echo "FAIL: Only $R1_COUNT reads generated (hung at 10k boundary?)"
    exit 1
fi

# Verify R1 and R2 counts match
if [ "$R1_COUNT" -ne "$R2_COUNT" ]; then
    echo "FAIL: R1 ($R1_COUNT) and R2 ($R2_COUNT) read counts don't match"
    exit 1
fi

# Verify we got approximately the requested number
MIN_READS=$((NUM_READS * 80 / 100))
if [ "$R1_COUNT" -lt "$MIN_READS" ]; then
    echo "FAIL: Only $R1_COUNT reads generated (expected ~$NUM_READS)"
    exit 1
fi

echo "PASS: Generated $R1_COUNT read pairs (requested $NUM_READS)"
