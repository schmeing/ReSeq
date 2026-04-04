#!/usr/bin/env bash
# Downloads large test data files from Zenodo for regression testing.
# Files are cached in test/data/ and verified by MD5 checksum.
# Usage: ./test/download_test_data.sh
#
# Zenodo record: https://zenodo.org/records/19383555
# DOI: 10.5281/zenodo.19383555

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
ZENODO_RECORD="https://zenodo.org/api/records/19383555/files"

mkdir -p "$DATA_DIR"

download_and_verify() {
    local filename="$1"
    local expected_md5="$2"
    local dest="$DATA_DIR/$filename"

    if [ -f "$dest" ]; then
        local actual_md5
        actual_md5=$(md5sum "$dest" | awk '{print $1}')
        if [ "$actual_md5" = "$expected_md5" ]; then
            echo "OK: $filename already exists and checksum matches"
            return 0
        else
            echo "WARN: $filename exists but checksum mismatch, re-downloading"
            rm -f "$dest"
        fi
    fi

    echo "Downloading $filename..."
    curl -L -o "$dest" "$ZENODO_RECORD/$filename/content"

    local actual_md5
    actual_md5=$(md5sum "$dest" | awk '{print $1}')
    if [ "$actual_md5" != "$expected_md5" ]; then
        echo "ERROR: MD5 mismatch for $filename"
        echo "  Expected: $expected_md5"
        echo "  Actual:   $actual_md5"
        rm -f "$dest"
        return 1
    fi
    echo "OK: $filename downloaded and verified"
}

download_and_verify "Hs-Nova-TruSeq.reseq" "c374ef7198effa8ad6ed7fc9e6d9431e"
download_and_verify "Hs-Nova-TruSeq.reseq.ipf" "fbea201bcae34ea025d6408e5e9e91d5"

echo ""
echo "All test data files ready in $DATA_DIR"
