#!/bin/sh -e

for chr in $(seq 1 19); do
    awk -v chr=$chr '$1 == chr { print $0 }' p10-CCA-501-merged.bed | head -500
done
