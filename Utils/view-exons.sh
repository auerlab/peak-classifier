#!/bin/sh -e

gzcat $(Utils/gff-name.sh) | \
    awk '$3 == "gene" || $3 == "exon" { print $1, $4, $5, $3, $7 }'

