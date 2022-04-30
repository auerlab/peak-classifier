#!/bin/sh -e

if [ $0 != ./test.sh ]; then
    printf "Must be run as ./test.sh"
    exit 1
fi

gff=../$(../Utils/gff-name.sh)
if [ ! -e $gff ]; then
    (cd .. && Utils/get-gff.sh)
fi

cd ..
make clean all
cd Test

# Use cave-man installed libs if available
export LD_LIBRARY_PATH=../../local/lib:/usr/lib64:/usr/lib

rm -f ../Mus_musculus.GRCm38.100-augmented*

printf "\n1-base overlaps:\n\n"
../peak-classifier test.bed.xz $gff \
    test-overlaps.tsv
../filter-overlaps test-overlaps.tsv test-filtered.tsv \
    five_prime_utr three_prime_utr intron exon \
    upstream1000 upstream10000 upstream100000 upstream200000 upstream300000 \
    upstream400000 upstream500000 upstream600000 upstream700000 upstream800000 upstream-beyond

printf "\n20%% peak overlaps:\n\n"
../peak-classifier --min-peak-overlap 0.2 test.bed.xz \
    $gff test-peak-20-overlaps.tsv
../filter-overlaps test-peak-20-overlaps.tsv test-peak-20-filtered.tsv \
    five_prime_utr three_prime_utr intron exon \
    upstream1000 upstream10000 upstream100000 upstream200000 upstream300000 \
    upstream400000 upstream500000 upstream600000 upstream700000 upstream800000 upstream-beyond

printf "\n20%% GFF feature overlaps:\n\n"
../peak-classifier --min-gff-overlap 0.2 test.bed.xz \
    $gff test-gff-20-overlaps.tsv
../filter-overlaps test-gff-20-overlaps.tsv test-gff-20-filtered.tsv \
    five_prime_utr three_prime_utr intron exon \
    upstream1000 upstream10000 upstream100000 upstream200000 upstream300000 \
    upstream400000 upstream500000 upstream600000 upstream700000 upstream800000 upstream-beyond

printf "\n20%% either peak or GFF feature overlaps:\n\n"
../peak-classifier --min-gff-overlap 0.2 --min-gff-overlap 0.2 \
    --min-either-overlap test.bed.xz \
    $gff test-either-20-overlaps.tsv
../filter-overlaps test-either-20-overlaps.tsv test-either-20-filtered.tsv \
    five_prime_utr three_prime_utr intron exon \
    upstream1000 upstream10000 upstream100000 upstream200000 upstream300000 \
    upstream400000 upstream500000 upstream600000 upstream700000 upstream800000 upstream-beyond

printf "\nMidpoints only:\n\n"
../peak-classifier --midpoints test.bed.xz $gff \
    test-midpoint-overlaps.tsv
../filter-overlaps test-midpoint-overlaps.tsv test-filtered.tsv \
    five_prime_utr three_prime_utr intron exon \
    upstream1000 upstream10000 upstream100000 upstream200000 upstream300000 \
    upstream400000 upstream500000 upstream600000 upstream700000 upstream800000 upstream-beyond
