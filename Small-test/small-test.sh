#!/bin/sh -e

##########################################################################
#   Function description:
#       Pause until user presses return
##########################################################################

pause()
{
    local junk
    
    printf "Press return to continue..."
    read junk
}

if [ $0 != ./small-test.sh ]; then
    printf "Must be run as ./small-test.sh"
    exit
fi

rm -f small-test-augmented*.bed small-test-overlaps.tsv
cd ..
make clean all
cd Small-test

printf "Viewing pruned GFF input...\n"
pause
awk '{ if ( $0 ~ "^#" ) print $0; else printf("%s\t%u\t%u\t%u\t%s(%s)\n",
    $1, $4, $5, $5 - $4, $3, $7); }' small-test.gff3 | more

printf "Running peak-classifier...\n"
pause
../peak-classifier small-test.bed small-test.gff3 small-test-overlaps.tsv

printf "Viewing augmented GFF data...\n"
pause
more small-test-augmented.bed

printf "Viewing sorted GFF data...\n"
pause
more small-test-augmented+sorted.bed

printf "Viewing overlaps...\n"
pause
more small-test-overlaps.tsv

printf "Filtering overlaps...\n"
pause
../filter-overlaps small-test-overlaps.tsv small-test-filtered.tsv intron exon
more small-test-filtered.tsv
