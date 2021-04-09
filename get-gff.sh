#!/bin/sh -e

if [ -z $(ls Mus_musculus.GRC*.gff3.gz) ]; then
    printf "No GFF found.\n"
    release=103
    file=Mus_musculus.GRCm39.$release.gff3.gz
    url=http://ftp.ensembl.org/pub/release-103/gff3/mus_musculus/$file
    fetch $url
fi
