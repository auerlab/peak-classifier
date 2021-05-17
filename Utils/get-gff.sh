#!/bin/sh -e

gff=$(Utils/gff-name)
release=$(echo $gff | cut -d . -f 3)
if [ ! -e $gff ]; then
    url=http://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus/$gff
    fetch $url
else
    printf "$gff already exists.  Remove it and rerun $0 to force download.\n"
fi
