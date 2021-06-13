#!/bin/sh -e

gff=$(Utils/gff-name.sh)
release=$(echo $gff | cut -d . -f 3)

if which fetch; then
    fetch=fetch;
elif which curl; then
    fetch='curl -O';
elif which wget; then
    fetch=wget;
else
    printf "$0: No fetch program found, aborting.\n"
    exit 1
fi
if [ ! -e $gff ]; then
    url=http://ftp.ensembl.org/pub/release-$release/gff3/mus_musculus/$gff
    $fetch $url
else
    printf "$gff already exists.  Remove it and rerun $0 to force download.\n"
fi
