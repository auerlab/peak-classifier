#!/bin/sh -e

##########################################################################
#   Script description:
#       List all feature types in a GFF
#       
#   History:
#   Date        Name        Modification
#   2021-04-10  Jason Bacon Begin
##########################################################################

usage()
{
    printf "Usage: $0 file.gff[.gz|bz2|xz]\n"
    exit 1
}


##########################################################################
#   Main
##########################################################################

if [ $# != 1 ]; then
    usage
fi
file=$1
ext=${file#*.gff3.}
tool=$(echo $ext | tr -d 2)

${tool}cat $file | awk '$1 !~ "^#" { print $3 }' | sort -u
