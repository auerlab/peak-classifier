#!/bin/sh -e

printf "Extracting...\n"
awk '$3 ~ "gene" && $1 ~ /^[0-9]+$/ { printf("%s\t%s\t%s\t%s\t%s\n",
		    $1, $4, $5, $3, $7); }' Mus*.gff3 > mus.bed

printf "Sorting...\n"
bedtools sort -i mus.bed > mus-sorted.bed

printf "Merging...\n"
bedtools merge -i mus-sorted.bed > mus-merged.bed

printf "Diffing...\n"
wc mus-sorted.bed mus-merged.bed
diff mus-sorted.bed mus-merged.bed | more
