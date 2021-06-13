#!/bin/sh -e

gzcat $(Utils/gff-name.sh) \
    | awk '{ if ( $0 ~ "^#" ) print $0; else printf("%s\t%u\t%u\t%u\t%s(%s)\n",
						$1, $4, $5, $5 - $4, $3, $7); }'

