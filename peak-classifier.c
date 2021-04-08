/***************************************************************************
 *  Description:
 *      Classify peaks in a bed file according to features found in a GFF.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "peak-classifier.h"

int     main(int argc,char *argv[])

{
    int     c;
    FILE    *bed_stream,
	    *gff_stream;
    
    if ( argc < 3 )
	usage(argv);
    
    /* Process flags */
    for (c = 1; c < argc; ++c)
    {
	if ( (*argv[c] == '-') && (strcmp(argv[c],"-") != 0) )
	    ;
    }
    
    if ( strcmp(argv[c], "-") == 0 )
	bed_stream = stdin;
    else
	if ( (bed_stream = fopen(argv[c], "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    
    if ( strcmp(argv[++c], "-") == 0 )
	gff_stream = stdin;
    else
	if ( (gff_stream = fopen(argv[c], "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    
    return classify(bed_stream, gff_stream);
}


/***************************************************************************
 *  Description:
 *      Read through BED and GFF files, classifying each peak in the BED
 *      file according to overlapping features in the GFF.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     classify(FILE *bed_stream, FILE *gff_stream)

{
    bed_feature_t   bed_feature;

    bed_skip_header(bed_stream);
    while ( bed_read_feature(bed_stream, &bed_feature) == BIO_READ_OK )
    {
	printf("%s\n", bed_feature.chromosome);
    }
    return EX_OK;
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s BED-file GFF-file\n", argv[0]);
    exit(EX_USAGE);
}

