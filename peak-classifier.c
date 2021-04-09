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
#include <biostring.h>
#include <stdbool.h>
#include "peak-classifier.h"

int     main(int argc,char *argv[])

{
    int     c,
	    status;
    FILE    *bed_stream,
	    *gff_stream;
    bool    bed_is_pipe = false,
	    gff_is_pipe = false;
    
    if ( argc < 3 )
	usage(argv);
    
    /* Process flags */
    for (c = 1; (c < argc) && (*argv[c] == '-') && (strcmp(argv[c],"-") != 0);
	 ++c)
	;
    
    if ( strcmp(argv[c], "-") == 0 )
	bed_stream = stdin;
    else
	if ( (bed_stream = bio_open(argv[c], "r", &bed_is_pipe)) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    
    if ( strcmp(argv[++c], "-") == 0 )
	gff_stream = stdin;
    else
	if ( (gff_stream = bio_open(argv[c], "r", &gff_is_pipe)) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    
    status = classify(bed_stream, gff_stream);
    if ( bed_is_pipe )
	pclose(bed_stream);
    else
	fclose(bed_stream);
    if ( gff_is_pipe )
	pclose(gff_stream);
    else
	fclose(gff_stream);
    return status;
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
    uint64_t    last_pos = 0;
    char        last_chrom[BIO_CHROMOSOME_MAX_CHARS + 1];
    bed_feature_t   bed_feature;

    bed_skip_header(bed_stream);
    while ( bed_read_feature(bed_stream, &bed_feature) == BIO_READ_OK )
    {
	bed_check_order(&bed_feature, last_chrom, last_pos);
	bed_write_feature(stdout, &bed_feature, BED_FIELD_ALL);
    }
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Make sure the BED input is sorted
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-08  Jason Bacon Begin
 ***************************************************************************/

void    bed_check_order(bed_feature_t *bed_feature, char last_chrom[],
			uint64_t last_pos)

{
    if ( chromosome_name_cmp(bed_feature->chromosome, last_chrom) == 0 )
    {
	if ( bed_feature->start_pos < last_pos )
	{
	    fprintf(stderr, "peak-classifier: BED file not sorted by start position.\n");
	    exit(EX_DATAERR);
	}
    }
    else if ( chromosome_name_cmp(bed_feature->chromosome, last_chrom) < 0 )
    {
	fprintf(stderr, "peak-classifier: BED file not sorted by start chromosome.\n");
	exit(EX_DATAERR);
    }
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s BED-file GFF-file\n", argv[0]);
    exit(EX_USAGE);
}


FILE    *bio_open(char *filename, char *mode, bool *is_pipe)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[CMD_MAX + 1];
    
    if ( (strcmp(mode, "r") != 0 ) && (strcmp(mode, "w") != 0) )
    {
	fprintf(stderr, "bio_open(): Only \"r\" and \"w\" modes supported.\n");
	return NULL;
    }
    
    if ( ext == NULL )
    {
	fprintf(stderr, "bio_open(): No filename extension on %s.\n", filename);
	return NULL;
    }

    if ( *mode == 'r' )
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, CMD_MAX, "gzcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, CMD_MAX, "bzcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, CMD_MAX, "xzcat %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
    else    // "w"
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, CMD_MAX, "gzip -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, CMD_MAX, "bzip2 -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, CMD_MAX, "xz -c > %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
}
