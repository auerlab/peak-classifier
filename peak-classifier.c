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
#include <sys/param.h>
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
    uint64_t        last_pos = 0,
		    peak_center;
    char            last_chrom[BIO_CHROMOSOME_MAX_CHARS + 1];
    bed_feature_t   bed_feature;
    gff_feature_t   gff_feature;
    int             dist;
    unsigned long   total_bed_features = 0,
		    exon_overlaps = 0,
		    inside_exon = 0,
		    outside_exon = 0;
    bio_overlap_t   overlap;

    bed_skip_header(bed_stream);
    gff_skip_header(gff_stream);
    while ( bed_read_feature(bed_stream, &bed_feature) == BIO_READ_OK )
    {
	bed_check_order(&bed_feature, last_chrom, last_pos);
	
	/* Skip GFF features before this bed feature */
	while ( (gff_read_feature(gff_stream, &gff_feature) == BIO_READ_OK) 
		&& ((strcmp(gff_feature.feature,"exon") != 0) ||
		    ( (dist=bed_gff_cmp(&bed_feature, &gff_feature,
					&overlap)) > 0)) )
	    ;
	if ( dist == 0 )
	{
	    puts("===");
	    bed_write_feature(stdout, &bed_feature, BED_FIELD_ALL);
	    gff_write_feature(stdout, &gff_feature, BED_FIELD_ALL);
	    bio_print_overlap(&overlap, "Peak", "Exon");
	    peak_center = (BED_END_POS(&bed_feature) +
			    BED_START_POS(&bed_feature)) / 2;
	    printf("Peak center     : %" PRIu64 , peak_center);
	    if ( (peak_center + 1 >= GFF_START_POS(&gff_feature)) &&
		 (peak_center <= GFF_END_POS(&gff_feature)) )
	    {
		puts(" (inside exon)");
		++inside_exon;
	    }
	    else
	    {
		puts(" (outside exon)");
		++outside_exon;
	    }
	    ++exon_overlaps;
	}
	++total_bed_features;
    }
    printf("\nTotal peaks:    %lu\n"
	   "Exon overlaps:  %lu\n"
	   "Center inside:  %lu\n"
	   "Center outside: %lu\n",
	   total_bed_features, exon_overlaps, inside_exon, outside_exon);
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
    if ( chromosome_name_cmp(BED_CHROMOSOME(bed_feature), last_chrom) == 0 )
    {
	if ( BED_START_POS(bed_feature) < last_pos )
	{
	    fprintf(stderr, "peak-classifier: BED file not sorted by start position.\n");
	    exit(EX_DATAERR);
	}
    }
    else if ( chromosome_name_cmp(BED_CHROMOSOME(bed_feature), last_chrom) < 0 )
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


int     bed_gff_cmp(bed_feature_t *bed_feature, gff_feature_t *gff_feature,
		    bio_overlap_t *overlap)

{
    int         chromosome_cmp;
    uint64_t    bed_start, bed_end, bed_len,
		gff_start, gff_end, gff_len;
    
    chromosome_cmp = chromosome_name_cmp(BED_CHROMOSOME(bed_feature),
					 GFF_SEQUENCE(gff_feature));
    if ( chromosome_cmp == 0 )
    {
	/*
	 *  BED positions are 0-based, with end non-inclusive, which can
	 *  also be viewed as an inclusive 1-based coordinate
	 *  GFF is 1-based, both ends inclusive
	 */
	
	if ( BED_END_POS(bed_feature) < GFF_START_POS(gff_feature) )
	{
	    bio_set_overlap(overlap, 0, 0, 0, 0);
	    return -1;
	}
	else if ( BED_START_POS(bed_feature) + 1 > GFF_END_POS(gff_feature) )
	{
	    bio_set_overlap(overlap, 0, 0, 0, 0);
	    return 1;
	}
	else
	{
	    bed_start = BED_START_POS(bed_feature);
	    bed_end = BED_END_POS(bed_feature);
	    gff_start = GFF_START_POS(gff_feature);
	    gff_end = GFF_END_POS(gff_feature);
	    bed_len = bed_end - bed_start;
	    gff_len = gff_end - gff_start + 1;
	    bio_set_overlap(overlap, bed_len, gff_len,
			    MAX(bed_start+1, gff_start),
			    MIN(bed_end, gff_end));
	    return 0;
	}
    }
    return chromosome_cmp;
}
