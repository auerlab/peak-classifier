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
#include <bedio.h>
#include <gffio.h>
#include "peak-classifier.h"

int     main(int argc,char *argv[])

{
    int     c,
	    status;
    FILE    *bed_stream,
	    *gff_stream;
    
    if ( argc < 3 )
	usage(argv);
    
    /* Process flags */
    for (c = 1; (c < argc) && (*argv[c] == '-') && (strcmp(argv[c],"-") != 0);
	 ++c)
	;
    
    if ( strcmp(argv[c], "-") == 0 )
	bed_stream = stdin;
    else
	if ( (bed_stream = bio_fopen(argv[c], "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    
    if ( strcmp(argv[++c], "-") == 0 )
	gff_stream = stdin;
    else
	if ( (gff_stream = bio_fopen(argv[c], "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    
    status = classify(bed_stream, gff_stream);
    bio_fclose(bed_stream);
    bio_fclose(gff_stream);
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


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s BED-file GFF-file\n", argv[0]);
    exit(EX_USAGE);
}
