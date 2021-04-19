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
#include <stdbool.h>
#include <ctype.h>
#include <biostring.h>
#include <bedio.h>
#include <gffio.h>
#include <plist.h>
#include "peak-classifier.h"

int     main(int argc,char *argv[])

{
    int     c,
	    status;
    FILE    *peak_stream,
	    *gff_stream;
	    // Default, override with --upstream-boundaries
    char    *upstream_boundaries = "1000,10000,100000",
	    *p;
    
    if ( argc < 3 )
	usage(argv);
    
    /* Process flags */
    for (c = 1; (c < argc) && (memcmp(argv[c],"--",2) == 0);
	 ++c)
    {
	if ( strcmp(argv[c], "--upstream-boundaries") == 0 )
	{
	    upstream_boundaries = argv[++c];
	    for (p = upstream_boundaries; *p != '\0'; ++p)
		if ( !isdigit(*p) && (*p != ',') )
		{
		    fputs("peak-classifier: List should be comma-separated with no space.\n", stderr);
		    usage(argv);
		}
	}
	else
	    usage(argv);
    }

    fprintf(stderr, "bed file = %s\n", argv[c]);
    if ( strcmp(argv[c], "-") == 0 )
	peak_stream = stdin;
    else
	if ( (peak_stream = bio_fopen(argv[c], "r")) == NULL )
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

    if ( (status = filter_gff(gff_stream, upstream_boundaries)) == EX_OK )
    {
	status = EX_OK;
	//puts("Sorting...");
	//system("sort -n -k 1 -k 2 -k 3 gff-filtered.bed > gff-sorted.bed");
	//classify(peak_stream, feature_stream);
    }
    else
	fprintf(stderr, "peak-classifier: Error filtering GFF: %s\n",
		strerror(errno));
    bio_fclose(peak_stream);
    bio_fclose(gff_stream);
    return status;
}


/***************************************************************************
 *  Description:
 *      Filter the GFF file and insert explicit intron and upstream
 *      (promoter) regions, returning a FILE pointer to a BED file
 *      containing all features of interest.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-15  Jason Bacon Begin
 ***************************************************************************/

int     filter_gff(FILE *gff_stream, const char *upstream_boundaries)

{
    FILE            *bed_stream;
    bed_feature_t   bed_feature;
    gff_feature_t   gff_feature;
    char            *feature,
		    *strand;
    plist_t         plist = PLIST_INIT;
    //size_t          c;
    
    if ( (bed_stream = fopen("gff-filtered.bed", "w")) == NULL )
    {
	fprintf(stderr, "peak-classifier: Cannot write temp GFF: %s\n",
		strerror(errno));
	return EX_CANTCREAT;
    }
    fprintf(bed_stream, "#CHROM\tFirst\tLast+1\tStrand+Feature\n");
    
    //fprintf(stderr, "Parsing %s\n", upstream_boundaries);
    plist_from_csv(&plist, upstream_boundaries, MAX_UPSTREAM_BOUNDARIES);
    // Upstream features are 1 to first pos, first + 1 to second, etc.
    plist_add_position(&plist, 0);
    plist_sort(&plist, PLIST_ASCENDING);
    //for (c = 0; c < PLIST_COUNT(&plist); ++c)
    //    printf("%" PRIu64 "\n", PLIST_POSITIONS(&plist, c));

    // Write all of the first 4 fields to the feature file
    bed_set_fields(&bed_feature, 4);
    
    puts("Filtering...");
    gff_skip_header(gff_stream);
    while ( gff_read_feature(gff_stream, &gff_feature) == BIO_READ_OK )
    {
	feature = GFF_FEATURE(&gff_feature);
	if ( strcmp(feature, "###") == 0 )
	    fputs("###\n", bed_stream);
	else if ( strstr(feature, "gene") != NULL )
	{
	    // gff_plot_exons(gff_stream, &gff_feature);
	    // Write out upstream regions for likely regulatory elements
	    strand = GFF_STRAND(&gff_feature);
	    gff_to_bed(&bed_feature, &gff_feature);
	    
	    if ( *strand == '+' )
	    {
		generate_upstream_features(bed_stream, &gff_feature, &plist);
		bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
	    }
	    gff_process_subfeatures(gff_stream, bed_stream, &gff_feature);
	    if ( *strand == '-' )
	    {
		bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
		generate_upstream_features(bed_stream, &gff_feature, &plist);
	    }
	    fputs("###\n", bed_stream);
	}
	else if ( strcmp(feature, "chromosome") != 0 )
	{
	    gff_to_bed(&bed_feature, &gff_feature);
	    bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
	    fputs("###\n", bed_stream);
	}
    }
    fclose(bed_stream);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Process sub-features of a gene
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-19  Jason Bacon Begin
 ***************************************************************************/

void    gff_process_subfeatures(FILE *gff_stream, FILE *bed_stream,
				gff_feature_t *gene_feature)

{
    gff_feature_t   sub_feature;
    bed_feature_t   bed_feature;
    bool            first_exon = true,
		    exon,
		    utr5;
    uint64_t        intron_start,
		    intron_end;
    char            *feature,
		    *strand,
		    name[BED_NAME_MAX_CHARS + 1];
    
    bed_set_fields(&bed_feature, 4);
    while ( (gff_read_feature(gff_stream, &sub_feature) == BIO_READ_OK) &&
	    (strcmp(sub_feature.feature, "###") != 0) )
    {
	feature = GFF_FEATURE(&sub_feature);
	exon = (strcmp(feature, "exon") == 0);
	utr5 = (strcmp(feature, "five_prime_UTR") == 0);
	strand = GFF_STRAND(&sub_feature);
	
	// Generate introns between exons
	if ( exon )
	{
	    if ( !first_exon )
	    {
		intron_end = GFF_START_POS(&sub_feature) - 1;
		bed_set_chromosome(&bed_feature, GFF_SEQUENCE(&sub_feature));
		/*
		 *  BED start is 0-based and inclusive
		 *  GFF is 1-based and inclusive
		 */
		bed_set_start_pos(&bed_feature, intron_start);
		/*
		 *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
		 *  GFF is the same
		 */
		bed_set_end_pos(&bed_feature, intron_end);
		snprintf(name, BED_NAME_MAX_CHARS, "%sintron", strand);
		bed_set_name(&bed_feature, name);
		bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
	    }
	    
	    intron_start = GFF_END_POS(&sub_feature);
	    first_exon = false;
	}
	if ( exon || utr5 )
	{
	    gff_to_bed(&bed_feature, &sub_feature);
	    bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
	}
    }
}


/***************************************************************************
 *  Description:
 *      Copy GFF fields to a BED structure
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-19  Jason Bacon Begin
 ***************************************************************************/

void    gff_to_bed(bed_feature_t *bed_feature, gff_feature_t *gff_feature)

{
    char    name[BED_NAME_MAX_CHARS + 1],
	    *strand = GFF_STRAND(gff_feature);
    
    bed_set_chromosome(bed_feature, GFF_SEQUENCE(gff_feature));
    /*
     *  BED start is 0-based and inclusive
     *  GFF is 1-based and inclusive
     */
    bed_set_start_pos(bed_feature, GFF_START_POS(gff_feature) - 1);
    /*
     *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
     *  GFF is the same
     */
    bed_set_end_pos(bed_feature, GFF_END_POS(gff_feature));
    snprintf(name, BED_NAME_MAX_CHARS, "%s%s",
	     strand, GFF_FEATURE(gff_feature));
    bed_set_name(bed_feature, name);
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

int     classify(FILE *peak_stream, FILE *gff_stream)

{
    uint64_t        last_pos = 0,
		    peak_center;
    char            last_chrom[BIO_CHROMOSOME_MAX_CHARS + 1] = "",
		    *feature;
    bed_feature_t   bed_feature;
    gff_feature_t   gff_feature;
    int             dist,
		    gff_status;
    unsigned long   total_bed_features = 0,
		    feature_overlaps = 0,
		    inside_feature = 0,
		    outside_feature = 0;
    bio_overlap_t   overlap;

    /*
     *  Convert GFF to a bed file including potential promoter regions and
     *  introns, which are not listed in the GFF.
     */
    
    bed_skip_header(peak_stream);
    while ( bed_read_feature(peak_stream, &bed_feature) == BIO_READ_OK )
    {
	bed_check_order(&bed_feature, last_chrom, last_pos);
	
	/* Skip GFF features before this bed feature */
	while ( ((gff_status = gff_read_feature(gff_stream, &gff_feature))
		    == BIO_READ_OK) &&
		(dist=bed_gff_cmp(&bed_feature, &gff_feature, &overlap)) > 0 )
	{
	}
	
	/*
	 *  Inject possible promoter regions in ahead of "gene" features
	 *  in the GFF
	 */
	if ( strcmp(GFF_FEATURE(&gff_feature), "gene") == 0 )
	{
	}
	
	/*
	 *  Check for explicitly listed features
	 */
	if ( dist == 0 )
	{
	    feature = GFF_FEATURE(&gff_feature);
	    if ( (strcmp(feature, "five_prime_UTR") == 0) ||
		 (strcmp(feature, "lnc_RNA") == 0) ||
		 (strcmp(feature, "miRNA") == 0) ||
		 (strcmp(feature, "ncRNA") == 0) ||
		 (strcmp(feature, "exon") == 0) )
	    {
		puts("===");
		bed_write_feature(stdout, &bed_feature, BED_FIELD_ALL);
		gff_write_feature(stdout, &gff_feature, BED_FIELD_ALL);
		bio_print_overlap(&overlap, "peak", feature);
		peak_center = (BED_END_POS(&bed_feature) +
				BED_START_POS(&bed_feature)) / 2;
		printf("Peak center     : %" PRIu64 , peak_center);
		if ( (peak_center + 1 >= GFF_START_POS(&gff_feature)) &&
		     (peak_center <= GFF_END_POS(&gff_feature)) )
		{
		    printf(" (inside %s)\n", feature);
		    ++inside_feature;
		}
		else
		{
		    printf(" (outside %s)\n", feature);
		    ++outside_feature;
		}
		++feature_overlaps;
	    }
	}
	++total_bed_features;
    }
    printf("\nTotal peaks:      %lu\n"
	   "Feature overlaps: %lu\n"
	   "Center inside:    %lu\n"
	   "Center outside:   %lu\n",
	   total_bed_features, feature_overlaps, inside_feature, outside_feature);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Check regions upstream of GFF gene TSS for overlap with BED feature
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-10  Jason Bacon Begin
 ***************************************************************************/

void    check_promoter(bed_feature_t *bed_feature, gff_feature_t *gff_feature,
			uint64_t upstream_dist)

{
    uint64_t        promoter_start;
    gff_feature_t   gff_promoter;
    char            promoter_name[GFF_FEATURE_MAX_CHARS + 1];
    bio_overlap_t   overlap;
    
    /* Build a fake GFF feature to represent the upstream region */
    promoter_start = GFF_START_POS(gff_feature) - upstream_dist;
    GFF_SET_START_POS(&gff_promoter, promoter_start);
    GFF_SET_END_POS(&gff_promoter, GFF_START_POS(gff_feature) - 1);
    puts("======");
    printf("Peak at: %" PRIu64 " to %" PRIu64 "\n",
	    BED_START_POS(bed_feature), BED_END_POS(bed_feature));
    printf("Gene at: %" PRIu64 " to %" PRIu64 "\n",
	    GFF_START_POS(gff_feature), GFF_END_POS(gff_feature));
    printf("Trying promoter %" PRIu64 " to %" PRIu64 "\n",
	    GFF_START_POS(&gff_promoter), GFF_END_POS(&gff_promoter));
    if ( bed_gff_cmp(bed_feature, &gff_promoter, &overlap) == 0 )
    {
	GFF_SET_FEATURE(&gff_promoter, promoter_name);
	snprintf(promoter_name, GFF_FEATURE_MAX_CHARS,
		"promoter%" PRIu64, upstream_dist);
	puts("======");
	bed_write_feature(stdout, bed_feature, BED_FIELD_ALL);
	gff_write_feature(stdout, &gff_promoter, BED_FIELD_ALL);
    }
}


/***************************************************************************
 *  Description:
 *      Generate upstream region features from a GFF feature and a list
 *      of upstream distances
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    generate_upstream_features(FILE *feature_stream,
				   gff_feature_t *gff_feature, plist_t *plist)

{
    bed_feature_t   bed_feature[MAX_UPSTREAM_BOUNDARIES];
    char            *strand,
		    name[BED_NAME_MAX_CHARS + 1];
    int             c;
    
    strand = GFF_STRAND(gff_feature);

    for (c = 0; c < PLIST_COUNT(plist) - 1; ++c)
    {
	bed_set_fields(&bed_feature[c], 4);
	bed_set_chromosome(&bed_feature[c], GFF_SEQUENCE(gff_feature));
	/*
	 *  BED start is 0-based and inclusive
	 *  GFF is 1-based and inclusive
	 *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
	 *  GFF is the same
	 */
	if ( *strand == '+' )
	{
	    bed_set_start_pos(&bed_feature[c],
			      GFF_START_POS(gff_feature) - 
			      PLIST_POSITIONS(plist, c + 1) - 1);
	    bed_set_end_pos(&bed_feature[c],
			    GFF_START_POS(gff_feature) -
			    PLIST_POSITIONS(plist, c) - 1);
	}
	else
	{
	    bed_set_start_pos(&bed_feature[c],
			      GFF_END_POS(gff_feature) +
			      PLIST_POSITIONS(plist, c));
	    bed_set_end_pos(&bed_feature[c],
			    GFF_END_POS(gff_feature) + 
			    PLIST_POSITIONS(plist, c + 1));
	}
	
	snprintf(name, BED_NAME_MAX_CHARS, "%supstream%" PRIu64,
		 strand, PLIST_POSITIONS(plist, c + 1));
	bed_set_name(&bed_feature[c], name);
    }
    
    if ( *strand == '-' )
    {
	for (c = 0; c < PLIST_COUNT(plist) - 1; ++c)
	    bed_write_feature(feature_stream, &bed_feature[c], BED_FIELD_ALL);
    }
    else
    {
	for (c = PLIST_COUNT(plist) - 2; c >= 0; --c)
	    bed_write_feature(feature_stream, &bed_feature[c], BED_FIELD_ALL);
    }
}


/***************************************************************************
 *  Description:
 *      Plot exons in the given gene
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-18  Jason Bacon Begin
 ***************************************************************************/

void    gff_plot_exons(FILE *gff_stream, gff_feature_t *gene)

{
    gff_feature_t   sub_feature;
    uint64_t        gene_len = gene->end_pos - gene->start_pos,
		    start,
		    end;
    int             line_ch,
		    c;
    long            file_pos;
    
    file_pos = ftell(gff_stream);
    printf("Gene: %" PRIu64 ", %" PRIu64 "\n", gene->start_pos, gene->end_pos);
    while ( (gff_read_feature(gff_stream, &sub_feature) != BIO_READ_EOF) &&
	    (sub_feature.start_pos <= gene->end_pos) )
    {
	if ( strcmp(sub_feature.feature, "exon") == 0 )
	{
	    /*printf("Exon: %" PRIu64 ", %" PRIu64 "\n",
		    sub_feature.start_pos, sub_feature.end_pos);*/
	    start = (sub_feature.start_pos - gene->start_pos) * 78 / gene_len;
	    end = (sub_feature.end_pos - gene->start_pos) * 78 / gene_len + 1;
	    line_ch = *sub_feature.strand == '+'? '>' : '<';
	    //printf("%" PRIu64 ", %" PRIu64 "\n", start, end);
	    for (c = 0; c < start; ++c)
		putchar(line_ch);
	    while ( c++ < end )
		putchar('E');
	    while ( c++ < 78 )
		putchar(line_ch);
	    putchar('\n');
	}
    }
    fseek(gff_stream, file_pos, SEEK_SET);
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s [--upstream-boundaries pos[,pos ...]] BED-file GFF-file\n", argv[0]);
    exit(EX_USAGE);
}
