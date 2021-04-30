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
#include <limits.h>
#include <sys/stat.h>
#include <biostring.h>
#include <bedio.h>
#include <gffio.h>
#include <plist.h>
#include "peak-classifier.h"

int     main(int argc,char *argv[])

{
    int     c,
	    ch,
	    status;
    FILE    *peak_stream,
	    *gff_stream,
	    *intersect_pipe;
	    // Default, override with --upstream-boundaries
    char    *upstream_boundaries = "1000,10000,100000",
	    *p,
	    cmd[CMD_MAX + 1],
	    *redirect_overwrite,
	    *redirect_append,
	    *overlaps_filename,
	    *sort;
    
    if ( argc < 4 )
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

    if ( strcmp(argv[++c], "-") == 0 )
    {
	overlaps_filename = "";
	redirect_overwrite = "";
	redirect_append = "";
    }
    else
    {
	overlaps_filename = argv[c];
	redirect_overwrite = " > ";
	redirect_append = " >> ";
    }

    if ( (status = gff_augment(gff_stream, upstream_boundaries)) == EX_OK )
    {
	// GNU sort is much faster and can use threads
	if ( system("which gsort") == 0 )
	    sort = "gsort";
	else
	    sort = "sort";
	snprintf(cmd, CMD_MAX, "grep -v '^#' pc-gff-augmented.bed | "
		"%s -S 1G -n -k 1 -k 2 -k 3 > pc-gff-sorted.bed", sort);
	fprintf(stderr, "Sorting with %s...\n", sort);
	system(cmd);
    
	fputs("Finding intersects...\n", stderr);
	snprintf(cmd, CMD_MAX,
		"printf '#Chr\tP-start\tP-end\tF-start\tF-end\tF-name\tStrand\tOverlap\n'%s%s",
		redirect_overwrite, overlaps_filename);
	system(cmd);
	
	snprintf(cmd, CMD_MAX,
		 "bedtools intersect -a - -b pc-gff-sorted.bed -wao"
		 "| cut -f 1,2,3,7,8,9,11,12%s%s",
		 redirect_append, overlaps_filename);

	// Alternative to bedtools intersect:
	// classify(peak_stream, feature_stream);
	if ( (intersect_pipe = popen(cmd, "w")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot pipe data to bedtools intersect.\n",
		    argv[0]);
	    return EX_CANTCREAT;
	}
	while ( (ch = getc(peak_stream)) != EOF )
	    putc(ch, intersect_pipe);
	pclose(intersect_pipe);
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

int     gff_augment(FILE *gff_stream, const char *upstream_boundaries)

{
    FILE            *bed_stream;
    bed_feature_t   bed_feature = BED_INIT;
    gff_feature_t   gff_feature = GFF_INIT;
    char            *feature,
		    strand;
    plist_t         plist = PLIST_INIT;
    
    if ( (bed_stream = fopen("pc-gff-augmented.bed", "w")) == NULL )
    {
	fprintf(stderr, "peak-classifier: Cannot write temp GFF: %s\n",
		strerror(errno));
	return EX_CANTCREAT;
    }
    fprintf(bed_stream, "#CHROM\tFirst\tLast+1\tStrand+Feature\n");
    
    plist_from_csv(&plist, upstream_boundaries, MAX_UPSTREAM_BOUNDARIES);
    // Upstream features are 1 to first pos, first + 1 to second, etc.
    plist_add_position(&plist, 0);
    plist_sort(&plist, PLIST_ASCENDING);

    // Write all of the first 4 fields to the feature file
    bed_set_fields(&bed_feature, 6);
    bed_set_score(&bed_feature, 0);
    
    fputs("Filtering...\n", stderr);
    gff_skip_header(gff_stream);
    while ( gff_read_feature(gff_stream, &gff_feature) == BIO_READ_OK )
    {
	// FIXME: Create a --autosomes-only flag to activate this check
	if ( strisnum(GFF_SEQUENCE(&gff_feature)) )
	{
	    feature = GFF_NAME(&gff_feature);
	    if ( strcmp(feature, "###") == 0 )
		fputs("###\n", bed_stream);
	    else if ( strstr(feature, "gene") != NULL )
	    {
		// Write out upstream regions for likely regulatory elements
		strand = GFF_STRAND(&gff_feature);
		gff_to_bed(&bed_feature, &gff_feature);
		bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
		
		if ( strand == '+' )
		    generate_upstream_features(bed_stream, &gff_feature, &plist);
		gff_process_subfeatures(gff_stream, bed_stream, &gff_feature);
		if ( strand == '-' )
		    generate_upstream_features(bed_stream, &gff_feature, &plist);
		fputs("###\n", bed_stream);
	    }
	    else if ( strcmp(feature, "chromosome") != 0 )
	    {
		gff_to_bed(&bed_feature, &gff_feature);
		bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
		fputs("###\n", bed_stream);
	    }
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
    gff_feature_t   subfeature = GFF_INIT;
    bed_feature_t   bed_feature = BED_INIT;
    bool            first_exon = true,
		    exon,
		    utr;
    uint64_t        intron_start,
		    intron_end;
    char            *feature,
		    strand,
		    name[BED_NAME_MAX_CHARS + 1];

    bed_set_fields(&bed_feature, 6);
    strand = GFF_STRAND(gene_feature);
    if ( bed_set_strand(&bed_feature, strand) != BIO_DATA_OK )
    {
	fputs("gff_process_subfeatures().\n", stderr);
	exit(EX_DATAERR);
    }
    
    while ( (gff_read_feature(gff_stream, &subfeature) == BIO_READ_OK) &&
	    (strcmp(subfeature.name, "###") != 0) )
    {
	feature = GFF_NAME(&subfeature);
	exon = (strcmp(feature, "exon") == 0);
	utr = (strstr(feature, "UTR") != NULL);

	// mRNA or lnc_RNA mark the start of a new set of exons
	if ( (strstr(subfeature.name, "RNA") != NULL) ||
	     (strstr(subfeature.name, "_transcript") != NULL) ||
	     (strstr(subfeature.name, "gene_segment") != NULL) )
	    first_exon = true;
	
	// Generate introns between exons
	if ( exon )
	{
	    if ( !first_exon )
	    {
		intron_end = GFF_START_POS(&subfeature) - 1;
		bed_set_chromosome(&bed_feature, GFF_SEQUENCE(&subfeature));
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
		snprintf(name, BED_NAME_MAX_CHARS, "intron");
		bed_set_name(&bed_feature, name);
		bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
	    }
	    
	    intron_start = GFF_END_POS(&subfeature);
	    first_exon = false;
	}
	
	gff_to_bed(&bed_feature, &subfeature);
	bed_write_feature(bed_stream, &bed_feature, BED_FIELD_ALL);
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
    char            strand,
		    name[BED_NAME_MAX_CHARS + 1];
    int             c;
    
    strand = GFF_STRAND(gff_feature);

    for (c = 0; c < PLIST_COUNT(plist) - 1; ++c)
    {
	bed_set_fields(&bed_feature[c], 6);
	bed_set_strand(&bed_feature[c], strand);
	bed_set_chromosome(&bed_feature[c], GFF_SEQUENCE(gff_feature));
	/*
	 *  BED start is 0-based and inclusive
	 *  GFF is 1-based and inclusive
	 *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
	 *  GFF is the same
	 */
	if ( strand == '+' )
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
	
	snprintf(name, BED_NAME_MAX_CHARS, "upstream%" PRIu64,
		 PLIST_POSITIONS(plist, c + 1));
	bed_set_name(&bed_feature[c], name);
    }
    
    if ( strand == '-' )
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


void    usage(char *argv[])

{
    fprintf(stderr,
	    "Usage: %s [--upstream-boundaries pos[,pos ...]] "
	    "peaks.bed features.gff overlaps.tsv\n", argv[0]);
    exit(EX_USAGE);
}
