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
#include <unistd.h>
#include <assert.h>
#include <xtend.h>
#include <biolibc/bed.h>
#include <biolibc/gff.h>
#include <biolibc/plist.h>
#include "peak-classifier.h"

int     main(int argc,char *argv[])

{
    int     c,
	    status;
    double  min_peak_overlap = 1.0e-9,
	    min_gff_overlap = 1.0e-9;
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
	    *min_overlap_flags = "",
	    *end,
	    *gff_stem,
	    augmented_filename[PATH_MAX + 1],
	    sorted_filename[PATH_MAX + 1];
    bool    midpoints_only = false;
    bed_feature_t   bed_feature;
    struct stat     file_info;
    
    if ( argc < 4 )
	usage(argv);
    
    /* Process flags */
    for (c = 1; (c < argc) && (memcmp(argv[c],"--",2) == 0); ++c)
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
	else if ( strcmp(argv[c], "--min-peak-overlap") == 0 )
	{
	    min_peak_overlap = strtod(argv[++c], &end);
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[c], "--min-gff-overlap") == 0 )
	{
	    min_gff_overlap = strtod(argv[++c], &end);
	    if ( *end != '\0' )
		usage(argv);
	}
	else if ( strcmp(argv[c], "--min-either-overlap") == 0 )
	    min_overlap_flags = "-e";
	else if ( strcmp(argv[c], "--midpoints") == 0 )
	    midpoints_only = true;
	else
	    usage(argv);
    }

    if ( strcmp(argv[c], "-") == 0 )
	peak_stream = stdin;
    else
    {
	assert(valid_extension(argv[c], ".bed"));
	if ( (peak_stream = xt_fopen(argv[c], "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
    }
    
    if ( strcmp(argv[++c], "-") == 0 )
    {
	gff_stream = stdin;
	gff_stem = "unknown-stdin-gff";
    }
    else
    {
	assert(valid_extension(argv[c], ".gff3"));
	if ( (gff_stream = xt_fopen(argv[c], "r")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[c],
		    strerror(errno));
	    exit(EX_NOINPUT);
	}
	gff_stem = argv[c];
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
	assert(valid_extension(overlaps_filename, ".tsv"));
	redirect_overwrite = " > ";
	redirect_append = " >> ";
    }

    // Already verified .gff3[.*z] extension above
    *strstr(gff_stem, ".gff3") = '\0';
    snprintf(augmented_filename, PATH_MAX, "%s-augmented.bed", gff_stem);
    if ( stat(augmented_filename, &file_info) == 0 )
	fprintf(stderr, "Using existing %s...\n", augmented_filename);
    else if ( gff_augment(gff_stream, upstream_boundaries, augmented_filename) != EX_OK )
    {
	fprintf(stderr, "gff_augment() failed.  Removing %s...\n",
		augmented_filename);
	unlink(augmented_filename);
	exit(EX_DATAERR);
    }
    
    snprintf(sorted_filename, PATH_MAX, "%s-augmented+sorted.bed", gff_stem);
    if ( stat(sorted_filename, &file_info) == 0 )
	fprintf(stderr, "Using existing %s...\n", sorted_filename);
    else
    {
	// LC_ALL=C makes sort assume 1 byte/char, which improves speed
	snprintf(cmd, CMD_MAX, "env LC_ALL=C grep -v '^#' %s | "
		"sort -n -k 1 -k 2 -k 3 > %s\n",
		augmented_filename, sorted_filename);
	fputs("Sorting...\n", stderr);
	if ( (status = system(cmd)) != 0 )
	{
	    fprintf(stderr, "Sort failed.  Removing %s...\n", sorted_filename);
	    unlink(sorted_filename);
	    exit(EX_DATAERR);
	}
    }
    
    fputs("Finding intersects...\n", stderr);
    snprintf(cmd, CMD_MAX,
	    "printf '#Chr\tP-start\tP-end\tF-start\tF-end\tF-name\tStrand\tOverlap\n'%s%s",
	    redirect_overwrite, overlaps_filename);
    if ( (status = system(cmd)) == 0 )
    {
	/*
	 *  Peaks not overlapping anything else are labeled
	 *  upstream-beyond.  The entire peak length must overlap the
	 *  beyond region since none of it overlaps anything else.
	 */
	 
	// Insert "set -x; " for debugging
	snprintf(cmd, CMD_MAX,
		 "bedtools intersect -a - -b %s -f %g -F %g %s -wao"
		 "| awk 'BEGIN { OFS=IFS; } { if ( $8 == -1 ) "
		    "$9 = \"upstream-beyond\"; $12 = $3 - $2; "
		    "printf(\"%%s\\t%%d\\t%%d\\t%%d\\t%%d\\t"
		    "%%s\\t%%s\\t%%s\\n\", "
		    "$1, $2, $3, $7, $8, $9, $11, $12); }' %s%s\n",
		sorted_filename, min_peak_overlap, min_gff_overlap,
		min_overlap_flags, redirect_append, overlaps_filename);

	if ( (intersect_pipe = popen(cmd, "w")) == NULL )
	{
	    fprintf(stderr, "%s: Cannot pipe data to bedtools intersect.\n",
		    argv[0]);
	    return EX_CANTCREAT;
	}
	while ( bed_read_feature(peak_stream, &bed_feature, BED_FIELD_ALL) != EOF )
	{
	    if ( midpoints_only )
	    {
		// Replace peak start/end with midpoint coordinates
		bed_set_start_pos(&bed_feature,
		    (BED_START_POS(&bed_feature) + BED_END_POS(&bed_feature))
		    / 2);
		bed_set_end_pos(&bed_feature, BED_START_POS(&bed_feature) + 1);
	    }
	    bed_write_feature(intersect_pipe, &bed_feature,
			      BED_FIELD_ALL);
	}
	pclose(intersect_pipe);
    }
    xt_fclose(peak_stream);
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

int     gff_augment(FILE *gff_stream, const char *upstream_boundaries,
		    const char *augmented_filename)

{
    FILE            *bed_stream;
    bed_feature_t   bed_feature = BED_INIT;
    gff_feature_t   gff_feature = GFF_INIT;
    char            *feature,
		    strand;
    plist_t         plist = PLIST_INIT;
    
    if ( (bed_stream = fopen(augmented_filename, "w")) == NULL )
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
    
    fputs("Augmenting GFF3 data...\n", stderr);
    gff_skip_header(gff_stream);
    while ( gff_read_feature(gff_stream, &gff_feature, GFF_FIELD_ALL) == BIO_READ_OK )
    {
	// FIXME: Create a --autosomes-only flag to activate this check
	if ( strisint(GFF_SEQUENCE(&gff_feature), 10) )
	{
	    feature = GFF_NAME(&gff_feature);
	    // FIXME: Rely on parent IDs instead of ###?
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
    xt_fclose(gff_stream);
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
		    exon;
    uint64_t        intron_start = 0,   // Silence bogus warning from GCC
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
    
    while ( (gff_read_feature(gff_stream, &subfeature, GFF_FIELD_ALL) == BIO_READ_OK) &&
	    (strcmp(subfeature.name, "###") != 0) )
    {
	feature = GFF_NAME(&subfeature);
	exon = (strcmp(feature, "exon") == 0);

	// mRNA or lnc_RNA mark the start of a new set of exons
	if ( (strstr(subfeature.name, "RNA") != NULL) ||
	     (strstr(subfeature.name, "transcript") != NULL) ||
	     (strstr(subfeature.name, "gene_segment") != NULL) ||
	     (strstr(subfeature.name, "_overlapping_ncrna") != NULL) )
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
	    "\nUsage: %s [--upstream-boundaries pos[,pos ...]] "
	    "[--min-peak-overlap x.y] [--min-gff-overlap x.y] [--midpoints] "
	    "peaks.bed features.gff3 overlaps.tsv\n\n", argv[0]);
    fputs("Upstream boundaries are distances upstream from TSS, for which we want\n"
	  "overlaps reported.  The default is 1000,10000,100000, which means features\n"
	  "are generated for 1 to 1000, 1001 to 10000, and 10001 to 100000 bases\n"
	  "upstream.  Peaks that do not overlap any of these or other features are\n"
	  "reported as 'upstream-beyond.\n\n"
	  "The minimum peak/gff overlap must range from 1.0e-9 (the default, which\n"
	  "corresponds to a single base) to 1.0. These values are passed directlry to\n"
	  "bedtools intersect -f/-F.\n"
	  "They must be used with great caution since the size of peaks and GFF\n"
	  "features varies greatly.\n\n"
	  "--min-either-overlap indicates that either the minimum peak or the minimum\n"
	  "GFF feature overlap satisfies the overlap requirement.  Otherwise, both\n"
	  "overlap requirements must be met.\n\n"
	  "--midpoints indicates that we are only interested in which feature contains\n"
	  "the midpoint of each peak.  This is the same as --min-peak-overlap 0.5\n"
	  "in cases where half the peak is contained in a feature, but can also report\n"
	  "overlaps with features too small to contain this much overlap.\n\n", stderr);
    exit(EX_USAGE);
}
