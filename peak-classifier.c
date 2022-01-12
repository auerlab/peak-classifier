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
#include <xtend/string.h>
#include <xtend/file.h>
#include <biolibc/bed.h>
#include <biolibc/gff.h>
#include <biolibc/pos-list.h>
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
	    cmd[PEAK_CMD_MAX + 1],
	    *redirect_overwrite,
	    *redirect_append,
	    *overlaps_filename,
	    *min_overlap_flags = "",
	    *end,
	    *gff_stem,
	    augmented_filename[PATH_MAX + 1],
	    sorted_filename[PATH_MAX + 1];
    bool    midpoints_only = false;
    bl_bed_t   bed_feature;
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
	snprintf(cmd, PEAK_CMD_MAX, "env LC_ALL=C grep -v '^#' %s | "
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
    snprintf(cmd, PEAK_CMD_MAX,
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
	snprintf(cmd, PEAK_CMD_MAX,
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
	while ( bl_bed_read(&bed_feature, peak_stream, BL_BED_FIELD_ALL) != EOF )
	{
	    if ( midpoints_only )
	    {
		// Replace peak start/end with midpoint coordinates
		bl_bed_set_chrom_start(&bed_feature,
		    (BL_BED_CHROM_START(&bed_feature) + BL_BED_CHROM_END(&bed_feature))
		    / 2);
		bl_bed_set_chrom_end(&bed_feature, BL_BED_CHROM_START(&bed_feature) + 1);
	    }
	    bl_bed_write(&bed_feature, intersect_pipe, BL_BED_FIELD_ALL);
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
    FILE        *bed_stream;
    bl_bed_t    bed_feature = BL_BED_INIT;
    bl_gff_t    gff_feature = BL_GFF_INIT;
    char        *feature,
		strand;
    bl_pos_list_t      pos_list = BL_POS_LIST_INIT;
    
    if ( (bed_stream = fopen(augmented_filename, "w")) == NULL )
    {
	fprintf(stderr, "peak-classifier: Cannot write temp GFF: %s\n",
		strerror(errno));
	return EX_CANTCREAT;
    }
    fprintf(bed_stream, "#CHROM\tFirst\tLast+1\tStrand+Feature\n");
    
    bl_pos_list_from_csv(&pos_list, upstream_boundaries, MAX_UPSTREAM_BOUNDARIES);
    // Upstream features are 1 to first pos, first + 1 to second, etc.
    bl_pos_list_add_position(&pos_list, 0);
    bl_pos_list_sort(&pos_list, BL_POS_LIST_ASCENDING);

    // Write all of the first 4 fields to the feature file
    // Done within bl_gff_to_bed() now
    //bl_bed_set_fields(&bed_feature, 6);
    //bl_bed_set_score(&bed_feature, 0);
    
    fputs("Augmenting GFF3 data...\n", stderr);
    bl_gff_skip_header(gff_stream);
    while ( bl_gff_read(&gff_feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK )
    {
	// FIXME: Create a --autosomes-only flag to activate this check
	if ( strisint(BL_GFF_SEQUENCE(&gff_feature), 10) )
	{
	    feature = BL_GFF_FEATURE(&gff_feature);
	    // FIXME: Rely on parent IDs instead of ###?
	    if ( strcmp(feature, "###") == 0 )
		fputs("###\n", bed_stream);
	    else if ( strstr(feature, "gene") != NULL )
	    {
		// Write out upstream regions for likely regulatory elements
		strand = BL_GFF_STRAND(&gff_feature);
		bl_gff_to_bed(&gff_feature, &bed_feature);
		bl_bed_write(&bed_feature, bed_stream, BL_BED_FIELD_ALL);
		
		if ( strand == '+' )
		    generate_upstream_features(bed_stream, &gff_feature, &pos_list);
		gff_process_subfeatures(gff_stream, bed_stream, &gff_feature);
		if ( strand == '-' )
		    generate_upstream_features(bed_stream, &gff_feature, &pos_list);
		fputs("###\n", bed_stream);
	    }
	    else if ( strcmp(feature, "chromosome") != 0 )
	    {
		bl_gff_to_bed(&gff_feature, &bed_feature);
		bl_bed_write(&bed_feature, bed_stream, BL_BED_FIELD_ALL);
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
				bl_gff_t *gene_feature)

{
    bl_gff_t   subfeature = BL_GFF_INIT;
    bl_bed_t   bed_feature = BL_BED_INIT;
    bool            first_exon = true,
		    exon;
    uint64_t        intron_start = 0,   // Silence bogus warning from GCC
		    intron_end;
    char            *feature,
		    strand,
		    name[BL_BED_NAME_MAX_CHARS + 1];

    bl_bed_set_fields(&bed_feature, 6);
    strand = BL_GFF_STRAND(gene_feature);
    if ( bl_bed_set_strand(&bed_feature, strand) != BL_DATA_OK )
    {
	fputs("gff_process_subfeatures(): bl_bed_set_strand() failed..\n", stderr);
	exit(EX_DATAERR);
    }
    
    while ( (bl_gff_read(&subfeature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK) &&
	    (strcmp(BL_GFF_FEATURE(&subfeature), "###") != 0) )
    {
	feature = BL_GFF_FEATURE(&subfeature);
	exon = (strcmp(feature, "exon") == 0);

	// mRNA or lnc_RNA mark the start of a new set of exons
	if ( (strstr(BL_GFF_FEATURE(&subfeature), "RNA") != NULL) ||
	     (strstr(BL_GFF_FEATURE(&subfeature), "transcript") != NULL) ||
	     (strstr(BL_GFF_FEATURE(&subfeature), "gene_segment") != NULL) ||
	     (strstr(BL_GFF_FEATURE(&subfeature), "_overlapping_ncrna") != NULL) )
	    first_exon = true;
	
	// Generate introns between exons
	if ( exon )
	{
	    if ( !first_exon )
	    {
		intron_end = BL_GFF_START(&subfeature) - 1;
		bl_bed_set_chrom_cpy(&bed_feature, BL_GFF_SEQUENCE(&subfeature),
				 BL_CHROM_MAX_CHARS + 1);
		/*
		 *  BED start is 0-based and inclusive
		 *  GFF is 1-based and inclusive
		 */
		bl_bed_set_chrom_start(&bed_feature, intron_start);
		/*
		 *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
		 *  GFF is the same
		 */
		bl_bed_set_chrom_end(&bed_feature, intron_end);
		snprintf(name, BL_BED_NAME_MAX_CHARS, "intron");
		bl_bed_set_name_cpy(&bed_feature, name, BL_BED_NAME_MAX_CHARS + 1);
		bl_bed_write(&bed_feature, bed_stream, BL_BED_FIELD_ALL);
	    }
	    
	    intron_start = BL_GFF_END(&subfeature);
	    first_exon = false;
	}
	
	bl_gff_to_bed(&subfeature, &bed_feature);
	bl_bed_write(&bed_feature, bed_stream, BL_BED_FIELD_ALL);
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
				   bl_gff_t *gff_feature, bl_pos_list_t *pos_list)

{
    bl_bed_t   bed_feature[MAX_UPSTREAM_BOUNDARIES];
    char            strand,
		    name[BL_BED_NAME_MAX_CHARS + 1];
    int             c;
    
    strand = BL_GFF_STRAND(gff_feature);

    for (c = 0; c < BL_POS_LIST_COUNT(pos_list) - 1; ++c)
    {
	bl_bed_set_fields(&bed_feature[c], 6);
	bl_bed_set_strand(&bed_feature[c], strand);
	bl_bed_set_chrom_cpy(&bed_feature[c], BL_GFF_SEQUENCE(gff_feature),
			     BL_CHROM_MAX_CHARS + 1);
	/*
	 *  BED start is 0-based and inclusive
	 *  GFF is 1-based and inclusive
	 *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
	 *  GFF is the same
	 */
	if ( strand == '+' )
	{
	    bl_bed_set_chrom_start(&bed_feature[c],
			      BL_GFF_START(gff_feature) - 
			      BL_POS_LIST_POSITIONS_AE(pos_list, c + 1) - 1);
	    bl_bed_set_chrom_end(&bed_feature[c],
			    BL_GFF_START(gff_feature) -
			    BL_POS_LIST_POSITIONS_AE(pos_list, c) - 1);
	}
	else
	{
	    bl_bed_set_chrom_start(&bed_feature[c],
			      BL_GFF_END(gff_feature) +
			      BL_POS_LIST_POSITIONS_AE(pos_list, c));
	    bl_bed_set_chrom_end(&bed_feature[c],
			    BL_GFF_END(gff_feature) + 
			    BL_POS_LIST_POSITIONS_AE(pos_list, c + 1));
	}
	
	snprintf(name, BL_BED_NAME_MAX_CHARS, "upstream%" PRIu64,
		 BL_POS_LIST_POSITIONS_AE(pos_list, c + 1));
	bl_bed_set_name_cpy(&bed_feature[c], name, BL_BED_NAME_MAX_CHARS + 1);
    }
    
    if ( strand == '-' )
    {
	for (c = 0; c < BL_POS_LIST_COUNT(pos_list) - 1; ++c)
	    bl_bed_write(&bed_feature[c], feature_stream, BL_BED_FIELD_ALL);
    }
    else
    {
	for (c = BL_POS_LIST_COUNT(pos_list) - 2; c >= 0; --c)
	    bl_bed_write(&bed_feature[c], feature_stream, BL_BED_FIELD_ALL);
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
