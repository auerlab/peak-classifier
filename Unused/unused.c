

/***************************************************************************
 *  Description:
 *      Read through BED and GFF files, classifying each peak in the BED
 *      file according to overlapping features in the GFF.
 *      This function is incomplete as this task is being outsourced to
 *      bedtools intersect for now. Should needs arise that bedtools cannot
 *      meet, development will resume.
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
    bed_feature_t   bed_feature = BED_INIT;
    gff_feature_t   gff_feature = GFF_INIT;
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
	if ( strcmp(GFF_NAME(&gff_feature), "gene") == 0 )
	{
	}
	
	/*
	 *  Check for explicitly listed features
	 */
	if ( dist == 0 )
	{
	    feature = GFF_NAME(&gff_feature);
	    if ( (strcmp(feature, "five_prime_UTR") == 0) ||
		 (strcmp(feature, "lnc_RNA") == 0) ||
		 (strcmp(feature, "miRNA") == 0) ||
		 (strcmp(feature, "ncRNA") == 0) ||
		 (strcmp(feature, "exon") == 0) )
	    {
		fputs("===\n", stderr);
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
 *      Not currently used: Task handled by bedtools intersect
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-10  Jason Bacon Begin
 ***************************************************************************/

void    check_promoter(bed_feature_t *bed_feature, gff_feature_t *gff_feature,
			uint64_t upstream_dist)

{
    uint64_t        promoter_start;
    gff_feature_t   gff_promoter = GFF_INIT;
    char            promoter_name[GFF_NAME_MAX_CHARS + 1];
    bio_overlap_t   overlap;
    
    /* Build a fake GFF feature to represent the upstream region */
    promoter_start = GFF_START_POS(gff_feature) - upstream_dist;
    GFF_SET_START_POS(&gff_promoter, promoter_start);
    GFF_SET_END_POS(&gff_promoter, GFF_START_POS(gff_feature) - 1);
    fputs("======\n", stderr);
    printf("Peak at: %" PRIu64 " to %" PRIu64 "\n",
	    BED_START_POS(bed_feature), BED_END_POS(bed_feature));
    printf("Gene at: %" PRIu64 " to %" PRIu64 "\n",
	    GFF_START_POS(gff_feature), GFF_END_POS(gff_feature));
    printf("Trying promoter %" PRIu64 " to %" PRIu64 "\n",
	    GFF_START_POS(&gff_promoter), GFF_END_POS(&gff_promoter));
    if ( bed_gff_cmp(bed_feature, &gff_promoter, &overlap) == 0 )
    {
	GFF_SET_NAME(&gff_promoter, promoter_name);
	snprintf(promoter_name, GFF_NAME_MAX_CHARS,
		"promoter%" PRIu64, upstream_dist);
	fputs("======\n", stderr);
	bed_write_feature(stdout, bed_feature, BED_FIELD_ALL);
	gff_write_feature(stdout, &gff_promoter, BED_FIELD_ALL);
    }
}

