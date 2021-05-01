/***************************************************************************
 *  Description:
 *      Filter overlaps file produced by peak-classifier
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-30  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>
#include <dsvio.h>
#include "filter-overlaps.h"

int     main(int argc,char *argv[])

{
    char    *overlaps_file,
	    *output_file,
	    **features;

    switch(argc)
    {
	case    1:
	case    2:
	case    3:
	    usage(argv);

	default:
	    overlaps_file = argv[1];
	    output_file = argv[2];
	    features = argv + 3;
	    break;
    }
    return filter_overlaps(overlaps_file, output_file, features);
}


/***************************************************************************
 *  Description:
 *      Process overlaps
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-30  Jason Bacon Begin
 ***************************************************************************/

int     filter_overlaps(const char *overlaps_file, const char *output_file,
			char *features[])

{
    FILE        *infile,
		*outfile;
    dsv_line_t  dsv_line = DSV_INIT,
		keeper = DSV_INIT,
		last_line = DSV_INIT;
    int         delim;
    size_t      keeper_rank,
		new_rank,
		c;
    unsigned long   unique_peaks = 0,
		    feature_overlaps[MAX_OVERLAP_FEATURES];
    
    if ( strcmp(overlaps_file, "-") == 0 )
	infile = stdin;
    else if ( (infile = bio_fopen(overlaps_file, "r")) == NULL )
    {
	fprintf(stderr, "filter-overlaps: Cannot open %s: %s\n",
		overlaps_file, strerror(errno));
	return EX_NOINPUT;
    }
    
    if ( strcmp(output_file, "-") == 0 )
	outfile = stdout;
    else if ( (outfile = bio_fopen(output_file, "w")) == NULL )
    {
	fprintf(stderr, "filter-overlaps: Cannot open %s: %s\n",
		overlaps_file, strerror(errno));
	return EX_CANTCREAT;
    }
    
    for (c = 0; c < MAX_OVERLAP_FEATURES; ++c)
	feature_overlaps[c]= 0;
    
    delim = dsv_read_line(infile, &dsv_line, "\t");
    while ( delim != EOF )
    {
	dsv_free_line(&last_line);
	dsv_copy_line(&last_line, &dsv_line);
	/*
	 *  If this is a keeper (in the features list), check subsequent
	 *  lines with the same peak for higher ranking features.  Input
	 *  is sorted by peak position, so lines with the same peak should
	 *  be contiguous.
	 */
	if ( (keeper_rank = feature_rank(&dsv_line, features)) != 0 )
	{
	    //fprintf(stderr, "%s %zu\n", DSV_FIELD(&dsv_line, 6), keeper_rank);
	    dsv_copy_line(&keeper, &dsv_line);
	    dsv_free_line(&dsv_line);
	    while ( ((delim = dsv_read_line(infile, &dsv_line, "\t")) != EOF)
		    && same_peak(&dsv_line, &keeper) )
	    {
		new_rank = feature_rank(&dsv_line, features);
		// If new feature has a higher rank, replace the old one
		if ( (new_rank != 0) && (new_rank < keeper_rank) )
		{
		    /*fprintf(stderr, "%s:%zu outranks %s:%zu.\n",
			    DSV_FIELD(&dsv_line, 6), new_rank,
			    DSV_FIELD(&keeper, 6), keeper_rank);*/
		    dsv_free_line(&keeper);
		    dsv_copy_line(&keeper, &dsv_line);
		    keeper_rank = new_rank;
		}
	    }
	    ++feature_overlaps[keeper_rank - 1];
	    dsv_write_line(outfile, &keeper);
	    dsv_free_line(&keeper);
	}
	else
	{
	    // Not an interesting feature, toss it
	    dsv_free_line(&dsv_line);
	    delim = dsv_read_line(infile, &dsv_line, "\t");
	}
	if ( (delim != EOF) && !same_peak(&dsv_line, &last_line) )
	    ++unique_peaks;
    }
    fclose(infile);
    fclose(outfile);
    
    printf("Total unique peaks: %zu\n", unique_peaks);
    for (c = 0; features[c] != NULL; ++c)
	printf("Overlaps with %-20s: %7zu (%2zu%%)\n", features[c],
		feature_overlaps[c], 100 * feature_overlaps[c] / unique_peaks);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      See if a DSV line has one of the features for which we're filtering.
 *      Return the 1-based position in the list if so, 0 if not.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

size_t  feature_rank(dsv_line_t *line, char *features[])

{
    int     c;
    
    for (c = 0; features[c] != NULL; ++c)
	if ( strcasecmp(DSV_FIELD(line, 6), features[c]) == 0 )
	    return c+1;
    return 0;
}


/***************************************************************************
 *  Description:
 *      Return true if two features represent the same peak, as evidenced
 *      by the start and end positions in columns 2 and 3.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

bool    same_peak(dsv_line_t *line1, dsv_line_t *line2)

{
    char        *start1,
		*end1,
		*start2,
		*end2;
    
    start1 = DSV_FIELD(line1, 2);
    end1 = DSV_FIELD(line1, 3);
    start2 = DSV_FIELD(line2, 2);
    end2 = DSV_FIELD(line2, 3);
    return (strcmp(start1, start2) == 0) && (strcmp(end1, end2) == 0);
}


void    usage(char *argv[])

{
    fprintf(stderr, "Usage: %s overlap-file.tsv feature [feature ...]\n", argv[0]);
    fprintf(stderr, "Example: %s overlaps.tsv exon intron upstream\n", argv[0]);
    exit(EX_USAGE);
}
