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
#include "overlap-filter.h"

int     main(int argc,char *argv[])

{
    char    *overlaps_file,
	    **features;

    switch(argc)
    {
	case    1:
	case    2:
	    usage(argv);

	default:
	    overlaps_file = argv[1];
	    features = argv + 2;
	    break;
    }
    return filter_overlaps(overlaps_file, features);
}


/***************************************************************************
 *  Description:
 *      Process overlaps
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-30  Jason Bacon Begin
 ***************************************************************************/

int     filter_overlaps(char *overlaps_file, char *features[])

{
    FILE        *infile;
    dsv_line_t  dsv_line = DSV_INIT,
		keeper = DSV_INIT;
    int         delim;
    size_t      keeper_rank,
		new_rank;
    unsigned long   total_overlaps = 0;
    
    if ( strcmp(overlaps_file, "-") == 0 )
	infile = stdin;
    else if ( (infile = bio_fopen(overlaps_file, "r")) == NULL )
    {
	fprintf(stderr, "filter-overlaps: Cannot open %s: %s\n",
		overlaps_file, strerror(errno));
	return EX_NOINPUT;
    }
    
    delim = dsv_read_line(infile, &dsv_line, "\t");
    while ( delim != EOF )
    {
	/*
	 *  If this is a keeper (in the features list), check subsequent
	 *  lines with the same peak for higher ranking features.
	 */
	if ( (keeper_rank = feature_rank(&dsv_line, features)) != 0 )
	{
	    dsv_copy_line(&keeper, &dsv_line);
	    dsv_free_line(&dsv_line);
	    while ( ((delim = dsv_read_line(infile, &dsv_line, "\t")) != EOF)
		    && same_peak(&dsv_line, &keeper) )
	    {
		new_rank = feature_rank(&dsv_line, features);
		// If new feature has a higher rank, replace the old one
		if ( (new_rank != 0) && (new_rank < keeper_rank) )
		{
		    dsv_free_line(&keeper);
		    dsv_copy_line(&keeper, &dsv_line);
		}
	    }
	    dsv_write_line(stdout, &keeper);
	    dsv_free_line(&keeper);
	}
	else
	{
	    // Not an interesting feature, toss it
	    dsv_free_line(&dsv_line);
	    delim = dsv_read_line(infile, &dsv_line, "\t");
	}
	++total_overlaps;   // Don't count multiple overlaps of same peak
    }
    fclose(infile);
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
	if ( strstr(DSV_FIELD(line, 6), features[c]) != NULL )
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
