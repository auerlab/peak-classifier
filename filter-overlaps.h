
#define MAX_OVERLAP_FEATURES    64

void    usage(char *argv[]);
int     filter_overlaps(const char *overlaps_file, const char *output_file,
	char *features[]);
size_t  feature_rank(xt_dsv_line_t *line, char *features[]);
bool    same_peak(xt_dsv_line_t *line1, xt_dsv_line_t *line2);

