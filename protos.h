/* peak-classifier.c */
int main(int argc, char *argv[]);
int filter_gff(FILE *gff_stream, FILE **feature_stream);
int classify(FILE *peak_stream, FILE *gff_stream);
void check_promoter(bed_feature_t *bed_feature, gff_feature_t *gff_feature, uint64_t upstream_dist);
void usage(char *argv[]);
