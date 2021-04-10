/* peak-classifier.c */
int main(int argc, char *argv[]);
int classify(FILE *bed_stream, FILE *gff_stream);
void check_promoter(bed_feature_t *bed_feature, gff_feature_t *gff_feature, uint64_t upstream_dist);
void usage(char *argv[]);
