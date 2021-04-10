/* peak-classifier.c */
int main(int argc, char *argv[]);
int classify(FILE *bed_stream, FILE *gff_stream);
void bed_check_order(bed_feature_t *bed_feature, char last_chrom[], uint64_t last_pos);
void usage(char *argv[]);
FILE *bio_open(char *filename, char *mode, _Bool *is_pipe);
int bed_gff_cmp(bed_feature_t *bed_feature, gff_feature_t *gff_feature, bio_overlap_t *overlap);
