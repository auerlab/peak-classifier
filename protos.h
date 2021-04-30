/* gff2bed.c */
/* peak-classifier.c */
int main(int argc, char *argv[]);
int gff_augment(FILE *gff_stream, const char *upstream_boundaries);
void gff_process_subfeatures(FILE *gff_stream, FILE *bed_stream, gff_feature_t *gene_feature);
int classify(FILE *peak_stream, FILE *gff_stream);
void check_promoter(bed_feature_t *bed_feature, gff_feature_t *gff_feature, uint64_t upstream_dist);
void generate_upstream_features(FILE *feature_stream, gff_feature_t *gff_feature, plist_t *plist);
void usage(char *argv[]);
