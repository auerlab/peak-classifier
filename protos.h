/* peak-classifier.c */
int main(int argc, char *argv[]);
int filter_gff(FILE *gff_stream, const char *upstream_boundaries);
void gff_process_subfeatures(FILE *gff_stream, FILE *bed_stream, gff_feature_t *gene_feature);
void gff_to_bed(bed_feature_t *bed_feature, gff_feature_t *gff_feature);
int classify(FILE *peak_stream, FILE *gff_stream);
void check_promoter(bed_feature_t *bed_feature, gff_feature_t *gff_feature, uint64_t upstream_dist);
void generate_upstream_features(FILE *feature_stream, gff_feature_t *gff_feature, plist_t *plist);
void gff_plot_subfeature(FILE *stream, gff_feature_t *gene, gff_feature_t *subfeature);
void usage(char *argv[]);
