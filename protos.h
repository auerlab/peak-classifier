/* peak-classifier.c */
int main(int argc, char *argv[]);
int gff3_augment(FILE *gff3_stream, const char *upstream_boundaries, const char *augmented_filename);
void gff3_process_subfeatures(FILE *gff3_stream, FILE *bed_stream, bl_gff3_t *gene_feature);
void generate_upstream_features(FILE *feature_stream, bl_gff3_t *gff3_feature, bl_pos_list_t *pos_list);
void usage(char *argv[]);
