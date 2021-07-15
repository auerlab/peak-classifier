/* peak-classifier.c */
int main(int argc, char *argv[]);
int gff_augment(FILE *gff_stream, const char *upstream_boundaries, const char *augmented_filename);
void gff_process_subfeatures(FILE *gff_stream, FILE *bed_stream, bl_gff_t *gene_feature);
void generate_upstream_features(FILE *feature_stream, bl_gff_t *gff_feature, bl_pos_list_t *pos_list);
void usage(char *argv[]);
