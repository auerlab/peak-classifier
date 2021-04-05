# peak-classifier
Classify ChIP/ATAC-Seq peaks based on features provided in a GFF

Peaks are provided in a BED file sorted by chromosome and position.  The GFF
must also be sorted by chromosome and position, which is the default for most
data sources.

Peak-classifier efficiently identifies overlapping features in the GFF in a
single pass, outputting an annotated BED file with additional columns to
describe the feature.  If a peak overlaps multiple features, a separate line
is output for each.
