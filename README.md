# peak-classifier
Classify ChIP/ATAC-Seq peaks based on features provided in a GFF

Peaks are provided in a BED file sorted by chromosome and position.  The GFF
must also be sorted by chromosome and position, which is the default for common
data sources.

Peak-classifier efficiently identifies overlapping features in the GFF in a
single pass, outputting an annotated BED-like TSV file with additional columns
to describe the feature.  If a peak overlaps multiple features, a separate
line is output for each.

Alternative approaches to this problem include R scripting with a tool such
as ChIPPeakAnno or multistep processing of the GFF using awk and bedtools to
extract features that are not explicitly listed in a GFF, such as introns and
intergenic regions.

In contrast, peak-classifier is a simple Unix command that takes a BED file
and a GFF file as inputs and reports all peak classifications in a matter of
seconds.

Admittedly, this particular problem doesn't really justify writing optimal C
code, since the crappiest implementation I can imagine would not take more
than hours to run for a typical ATAC-Seq peak set.  However:

    * It's an opportunity to develop and test biolibc code that will be
      useful for bigger data
    * This is more about making peak classification convenient than making it
      fast
    * It never hurts to hone your C skills
    * There's no such thing as a program that's too fast

