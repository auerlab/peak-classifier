#!/usr/bin/env python3

#############################################################################
#   Main Program
#
#   Description:
#       Generate a bar graph representing the features in a small section
#       of a BED file, such as a gene
#   
#   History: 
#   Date        Name        Modification
#   2021-04-21  Jason Bacon Begin

import sys,os
# Pkgsrc
# Requires math/py-matplotlib and x11/py-Tk
# import matplotlib
# matplotlib.use('TkAgg')   # Seems no longer necessary
# End pkgsrc
import matplotlib.pyplot as pyplot
import numpy

#############################################################################
#   Process command line args

if len(sys.argv) != 2:
    sys.exit("Usage: sys.argv[0] bed-file\n")
bed_file=sys.argv[1]

#############################################################################
#   Read BED file into lists of positions and feature names for bar graph

f=open(bed_file,"r")
features=[]
lengths=[]
fstarts=[]
for line in f.readlines():
    a=line.split()
    features.append(a[3])
    lengths.append(int(a[2])-int(a[1]))
    fstarts.append(int(a[1]))
f.close()

#############################################################################
#   Set up plot arguments

pyplot.rcdefaults()
fig, ax = pyplot.subplots()
y_pos = numpy.arange(len(features))
ax.barh(y_pos, lengths, 0.3, left=fstarts, align='center')
ax.set_yticks(y_pos)
ax.set_yticklabels(features)
ax.invert_yaxis()  # labels read top-to-bottom
ax.set_xlabel('Position')
ax.set_title('Gene subfeatures')

pyplot.show()
