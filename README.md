DETECT
======

#DETECT - a Density Estimation Tool for Enzyme CassificaTion#

Original publication by <a href="http://bioinformatics.oxfordjournals.org/content/26/14/1690">Hung et al., 2010</a>

DETECT implements a density estimation determination for enzyme function based on the sequence. Given a FASTA file, it produces a number of EC number predictions.

This is the python implementation that includes several optimizations compared to the original version.

Improvements mainly include lowering the nubmer of file system operations and using an SQLite database for density estimation tables.
