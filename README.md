DETECT - a Density Estimation Tool for Enzyme CassificaTion
======

Original publication by <a href="http://bioinformatics.oxfordjournals.org/content/26/14/1690">Hung, Wasmuth, Sanford and Parkinson (2010)</a>

DETECT implements a density estimation determination for enzyme function based on the sequence. Given a FASTA file, it produces a number of EC number predictions.

This is the python implementation that includes several optimizations compared to the original version.

Improvements mainly include lowering the nubmer of file system operations and using an SQLite database for density estimation tables.

##Examples##
Sample run
<pre><code>python detect.py --verbose target_sequence.fasta</code></pre>
Running DETECT for a large sequence archive
<pre><code>python detect.py target_seqeunces.fasta --tabular_output --output_file predictions.tsv</code></pre>
Built-in help available
<pre><code>python detect.py --help</code></pre>

##Dependencies###
+ <a href="http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download">NCBI BLAST</a>
+ <a href="http://emboss.sourceforge.net/download/">EMBOSS suite</a>
+ <a href="http://www.sqlite.org/download.html">SQLite3</a>
+ <a href="http://www.python.org/">python2.7.3</a>
