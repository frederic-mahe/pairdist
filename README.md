pairdist
========

PairDist creates a neighbor-joining tree from a set of sequences,
using maximum likelihood distances based on pairwise sequence
alignment.

## PairDist's logic ##

The classical Neighbor-joining approach operates on a matrix of
pairwise distances calculated from a multiple sequence alignment. In
such an alignment, sequences are aligned, usually by a software
package, in such a way so that the overall mismatches among all
sequences are minimized according to an optimization criterion. From a
computational point of view, multiple sequence alignments are a
*n*-complete problem for which an optimal solution cannot be achieved in
a reasonable time frame. Many software packages with a variety of
alignment strategies, optimization criteria, and numerous
user-definable parameters exist.

A key problem in alignments is the question of positional
homology—which nucleotide or amino acid positions are homologous to
each other, and thus suitable for comparison, e.g. when calculating a
maximum likelihood distance of two sequences. For regions in the
alignment with a very high level of substitutions and
insertions-deletions, such as seen in hyper-variable gene regions,
statements of positional homology can be highly speculative, and many
equal or nearly equal solutions may exist. As a consequence, such
regions are often excluded from analyses, because wrong assumptions
about positional homology can have deteriorating effects on both
precision and accuracy of the results. On the other hand, excluding
data loses information, as inability to generate a multiple sequence
alignment for certain regions does not necessarily mean that those
regions do not contain valuable information.

PairDist is an attempt to overcome this problem for nucleotide
sequences that were previously problematic, by reducing the data to
smaller taxonomic sets so they are more easily alignable. Pairdist.py
is a python script that connects the commands *clustalw2* from the
[Clustal package](http://www.clustal.org/clustal2 "Clustal website")
with *dnadist* and *neighbor* from the
[Phylip package](http://evolution.genetics.washington.edu/phylip.html
"Phylip homepage"). Rather than calculating sequence distance from a
full multiple sequence alignment, each sequence pair is aligned
independently (with clustal) prior to the calculation of the Maximum
Likelihood distance with *dnadist*. From the resulting pairwise
distances, a matrix is generated which serves as input for *neighbor*
which finally calculates a neighbor-joining tree. A bootstrap option
is available, where the two-sequence alignment is bootstrapped prior
to distance calculation. PairDist can handle protein sequences, but
this option is considered as experimental.

## Requirements ##

Besides the Phylip and Clustal packages, the python module
[Biopython](http://www.biopython.org) is required to run
pairdist.py. Pairdist has been tested with Biopython version 1.59;
newer versions are likely to work as well. An alternative for
*clustalw2* is available in the Biopython package and integrated in
pairdist.py, but execution time is greatly decreased when *clustalw2*
is not available. After a standard installation of Biopython, Phylip
and Clustal packages, pairdist.py should run without changes, assuming
that *clustalw2*, *dnadist*, and *neighbor* are in your path and
available system-wide. For installation details of the prerequisite
software packages please consult their respective manuals.

## Installation and use ##

The program pairdist.py is written in the
[python programming language](http://www.python.org) and available for
download at: https://github.com/frederic-mahe/pairdist. Unzip the
package and copy the executable pairdist.py either into the folder
where your data lives, or in any location of your system path,
e.g. `/bin`, `/usr/bin`, or `/usr/local/bin` or in most GNU/Linux or
Mac systems. The details may vary according to the specific set-up of
your computer.

Pairdist.py is a simple command line tool. A short help describing the
usage and options can be printed as such:

```
python pairdist.py -h
```

Given an input file in FASTA format (aligned sequences or not), the
program is called as:

```
python pairdist.py -i sequences.fas
```

where `sequences.fas` is to be replaced by the file name of your input
FASTA file.

Among various intermediate files produced by the script, the resulting
tree in Newick-format is written to a file which has the same name as
your input file with the suffix '.tree' added. It can be read as input
and displayed by most applications for visualization of phylogenetic
trees, e.g. [FigTree](http://tree.bio.ed.ac.uk/software/figtree/
"FigTree: viewer of phylogenetic trees").

With the options `-b` and `-n`, a bootstrap run is performed with a
number of replicates specified with `-n`, e.g.

```
python pairdist.py -i sequences.fas -b -n 100
```

will, in addition to the NJ-tree, also calculate 100 bootstrap
replicates, written to the file `pairdist_bootstrap.trees`. Bootstrap
trees and NJ-tree are then merged into a single output file, named as
above. This file contains the NJ-tree with branch lengths together
with the bootstrap frequencies, and can be displayed using
[FigTree](http://tree.bio.ed.ac.uk/software/figtree/ "FigTree: viewer
of phylogenetic trees") or other programs.
