pairdist
========

Creates a neighbor-joining tree from a set of sequences, using maximum
likelihood distances based on pairwise sequence alignment.

The classical Neighbor-joining approach operates on a matrix of
pairwise distances calculated from a multiple sequence alignment. In
such an alignment, sequences are aligned, usually by a software
package, in such a way so that the overall mismatches among all
sequences are minimized according to an optimization criterion. From a
computational point of view, multiple sequence alignments are a
n-complete problem for which an optimal solution cannot be achieved in
a reasonable time frame. Many software packages with a variety of
alignment strategies, optimization criteria, and numerous
user-definable parameters exist.