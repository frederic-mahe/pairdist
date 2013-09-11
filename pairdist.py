#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Pairdist creates a neighbor-joining tree from a set of sequences,
using maximum likelihood distances based on pairwise sequence
alignments.
"""

__author__ = "Frank Kauff <frank.kauff@gmx.de>"

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This version of pairdist has been modified by Frédéric Mahé
# <mahe@rhrk.uni-kl.fr>

import os
import sys
import random
import itertools
import subprocess
from optparse import OptionParser
from Bio import pairwise2, AlignIO
from Bio.Align.Applications import ClustalwCommandline

# Needed for proper execution:
# - PHYLIP <http://evolution.genetics.washington.edu/phylip.html>
# - CLUSTAL <http://www.clustal.org>
# - Biopython <http://www.biopython.org>

# With recent versions of phylip, the commands "dnadist" and
# "neighbor" has to be preceeded by the command "phylip": "phylip
# dnadist" for instance. Older versions can call the sub-programs
# directly.
CLUSTALWCOMMAND = 'clustalw'
PROTDIST = "phylip dnadist"
NEIGHBOR = 'phylip neighbor'

# filenames for intermediate files.
CLUSTALFASTA = 'clustal.fas'
CLUSTALALIGNMENT = 'clustal.aln'
INFILE_PROTDIST = 'infile'
OUTFILE_PROTDIST = 'outfile'
INFILE_NEIGHBOR = 'infile'
OUTFILE_NEIGHBOR = 'outfile'
PROTDIST_COMMANDFILE = 'protcommand'
NEIGHBOR_COMMANDFILE = 'njcommand'

# input commands for protdist and neighbor
NJ_TREE = 'outtree'

# If no clustalw is available, Bio.pairwise (slower) can be used as a
# substitute. Change ALIGNMENT_METHOD accordingly
PAIRWISE2 = 1
CLUSTALW = 2
ALIGNMENT_METHOD = CLUSTALW

PROTDIST_COMMANDS = """2
m
d
%d
y
"""

NEIGHBOR_COMMANDS = """2
3
y
"""

INFILE_FORMAT = '%d %d\n%-10s %s\n%-10s %s\n'

#******************************************************************************#
#                                                                              #
#                                  Functions                                   #
#                                                                              #
#******************************************************************************#

def option_parser():
    """
    Parse arguments from command line.
    """
    desc = """Pairdist creates a neighbor-joining tree from a set of
    sequences, using maximum likelihood distances based on pairwise
    sequence alignments."""
    
    parser = OptionParser(usage="usage: %prog --input_file fasta_file [--strip --bootstrap --nreps integer]",
                          description = desc,
                          version = "%prog version 1.0")

    parser.add_option("-i", "--input_file",
                      metavar = "<FILENAME>",
                      action = "store",
                      dest = "input_file",
                      help = "set <FILENAME> as input fasta file.")

    parser.add_option("-s", "--strip",
                      dest = "strip",
                      action = "store_true",
                      default = False,
                      help = "strip alignment from unmatched ends")

    parser.add_option("-b", "--bootstrap",
                      dest = "bootstrap",
                      action = "store_true",
                      default = False,
                      help = "calculate bootstraps")

    parser.add_option("-n", "--nreps",
                      dest = "nreps",
                      default = 100,
                      type = "int",
                      help = "number of bootstrap replicates")

    (options, args) = parser.parse_args()
    
    return options.input_file, options.strip, options.bootstrap, options.nreps


def sanity_check():
    """Check for the presence of required third-party softwares
    (PHYLIP, CLUSTAL and BIOPYTHON)."""
    protdist = "dnadist"
    neighbor = "neighbor"
    clustalw = "clustalw"
    programs = ("phylip", "dnadist", "protdist", "neighbor", "clustalw", "clustalw2")
    programs_status = dict(zip(programs, [False] * len(programs)))
    for program in programs:
        try:
            returncode = subprocess.call("which" + " " + program, shell=True)
            if returncode is 0:
                programs_status[program] = True
        except OSError as e:
            print >>sys.stderr, "Execution failed:", e
    # Abort if anything is missing
    if programs_status["phylip"] is False and (programs_status["dnadist"] or programs_status["neighbor"] or programs_status["protdist"]) is False:
        print >>sys.stderr, "Phylip package is missing!"
        sys.exit(-1)
    if programs_status["clustalw"] is False and programs_status["clustalw2"] is False:
        print >>sys.stderr, "Clustal package is missing!"
        sys.exit(-1)
    # Target the correct commands
    if programs_status["phylip"]:
        protdist = "phylip dnadist"
        neighbor = "phylip neighbor"
    if programs_status["clustalw2"]:
        clustalw = "clustalw2"
    return clustalw, protdist, neighbor


def replace_safenames(tree, sndict):
    """Replaces the safe OTU descriptions with original OTUs."""
    for sn in sndict.keys():
        tree = tree.replace(sn, sndict[sn])
    return tree


def get_all_possible_pairs(sequences):
    """Returns a list of all possible unique sequence pairs."""
    pairs = [pair for pair in itertools.combinations(sequences, 2)]
    return pairs


def write_fasta(filename, seqs):
    """Writes a file of sequences in FASTA format."""
    f = open(filename, 'w')
    for s in seqs:
        f.write('>%s\n%s\n' % (s[0], s[1]))
    f.close()


def align(method, seq1, seq2):
    """Calls alignment procedures depending on alignment method."""
    if method == PAIRWISE2:
        return align_pairwise(seq1, seq2)
    elif method == CLUSTALW:
        return align_clustalw(seq1, seq2)
    else:
        sys.exit('Unknown aligment method')


def align_clustalw(seq1, seq2):
    """Executes clustalw alignment using the ClustalW module of Biopython."""
    write_fasta(CLUSTALFASTA, (seq1, seq2))
    cline = ClustalwCommandline(CLUSTALWCOMMAND,
                                infile=CLUSTALFASTA, outfile=CLUSTALALIGNMENT)
    stdout, stderr = cline()
    clustalalignment = AlignIO.read(CLUSTALALIGNMENT, "clustal")
    length = clustalalignment.get_alignment_length()
    allseqs = list(clustalalignment)
    bestseq1, bestseq2 = allseqs[0].seq.tostring(), allseqs[1].seq.tostring()
    # clean clustal files
    os.remove(CLUSTALFASTA)
    os.remove(CLUSTALALIGNMENT)
    os.remove("clustal.dnd")
    return bestseq1, bestseq2, length


def align_pairwise(seq1, seq2):
    """Executes alignment using the Bio.pairwise2 module of Biopython."""
    seqname1, seqname2 = seq1[0], seq2[0]
    sequence1, sequence2 = seq1[1], seq2[1]
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    bestseq1, bestseq2 = alignments[0][0], alignments[0][1]
    length = alignments[0][-1]
    return bestseq1, bestseq2, length


def all_pairwise_alignments(sequences, stripgaps=False):
    """
    Gets a list of sequences [(name, sequence), (name,sequence),...].
    Returns a list of pairs [((name,seq),(name,seq)),...] with all
    pairwise alignments. Positions with only gaps in both aligned
    sequences can be stripped off.
    """
    allpairs = get_all_possible_pairs(sequences)
    all_alignments = []
    for i, p in enumerate(allpairs):
        print '%d of %d ' % (i+1, len(allpairs))
        bestseq1, bestseq2, length = align(ALIGNMENT_METHOD, p[0], p[1])
        if stripgaps:
            bs1 = ''
            bs2 = ''
            for i in range(len(bestseq1)):
                if bestseq1[i] != '-' and bestseq2[i] != '-':
                    bs1 += bestseq1[i]
                    bs2 += bestseq2[i]
            all_alignments.append(((p[0][0], bs1), (p[1][0], bs2)))
        else:
            all_alignments.append(((p[0][0], bestseq1), (p[1][0], bestseq2)))
    print '%d pairwise alignments calculated.' % len(all_alignments)
    return all_alignments


def bootstrap_record(sequences):
    """Returns a bootstrap record for the sequences"""
    sitesm = zip(*(zip(*sequences)[1]))
    bootstrapsitesm = [sitesm[random.randint(0, len(sitesm) - 1)]
                       for i in range(len(sitesm))]
    bootstrapseqs = map(''.join, zip(*bootstrapsitesm))
    bootstrap_record = zip(zip(*sequences)[0], bootstrapseqs)
    return bootstrap_record


def pairdisttree(alignments, seqnames, bootstrap=False):
    """Returns a tree (plain text newick) from a bunch of sequences
    [(name,Seq()),...]"""
    # prepare input file for PROTDIST and call PROTDIST
    infile = ''
    npairs = len(alignments)
    for al in alignments:
        length = len(al[0][1])
        if bootstrap:
            brec = bootstrap_record(al)
            infile += INFILE_FORMAT % (2, length, brec[0][0],
                                       brec[0][1], brec[1][0], brec[1][1])
        else:
            infile += INFILE_FORMAT % (2, length, al[0][0],
                                       al[0][1], al[1][0], al[1][1])

    if os.path.exists(OUTFILE_PROTDIST):
        os.remove(OUTFILE_PROTDIST)
    inf = open(INFILE_PROTDIST, 'w')
    inf.write(infile)
    inf.close()

    pc = open(PROTDIST_COMMANDFILE, 'w')
    pc.write(PROTDIST_COMMANDS % npairs)
    pc.close()

    # Execute PROTDIST
    os.system('%s < %s > log' % (PROTDIST, PROTDIST_COMMANDFILE))

    # Read output of PROTDIST and verify length
    distdict = {}
    outfile = open(OUTFILE_PROTDIST).readlines()
    if len(outfile) != npairs * 3:
        sys.exit('Error in PROTDIST outfile\n' + '\n'.join(outfile))

    # format PROTDIST output and generate input file for neighbor
    while outfile:
        l1, l2 = outfile[1].split(), outfile[2].split()
        seq1, seq2 = l1[0], l2[0]
        distance = float(l1[-1])
        # fill both triangles of matrix so we know how to look them up
        distdict[(seq1, seq2)] = distance
        distdict[(seq2, seq1)] = distance
        outfile = outfile[3:]
        if os.path.exists(INFILE_NEIGHBOR):
            os.remove(INFILE_NEIGHBOR)
        if os.path.exists(OUTFILE_NEIGHBOR):
            os.remove(OUTFILE_NEIGHBOR)
        if os.path.exists(NJ_TREE):
            os.remove(NJ_TREE)

    # now we generate a matrix for input to neighbor
    dm = open(INFILE_NEIGHBOR, 'w')
    dm.write('    %d\n' % len(seqnames))
    for i, s1 in enumerate(seqnames):
        dm.write('%-10s' % s1)
        for j, s2 in enumerate(seqnames):
            if i != j:
                dm.write('  %1.6f' % distdict[(s1, s2)])
            else:
                dm.write('  0.000000')
        dm.write('\n')
    dm.close()

    nj = open(NEIGHBOR_COMMANDFILE, 'w')
    nj.write(NEIGHBOR_COMMANDS)
    nj.close()
    os.system('%s <%s >neighbor.log' % (NEIGHBOR, NEIGHBOR_COMMANDFILE))

    ot = open(NJ_TREE).read().replace('\n', '')
    return ot


if __name__ == '__main__':

    input_file, strip, bootstrap, nreps = option_parser()
    CLUSTALWCOMMAND, PROTDIST, NEIGHBOR = sanity_check()

    # create tmp files and pointers

    # parse the fasta file (normal IO, degap reads if needed)
    seqs = [(s.id, s.seq.tostring())
            for s in list(AlignIO.read(input_file, 'fasta'))]
    
    # use safenames
    safenames = {}
    seqsnew = []
    # PHYLIP software needs strict phylip format, i.e. sequence names
    # of exactly 10 characters. We substitute original sequence names
    # during script execution
    for i, (seqname, seq) in enumerate(seqs):
        safename = str(i).zfill(3) + seqname[:4] + str(i).zfill(3)
        safename = safename[:10]  # 10 characters max
        safename = ''.join([c for c in safename if c.isalnum()])
        safename = safename.zfill(10)
        safenames[safename] = seqname
        seqsnew.append((safename, seq.replace('-', '')))
        # print seqname,' -> ',safename

    # do the pairwise alignments
    all_alignments = all_pairwise_alignments(seqsnew, stripgaps=strip)
    safenamesorder = zip(*seqsnew)[0]

    trees = []
    if bootstrap:
        log = open('pairdist_bootstrap.log', 'w')
        for i in range(nreps):
            print 'Bootstrap %d / %d' % (i+1, nreps)
            bs_tree = pairdisttree(all_alignments,
                                   safenamesorder, bootstrap=True)
            trees.append(bs_tree)
            log.write(bs_tree+'\n')
            log.flush()
        log.close()
    else:
        trees = [pairdisttree(all_alignments, safenamesorder)]

    otn = open(input_file + '.tree', 'w')
    otn.write('begin trees;\n')
    for i, t in enumerate(trees):
        otn.write('tree pdtree%d = %s;\n' %
                  (i, replace_safenames(t, safenames)))
    otn.write('end;\n')
    otn.close()

sys.exit(0)

# python pairdist.py test.fas_aln
# rm log clustal.* infile outfile protcommand neighbor.log njcommand *tree
