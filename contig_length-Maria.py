#!/usr/bin/env python
###############################################################################
#
# __contig_length__.py - Grab sequence information from a directory of fna files.
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2014"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = ""
__status__ = "Development"

###############################################################################

import argparse
import sys
import datetime

from multiprocessing import Pool
from subprocess import Popen, PIPE
import glob 
from Bio import SeqIO
from Bio.Seq import Seq

#import os
#import errno

#import numpy as np
#np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

  # classes here

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    p = Popen(cmd.split(' '), stdout=PIPE)
    return p.communicate()

def printHeader():
    print "\t".join(["ORF","ORF_length","ORF_GC_content","contig","contig_length"])

def printOUT(orf_len,orf_gc,contig_len):
    for orf in orf_len.keys():
        contig = orf.split("_")[2] + orf.split("_")[3]
        print "\t".join([orf, orf_len[orf], orf_gc[orf], contig, contig_len[contig]])

def doWork( args ):
    """ Main wrapper"""
    """Read in a directory containing multiple fna files, then file by file count the length of contigs.
    Then write this information to a file""" 
    #global variables    
    #listing = glob.glob('%s/*.fna' % args.genomes_directory)
    ORF_dict = {} #dictionary to store index information of genomes
    GC_dict = {}
    contigs_dict = {} #store unique contig names
    count = 0
    
    # Read through ORF file
    with open(args.orf_file,"r") as fh:
        for l in fh:
            if ">" in l:
                # capture accession information
                accession   = l.split()[0][1:]
                start       = int(l.split()[2])
                stop        = int(l.split()[4])
                GC_cont     = l.split()[-1].split(";")[-1].split("=")[-1]
                contig      = l.split()[0].split("_")[2] + l.split()[0].split("_")[3]
                length = 0 
                if stop > start:
                    length = (stop - start) + 1
                else: 
                    length = (start - stop) + 1
                contigs_dict[contig] = None # add contig to dict
                ORF_dict[accession] = length
                GC_dict[accession] = GC_cont
    
    # read through contig fasta file
    for accession,sequence in SeqIO.to_dict(SeqIO.parse(args.contig_file,"fasta")).items():
        contig = accession.split("_")[2] + accession.split("_")[3]
        if contig in contigs_dict:
            contigs_dict[contig] = len(sequence)
    
    printHeader()
    printOUT(ORF_dict, GC_dict, contigs_dict)
 
    """
# run somethign external in threads
pool = Pool(6)
cmds = ['ls -l', 'ls -alh', 'ps -ef']
print pool.map(runCommand, cmds)
"""

    """
# parse a file
try:
with open(filename, "r") as fh:
for line in fh:
print line
except:
print "Error opening file:", filename, exc_info()[0]
raise
"""

    """
fig = plt.figure()

#-----
# make a 3d plot
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0],
points[:,1],
points[:,2],
#edgecolors='none',
#c=colors,
#s=2,
#marker='.'
)

#-----
# make a 2d plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(points[:,0],
points[:,1],
'*g')

#-----
# show figure
plt.show()
# or save figure
plt.savefig(filename,dpi=300,format='png')

#-----
# clean up!
plt.close(fig)
del fig
"""

    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-orf','--orf_file', help="File containing open reading frames in fasta format")
    parser.add_argument('-contig','--contig_file', help="File containing contigs in fasta format")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
