#!/usr/bin/env python
###############################################################################
#
# __edit_phyla_file__.py - description!
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

from multiprocessing import Pool
from subprocess import Popen, PIPE

import math
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

def doWork( args ):
    """ Main wrapper"""
    # objects
    ids_dict = {}
    listing = glob.glob('%s/*.fna' % args.genome_dir)
    genome_length_dir = {}
    bodysite_dir = {}
    
    # read in ordered ids file
    """
    genome_tree_id  img_id  order   tax_string
    """
    with open(args.tax_file,"r") as fh:
        header = fh.readline() #capture header
        for l in fh: #line by line
            tabs= l.split("\t")
            img= tabs[1]
            gt= tabs[0]
            ids_dict[gt]=img
    #-----
    """
    img -> bodysite
    """
    #read in bodysite file
    with open(args.bodysite,"r") as fh:
        # no header
        for l in fh: #line by line
            tabs = l.split("\t")
            gt_id= tabs[0]
            img_id= tabs[1]
            bodysite= tabs[2].rstrip()
            bodysite_dir[gt]= bodysite
    #-----
    """ 
    Calculate the total genome length
    """
    #read in genomes directory
    for c_file in listing:
        #img_id = c_file.split("/")[2].split(".")[0] 
        #if count <1:
        img_id = c_file.split("/")[-1].split(".")[0]
        #print img_id
        for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
            try: 
                genome_length_dir[img_id]+= len(sequence)
            except KeyError:
                genome_length_dir[img_id]= len(sequence)
    
    #-----        
    # header
    print "\t".join(["genome_tree_id",
                     "img_id",
                     "phylum",
                     "genome_length",
                     "bodysite"
                     ])
    # read in phylum file
    with open(args.phyla_file,"r") as fh:
        for l in fh:
            tabs= l.split("\t")
            gt= tabs[0]
            phyla= tabs[1].rstrip()
            print "\t".join([gt,
                             ids_dict[gt],
                             phyla,
                             str(genome_length_dir[ids_dict[gt]]),
                             bodysite_dir[gt]
                             ])
     
            
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
    parser.add_argument('-phyla_file','--phyla_file', help="File containing phylum information")
    parser.add_argument('-tax_file','--tax_file', help="File containing img ids and genome tree ids")
    parser.add_argument('-g','--genome_dir', help="Directory containing genomes")
    parser.add_argument('-bodysite','--bodysite', help="File containing bodysite information")
    #parser.add_argument('input_file2', help="gut_img_ids")
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
