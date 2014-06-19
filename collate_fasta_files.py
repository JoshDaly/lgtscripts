#!/usr/bin/env python
###############################################################################
#
# __collate_fasta_file__.py - Collate multiple fasta files into a single file. 
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
import glob

from multiprocessing import Pool
from subprocess import Popen, PIPE

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

class metadataParser(object):
    def __init__(self,l):
        self.readMETA(l)
        
    def readMETA(self,l):
        tabs        = l.rstrip().split("\t")
        self.img_id      = tabs[0]
        self.genome_name = tabs[4]
    
class metadataDB(object):
    def __init__(self):
        self.plasmid_db = {}
        self.virus_db   = {}
    
    def addPlasmid(self, l):
        line = metadataParser(l)
        self.plasmid_db[line.img_id] =  line.genome_name
        
    def addVirus(self, l):
        line = metadataParser(l) 
        self.virus_db[line.img_id] =  line.genome_name
        
    def checkIDvirus(self, id):
        if id in self.virus_db:
            return True
    
    def checkIDplasmid(self, id):
        if id in self.plasmid_db:
            return True
    
    def printDict(self):
        for id in self.virus_db.keys():
            print id
    

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
    listing = glob.glob('%s/*/*[0-9].fna' % args.genomes_directory)
    META    = metadataDB() # call class
    
    # read in viral/plasmid metadata files
    with open(args.plasmid_metadata,"r") as fh:
        header = fh.readline() # capture header
        for l in fh:
            META.addPlasmid(l)
        
    with open(args.virus_metadata,"r") as fh:
        header = fh.readline() # capture header
        for l in fh:
            META.addVirus(l)

    for c_file in listing:
        img_id = c_file.split("/")[-1].split(".")[0]
        #print img_id
        # parse fasta file using biopython
        
        if META.checkIDplasmid(img_id):
            print img_id
            for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
                pass
        if META.checkIDvirus(img_id):
            print img_id
            for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
                pass
        #    print accession
            #if META.checkIDplasmid(img_id):
            #    print img_id
            #if META.checkIDvirus(img_id):
            #    print img_id
        
            

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
    parser.add_argument('-g','--genomes_directory', help="...")
    parser.add_argument('-p','--plasmid_metadata', help="...")
    parser.add_argument('-v','--virus_metadata', help="...")
    parser.add_argument('-o','--output_file', help="...")
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
