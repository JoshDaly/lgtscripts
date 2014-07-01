#!/usr/bin/env python
###############################################################################
#
# __DB_paired_table__.py - DATABASE FILE: create paired genome table for database!
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

#from Bio import SeqIO
#from Bio.Seq import Seq

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
        
class METAparser(object):
    def __init__(self,l):
        self.readMetadata(l)
        
    def readMetadata(self,l):
        """ taxon  = 0, genome_name = 4"""
        tabs = l.rstrip().split("\t")
        self._img_id        = tabs[0]
        self._genome_name   = tabs[4]
        self._body_site     = tabs[27]
        
class ANIparser(object):
    def __init__(self,l):
        self.readANI(l)
        
    def readANI(self,l):
        """IMG_Taxon1      Species IMG_taxon2      Species2        ANI1    ANI2    AF1     AF2"""
        tabs = l.rstrip().split("\t")
        self._img_id_1  = tabs[0]
        self._img_id_2  = tabs[2]
        self._ANI_1     = tabs[4]
        self._ANI_2     = tabs[5]
        
class genomeTreeParser(object):
    def __init__(self,l):
        self.readGT(l)
    
    def readGT(self,l):
        """genome_tree_id  img_id  order   tax_string"""
        tabs = l.rstrip().split("\t")
        self._gt_id  = tabs[0]  
        self._img_id = tabs[1]

class paired_data(object):
    def __init__(self):
        self.img_to_gt_dict = {} # dict to store img-> genome tree ids
        self.img_metadata_dict = {}
        self.ANI_with_gt = {}
        
    def addGT(self,l):
        line = genomeTreeParser(l)
        self.img_to_gt_dict[line._img_id] = line._gt_id
        
    def printGTs(self):
        for key in self.img_to_gt_dict.keys():
            print "\t".join([key,self.img_to_gt_dict[key]])
        
    def addMETA(self,l):
        line = METAparser(l)
        #self.img_metadata_dict[line._img_id] = 
    
    def checkID(self,l,type): # check if it has genome tree id
        if type == "ANI":
            line =  ANIparser(l)
            if line._img_id_1 in self.img_to_gt_dict and line._img_id_2 in self.img_to_gt_dict:
                return True
        if type == "META"
            line = METAparser(l)
            if line._img_id in self.img_to_gt_dict:
                return True


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
    PD = paired_data() # call class
    
    # read in genome tree file
    with open(args.genome_tree_file,"r") as fh:
        header = fh.readline()
        for l in fh:
            PD.addGT(l)
    printGTs()
    """
    # read in ANI file
    with open(args.ANI_file,"r") as fh:
        header = fh.readline()
        for l in fh:
            if PD.checkID(l, "ANI"): # ID in genome tree list
                
        
    # read in IMG metadata
    with open(args.img_metadata,"r") as fh:
        header = fh.readline()
        for l in fh:
            if PD.checkID(l, "META"): # ID in genome tree list
    """            
                    
            
            
            
            
    """
# parse fasta file using biopython
for accession,sequence in SeqIO.to_dict(SeqIO.parse(c_file,"fasta")).items():
if accession in genomes_dict:
pass
else:
#print accession
genomes_dict[accession] = [len(sequence),img_id, sequence.seq
"""  
    

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
    parser.add_argument('-gt','--genome_tree_file', help="...")
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
