#!/usr/bin/env python
###############################################################################
#
# __parse_methylase_genes__.py - Description!
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

class annotationParser(object):
    def __init__(self,l):
        self.readAnnotationFile(l)
    
    def readAnnotationFile(self,l):
        tabs = l.rstrip().split("\t")
        self.uid            = tabs[0].split("-")[-1] 
        self.cog_id         = tabs[5]
        self.cog_annotation = tabs[6]
        self.annotation     = tabs[8]
        
class annotationDB(object):
    def __init__(self):
        self.anno_db = {}
    
    def addMethylaseGene(self,l):
        line = annotationParser(l)
        #self.anno_db[line.uid] = [line.annotation,line.cog_annotation] # add uid to dictionary
        print line.uid
    
    def returnUIDs(self):
        for uid in self.anno_db.keys():
            print "\t".join([uid,self.anno_db[uid][0],self.anno_db[uid][1]])

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
    # 1. read in prodigalled.annoated file
    # 2. Capture methylase gene unique names
    # 3. read in prodigalled file and grab fastaa seqs
    # 4. Run promer/blast of gut_oral_methylase genes vs methylase gene DB
    
    # objects
    ANNO = annotationDB()
    
    # read in annotation file
    with open(args.anno_file,"r") as fh:
        for l in fh:
            #print l.rstrip()
            ANNO.addMethylaseGene(l)
    #ANNO.returnUIDs() # print uids containing methylase or restriction
        

    
    
    # parse fasta file using biopython
    #for accession,sequence in SeqIO.to_dict(SeqIO.parse(args.fasta_file,"fasta")).items():
    #    print ">"+accession
    #    print sequence.seq
            
    

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
    parser.add_argument('-a','--anno_file', help="...")
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
