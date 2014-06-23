#!/usr/bin/env python
###############################################################################
#
# __ANI_stats__.py - Description!
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
import math

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

class ANIparser(object):
    def __init__(self,l):
        self.readANI(l)
        
    def readANI(self,l):
        tabs = l.rstrip().split("\t")
        self._img_id_1  = tabs[0]
        self._img_id_2  = tabs[2]
        self._species_1 = tabs[1]
        self._species_2 = tabs[3]
        self._ANI_1     = tabs[4]
        self._ANI_2     = tabs[5]
        self._AFI_1     = tabs[6]
        self._AFI_2     = tabs[7]

class ANIDB(object):
    def __init__(self):
        self.pairs_db = {}
        self.ANI_scores = {}

    def addPAIR(self,l,id):
        line = ANIparser(l)
        self.pairs_db[id] = [line._img_id_1, line._img_id_2]
        
    def checkPAIR(self,l):
        line = ANIparser(l)
        try:
            for key in self.pairs_db.keys():
                if self.pairs_db[key][0] == line._img_id_1 and self.pairs_db[key][1] == line._img_id_2:
                    return False
                else: 
                    return True
        except KeyError:
            pass

    def addScores(self,l,uid):
        line = ANIparser(l)
        self.ANI_scores[uid] = [line._ANI_1_, line._ANI_2]
    
    def getAverageDIFF(self):
        # difference between ANI scores
        counter = 0 
        total = 0
        for key in  self.ANI_scores.keys():
            total = total + math.fabs(self.ANI_scores[key][0]-self.ANI_scores[key][1])
            counter+=1
        print str(total/counter)     
    

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
    uid = 1
    ANI = ANIDB()
    count_break = 0
    
    # open ANI file
    with open(args.ANI_file,"r") as fh:
        header = fh.readline()  # capture header
        for l in fh:
            if ANI.checkPAIR(l):
                ANI.addPAIR(l, uid)
                ANI.addScores(l, uid)
                uid +=1
            if count_break >= 100:
                break
            count_break +=1  
    ANI.getAverageDIFF() 
                        
            
            
            
            
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
    parser.add_argument('-ani','--ANI_file', help="...")
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
