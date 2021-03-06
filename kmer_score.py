#!/usr/bin/env python
###############################################################################
#
# __script_name__.py - description!
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
from scipy.spatial.distance import pdist
import glob 

from multiprocessing import Pool
from subprocess import Popen, PIPE

#import os
#import errno

import numpy as np
np.seterr(all='raise')

#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from pylab import plot,subplot,axis,stem,show,figure


###############################################################################
###############################################################################
###############################################################################
###############################################################################

class LGTInfoStore(object):
    """Helper class for storing information about directionality of LGT events"""
    # __init__ links to object, and must be present in every Class! 
    def __init__(self, lgtTmerArray):
        self.lgtTmer = lgtTmerArray
        self.genomeTmers = {}
        self.lgtGenomes = {}
        self.what = None

    def addWhat(self, what):
        self.what = what
    
    def addLGT(self,lgt_id, GID1, GID2):
        self.lgtGenomes[lgt_id] = [GID1,GID2]
        
    def addGenome(self, GID):
        """add a holder for a genome but don't worry about the the tmer yet"""
        self.genomeTmers[GID] = None
        
    def addGenomeTmer(self, GID, tmer):
        """add a tmer for a genome"""
        self.genomeTmers[GID] = tmer

    def getClosestGID(self):
        """Calculate the kmer score
        
        returns (score, (closestGID, dist), (furthestGID, dist))
        """
        LGTs = self.lgtGenomes.keys()
        dgs = []
        for lgt in LGTs:
            dg1 = pdist([self.genomeTmers[self.lgtGenomes[lgt][0]], self.lgtTmer ])
            dg2 = pdist([self.genomeTmers[self.lgtGenomes[lgt][1]], self.lgtTmer ]) 
        score = dg1/(dg1+dg2)
        return score
        #if score >= 0.5:
        #    return (score, (self.lgtGenomes[lgt][1], dg2), (self.lgtGenomes[lgt][0], dg1))
        #return (score, (self.lgtGenomes[lgt][0], dg1), (self.lgtGenomes[lgt][1], dg2))


    def __str__(self):
        """print function"""
        return "GID1: %s GID2: %s" % (i for i in self.genomeTmers.keys())                                              

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

def printDict(dict):
    for key in dict.keys():
        print "\t".join([str(key),str(dict[key])])

def printHeader():
    print "\t".join(["score","instances"])

def getIDs(header):
    genome_1= header.rstrip().split("-")[2].split(":")[1] 
    genome_2= header.rstrip().split("-")[6].split(":")[1]
    lgt_id= header.rstrip().split("-")[0]
    return (lgt_id,genome_1,genome_2)


def doWork( args ):
    """ Main wrapper"""
    LGT_dict = {}       # hash of LGTInfoStore objects
    G2L_dict = {}       # GID to LGT_dict keys
    Dist_dict = {}      # hash of rounded dist scores
    
    listing = glob.glob('%s/*' % args.kmers_directory) # list of directories
    #print listing
    
    with open(args.lgts, 'r') as lgt_fh:
        tmp_array = []
        (LGT_id, GID1, GID2) = (None, None, None) 
        for l in lgt_fh:
            fields = l.rstrip().split("\t")
            # get IDs
            if l[0:2] == "ID":
                pass
            else:
                (LGT_id, GID1, GID2) = getIDs(fields[0])
                tmp_array.append([float(i) for i in fields[1:]])
        
        lgt_tmer = np.mean(tmp_array, axis=0) # get LGT tmer
        LGT_dict[LGT_id] = LGTInfoStore(lgt_tmer) # Call class functions
        LGT_dict[LGT_id].addLGT(LGT_id, GID1, GID2) # Add LGT
        LGT_dict[LGT_id].addGenome(GID1) # add genome
        LGT_dict[LGT_id].addGenome(GID2) # add genome
        
        #-----
        """genome1"""
        with open(args.genome1,'r') as g1_fh:
            g1_tmp_array = []
            for l in g1_fh:
                fields = l.rstrip().split("\t")
                if l[0:2] == "ID":
                    pass
                else:
                    g1_tmp_array.append([float(i) for i in fields[1:]]) # misses first element 
            g1_tmer = np.mean(g1_tmp_array, axis=0)
            LGT_dict[LGT_id].addGenomeTmer(GID1,g1_tmer)
        
        #-----
        """genome2"""
        with open(args.genome2,'r') as g2_fh:
            g2_tmp_array = []
            for l in g2_fh:
                fields = l.rstrip().split("\t")
                if l[0:2] == "ID":
                    pass
                else:
                    g2_tmp_array.append([float(i) for i in fields[1:]]) # misses first element 
            g2_tmer = np.mean(g2_tmp_array, axis=0)
            LGT_dict[LGT_id].addGenomeTmer(GID2,g2_tmer)
        
        """round scores and add to dict"""
        rounded_score = float(np.round(LGT_dict[LGT_id].getClosestGID(),decimals=2))
        try:
            Dist_dict[rounded_score]+=1
        except KeyError:
            Dist_dict[rounded_score]=1
    printHeader()
    printDict(Dist_dict)
         
    """
    
    print LGT_dict[LGT_id]
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
    parser.add_argument('-lgts','--lgts', help="")
    parser.add_argument('-genome1','--genome1', help="")
    parser.add_argument('-genome2','--genome2', help="")
    parser.add_argument('-kmers_directory','--kmers_directory', help="")
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
