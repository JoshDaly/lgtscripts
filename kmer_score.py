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
from scipy.spatial import cdist

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
    def __init__(self, lgtTmerArray):
        self.lgtTmer = lgtTmerArray
        self.genomeTmers = {}
        self.what = None

    def addWhat(self, what):
        self.what = what

    def addGenome(self, GID):
        """add a holder for a genome but don't worry about the the tmer yet"""
        self.genomeTmers[GID] = None

    def addGenomeTmer(self, GID, tmer):
        """add a tmer for a genome"""
        self.genomeTmers[GID] = tmer

    def getClosestGID(self):
        """Calculate the score for lkdhlkdsjfhg"""
        GIDs = self.genomeTmers.keys()
        dgs = []
        for GID in GIDs:
            dgs.append(cdist(self.genomeTmers[GID], self.lgtTmer))
        dgg = cdist(self.genomeTmers[GIDs[0]], self.genomeTmers[GIDs[1]])
        
        # do magic math...
        score = whatever
        if score < ??:
            return (GIDs[0], score)
        else:
            return (GIDs[1], score)

    def __str__(self):
        """print function"""
        return "GID1: %s GID2: %s what: %d" % (self.genomeTmers.keys(), self.what)                                              

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

def getIDs(header):
    return (,,)

def doWork( args ):
    """ Main wrapper"""
    LGT_dict = {}       # hash of LGTInfoStore objects
    G2L_dict = {}       # GID to LGT_dict keys
    
    with open(args.lgts, 'r') as lgt_fh:
        tmp_array = []
        (LGT_id, GID1, GI2) = (None, None, None)
        for line in lgt_fh:
            fields = line.rstrip().split("\t")
            # get IDs
            (LGT_id, GID1, GI2) = getIDs(fields[0])
            tmp_array.append([float(i) for i in fields[1:]])
            
        # get LGT tmer
        tmer = np.mean(tmp_array, axis=0)
        # we have a new LGT_id and corresponding tmer
        LGT_dict[LGT_id] = LGTInfoStore(tmer)
        LGT_dict[LGT_id].addGenome(GID1)
        LGT_dict[LGT_id].addGenome(GID2)
            
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
    parser.add_argument('-i','--input_file', help="...")
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
