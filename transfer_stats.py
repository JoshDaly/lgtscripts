#!/usr/bin/env python
###############################################################################
#
# __transfer_stats__.py - Statistics of lgt transfer events.
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
class TransferParser:
    """Wrapper class for parsing transfer files"""
    """img_id_a        genome_tree_id_a        contig_a        contig_length_a start_a stop_a  length_a        img_id_b        genome_tree_id_b        contig_b        contig_length_b start_b stop_b  length_b"""
    #constants to make the code readable
    _IMG_ID_1       = 0
    _GT_ID_1        = 1
    _CONTIG_1       = 2
    _CONTIG_LEN_1   = 3
    _START_1        = 4
    _STOP_1         = 5
    _LEN_1          = 6
    _IMG_ID_2       = 7
    _GT_ID_2        = 8
    _CONTIG_2       = 9
    _CONTIG_LEN_2   = 10
    _START_2        = 11
    _STOP_2         = 12
    _LEN_2          = 13
    
    def __init__(self):
        self.prepped = False
    
    def readTrans(self,fh):
        
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fh: # search for the first record
                        if l[0:3] =="img": # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fh:
                fields = l.split("\t")
                yield [fields[0],
                       fields[1],
                       fields[2],
                       int(fields[3]),
                       int(fields[4]),
                       int(fields[5]),
                       int(fields[6]),
                       fields[7],
                       fields[8],
                       fields[9],
                       int(fields[10]),
                       int(fields[11]),
                       int(fields[12]),
                       int(fields[13])]
            break # done!
#class HitsLengthParser:
    
    #constants to make the code more readable
    
    
                        
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
    #objects
    genomes = {} # dictionary to store all things genome
    
    #parse transfer file
    TP=  TransferParser()
    with open(args.transfer_file,"r") as fh:
        for hit in TP.readTrans(fh):
            # add genome to dictionary
            try:
               genomes[hit[TP._IMG_ID_1]][hit[TP._IMG_ID_2]]["transfer_length"]+=[hit[TP._LEN_1]]
            except KeyError:
                try:
                    genomes[hit[TP._IMG_ID_1]][hit[TP._IMG_ID_2]]={"transfer_length":hit[TP._LEN_1]}
                except KeyError:
                    pass
    print genomes
            
            
            
            
            
    
    

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
    parser.add_argument('-transfer_file','--transfer_file', help="File containing transfer events")
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
