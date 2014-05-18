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
import glob

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

class NucMerParser:
    """Wrapper class for parsing nucmer output"""
    # constants to make the code more readable
    _START_1  = 0
    _END_1    = 1
    _START_2  = 2
    _END_2    = 3
    _LEN_1    = 4
    _LEN_2    = 5
    _IDENTITY = 6
    _ID_1     = 7
    _ID_2     = 8

    def __init__(self):
        self.prepped = False

    def reset(self):
        self.prepped = False

    def readNuc(self, fp):
        """Read through a nucmer coords file

        this is a generator function
        """
        line = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not self.prepped:
                # we still need to strip out the header
                    for l in fp: # search for the first record
                        
                        if l[0] == '=': # next line is good
                            self.prepped = True
                            break
            # file should be prepped now
            for l in fp:
                fields = l.split('|')
                yield ([int(i) for i in fields[0].split()] +
                       [int(i) for i in fields[1].split()] +
                       [int(i) for i in fields[2].split()] +
                       [float(i) for i in fields[3].split()] +
                       fields[4].split())
            break # done!

###############################################################################
###############################################################################
###############################################################################
###############################################################################

  
class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readFasta(self, fp): # this is a generator function
        header = None
        seq = None
        while True:
            for l in fp:
                if l[0] == '>': # fasta header line
                    if header is not None:
                        # we have reached a new sequence
                        yield header, "".join(seq)
                    header = l.rstrip()[1:].partition(" ")[0] # save the header we just saw
                    seq = []
                else:
                    seq.append(l.rstrip())
            # anything left in the barrel?
            if header is not None:
                yield header, "".join(seq)
            break  
  
  
  

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

    NP = NucMerParser()
    
    #Need to include glob command to parse MANY coords files
    listing = glob.glob('%s/*.coords' % args.coords_dir)
    # tracking variables
    hits = 0 # total number of hits between the genome
    length = 0 # total length of DNA shared between genomes
    interaction = 0 # variable to indicate whether it has data
       
    #print file header
    print "\t".join(["isolate_id_a",
                     "isolate_id_b",
                     "hits", 
                     "length",
                     ])   
    
    """    Script for calculating the number of hits between genomes
           and the total length of shared genes*
    """
    for c_file in listing:
        #print c_file 
        # split the img IDs from the file name in the form ../gut_oral_all_v_all/2500069000v637000103.gen.coords
        isolate_id_a = c_file.split('/')[1].split('v')[0]
        isolate_id_b = c_file.split('/')[1].split('v')[1]

        # open the file as fh
        with open(c_file, "r") as fh:
            # apply NucmerParers (NP) function readNuc on fh
            for hit in NP.readNuc(fh):
                # does it satisfy the criteria?
                try:
                    if (hit[NP._IDENTITY] >= 99 
                        and hit[NP._LEN_1] >= 500 
                        and hit[NP._LEN_2] >= 500
                        ):
                                hits += 1
                                length = (hit[NP._LEN_1] + hit[NP._LEN_2]) / 2 + length\
                                #print "Yar, she be dirty"
                            #else:
                                #print "Yep, theys be dirty!" + [hit[NP._ID_1]] + "\t" + [hit[NP._ID_2]]                           
                except KeyError:
                    pass
            
               
            print "\t".join([isolate_id_a,
                             isolate_id_b,
                             str(hits),
                             str(length),
                             ])
            
            # reset hits and length to zero
            hits = 0
            length = 0
            # Nucmer parser needs to be reset to read in new nucmer coords file
            NP.reset()
    
    

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
    parser.add_argument('-cd','--coords_dir', help="directory containing nucmer coords files")
    #parser.add_argument('-cf','--contigs_file', help="file containing contig length information")
    #parser.add_argument('-d','--dirty_contigs', help="file containing dirty contig IDs")
    #parser.add_argument('-g','--anno_file',help="File containing annotation information")
    #parser.add_argument('-id','--id_file',help="File containing matching genome tree and img IDs")
    #parser.add_argument('-gi','--gut_ids', help="gut_img_ids")
    #parser.add_argument('-oi','--oral_ids', help="oral_img_ids")
    #parser.add_argument('-bi','--both_ids', help="ids_present_gut_and_oral.csv")
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
