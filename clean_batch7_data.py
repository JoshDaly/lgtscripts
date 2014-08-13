#!/usr/bin/env python
###############################################################################
#
# __clean_batch7_data__.py - Description!
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

# put classes here 

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
    aniLookUp       = {}
    cleanData       = {}
    totalHitData    = {}
    metaData        = {}
    aniGroup        = {}
    group           = 0 
    
    # read in metadata file
    with open(args.metadata, 'r') as fh:
        header = fh.readline() # capture header
        for l in fh:
            tabs = l.rstrip().split("\t")
            gid                 = tabs[0] 
            bodySite            = tabs[57]
            genus               = tabs[11]
            phylum              = tabs[7]
            genomeLength        = tabs[71]
            sequencingCentre    = tabs[6]
            sequencingPlatform  = tabs[-1]
            status              = tabs[3]
            
            metaData[gid] = [bodySite, genus,phylum, genomeLength, sequencingCentre, sequencingPlatform, status]
            
    # read in ani file 
    with open(args.ani, 'r') as fh:
        for l in fh:
            tabs = l.rstrip().split("\t")
            gid1 = tabs[0]
            gid2 = tabs[1]
            ani  = tabs[2]
            try:
                aniLookUp[gid1][gid2] = ani 
                   
            except KeyError: 
                aniLookUp[gid1] = {gid2 : ani}
    
    # read in clean hit data
    with open(args.hitData, 'r') as fh:
        for l in fh:
            tabs = l.rstrip().split("\t")
            gid1 = tabs[0]
            gid2 = tabs[1]
            hits = tabs[2]
            try: 
                cleanData[gid1][gid2] = hits
            except KeyError:
                cleanData[gid1] = {gid2:hits}
                
    # read in dirty hit data
    with open(args.totalHitData, 'r') as fh:
        for l in fh:
            tabs = l.rstrip().split("\t")
            gid1 = tabs[0]
            gid2 = tabs[1]
            hits = int(tabs[2])
            if hits>0:
                try: 
                    print "\t".join([gid1,
                                     metaData[gid1][0],
                                     metaData[gid1][1],
                                     metaData[gid1][2],
                                     metaData[gid1][3],
                                     metaData[gid1][4],
                                     metaData[gid1][5],
                                     metaData[gid1][6],
                                     gid2, 
                                     metaData[gid2][0],
                                     metaData[gid2][1],
                                     metaData[gid2][2],
                                     metaData[gid2][3],
                                     metaData[gid2][4],
                                     metaData[gid2][5],
                                     metaData[gid2][6],
                                     str(hits),
                                     str(cleanData[gid1][gid2])
                                     ])
                except KeyError:
                    print "\t".join([gid1,
                                     metaData[gid1][0],
                                     metaData[gid1][1],
                                     metaData[gid1][2],
                                     metaData[gid1][3],
                                     metaData[gid1][4],
                                     metaData[gid1][5],
                                     metaData[gid1][6],
                                     gid2, 
                                     metaData[gid2][0],
                                     metaData[gid2][1],
                                     metaData[gid2][2],
                                     metaData[gid2][3],
                                     metaData[gid2][4],
                                     metaData[gid2][5],
                                     metaData[gid2][6],
                                     str(hits),
                                     str(0)
                                     ])
            
    return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-a','--ani', help="...")
    parser.add_argument('-hit','--hitData', help="...")
    parser.add_argument('-m','--metadata', help="...")
    parser.add_argument('-t','--totalHitData', help="...")
    #parser.add_argument('input_file2', help="gut_img_ids")
    #parser.add_argument('input_file3', help="oral_img_ids")
    #parser.add_argument('input_file4', help="ids_present_gut_and_oral.csv")
    #parser.add_argument('output_file', help="output file")
    #parser.add_argument('positional_arg3', nargs='+', help="Multiple values")
    #parser.add_argument('-X', '--optional_X', action="store_true", default=False, help="flag")
    #parser.add_argument('-X', '--optional_X', action="store_true", type=int,default=False, help="flag")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
